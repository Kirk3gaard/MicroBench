"""Sample read-length/quality stats from a remote unaligned BGZF BAM via HTTP range requests."""
import zlib, struct, re, urllib.request, math, random

def get_range(url, start, end, tries=4):
    import time
    for k in range(tries):
        try:
            return _get_range(url, start, end)
        except Exception:
            if k == tries - 1: raise
            time.sleep(3 * (k + 1))

def _get_range(url, start, end):
    req = urllib.request.Request(url, headers={"Range": f"bytes={start}-{end}"})
    with urllib.request.urlopen(req, timeout=120) as r:
        return r.read()

def file_size(url):
    req = urllib.request.Request(url, method="HEAD")
    with urllib.request.urlopen(req, timeout=60) as r:
        return int(r.headers["Content-Length"])

def inflate_bgzf(buf):
    """Decompress consecutive BGZF blocks starting at the first block found in buf."""
    m = re.search(b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00", buf)
    if not m: return b""
    i = m.start(); out = []
    while i + 18 <= len(buf):
        if buf[i:i+4] != b"\x1f\x8b\x08\x04": break
        bsize = struct.unpack("<H", buf[i+16:i+18])[0] + 1
        if i + bsize > len(buf): break
        try: out.append(zlib.decompress(buf[i+18:i+bsize-8], -15))
        except zlib.error: break
        i += bsize
    return b"".join(out)

def header(url):
    for size in (262144, 2097152, 8388608):
        raw = inflate_bgzf(get_range(url, 0, size - 1))
        if len(raw) >= 8:
            l = struct.unpack("<i", raw[4:8])[0]
            if len(raw) >= 8 + l:
                return raw[8:8+l].decode(errors="replace")
    raise ValueError("could not read BAM header")

UUID = re.compile(rb"[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}\x00")

def parse_records(data):
    """Find first record by heuristic, then walk records. Yields (name, length, meanQ)."""
    for m in UUID.finditer(data):
        s = m.start() - 36  # record start (block_size field)
        if s < 0: continue
        if data[s+12] != 37 or struct.unpack("<i", data[s+4:s+8])[0] != -1: continue
        if struct.unpack("<i", data[s+24:s+28])[0] != -1: continue
        break
    else:
        return
    i = s
    while i + 36 <= len(data):
        bs, ref, pos, lrn = struct.unpack("<iiiB", data[i:i+13])
        ncig = struct.unpack("<H", data[i+16:i+18])[0]
        lseq = struct.unpack("<i", data[i+20:i+24])[0]
        if ref != -1 or lrn > 255 or bs <= 0: return
        end = i + 4 + bs
        if end > len(data): return
        q0 = i + 36 + lrn + 4*ncig + (lseq+1)//2
        qual = data[q0:q0+lseq]
        if lseq > 0 and qual[0] != 255:
            # accurate mean Q via error probability
            pe = sum(PE[q] for q in qual) / lseq
            mq = -10*math.log10(pe)
        else:
            mq = None
        yield data[i+36:i+36+lrn-1], lseq, mq
        i = end

PE = [10**(-q/10) for q in range(256)]

def sample(url, n_chunks=40, chunk=1_000_000, seed=1):
    size = file_size(url)
    lens, qs = [], []
    starts = [int(size * (k + 0.5) / n_chunks) for k in range(n_chunks)]
    for st in starts:
        data = inflate_bgzf(get_range(url, st, min(size-1, st + chunk)))
        for _, l, q in parse_records(data):
            if l > 0:
                lens.append(l)
                if q is not None: qs.append(q)
    return size, lens, qs

def n50(lens):
    s = sorted(lens, reverse=True); half = sum(s)/2; c = 0
    for l in s:
        c += l
        if c >= half: return l
    return 0


def full_pass(url, size, step=8_000_000, overlap=2_000_000, threads=4):
    """Read the whole BAM (for small files) and return exact per-read lengths and qualities."""
    from concurrent.futures import ThreadPoolExecutor
    def job(s):
        return list(parse_records(inflate_bgzf(get_range(url, s, min(size - 1, s + step + overlap)))))
    seen = {}
    with ThreadPoolExecutor(threads) as ex:
        for recs in ex.map(job, range(0, size, step)):
            for name, l, q in recs:
                seen[name] = (l, q)
    lens = [l for l, q in seen.values() if l > 0]
    qs = [q for l, q in seen.values() if l > 0 and q is not None]
    return lens, qs
