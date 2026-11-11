"""Sample complete BAM records from a remote unaligned BAM and write a local BAM (pysam-readable)."""
import sys, struct, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "scripts"))
import bamsample as bs
from concurrent.futures import ThreadPoolExecutor

def records_raw(data):
    """Yield raw record bytes (incl. block_size) starting from first plausible record."""
    import re
    for m in bs.UUID.finditer(data):
        s = m.start() - 36
        if s < 0: continue
        if data[s+12] != 37 or struct.unpack("<i", data[s+4:s+8])[0] != -1: continue
        if struct.unpack("<i", data[s+24:s+28])[0] != -1: continue
        break
    else:
        return
    i = s
    while i + 36 <= len(data):
        bs_, ref = struct.unpack("<ii", data[i:i+8])
        if ref != -1 or bs_ <= 0: return
        end = i + 4 + bs_
        if end > len(data): return
        yield data[i:end]
        i = end

def main(url, size, n_chunks, out, chunk=1_000_000, threads=4):
    raw = bs.inflate_bgzf(bs.get_range(url, 0, 4_000_000))
    l_text = struct.unpack("<i", raw[4:8])[0]
    p = 8 + l_text
    n_ref = struct.unpack("<i", raw[p:p+4])[0]; p += 4
    assert n_ref == 0
    header = raw[:p]
    starts = [int(size * (k + 0.5) / n_chunks) for k in range(n_chunks)]
    def job(k_st):
        k, st = k_st
        try:
            d = bs.inflate_bgzf(bs.get_range(url, st, min(size - 1, st + chunk)))
            return k, list(records_raw(d))
        except Exception as e:
            print("fail", k, e, file=sys.stderr); return k, []
    seen = set(); n = 0
    with open(out + ".raw", "wb") as fh, ThreadPoolExecutor(threads) as ex:
        fh.write(header)
        for k, recs in ex.map(job, enumerate(starts)):
            for r in recs:
                name = r[36:36+r[12]-1]
                if name in seen: continue
                seen.add(name); fh.write(r); n += 1
            print(f"\rchunk {k+1}/{n_chunks} reads {n}", end="", file=sys.stderr)
    print(file=sys.stderr)
    import pysam
    pysam.tabix_compress(out + ".raw", out, force=True); os.remove(out + ".raw")

if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4])
