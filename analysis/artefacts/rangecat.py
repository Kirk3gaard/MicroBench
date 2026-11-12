#!/usr/bin/env python3
"""Stream a large remote file to stdout using parallel HTTP range requests, in order.
Lets a BAM of hundreds of GB be processed without storing it:  rangecat.py URL SIZE | samtools ...

usage: rangecat.py URL SIZE [threads=8] [chunk_mb=32]
"""
import sys, os, subprocess, time, hashlib
from concurrent.futures import ThreadPoolExecutor

url, size = sys.argv[1], int(sys.argv[2])
threads = int(sys.argv[3]) if len(sys.argv) > 3 else 8
chunk = (int(sys.argv[4]) if len(sys.argv) > 4 else 32) * 1024 * 1024

def get(start):
    """One byte range via curl (robust retries and timeouts); retried until complete."""
    end = min(size, start + chunk) - 1
    for attempt in range(30):
        try:
            # give up quickly on stalled transfers (< 200 kB/s for 15 s) and retry instead
            data = subprocess.run(["curl", "-sf", "--max-time", "180", "--speed-limit", "200000", "--speed-time", "15",
                                   "-r", f"{start}-{end}", url], capture_output=True, timeout=240).stdout
            if len(data) == end - start + 1:
                return data
        except Exception:
            pass
        time.sleep(1 + 2 * attempt)
    raise SystemExit(f"failed at byte {start}")

out = sys.stdout.buffer
md5 = hashlib.md5() if os.environ.get("MD5_OUT") else None
starts = list(range(0, size, chunk))
with ThreadPoolExecutor(threads) as ex:
    # keep at most 4*threads chunks in flight so memory stays bounded
    pending, i = [], 0
    while i < len(starts) or pending:
        while i < len(starts) and len(pending) < 4 * threads:
            pending.append(ex.submit(get, starts[i])); i += 1
        data = pending.pop(0).result()
        if md5: md5.update(data)
        try:
            out.write(data)
        except BrokenPipeError:
            break
try: out.flush()
except BrokenPipeError: pass
if md5: open(os.environ["MD5_OUT"], "w").write(md5.hexdigest() + "\n")
