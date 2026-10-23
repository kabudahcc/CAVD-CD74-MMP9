# -*- coding: utf-8 -*-
"""纯 Python 远程 tabix 客户端 v2：tbi 本地缓存 + 多线程分块抓取 + 行完整性过滤"""
import io, os, struct, zlib, urllib.request
from concurrent.futures import ThreadPoolExecutor

CACHE = r'C:\Users\98757\Desktop\生信分析\data\mr\tbi_cache'
os.makedirs(CACHE, exist_ok=True)

def _bgzf_decompress(data):
    out = []
    d = zlib.decompressobj(wbits=31)
    while data:
        try:
            out.append(d.decompress(data))
        except zlib.error:
            break
        out.append(d.flush())
        data = d.unused_data
        if not data:
            break
        d = zlib.decompressobj(wbits=31)
    return b''.join(out)

def _http_range(url, start, end, retry=3):
    for i in range(retry):
        try:
            req = urllib.request.Request(url, headers={'Range': f'bytes={start}-{end}', 'User-Agent': 'Mozilla/5.0'})
            return urllib.request.urlopen(req, timeout=60).read()
        except Exception:
            if i == retry - 1:
                raise
    return b''

def _reg2bins(beg, end):
    bins = [0]
    for start in (1, 9, 73, 585, 4681):
        shift = {1: 26, 9: 23, 73: 20, 585: 17, 4681: 14}[start]
        for k in range(start + (beg >> shift), start + ((end - 1) >> shift) + 1):
            bins.append(k)
    return bins

class TBI:
    def __init__(self, tbi_bytes):
        raw = _bgzf_decompress(tbi_bytes)
        buf = io.BytesIO(raw)
        assert buf.read(4) == b'TBI\x01', 'not a tabix index'
        rd = lambda f: struct.unpack(f, buf.read(struct.calcsize(f)))[0]
        self.n_ref = rd('<i'); self.fmt = rd('<i'); self.col_seq = rd('<i')
        self.col_beg = rd('<i'); self.col_end = rd('<i'); self.meta = rd('<i')
        self.skip = rd('<i'); l_nm = rd('<i')
        self.names = buf.read(l_nm).rstrip(b'\x00').decode().split('\x00')
        self.refs = {}
        for nm in self.names:
            n_bin = rd('<i')
            binmap = {}
            for _ in range(n_bin):
                b = rd('<I'); n_chunk = rd('<i')
                binmap[b] = [(rd('<Q'), rd('<Q')) for _ in range(n_chunk)]
            n_intv = rd('<i')
            self.refs[nm] = (binmap, [rd('<Q') for _ in range(n_intv)])

    def spans(self, chrom, beg, end):
        ref = None
        for cand in (chrom, 'chr' + chrom, chrom.replace('chr', '')):
            if cand in self.refs:
                ref = cand; break
        if ref is None:
            raise KeyError(f'{chrom} not in index; refs sample: {self.names[:8]}')
        binmap, intv = self.refs[ref]
        min_off = intv[beg >> 14] if (beg >> 14) < len(intv) else 0
        chunks = []
        for b in _reg2bins(beg, end):
            for cb, ce in binmap.get(b, []):
                if ce > min_off:
                    chunks.append((cb >> 16, (ce >> 16) + 65536))
        chunks.sort()
        merged = []
        for s, e in chunks:
            if merged and s <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        return ref, merged

def _load_tbi(url):
    fn = os.path.join(CACHE, url.rsplit('/', 1)[-1] + '.tbi')
    if not os.path.exists(fn):
        data = urllib.request.urlopen(urllib.request.Request(url + '.tbi', headers={'User-Agent': 'Mozilla/5.0'}), timeout=300).read()
        open(fn, 'wb').write(data)
    return TBI(open(fn, 'rb').read())

def tabix_query(url, chrom, beg, end, nfields=None, threads=4):
    tbi = _load_tbi(url)
    ref, spans = tbi.spans(chrom, beg, end)
    parts = [None] * len(spans)
    def work(i):
        s, e = spans[i]
        return _bgzf_decompress(_http_range(url, s, e))
    with ThreadPoolExecutor(max_workers=threads) as ex:
        for i, data in zip(range(len(spans)), ex.map(work, range(len(spans)))):
            parts[i] = data
    blob = b'\n'.join(p.rstrip(b'\n') for p in parts if p)
    lines = blob.decode('utf-8', errors='ignore').split('\n')
    out = []
    for ln in lines:
        if not ln.strip():
            continue
        if nfields and ln.count('\t') + 1 != nfields:
            continue  # 丢弃分块边缘的截断行
        f = ln.split('\t')
        try:
            pos = int(f[1]) if nfields is None else None
        except Exception:
            pos = None
        out.append(ln)
    # 区域内二次过滤由调用方做（不同文件列序不同）
    return out

def fetch_header(url):
    txt = _bgzf_decompress(_http_range(url, 0, 524287)).decode('utf-8', errors='ignore')
    for ln in txt.split('\n'):
        if ln.startswith('#') and '\t' in ln:
            return ln.lstrip('#').split('\t')
    return txt.split('\n')[0].split('\t')
