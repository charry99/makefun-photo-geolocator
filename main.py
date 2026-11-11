import urllib.request
import urllib.error
import urllib.parse
import json
import hashlib
import math
import struct
import sys
import os

class JPEGError(Exception):
    pass

ZIGZAG = [0,1,8,16,9,2,3,10,17,24,32,25,18,11,4,5,12,19,26,33,40,48,41,34,27,20,13,6,7,14,21,28,35,42,49,56,57,50,43,36,29,22,15,23,30,37,44,51,58,59,52,45,38,31,39,46,53,60,61,54,47,55,62,63]

class EntropyReader:
    def __init__(self, data, pos):
        self.data, self.pos, self.bits, self.nbits = data, pos, 0, 0
        self.marker = None
    def byte(self):
        if self.pos >= len(self.data):
            raise JPEGError('truncated entropy data')
        b = self.data[self.pos]
        self.pos += 1
        if b != 255:
            return b
        if self.pos >= len(self.data):
            raise JPEGError('truncated marker')
        nxt = self.data[self.pos]
        if nxt == 0:
            self.pos += 1
            return 255
        self.marker = nxt
        raise JPEGError('marker encountered during entropy read')
    def bit(self):
        if self.nbits == 0:
            self.bits = self.byte()
            self.nbits = 8
        self.nbits -= 1
        return (self.bits >> self.nbits) & 1
    def value(self, n):
        v = 0
        for _ in range(n):
            v = (v << 1) | self.bit()
        return v
    def align(self):
        self.bits = self.nbits = 0
    def restart(self, expected):
        self.align()
        if self.pos >= len(self.data) or self.data[self.pos] != 255:
            raise JPEGError('expected restart marker RST%d' % expected)
        while self.pos < len(self.data) and self.data[self.pos] == 255:
            self.pos += 1
        if self.pos >= len(self.data):
            raise JPEGError('truncated restart marker')
        got = self.data[self.pos]
        self.pos += 1
        if got != 0xD0 + expected:
            raise JPEGError('expected restart marker RST%d, got FF%02X' % (expected, got))

def build_huffman(counts, symbols, cls, ident):
    table = {}
    code = 0
    k = 0
    for length, count in enumerate(counts, 1):
        for _ in range(count):
            if k >= len(symbols):
                raise JPEGError('short Huffman symbols')
            table[(length, code)] = symbols[k]
            code += 1
            k += 1
        code <<= 1
    if k != len(symbols):
        raise JPEGError('extra Huffman symbols')
    return (table, 'DC' if cls == 0 else 'AC', ident)

def huff_decode(br, ht):
    table, cls, ident = ht
    code = 0
    seen = ''
    for length in range(1, 17):
        b = br.bit()
        seen += str(b)
        code = (code << 1) | b
        if (length, code) in table:
            return table[(length, code)]
    raise JPEGError('invalid %s Huffman code table %d bits=%s' % (cls, ident, seen))

def receive_extend(br, n):
    if n == 0:
        return 0
    v = br.value(n)
    if v < (1 << (n - 1)):
        v -= (1 << n) - 1
    return v

def idct(block):
    # Matrix form of the standard separable 8x8 inverse DCT.
    basis = [[(1 / math.sqrt(2) if u == 0 else 1) * math.cos((2*x+1)*u*math.pi/16) / 2 for u in range(8)] for x in range(8)]
    out = []
    for y in range(8):
        row = []
        for x in range(8):
            s = 0.0
            for v in range(8):
                by = basis[y][v]
                for u in range(8):
                    s += by * basis[x][u] * block[v*8+u]
            row.append(max(0, min(255, int(math.floor(s + 128.5)))))
        out.append(row)
    return out

def decode_jpeg(data: bytes):
    if not data.startswith(b'\xff\xd8'):
        raise JPEGError('missing SOI')
    pos = 2
    qtables, htables = {}, {}
    frame = None
    restart_interval = 0
    while pos < len(data):
        while pos < len(data) and data[pos] != 255:
            pos += 1
        while pos < len(data) and data[pos] == 255:
            pos += 1
        if pos >= len(data):
            break
        marker = data[pos]
        pos += 1
        if marker == 0xD9:
            break
        if marker in range(0xD0, 0xD8) or marker == 0x01:
            continue
        if pos + 2 > len(data):
            raise JPEGError('truncated segment length')
        length = (data[pos] << 8) | data[pos+1]
        if length < 2 or pos + length > len(data):
            raise JPEGError('invalid segment length')
        seg = data[pos+2:pos+length]
        pos += length
        if marker == 0xC2:
            raise JPEGError('progressive SOF2 JPEG is unsupported')
        if marker == 0xDB:
            i = 0
            while i < len(seg):
                info = seg[i]; i += 1
                precision, tid = info >> 4, info & 15
                if precision not in (0, 1) or i + 64*(precision+1) > len(seg):
                    raise JPEGError('invalid quantization table')
                vals = []
                for _ in range(64):
                    if precision:
                        vals.append((seg[i] << 8) | seg[i+1]); i += 2
                    else:
                        vals.append(seg[i]); i += 1
                qtables[tid] = [0]*64
                for z, val in enumerate(vals):
                    qtables[tid][ZIGZAG[z]] = val
        elif marker == 0xC0:
            if len(seg) < 6: raise JPEGError('short SOF0')
            precision, height, width, n = seg[0], (seg[1]<<8)|seg[2], (seg[3]<<8)|seg[4], seg[5]
            if precision != 8 or n not in (1,3) or len(seg) != 6 + 3*n:
                raise JPEGError('unsupported SOF0 layout')
            comps = {}
            i = 6
            for _ in range(n):
                cid, samp, qid = seg[i:i+3]; i += 3
                h, v = samp >> 4, samp & 15
                if not (1 <= h <= 2 and 1 <= v <= 2): raise JPEGError('unsupported sampling factor')
                comps[cid] = {'h':h, 'v':v, 'q':qid}
            frame = {'w':width, 'h':height, 'comps':comps, 'n':n}
        elif marker == 0xC4:
            i = 0
            while i < len(seg):
                info = seg[i]; i += 1
                cls, tid = info >> 4, info & 15
                if cls not in (0,1) or i+16 > len(seg): raise JPEGError('invalid DHT')
                counts = list(seg[i:i+16]); i += 16
                total = sum(counts)
                if i+total > len(seg): raise JPEGError('short DHT symbols')
                syms = list(seg[i:i+total]); i += total
                htables[(cls,tid)] = build_huffman(counts, syms, cls, tid)
        elif marker == 0xDD:
            if len(seg) != 2: raise JPEGError('invalid DRI')
            restart_interval = (seg[0]<<8)|seg[1]
        elif marker == 0xDA:
            if frame is None: raise JPEGError('SOS before SOF0')
            ns = seg[0]
            if len(seg) != 1 + 2*ns + 3: raise JPEGError('invalid SOS')
            scan = {}
            i = 1
            for _ in range(ns):
                cid, sel = seg[i:i+2]; i += 2
                if cid not in frame['comps']: raise JPEGError('unknown scan component')
                scan[cid] = (htables.get((0,sel>>4)), htables.get((1,sel&15)))
                if scan[cid][0] is None or scan[cid][1] is None: raise JPEGError('missing Huffman table')
            if seg[i:i+3] != bytes((0,63,0)): raise JPEGError('non-baseline spectral scan')
            break
    else:
        raise JPEGError('missing SOS')
    if frame is None: raise JPEGError('missing SOF0')
    w,h,comps,n = frame['w'],frame['h'],frame['comps'],frame['n']
    if w <= 0 or h <= 0: raise JPEGError('invalid image size')
    maxh=max(c['h'] for c in comps.values()); maxv=max(c['v'] for c in comps.values())
    for c in comps.values():
        if c['q'] not in qtables: raise JPEGError('missing quantization table %d' % c['q'])
    mcu_w, mcu_h = maxh*8, maxv*8
    mcu_cols=(w+mcu_w-1)//mcu_w; mcu_rows=(h+mcu_h-1)//mcu_h
    br=EntropyReader(data,pos)
    pred={cid:0 for cid in comps}
    planes={cid:[[0]*(mcu_cols*c['h']*8) for _ in range(mcu_rows*c['v']*8)] for cid,c in comps.items()}
    rst=0; mcu_count=0
    for my in range(mcu_rows):
        for mx in range(mcu_cols):
            if restart_interval and mcu_count and mcu_count % restart_interval == 0:
                br.restart(rst); rst=(rst+1)&7
                pred={cid:0 for cid in comps}
            for cid,c in comps.items():
                if cid not in scan: continue
                dc_ht,ac_ht=scan[cid]
                for by in range(c['v']):
                    for bx in range(c['h']):
                        block=[0]*64
                        cat=huff_decode(br,dc_ht)
                        if cat>11: raise JPEGError('invalid DC category %d' % cat)
                        pred[cid]+=receive_extend(br,cat)
                        block[0]=pred[cid]*qtables[c['q']][0]
                        k=1
                        while k<64:
                            sym=huff_decode(br,ac_ht); run=sym>>4; size=sym&15
                            if sym==0: break
                            if sym==0xF0:
                                k+=16
                                if k>64: raise JPEGError('bad AC run: ZRL passes coefficient 63')
                                continue
                            if size==0: raise JPEGError('bad AC run: zero size symbol')
                            k+=run
                            if k>=64: raise JPEGError('bad AC run: run passes coefficient 63')
                            z=ZIGZAG[k]
                            block[z]=receive_extend(br,size)*qtables[c['q']][z]
                            k+=1
                        pixels=idct(block)
                        ox=(mx*c['h']+bx)*8; oy=(my*c['v']+by)*8
                        plane=planes[cid]
                        for yy in range(8): plane[oy+yy][ox:ox+8]=pixels[yy]
            mcu_count+=1
    ids=list(comps)
    out=[]
    for y in range(h):
        for x in range(w):
            vals=[]
            for cid,c in comps.items():
                sx=x*c['h']//maxh; sy=y*c['v']//maxv
                vals.append(planes[cid][sy][sx])
            if n==1:
                r=g=b=vals[0]
            else:
                yy,cb,cr=vals[0], vals[1]-128, vals[2]-128
                r=max(0,min(255,int(math.floor(yy+1.402*cr+0.5))))
                g=max(0,min(255,int(math.floor(yy-0.344136*cb-0.714136*cr+0.5))))
                b=max(0,min(255,int(math.floor(yy+1.772*cb+0.5))))
            out.append((r,g,b))
    return w,h,out
def features(w,h,pix):
    n=len(pix); means=[sum(p[i] for p in pix)/n for i in range(3)]
    bright=sum(sum(p)/3 for p in pix)/n
    sat=sum((max(p)-min(p))/max(p) if max(p) else 0 for p in pix)/n
    colorful=math.sqrt(sum((p[0]-p[1])**2 for p in pix)/n + sum(((p[0]+p[1])/2-p[2])**2 for p in pix)/n)
    top=max(1,h//3); bot=2*h//3
    top_pix=[pix[y*w+x] for y in range(top) for x in range(w)]
    sky=[p for p in top_pix if p[2]>p[0]+10 and p[2]>p[1]+5]
    veg=sum(1 for r,g,b in pix if g>r+8 and g>b+8)/n
    roadpix=[pix[y*w+x] for y in range(bot,h) for x in range(w)] or pix
    road=sum(1 for r,g,b in roadpix if max(r,g,b)-min(r,g,b)<25 and 60<=(r+g+b)/3<=200)/len(roadpix)
    soil=sum(1 for r,g,b in pix if r>g>b and r-b>30)/n
    snow=sum(1 for p in pix if sum(p)/3>225 and max(p)-min(p)<20)/n
    red=sum(1 for r,g,b in pix if r>120 and r-g>40 and g-b>10)/n
    dry=sum(1 for r,g,b in pix if g>r>b and 120<=(r+g+b)/3<=210)/n
    horiz=[]; vert=0; strong=0
    for y in range(h):
        for x in range(1,w):
            a=pix[y*w+x-1]; b=pix[y*w+x]
            d=abs(sum(a)-sum(b))/3
            horiz.append(d)
            if d>25: strong+=1; vert+=1
    vstrong=0
    for y in range(1,h):
        for x in range(w):
            a=pix[(y-1)*w+x]; b=pix[y*w+x]
            if abs(sum(a)-sum(b))/3>25: vstrong+=1
    verticality=strong/(strong+vstrong) if strong+vstrong else .5
    left=[pix[y*w+x] for y in range(top) for x in range(w//2)] or top_pix
    right=[pix[y*w+x] for y in range(top) for x in range(w//2,w)] or top_pix
    sun=sum(sum(p)/3 for p in left)/len(left)-sum(sum(p)/3 for p in right)/len(right)
    return {'mean':means,'brightness':bright,'sat':sat,'colorful':colorful,'sky':len(sky)/len(top_pix),'skyblue':sum(p[2]-p[0] for p in top_pix)/len(top_pix),'veg':veg,'road':road,'soil':soil,'snow':snow,'red':red,'edge':sum(horiz)/len(horiz) if horiz else 0,'verticality':verticality,'sun':sun,'dry':dry}

REGIONS={
 'west':(45,-120,900),'east':(40,-85,750),'mexico':(20,-100,900),'brazil':(-12,-52,1100),'andes':(-30,-65,900),
 'weurope':(48,5,500),'eeurope':(52,25,750),'russia':(58,70,1200),'scandinavia':(63,15,900),'middleeast':(28,45,900),
 'northafrica':(28,10,1000),'africa':(-2,28,1100),'southasia':(22,79,900),'eastasia':(33,112,900),'seasia':(12,102,900),'oceania':(-28,140,1100),'japan':(37,132,600)}
def classify(f):
    g=f['veg']; warm=f['brightness']; dry=f['dry']; soil=f['soil']; red=f['red']; snow=f['snow']; road=f['road']
    if snow>.20:
        return 'scandinavia' if g>.12 and warm<175 else ('russia' if g<.10 else 'west')
    if red>.12 and dry<.10:
        return 'brazil' if g>.25 else ('africa' if g>.08 else 'oceania')
    if g>.30 and warm>90:
        return 'seasia' if warm>135 and f['sat']>.22 else ('brazil' if f['skyblue']<15 else 'africa')
    if soil>.16 and g<.18 and warm>125:
        return 'middleeast' if f['skyblue']>8 else 'northafrica'
    if g<.08 and warm<160:
        if snow>.05 or (road>.40 and f['skyblue']<25 and f['brightness']<130):
            return 'russia'
        if f['brightness']<100 and soil>.12:
            return 'andes'
    if g>.12 and road>.12:
        return 'weurope' if f['skyblue']>12 else ('east' if f['verticality']>.52 else 'eeurope')
    if f['sky']>.30 and g<.18:
        return 'mexico' if warm>145 else 'west'
    if f['sat']>.28 and g>.15:
        return 'southasia' if f['skyblue']<10 else 'eastasia'
    if f['verticality']>.56:
        return 'japan' if f['skyblue']>20 else 'eeurope'
    return 'eeurope'

def jitter(region, photo_id, seed):
    lat,lon,spread=REGIONS[region]
    digest=hashlib.sha256((seed+':'+photo_id).encode()).digest()
    a=int.from_bytes(digest[:4],'big')/0xffffffff*2-1
    b=int.from_bytes(digest[4:8],'big')/0xffffffff*2-1
    lat+=a*spread/111
    lon+=b*spread/(111*max(.01,math.cos(math.radians(lat))))
    return round(max(-89,min(89,lat)),6),round(max(-179,min(179,lon)),6)

def fallback(photo_id,seed):
    digest=hashlib.sha256((seed+':'+photo_id).encode()).digest()
    a=int.from_bytes(digest[:4],'big')/0xffffffff*2-1
    b=int.from_bytes(digest[4:8],'big')/0xffffffff*2-1
    lat=25+a*4000/111
    lon=20+b*4000/(111*math.cos(math.radians(lat)))
    return round(max(-89,min(89,lat)),6),round(max(-179,min(179,lon)),6)

def downsample(w,h,pix,limit=160):
    scale=max(1,math.ceil(max(w,h)/limit))
    nw=(w+scale-1)//scale; nh=(h+scale-1)//scale
    out=[]
    for y in range(nh):
        for x in range(nw):
            cells=[pix[yy*w+xx] for yy in range(y*scale,min(h,(y+1)*scale)) for xx in range(x*scale,min(w,(x+1)*scale))]
            out.append(tuple(sum(p[c] for p in cells)//len(cells) for c in range(3)))
    return nw,nh,out

def request(url, data=None):
    req=urllib.request.Request(url,data=data,headers={'Content-Type':'application/json'} if data is not None else {})
    last=None
    for attempt in range(3):
        try:
            with urllib.request.urlopen(req,timeout=60) as r: return r.read()
        except urllib.error.HTTPError as e:
            if 400<=e.code<500: raise
            last=e
        except (urllib.error.URLError,TimeoutError,ConnectionError) as e:
            last=e
        if attempt<2:
            import time
            time.sleep(2)
    raise last

def main():
    env=os.environ.get('MAKEFUN_ENV_URL','').rstrip('/')
    if not env:
        print('MAKEFUN_ENV_URL is required',file=sys.stderr); return 1
    prompt=os.environ.get('MAKEFUN_PROMPT_FILE')
    if prompt:
        try:
            with open(prompt,encoding='utf-8') as f: f.read()
        except OSError as e:
            print('prompt read failed: %s' % e,file=sys.stderr)
    seed=os.environ.get('MAKEFUN_TRIAL_SEED','')
    try:
        listing=json.loads(request(env+'/photos'))
    except Exception as e:
        print('photo list failed: %s' % e,file=sys.stderr); return 1
    answers=[]
    for item in listing['photos']:
        pid=item['id']
        try:
            url=urllib.parse.urljoin(env+'/',item['url'])
            raw=request(url)
            w,h,pix=decode_jpeg(raw)
            w,h,pix=downsample(w,h,pix)
            lat,lon=jitter(classify(features(w,h,pix)),pid,seed)
        except Exception as e:
            print('photo %s failed; using fallback: %s' % (pid,e),file=sys.stderr)
            lat,lon=fallback(pid,seed)
        answers.append({'id':pid,'lat':lat,'lon':lon})
    body=json.dumps({'answers':answers},separators=(',',':')).encode()
    try:
        result=json.loads(request(env+'/answers',body))
    except urllib.error.HTTPError as e:
        print('POST failed: %s %s' % (e.code,e.read().decode(errors='replace')),file=sys.stderr); return 1
    except Exception as e:
        print('POST failed: %s' % e,file=sys.stderr); return 1
    for a in answers: print('%s %s %s' % (a['id'],a['lat'],a['lon']))
    print('answered %d/%d' % (result.get('answered',len(answers)),len(answers)))
    return 0

if __name__=='__main__':
    sys.exit(main())
