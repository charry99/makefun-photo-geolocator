import hashlib
import json
import math
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request


ZIGZAG = [0,1,8,16,9,2,3,10,17,24,32,25,18,11,4,5,12,19,26,33,40,48,41,34,27,20,13,6,7,14,21,28,35,42,49,56,57,50,43,36,29,22,15,23,30,37,44,51,58,59,52,45,38,31,39,46,53,60,61,54,55,62,63]
REGIONS = {
    "west": (45, -120, 900), "east": (40, -85, 750), "mexico": (20, -100, 900),
    "brazil": (-12, -52, 1100), "andes": (-30, -65, 900), "weurope": (48, 5, 500),
    "eeurope": (52, 25, 750), "russia": (58, 70, 1200), "scandinavia": (63, 15, 900),
    "middleeast": (28, 45, 900), "northafrica": (28, 10, 1000), "africa": (-2, 28, 1100),
    "southasia": (22, 79, 900), "eastasia": (33, 112, 900), "seasia": (12, 102, 900),
    "oceania": (-28, 140, 1100), "japan": (37, 132, 600)
}


def http_request(url, data=None):
    req = urllib.request.Request(url, data=data, headers={"Content-Type": "application/json"} if data is not None else {})
    last = None
    for attempt in range(3):
        try:
            with urllib.request.urlopen(req, timeout=60) as response:
                return response.read()
        except urllib.error.HTTPError as exc:
            if 400 <= exc.code < 500:
                raise
            last = exc
        except Exception as exc:
            last = exc
        if attempt < 2:
            time.sleep(2)
    raise last


def receive_extend(bits, count):
    if count == 0:
        return 0
    threshold = 1 << (count - 1)
    return bits if bits >= threshold else bits - ((1 << count) - 1)


class BitReader:
    def __init__(self, data, pos):
        self.data = data
        self.pos = pos
        self.value = 0
        self.count = 0

    def bit(self):
        if self.count == 0:
            if self.pos >= len(self.data):
                raise ValueError("truncated JPEG entropy data")
            value = self.data[self.pos]
            self.pos += 1
            if value == 255:
                while self.pos < len(self.data) and self.data[self.pos] == 255:
                    self.pos += 1
                if self.pos >= len(self.data) or self.data[self.pos] != 0:
                    raise ValueError("unexpected marker in entropy data")
                self.pos += 1
            self.value = value
            self.count = 8
        self.count -= 1
        return (self.value >> self.count) & 1

    def bits(self, count):
        value = 0
        for _ in range(count):
            value = (value << 1) | self.bit()
        return value


def huffman_table(counts, symbols):
    table = {}
    code = 0
    index = 0
    for length, amount in enumerate(counts, 1):
        for _ in range(amount):
            table[(length, code)] = symbols[index]
            code += 1
            index += 1
        code <<= 1
    return table


def huffman_value(reader, table):
    code = 0
    for length in range(1, 17):
        code = (code << 1) | reader.bit()
        item = table.get((length, code))
        if item is not None:
            return item
    raise ValueError("invalid Huffman code")


def inverse_dct(coefficients):
    result = [[0] * 8 for _ in range(8)]
    for y in range(8):
        for x in range(8):
            total = 0.0
            for v in range(8):
                cv = 1 / math.sqrt(2) if v == 0 else 1.0
                for u in range(8):
                    cu = 1 / math.sqrt(2) if u == 0 else 1.0
                    total += cu * cv * coefficients[v * 8 + u] * math.cos((2 * x + 1) * u * math.pi / 16) * math.cos((2 * y + 1) * v * math.pi / 16)
            value = math.floor(total / 4 + 128.5)
            result[y][x] = max(0, min(255, value))
    return result


def decode_jpeg(data):
    raise NotImplementedError

# Correct scan order is fixed explicitly so every coefficient maps to its natural position.
ZIGZAG = [0,1,8,16,9,2,3,10,17,24,32,25,18,11,4,5,12,19,26,33,40,48,41,34,27,20,13,6,7,14,21,28,35,42,49,56,57,50,43,36,29,22,15,23,30,37,44,51,58,59,52,45,38,31,39,46,53,60,61,54,47,55,62,63]


def decode_jpeg(data):
    if data[:2] != b"\xff\xd8":
        raise ValueError("not a JPEG")
    pos = 2
    quant = {}
    dc_tables = {}
    ac_tables = {}
    components = {}
    width = height = None
    max_h = max_v = 1
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
            raise ValueError("truncated JPEG segment")
        length = (data[pos] << 8) | data[pos + 1]
        start = pos + 2
        end = pos + length
        if length < 2 or end > len(data):
            raise ValueError("invalid JPEG segment length")
        segment = data[start:end]
        pos = end
        if marker == 0xC2:
            raise ValueError("progressive JPEG is unsupported")
        if marker == 0xDB:
            i = 0
            while i < len(segment):
                info = segment[i]
                i += 1
                precision, table_id = info >> 4, info & 15
                size = 2 if precision else 1
                vals = []
                for _ in range(64):
                    vals.append(int.from_bytes(segment[i:i + size], "big"))
                    i += size
                natural = [0] * 64
                for k, value in enumerate(vals):
                    natural[ZIGZAG[k]] = value
                quant[table_id] = natural
        elif marker == 0xC4:
            i = 0
            while i < len(segment):
                info = segment[i]
                i += 1
                counts = list(segment[i:i + 16])
                i += 16
                symbols = list(segment[i:i + sum(counts)])
                i += sum(counts)
                table = huffman_table(counts, symbols)
                (dc_tables if info >> 4 == 0 else ac_tables)[info & 15] = table
        elif marker == 0xC0:
            if segment[0] != 8:
                raise ValueError("only 8-bit JPEG precision is supported")
            height = int.from_bytes(segment[1:3], "big")
            width = int.from_bytes(segment[3:5], "big")
            count = segment[5]
            if count not in (1, 3):
                raise ValueError("JPEG must have one or three components")
            i = 6
            for _ in range(count):
                cid, sampling, qid = segment[i:i + 3]
                i += 3
                hc, vc = sampling >> 4, sampling & 15
                if hc not in (1, 2) or vc not in (1, 2):
                    raise ValueError("unsupported sampling factor")
                components[cid] = {"h": hc, "v": vc, "q": qid, "dc": 0, "plane": None}
                max_h, max_v = max(max_h, hc), max(max_v, vc)
        elif marker == 0xDD:
            restart_interval = int.from_bytes(segment[:2], "big")
        elif marker == 0xDA:
            count = segment[0]
            selectors = {}
            i = 1
            for _ in range(count):
                cid, tables = segment[i:i + 2]
                i += 2
                selectors[cid] = (tables >> 4, tables & 15)
            if segment[i] != 0 or segment[i + 1] != 63 or segment[i + 2] != 0:
                raise ValueError("unsupported JPEG scan")
            if not width or not height or len(components) not in (1, 3):
                raise ValueError("missing JPEG frame")
            mcu_cols = (width + max_h * 8 - 1) // (max_h * 8)
            mcu_rows = (height + max_v * 8 - 1) // (max_v * 8)
            for comp in components.values():
                comp["plane"] = [[0] * (mcu_cols * comp["h"] * 8) for _ in range(mcu_rows * comp["v"] * 8)]
            reader = BitReader(data, pos)
            mcu_index = 0
            expected_rst = 0
            for my in range(mcu_rows):
                for mx in range(mcu_cols):
                    for cid, comp in components.items():
                        dcid, acid = selectors[cid]
                        for by in range(comp["v"]):
                            for bx in range(comp["h"]):
                                dc_size = huffman_value(reader, dc_tables[dcid])
                                comp["dc"] += receive_extend(reader.bits(dc_size), dc_size)
                                coeff = [0] * 64
                                coeff[0] = comp["dc"] * quant[comp["q"]][0]
                                k = 1
                                while k < 64:
                                    symbol = huffman_value(reader, ac_tables[acid])
                                    run, size = symbol >> 4, symbol & 15
                                    if size == 0:
                                        if run == 0:
                                            break
                                        if run != 15:
                                            raise ValueError("invalid AC run")
                                        k += 16
                                        continue
                                    k += run
                                    if k >= 64:
                                        raise ValueError("AC coefficient outside block")
                                    coeff[ZIGZAG[k]] = receive_extend(reader.bits(size), size) * quant[comp["q"]][ZIGZAG[k]]
                                    k += 1
                                block = inverse_dct(coeff)
                                ox, oy = (mx * comp["h"] + bx) * 8, (my * comp["v"] + by) * 8
                                plane = comp["plane"]
                                for yy in range(8):
                                    plane[oy + yy][ox:ox + 8] = block[yy]
                    mcu_index += 1
                    if restart_interval and mcu_index % restart_interval == 0 and mcu_index < mcu_cols * mcu_rows:
                        reader.count = 0
                        marker_pos = reader.pos
                        while marker_pos < len(data) and data[marker_pos] == 255:
                            marker_pos += 1
                        if marker_pos >= len(data) or data[marker_pos] != 0xD0 + expected_rst:
                            raise ValueError("missing restart marker")
                        reader.pos = marker_pos + 1
                        expected_rst = (expected_rst + 1) & 7
                        for comp in components.values():
                            comp["dc"] = 0
            pixels = []
            for y in range(height):
                for x in range(width):
                    samples = []
                    for comp in components.values():
                        sx = x * comp["h"] // max_h
                        sy = y * comp["v"] // max_v
                        samples.append(comp["plane"][sy][sx])
                    if len(samples) == 1:
                        rgb = (samples[0],) * 3
                    else:
                        yy, cb, cr = samples
                        red = math.floor(yy + 1.402 * (cr - 128) + 0.5)
                        green = math.floor(yy - 0.344136 * (cb - 128) - 0.714136 * (cr - 128) + 0.5)
                        blue = math.floor(yy + 1.772 * (cb - 128) + 0.5)
                        rgb = tuple(max(0, min(255, v)) for v in (red, green, blue))
                    pixels.append(rgb)
            return width, height, pixels
    raise ValueError("JPEG has no baseline scan")


def downsample(width, height, pixels):
    scale = max(1, (max(width, height) + 159) // 160)
    nw, nh = (width + scale - 1) // scale, (height + scale - 1) // scale
    result = []
    for oy in range(nh):
        for ox in range(nw):
            vals = [0, 0, 0]
            count = 0
            for y in range(oy * scale, min(height, (oy + 1) * scale)):
                for x in range(ox * scale, min(width, (ox + 1) * scale)):
                    p = pixels[y * width + x]
                    for c in range(3):
                        vals[c] += p[c]
                    count += 1
            result.append(tuple(v // count for v in vals))
    return nw, nh, result


def pixel_features(width, height, pixels):
    n = len(pixels)
    means = [sum(p[c] for p in pixels) / n for c in range(3)]
    brightnesses = [(r + g + b) / 3 for r, g, b in pixels]
    brightness = sum(brightnesses) / n
    sats = [(max(p) - min(p)) / max(p) if max(p) else 0 for p in pixels]
    sat = sum(sats) / n
    colorful = math.sqrt(sum((r - g) ** 2 for r, g, b in pixels) / n + sum(((r + g) / 2 - b) ** 2 for r, g, b in pixels) / n)
    top = [(y, x, pixels[y * width + x]) for y in range(max(1, (height + 2) // 3)) for x in range(width)]
    sky_count = sum(1 for _, _, (r, g, b) in top if b > r + 10 and b > g + 5)
    sky = sky_count / len(top)
    skyblue = sum(b - r for _, _, (r, g, b) in top) / len(top)
    veg = sum(1 for r, g, b in pixels if g > r + 8 and g > b + 8) / n
    bottom = [pixels[y * width + x] for y in range((2 * height) // 3, height) for x in range(width)]
    if not bottom:
        bottom = pixels
    road = sum(1 for p in bottom if max(p) - min(p) < 25 and 60 <= sum(p) / 3 <= 200) / len(bottom)
    soil = sum(1 for r, g, b in pixels if r > g > b and r - b > 30) / n
    snow = sum(1 for p in pixels if sum(p) / 3 > 225 and max(p) - min(p) < 20) / n
    red = sum(1 for r, g, b in pixels if r > 120 and r - g > 40 and g - b > 10) / n
    horiz = []
    vert = []
    for y in range(height):
        for x in range(width - 1):
            a, b = pixels[y * width + x], pixels[y * width + x + 1]
            horiz.append(abs(sum(a) - sum(b)) / 3)
    for y in range(height - 1):
        for x in range(width):
            a, b = pixels[y * width + x], pixels[(y + 1) * width + x]
            vert.append(abs(sum(a) - sum(b)) / 3)
    edge = sum(horiz) / len(horiz) if horiz else 0
    hc = sum(d > 25 for d in horiz)
    vc = sum(d > 25 for d in vert)
    verticality = hc / (hc + vc) if hc + vc else 0.5
    left = [brightnesses[y * width + x] for y in range(max(1, (height + 2) // 3)) for x in range(width // 2)]
    right = [brightnesses[y * width + x] for y in range(max(1, (height + 2) // 3)) for x in range(width // 2, width)]
    if not left or not right:
        left = brightnesses[:len(top)]
        right = left
    sun = sum(left) / len(left) - sum(right) / len(right)
    dry = sum(1 for r, g, b in pixels if g > r > b and 120 <= (r + g + b) / 3 <= 210) / n
    return {"means": means, "brightness": brightness, "sat": sat, "colorful": colorful, "sky": sky, "skyblue": skyblue, "veg": veg, "road": road, "soil": soil, "snow": snow, "red": red, "edge": edge, "verticality": verticality, "sun": sun, "dry": dry}


def classify(f):
    g, warm, dry, soil, red, snow, road = f["veg"], f["brightness"], f["dry"], f["soil"], f["red"], f["snow"], f["road"]
    if snow > 0.20:
        return "scandinavia" if g > 0.12 and warm < 175 else ("russia" if g < 0.10 else "west")
    if red > 0.12 and dry < 0.10:
        return "brazil" if g > 0.25 else ("africa" if g > 0.08 else "oceania")
    if g > 0.30 and warm > 90:
        return "seasia" if warm > 135 and f["sat"] > 0.22 else ("brazil" if f["skyblue"] < 15 else "africa")
    if soil > 0.16 and g < 0.18 and warm > 125:
        return "middleeast" if f["skyblue"] > 8 else "northafrica"
    if g < 0.08 and warm < 160:
        if snow > 0.05 or (road > 0.40 and f["skyblue"] < 25 and f["brightness"] < 130):
            return "russia"
        if f["brightness"] < 100 and soil > 0.12:
            return "andes"
    if g > 0.12 and road > 0.12:
        return "weurope" if f["skyblue"] > 12 else ("east" if f["verticality"] > 0.52 else "eeurope")
    if f["sky"] > 0.30 and g < 0.18:
        return "mexico" if warm > 145 else "west"
    if f["sat"] > 0.28 and g > 0.15:
        return "southasia" if f["skyblue"] < 10 else "eastasia"
    if f["verticality"] > 0.56:
        return "japan" if f["skyblue"] > 20 else "eeurope"
    return "eeurope"


def offsets(photo_id, seed):
    digest = hashlib.sha256((seed + ":" + str(photo_id)).encode()).hexdigest()
    return int(digest[:4], 16) / 0xffffffff * 2 - 1, int(digest[4:8], 16) / 0xffffffff * 2 - 1


def finish_coords(lat, lon):
    lat = max(-89, min(89, lat))
    lon = max(-179, min(179, lon))
    return round(lat, 6), round(lon, 6)


def jitter(region, photo_id, seed):
    lat, lon, spread = REGIONS[region]
    a, b = offsets(photo_id, seed)
    lat += a * spread / 111
    lon += b * spread / (111 * max(0.01, math.cos(math.radians(lat))))
    return finish_coords(lat, lon)


def fallback(photo_id, seed):
    a, b = offsets(photo_id, seed)
    lat = 25 + a * 4000 / 111
    lon = 20 + b * 4000 / (111 * math.cos(math.radians(lat)))
    return finish_coords(lat, lon)


def main():
    env = os.environ.get("MAKEFUN_ENV_URL", "").rstrip("/")
    if not env:
        print("MAKEFUN_ENV_URL is required", file=sys.stderr)
        return 1
    prompt_path = os.environ.get("MAKEFUN_PROMPT_FILE")
    if prompt_path:
        try:
            with open(prompt_path, "r", encoding="utf-8") as prompt_file:
                prompt_file.read()
        except Exception:
            pass
    seed = os.environ.get("MAKEFUN_TRIAL_SEED", "")
    try:
        listing = json.loads(http_request(env + "/photos"))
        items = listing["photos"]
    except Exception as exc:
        print(str(exc), file=sys.stderr)
        return 1
    answers = []
    for item in items:
        photo_id = item["id"]
        try:
            image_url = urllib.parse.urljoin(env + "/", item["url"])
            image_data = http_request(image_url)
            w, h, pixels = decode_jpeg(image_data)
            sw, sh, small = downsample(w, h, pixels)
            region = classify(pixel_features(sw, sh, small))
            lat, lon = jitter(region, photo_id, seed)
        except Exception:
            lat, lon = fallback(photo_id, seed)
        answers.append({"id": photo_id, "lat": lat, "lon": lon})
    payload = json.dumps({"answers": answers}, separators=(",", ":")).encode()
    try:
        http_request(env + "/answers", payload)
    except Exception as exc:
        print(str(exc), file=sys.stderr)
        return 1
    for answer in answers:
        print(f'{answer["id"]} {answer["lat"]} {answer["lon"]}')
    print(f"answered {len(answers)}/{len(items)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
