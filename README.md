# Photo geolocator

This Python 3.13 entry fetches the trial's 25 photos, estimates broad regions
from image pixels, and submits one coordinate answer for each photo.

- A pure Python baseline JPEG decoder reconstructs grayscale or RGB images.
- Box downsampling and simple color, terrain, sky, and edge features feed
  transparent hand-written region rules.
- SHA-256 seeded jitter makes each region guess deterministic; failed photos
  receive a deterministic global fallback.

Set `MAKEFUN_PROMPT_FILE`, `MAKEFUN_ENV_URL`, and optionally
`MAKEFUN_TRIAL_SEED`, then run `python3 main.py`.

The entry makes no model calls and needs no API key. Licensed under MIT;
see `LICENSE`.
