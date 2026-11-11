# Photo Geolocator

This entry reads the trial environment's photo list, fetches each JPEG, decodes baseline JPEG data in pure Python, and derives pixel statistics from a small box-averaged image. An ordered ruleset maps those features to one of the supplied world regions. Coordinates are jittered around each region with SHA-256 values derived from the trial seed and photo ID, so repeated runs produce the same answers. If one photo cannot be processed, a deterministic fallback coordinate is submitted for it.

The program uses only Python's standard library and makes no model calls. Its only network requests are to the environment's photo and answer endpoints.
