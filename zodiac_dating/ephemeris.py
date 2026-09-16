"""
Geocentric apparent ecliptic longitudes (J2000) for the seven classical bodies.

No fudges, no tolerance, no hidden offsets live here: this module only answers
"where was this body, geocentrically, in the J2000 ecliptic frame, at this
instant".  Everything else (boundaries, intervals, matching) is in the other
modules and is explicit in the input file.

Ephemeris: JPL DE441/DE440 (public domain, Park et al. 2021) through Skyfield
(Rhodes, MIT).  Skyfield also supplies the Delta-T (TT-UT1) model used when a
calendar instant is converted to a time argument, and that model is the one
source of uncertainty this engine cannot remove; see REFERENCES.md.
"""
import os
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np

BODIES = ("Sun", "Moon", "Mercury", "Venus", "Mars", "Jupiter", "Saturn")

# DE441 ships in two parts: part 1 covers -13200 to 1969, which is the span a
# historical dating engine needs; part 2 continues to the present.  DE440s and
# DE421 are the small all-in-one files, good for modern dates.  Whatever is
# actually loaded is printed, and the range is checked before any search runs.
DEFAULT_EPHEMERIS = "de441_part-1.bsp"

# Downloads land in this directory (override with ZODIAC_DATING_CACHE).
DEFAULT_CACHE = "~/.cache/zodiac_dating"
# Existing kernels are also looked for here(s) before any download is attempted.
KERNEL_DIRS_ENV = "ZODIAC_DATING_KERNEL_DIR"
DEFAULT_KERNEL_DIRS = ("~/.local/share/zodiac_dating", "/usr/local/share/zodiac_dating",
                       "/opt/data/fomenko_check/sescc")

EPHEMERIS_URLS = {
    "de441_part-1.bsp":
        "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de441_part-1.bsp",
    "de441_part-2.bsp":
        "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de441_part-2.bsp",
    "de440s.bsp": "https://ssd.jpl.nasa.gov/ftp/eph/planets/bsp/de440s.bsp",
    "de421.bsp": "https://ssd.jpl.nasa.gov/ftp/eph/planets/bsp/de421.bsp",
}


class EphemerisRangeError(RuntimeError):
    """Raised when a requested instant falls outside the loaded ephemeris."""


def cache_dir() -> Path:
    return Path(os.environ.get("ZODIAC_DATING_CACHE") or DEFAULT_CACHE).expanduser()


def kernel_dirs():
    dirs = [Path(d).expanduser() for d in DEFAULT_KERNEL_DIRS]
    extra = os.environ.get(KERNEL_DIRS_ENV, "")
    return [Path(d).expanduser() for d in extra.split(os.pathsep) if d] + dirs


def resolve_kernel(name_or_path: str) -> str:
    """Find an ephemeris kernel without downloading it if one is already here.

    Search order, and nothing else:

    1. the argument, if it is an existing path (absolute or relative);
    2. the cache directory (see ``cache_dir``);
    3. any directory in ``ZODIAC_DATING_KERNEL_DIR``, then the defaults in
       ``DEFAULT_KERNEL_DIRS``;
    4. otherwise the argument is returned unchanged and Skyfield will treat it
       as a URL to fetch into the cache.

    Nothing is silently substituted: the caller is told which file was loaded.
    """
    p = Path(name_or_path).expanduser()
    if p.exists():
        return str(p)
    name = p.name
    for d in [cache_dir()] + kernel_dirs():
        cand = d / name
        if cand.exists():
            return str(cand)
    return name_or_path


@dataclass
class Ephemeris:
    ts: object
    eph: object
    filename: str              # the name as requested
    path: str                  # where it was actually read from

    @classmethod
    def load(cls, filename: str = DEFAULT_EPHEMERIS):
        from skyfield.api import Loader
        ts = Loader(str(cache_dir())).timescale()
        ts.julian_calendar_cutoff = -10**9   # always the Julian calendar, no silent switch
        resolved = resolve_kernel(filename)
        loader = Loader(str(cache_dir()))
        if os.path.exists(resolved):
            eph = loader(resolved)
            path = os.path.abspath(resolved)
        else:
            url = EPHEMERIS_URLS.get(resolved)
            if url is None:
                raise FileNotFoundError(
                    f"ephemeris {filename!r} is not on disk and is not a known "
                    f"JPL kernel name; known names: {sorted(EPHEMERIS_URLS)}")
            eph = loader(url)                 # downloads into the cache directory
            path = str(cache_dir() / Path(url).name)
        return cls(ts=ts, eph=eph, filename=filename, path=path)

    def _body(self, name):
        key = {"Sun": "sun", "Moon": "moon", "Mercury": "mercury barycenter",
               "Venus": "venus barycenter", "Mars": "mars barycenter",
               "Jupiter": "jupiter barycenter", "Saturn": "saturn barycenter"}[name]
        if key not in self.eph:
            key = key.split()[0]
        return self.eph[key]

    def longitudes(self, jd_tt: np.ndarray) -> dict:
        """Geocentric apparent J2000 ecliptic longitude, degrees, per body.

        jd_tt: array of Julian dates on the TT scale.
        """
        from skyfield.framelib import ecliptic_J2000_frame
        t = self.ts.tt_jd(np.asarray(jd_tt, dtype=float))
        earth = self.eph["earth"]
        out = {}
        for name in BODIES:
            try:
                astrometric = earth.at(t).observe(self._body(name))
                _, lon, _ = astrometric.frame_latlon(ecliptic_J2000_frame)
                out[name] = np.asarray(lon.degrees) % 360.0
            except Exception as exc:                       # pragma: no cover
                raise EphemerisRangeError(
                    f"{name} could not be computed for the requested instants "
                    f"with {self.filename}: {exc}") from exc
        return out

    def range_years(self):
        return self.eph.spk.segments[0].start_jd, self.eph.spk.segments[0].end_jd

    def describe(self) -> str:
        lo, hi = self.range_years()
        return (f"{self.filename} (read from {self.path}); "
                f"JD(TT) {lo:.1f} .. {hi:.1f}")


@lru_cache(maxsize=4)
def load_ephemeris(filename: str = DEFAULT_EPHEMERIS) -> Ephemeris:
    return Ephemeris.load(filename)
