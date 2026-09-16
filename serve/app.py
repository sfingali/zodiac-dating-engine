#!/usr/bin/env python3
"""
Positions API for the zodiac dating engine's sky view.

The static page cannot compute anything: Skyfield is Python and the ephemeris is
1.5 GB.  So this serves positions and nothing else, over a tiny JSON surface:

    GET /api/health                      is it up, and which kernel
    GET /api/range                       the kernel's JD(TT) coverage
    GET /api/sky?jd=...                  seven bodies at one instant
    GET /api/sky?year=&month=&day=&hour=&minute=&second=
    GET /api/track?jd=&span_days=&n=     seven bodies on a grid, for drawing trails

Everything else the page draws - the constellation sectors - comes with /api/sky,
taken from the same derived IAU table the dating engine uses, so the wheel on the
page cannot drift from the numbers in the reports.

Time on the wire is Julian Date on the TT scale; a calendar date in the query is
read in the Julian calendar (which is what the engine prints and what a document
of the period means), converted with Skyfield's Delta T model.  Both the civil
input and the TT argument come back in the response, and so does Delta T, because
hiding that conversion is the one thing this project refuses to do.
"""
import os
import sys
import threading
from pathlib import Path

import numpy as np
from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware

REPO = Path(os.environ.get("ZODIAC_REPO", "/opt/data/repos/zodiac-dating-engine"))
sys.path.insert(0, str(REPO))

from zodiac_dating.boundaries import iau_set                      # noqa: E402
from zodiac_dating.engine import _jd_to_calendar                  # noqa: E402
from zodiac_dating.ephemeris import BODIES, load_ephemeris        # noqa: E402

MAX_TRACK_POINTS = 1200
_lock = threading.Lock()
_cache = {}

# A sky view that cannot show today is not a sky view.  DE441 part 1 covers
# -13200 to 1969; DE440s (32 MB) carries 1849 to 2150.  Both are kept loaded and
# the right one is chosen per request, so "any date" means any date either kernel
# covers.  Beyond both, the request fails loudly rather than extrapolating.
KERNEL_NAMES = ("de441_part-1.bsp", "de440s.bsp")


def _load(name):
    """Load once, keep forever: parsing these SPK files is not a per-request cost."""
    with _lock:
        if name not in _cache:
            cache = Path(os.environ.get("ZODIAC_DATING_CACHE")
                         or "~/.cache/zodiac_dating").expanduser()
            path = cache / name
            _cache[name] = load_ephemeris(str(path) if path.exists() else name)
            _cache.setdefault("_boundaries", iau_set())
        return _cache[name]


def eph_for(jd_tt):
    """The loaded kernel that covers this instant, preferring DE441 part 1 where
    the ranges overlap — it is the reference ephemeris for the epochs this project
    actually dates."""
    candidates = []
    for name in KERNEL_NAMES:
        e = _load(name)
        lo, hi = e.range_years()
        candidates.append((name, e, lo, hi))
    for want in ("de441_part-1.bsp", None):
        for name, e, lo, hi in candidates:
            if lo <= jd_tt <= hi and (want is None or name == want):
                return e, _cache["_boundaries"]
    spans = ", ".join(f"{n.split('.')[0]} {lo:.0f}..{hi:.0f}"
                      for n, _e, lo, hi in candidates)
    raise HTTPException(400, f"instant JD(TT) {jd_tt:.1f} is covered by no loaded "
                             f"kernel ({spans})")


def eph():
    """The default kernel, for parsing dates and reporting coverage."""
    return _load(KERNEL_NAMES[0]), _cache.setdefault("_boundaries", iau_set())


def calendar_string(jd_tt):
    y, m, d, hh, mm, ss = _jd_to_calendar(jd_tt)
    sign = "-" if y < 0 else "+"
    return f"{sign}{abs(y):04d}-{m:02d}-{d:02d} {hh:02d}:{mm:02d}:{ss:02d}"


def jd_from_calendar(year, month, day, hour, minute, second):
    """Civil (UT) calendar date -> TT Julian date, via Skyfield's Delta T model."""
    ts = eph()[0].ts
    t = ts.utc(year, month, day, hour, minute, second)
    return t


app = FastAPI(title="zodiac dating engine - positions", docs_url=None, redoc_url=None)
app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_methods=["GET"])


@app.get("/api/health")
def health():
    kernels = []
    for name in KERNEL_NAMES:
        e = _load(name)
        lo, hi = e.range_years()
        kernels.append({"name": name, "jd_start": lo, "jd_end": hi})
    _e, bs = eph()
    return {"ok": True, "kernels": kernels, "bodies": list(BODIES),
            "boundary_set": bs.name, "boundary_hash": bs.digest()}


@app.get("/api/range")
def range_():
    return {"kernels": [{"name": n, "jd_start": _load(n).range_years()[0],
                         "jd_end": _load(n).range_years()[1]} for n in KERNEL_NAMES]}


def _resolve_instant(jd, year, month, day, hour, minute, second):
    """Either an explicit JD(TT), or a civil calendar date to convert."""
    e, _ = eph()
    if jd is not None:
        jd_tt = float(jd)
        dt = float(e.ts.tt_jd(jd_tt).delta_t)
        return jd_tt, dt, {"given": "jd_tt", "value": jd_tt}
    if year is None:
        raise HTTPException(400, "give either jd, or year/month/day")
    month = 1 if month is None else month
    day = 1 if day is None else day
    hour = 0 if hour is None else hour
    minute = 0 if minute is None else minute
    second = 0.0 if second is None else second
    t = jd_from_calendar(year, month, day, hour, minute, second)
    jd_tt = float(t.tt)
    return jd_tt, float(t.delta_t), {"given": "calendar", "year": year, "month": month,
                                     "day": day, "hour": hour, "minute": minute,
                                     "second": second, "calendar": "julian"}


@app.get("/api/sky")
def sky(jd: float = Query(None), year: int = Query(None), month: int = Query(None),
        day: int = Query(None), hour: int = Query(None), minute: int = Query(None),
        second: float = Query(None), lat: int = Query(0)):
    e, bs = eph()
    jd_tt, delta_t, given = _resolve_instant(jd, year, month, day, hour, minute, second)
    e, bs = eph_for(jd_tt)
    try:
        lons = e.longitudes(np.array([jd_tt]))
    except Exception as exc:
        raise HTTPException(400, str(exc))
    t = e.ts.tt_jd(jd_tt)
    lats = {}
    if lat:
        lats = e.latitudes(np.array([jd_tt]))
    bodies = {}
    for b in BODIES:
        lon = float(lons[b][0])
        entry = {"lon": round(lon, 4), "sector": bs.sector_name(lon)}
        if lat:
            entry["lat"] = round(float(lats[b][0]), 4)
        bodies[b] = entry
    return {
        "jd_tt": jd_tt,
        "jd_ut1": jd_tt - delta_t / 86400.0,
        "delta_t_s": round(delta_t, 1),
        "tt_calendar": calendar_string(jd_tt),
        "given": given,
        "bodies": bodies,
        "boundaries": [{"lon": round(l, 6), "name": n} for l, n in bs.crossings],
        "boundary_set": bs.name, "boundary_hash": bs.digest(),
        "kernel": e.filename,
    }


@app.get("/api/track")
def track(jd: float = Query(None), span_days: float = Query(365.25),
          n: int = Query(240), lat: int = Query(0)):
    """A grid of instants around an instant, so the page can draw where each body
    has been and where it is going.  Retrograde arcs show up as loops."""
    e, _ = eph()
    if jd is None:
        raise HTTPException(400, "track needs jd")
    n = max(2, min(int(n), MAX_TRACK_POINTS))
    span = float(span_days)
    if not (0 < span <= 3652500):
        raise HTTPException(400, "span_days out of range")
    jds = float(jd) - span / 2.0 + (span * np.arange(n) / (n - 1))
    e, _ = eph_for(float(jds[len(jds) // 2]))
    try:
        lons = e.longitudes(jds)
    except Exception as exc:
        raise HTTPException(400, str(exc))
    out = {"jd_start": float(jds[0]), "jd_end": float(jds[-1]), "n": n,
           "jd_centre": float(jd), "span_days": span, "bodies": {}}
    for b in BODIES:
        out["bodies"][b] = [round(float(x), 3) for x in lons[b]]
    if lat:
        lats = e.latitudes(jds)
        for b in BODIES:
            out["bodies"][b + "_lat"] = [round(float(x), 3) for x in lats[b]]
    return out
