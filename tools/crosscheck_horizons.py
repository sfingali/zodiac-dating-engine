"""Cross-check the engine against JPL Horizons vectors.

Independence: Horizons serves its own SPK/ephemeris chain and its own light-time
and aberration handling.  We ask it for ICRF vectors and turn them into J2000
ecliptic longitudes ourselves, with the fixed mean obliquity of J2000, so the
comparison does not depend on Skyfield's frame code either.

Note on what the comparison shows: Horizons vectors are *geometric*, while the
engine reports *apparent* longitudes, so a residual of up to about 0.01 deg is
expected and is the known light-time and aberration difference, not an error.
Horizons may also decline to serve some of the outer bodies for the most ancient
dates; the tool reports that as "query failed" rather than hiding it.

Needs network access.
"""
import json, math, os, sys, urllib.parse, urllib.request
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from zodiac_dating.ephemeris import load_ephemeris, resolve_kernel, BODIES, DEFAULT_EPHEMERIS

KERNEL = resolve_kernel(os.environ.get("ZODIAC_DATING_KERNEL", DEFAULT_EPHEMERIS))
OBLIQUITY_J2000 = math.radians(23.439291111)
HORIZONS_IDS = {"Sun": "10", "Moon": "301", "Mercury": "199", "Venus": "299",
                "Mars": "499", "Jupiter": "599", "Saturn": "699"}


def horizons_vector(body, jd_tt):
    """ICRF position of body relative to geocentre, in AU, at JD(TT)."""
    start = f"JD{jd_tt:.6f}"
    stop = f"JD{jd_tt + 0.001:.6f}"
    params = {
        "format": "text", "COMMAND": f"'{HORIZONS_IDS[body]}'",
        "OBJ_DATA": "'NO'", "MAKE_EPHEM": "'YES'", "EPHEM_TYPE": "'VECTORS'",
        "CENTER": "'500@399'", "START_TIME": f"'{start}'", "STOP_TIME": f"'{stop}'",
        "STEP_SIZE": "'1 d'", "REF_SYSTEM": "'ICRF'", "REF_PLANE": "'FRAME'",
        "VEC_TABLE": "'2'", "OUT_UNITS": "'AU-D'", "CSV_FORMAT": "'YES'",
    }
    url = "https://ssd.jpl.nasa.gov/api/horizons.api?" + urllib.parse.urlencode(params)
    with urllib.request.urlopen(url, timeout=60) as r:
        text = r.read().decode()
    inside = text.split("$$SOE")[1].split("$$EOE")[0].strip().splitlines()
    fields = [f.strip() for f in inside[0].split(",")]
    # CSV layout: JDTDB, Calendar, X, Y, Z, VX, VY, VZ, ...
    x, y, z = float(fields[2]), float(fields[3]), float(fields[4])
    return np.array([x, y, z])


def j2000_ecliptic_longitude(v):
    x, y, z = v
    # rotate about the x axis by +obliquity to go ICRF -> ecliptic J2000
    y2 = y * math.cos(OBLIQUITY_J2000) + z * math.sin(OBLIQUITY_J2000)
    return math.degrees(math.atan2(y2, x)) % 360.0


def main():
    eph = load_ephemeris(KERNEL)
    for jd in (2086302.5, 2437000.5):
        print(f"\n=== JD(TT) {jd} ===")
        ours = eph.longitudes(np.array([jd]))
        print(f"{'body':<9}{'engine':>11}{'Horizons':>11}{'diff (deg)':>12}")
        diffs = []
        for body in BODIES:
            try:
                v = horizons_vector(body, jd)
                h = j2000_ecliptic_longitude(v)
            except Exception as exc:
                print(f"{body:<9}{'-':>11}{'-':>11}   query failed: {exc}")
                continue
            o = float(ours[body][0])
            d = (o - h + 180) % 360 - 180
            diffs.append(abs(d))
            print(f"{body:<9}{o:>11.4f}{h:>11.4f}{d:>12.5f}")
        if diffs:
            print(f"max |diff| = {max(diffs):.5f} deg")


if __name__ == "__main__":
    main()
