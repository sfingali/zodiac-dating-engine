#!/usr/bin/env python3
"""
Derive the true constellation boundaries along the ecliptic from the IAU
constellation boundary data bundled with Skyfield (which encodes the Delporte
1930 / IAU 1930 boundaries as tabulated by Roman 1987, CDS VI/42).

Method: walk the ecliptic (J2000, ecliptic latitude 0) in longitude, ask
Skyfield's constellation map which constellation each point is in, and bisect
every change of constellation to high precision.

Output: JSON table of the boundaries, in J2000 ecliptic longitude, written into
the package data directory (``--out`` overrides).  Sector names are the full
constellation names, spelled as the twelve signs are spelled in the published
tables this engine is meant to be compared against (Scorpio, Capricorn), so a
specification can be run unchanged under the IAU set or under a HOROS set; each
entry also carries the IAU three-letter abbreviation.

    python tools/derive_boundaries.py
    python tools/derive_boundaries.py --out /tmp/boundaries.json
"""
import argparse
import json
from pathlib import Path

import numpy as np
from skyfield.api import load_constellation_map, position_of_radec, load
from skyfield.framelib import ecliptic_J2000_frame

# IAU three-letter abbreviation -> full name.  The twelve are spelled as in the
# published zodiac tables (Scorpio, Capricorn) rather than the Latin forms.
FULL_NAMES = {
    "Ari": "Aries", "Tau": "Taurus", "Gem": "Gemini", "Cnc": "Cancer",
    "Leo": "Leo", "Vir": "Virgo", "Lib": "Libra", "Sco": "Scorpio",
    "Oph": "Ophiuchus", "Sgr": "Sagittarius", "Cap": "Capricorn",
    "Aqr": "Aquarius", "Psc": "Pisces",
}
DEFAULT_OUT = (Path(__file__).resolve().parents[1]
               / "zodiac_dating" / "data" / "ecliptic_constellations_iau_j2000.json")


def make_radec_fn(ts):
    """Return f(lon_deg, lat_deg) -> (ra_hours, dec_deg) in J2000."""
    R = ecliptic_J2000_frame.rotation_at(ts.tt_jd(2451545.0))   # ecliptic -> ICRS
    R = np.asarray(R)

    def f(lon_deg, lat_deg):
        lam = np.atleast_1d(np.radians(np.asarray(lon_deg, dtype=float)))
        bet = np.full_like(lam, np.radians(lat_deg))
        v = np.stack([np.cos(bet) * np.cos(lam),
                      np.cos(bet) * np.sin(lam),
                      np.sin(bet)], axis=0)
        u = R.T.dot(v)   # rotation_at maps ICRS -> ecliptic, so transpose
        ra = np.degrees(np.arctan2(u[1], u[0])) % 360.0
        dec = np.degrees(np.arcsin(np.clip(u[2], -1.0, 1.0)))
        return ra / 15.0, dec
    return f


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=str(DEFAULT_OUT))
    args = ap.parse_args()

    ts = load.timescale()
    to_radec = make_radec_fn(ts)
    constellation_at = load_constellation_map()

    def const_of(lon):
        ra, dec = to_radec(lon, 0.0)
        return np.atleast_1d(constellation_at(position_of_radec(ra, dec, epoch=ts.J2000)))

    step = 0.02
    lons = np.arange(0.0, 360.0, step)
    names = const_of(lons)

    crossings = []
    for i in range(len(lons) - 1):
        if names[i] != names[i + 1]:
            a, b = lons[i], lons[i + 1]
            na, nb = names[i], names[i + 1]
            for _ in range(60):                     # bisect to ~1e-15 deg
                m = 0.5 * (a + b)
                if const_of(m)[0] == na:
                    a = m
                else:
                    b = m
            crossings.append(((a + b) / 2.0, na, nb))

    name_at_0 = const_of(0.0)[0]
    print(f"ecliptic longitude 0 (J2000) lies in {name_at_0}")
    print(f"{len(crossings)} crossings found")
    for lon, na, nb in crossings:
        print(f"  {lon:10.4f}   {na:4s} -> {nb:4s} "
              f"({FULL_NAMES.get(na, na)} -> {FULL_NAMES.get(nb, nb)})")

    table = {
        "description": "Constellation boundaries along the ecliptic (J2000 ecliptic "
                       "longitude, latitude 0), derived from Skyfield's bundled IAU "
                       "constellation boundary map by tools/derive_boundaries.py. Each "
                       "entry marks where a sector begins: from this longitude until the "
                       "next entry the ecliptic lies in starts_sector.",
        "source": "Skyfield constellations.npz, which encodes the IAU (Delporte 1930) "
                  "boundaries as tabulated in Roman (1987), CDS VI/42, PASP 99, 695.",
        "epoch": "J2000.0 ecliptic longitude, degrees, geocentric direction",
        "crossings": [{"longitude_deg": round(l, 6),
                       "starts_sector": FULL_NAMES.get(nb, nb),
                       "ends_sector": FULL_NAMES.get(na, na),
                       "iau_abbreviation": nb}
                      for l, na, nb in crossings],
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as f:
        json.dump(table, f, indent=2)
        f.write("\n")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
