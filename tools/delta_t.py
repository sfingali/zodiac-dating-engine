#!/usr/bin/env python3
"""
Show the Delta T (TT - UT1) this engine uses, and what it costs in accuracy.

Every date you feed the engine is converted to a Terrestrial Time argument using
Skyfield's Delta T model, because the JPL ephemerides are functions of TT while
the drawings, chronicles and tables we date are in civil (UT1-ish) time.  That
conversion is the dominant uncertainty for anything before 1600, and it is not
something the ephemeris can fix: the planets are known to a milliarcsecond, the
Earth's rotation rate a thousand years ago is not.

    python tools/delta_t.py                 # sample the whole covered span
    python tools/delta_t.py --from-years -3000 --to-years 2000 --step 250

For scale: the Moon moves about 0.549 degrees per hour of time (13.18 deg/day),
Mercury in the mean about 0.0008 deg/h, Saturn about 0.00005 deg/h.  So a
disagreement between Delta T models of 1000 s - which is the order of the spread
between published models in the first millennium BC - moves the Moon by about
0.15 deg and the other bodies by nothing that matters.  For comparison the
ephemeris and frame handling in this engine agree with JPL Horizons to better
than 0.01 deg (tools/crosscheck_horizons.py).
"""
import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from zodiac_dating.ephemeris import load_ephemeris, resolve_kernel, DEFAULT_EPHEMERIS

MOON_DEG_PER_HOUR = 13.176358 / 24.0     # mean motion, degrees per hour


def jd_from_year(y):
    """Julian date at the start of the given Julian-calendar year, for sampling."""
    return 1721423.5 + (y - 1) * 365.25


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--from-years", type=int, default=-13000)
    ap.add_argument("--to-years", type=int, default=2000)
    ap.add_argument("--step", type=int, default=500)
    ap.add_argument("--kernel", default=DEFAULT_EPHEMERIS)
    args = ap.parse_args(argv)

    eph = load_ephemeris(resolve_kernel(args.kernel))
    lo, hi = eph.range_years()
    print(f"kernel: {eph.filename}\n  covers JD(TT) {lo:.1f} .. {hi:.1f}\n")
    print(f"{'year':>7}{'JD(TT)':>14}{'Delta T (s)':>14}{'= hours':>10}"
          f"{'Moon moves':>14}")
    for y in range(args.from_years, args.to_years + 1, args.step):
        jd = jd_from_year(y)
        if not (lo <= jd <= hi):
            continue
        dt = float(eph.ts.tt_jd(jd).delta_t)
        print(f"{y:>7}{jd:>14.1f}{dt:>14.1f}{dt/3600.0:>10.2f}"
              f"{dt/3600.0*MOON_DEG_PER_HOUR:>13.2f} deg")
    print("\nThis is the model value Skyfield applies; it is printed here so that it "
          "is never an unstated assumption.  It is used as a time argument, not as a "
          "geometric correction, and it is the same for every body.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
