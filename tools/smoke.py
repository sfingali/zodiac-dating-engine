#!/usr/bin/env python3
"""Offline self-test: ephemeris, boundary sets, and one anchored comparison.

Run it after installation, with no network access required:

    python tools/smoke.py                 # uses $ZODIAC_DATING_KERNEL or the default
    python tools/smoke.py --kernel /path/to/de441_part-1.bsp

What is checked:

1. the kernel loads and reports its own JD(TT) coverage;
2. the seven longitudes come out finite for an ancient instant and a medieval
   one inside that coverage;
3. the boundary sets load, and the IAU set really has thirteen sectors with
   Ophiuchus present;
4. for one instant, the longitudes sit within 0.01 deg of values taken
   independently from JPL Horizons.

On point 4: the Horizons numbers below were obtained from its geocentric ICRF
vectors, rotated to the J2000 ecliptic with the fixed mean obliquity - that is a
*geometric* position.  This engine reports *apparent* position (light-time and
aberration included), so the two differ by up to ~0.01 deg, mostly the annual
aberration of the observed body.  The tolerance encodes that known difference;
it is not a fudge factor, and it is not used anywhere in the dating code.

Exit status is 0 only if every check passes.
"""
import argparse
import os
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from zodiac_dating.boundaries import iau_set, get_set, MODERN_IAU
from zodiac_dating.ephemeris import BODIES, DEFAULT_EPHEMERIS, load_ephemeris, resolve_kernel

# JD(TT) 2437000.5 (~year 1961), inside DE441 part 1.
ANCHOR_JD = 2437000.5
ANCHOR_HORIZONS = {          # geometric, from JPL Horizons vectors (see docstring)
    "Sun": 346.9961, "Moon": 93.7209, "Mercury": 354.3042, "Venus": 319.3889,
    "Mars": 310.3616, "Jupiter": 271.2749, "Saturn": 286.9023,
}
ANCHOR_TOLERANCE_DEG = 0.01

ANCIENT_JD = 1500000.5       # ~year -900, near the start of DE441 part 1

failures = []


def check(label, ok, detail=""):
    print(f"  {'PASS' if ok else 'FAIL'}  {label}{('  ' + detail) if detail else ''}")
    if not ok:
        failures.append(label)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--kernel", default=os.environ.get("ZODIAC_DATING_KERNEL",
                                                       DEFAULT_EPHEMERIS))
    args = ap.parse_args(argv)

    path = resolve_kernel(args.kernel)
    eph = load_ephemeris(path)
    lo, hi = eph.range_years()
    print(f"kernel: {args.kernel}\n  read from: {eph.path}\n  covers JD(TT) {lo:.1f} .. {hi:.1f}")

    print("\n1. coverage")
    check("kernel covers the anchor instant", lo <= ANCHOR_JD <= hi)
    check("kernel covers the ancient instant", lo <= ANCIENT_JD <= hi)

    print("\n2. longitudes are finite")
    for jd, tag in ((ANCIENT_JD, "ancient"), (ANCHOR_JD, "anchor")):
        lons = eph.longitudes(np.array([jd]))
        bad = [b for b in BODIES if not np.isfinite(lons[b][0])]
        check(f"all seven bodies computed at the {tag} instant (JD {jd})", not bad,
              f"missing: {bad}" if bad else "")

    print("\n3. boundary sets")
    bs = iau_set()
    check("IAU set has thirteen sectors", len(bs.crossings) == 13,
          f"got {len(bs.crossings)}")
    names = [n for _, n in bs.crossings]
    check("Ophiuchus is present in the IAU set", "Ophiuchus" in names)
    for other in ("horos_cs_j2000", "horos_csn_j2000"):
        b = get_set(other)
        check(f"{other} loads (twelve sectors, Ophiuchus absent by construction)",
              len(b.crossings) == 12)

    print("\n4. anchored against JPL Horizons (geometric vs apparent)")
    lons = eph.longitudes(np.array([ANCHOR_JD]))
    for body in BODIES:
        got = float(lons[body][0])
        want = ANCHOR_HORIZONS[body]
        diff = abs((got - want + 180.0) % 360.0 - 180.0)
        check(f"{body:<8s} within {ANCHOR_TOLERANCE_DEG} deg",
              diff <= ANCHOR_TOLERANCE_DEG, f"engine {got:9.4f}  Horizons {want:9.4f}  "
                                            f"diff {diff:.5f}")

    print()
    if failures:
        print(f"{len(failures)} check(s) FAILED: " + ", ".join(failures))
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
