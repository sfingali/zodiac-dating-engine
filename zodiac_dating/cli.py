#!/usr/bin/env python3
"""Command-line interface.

    python -m zodiac_dating run examples/dendera_round_dr9.json
    python -m zodiac_dating boundaries            # print the boundary tables
    python -m zodiac_dating selftest              # ephemeris sanity check
"""
import argparse
import json
import sys

from .boundaries import get_set, iau_set, MODERN_IAU
from .engine import format_report, scan
from .specifications import Spec


def cmd_run(args):
    spec = Spec.from_file(args.spec)
    if args.step_hours:
        spec.step_hours = args.step_hours
    if args.from_jd is not None:
        spec.from_jd = args.from_jd
    if args.to_jd is not None:
        spec.to_jd = args.to_jd
    if args.tolerance is not None:
        spec.tolerance_deg = args.tolerance
    if args.ephemeris:
        spec.ephemeris = args.ephemeris

    def progress(n, k):
        print(f"  sampled {n} instants, {k} matched", file=sys.stderr)

    windows, stats = scan(spec, progress=progress)
    report = format_report(spec, windows, stats)
    print(report)
    if args.json_out:
        with open(args.json_out, "w") as f:
            json.dump({"spec": spec.note,
                       "specification": spec.describe().splitlines(),
                       "stats": stats,
                       "windows": [{"start_jd_tt": w.start_jd_tt,
                                    "end_jd_tt": w.end_jd_tt,
                                    "midpoint_jd_tt": w.midpoint_jd_tt,
                                    "width_days": w.width_days,
                                    "longitudes": w.longitudes,
                                    "deviations": w.deviations,
                                    "mean_abs_deviation": w.mean_abs_deviation,
                                    "order_ok": w.order_ok,
                                    "order_found": w.order_found}
                                   for w in windows]}, f, indent=2)
        print(f"wrote {args.json_out}", file=sys.stderr)
    return 0


def cmd_boundaries(args):
    for name in (MODERN_IAU, "horos_cs_j2000", "horos_csn_j2000"):
        b = get_set(name)
        print(f"== {b.name} ==   sha256 {b.digest()}")
        if b.provenance:
            print(f"   provenance: {b.provenance}")
        for lon, sector in b.crossings:
            print(f"   {lon:10.4f}  {sector}")
        spans = []
        for i, (lon, sector) in enumerate(b.crossings):
            nxt = b.crossings[(i + 1) % len(b.crossings)][0]
            spans.append(f"{sector} {(nxt - lon) % 360:6.2f}")
        print("   widths: " + ", ".join(spans))
        print()
    return 0


def cmd_selftest(args):
    """Compute the seven longitudes for a known instant and echo them."""
    from .ephemeris import load_ephemeris, BODIES
    import numpy as np
    eph = load_ephemeris(args.ephemeris or None)
    jd = args.jd
    lons = eph.longitudes(np.array([jd]))
    print(f"ephemeris: {eph.filename}")
    print(f"  read from: {eph.path}")
    print(f"  JD(TT) {jd}")
    for b in BODIES:
        print(f"  {b:<8s} {float(lons[b][0]):9.5f}")
    return 0


def main(argv=None):
    p = argparse.ArgumentParser(prog="zodiac_dating", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("run", help="scan a specification")
    r.add_argument("spec")
    r.add_argument("--step-hours", type=float)
    r.add_argument("--from-jd", type=float)
    r.add_argument("--to-jd", type=float)
    r.add_argument("--tolerance", type=float)
    r.add_argument("--ephemeris")
    r.add_argument("--json-out")
    r.set_defaults(func=cmd_run)

    b = sub.add_parser("boundaries", help="print the boundary tables in use")
    b.set_defaults(func=cmd_boundaries)

    s = sub.add_parser("selftest", help="print the seven longitudes for one instant")
    s.add_argument("--jd", type=float, required=True)
    s.add_argument("--ephemeris")
    s.set_defaults(func=cmd_selftest)

    args = p.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
