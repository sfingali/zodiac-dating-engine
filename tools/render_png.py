#!/usr/bin/env python3
"""Render a horoscope to PNG, from this engine's own positions.

There was no raster output anywhere in the project: the engine computes numbers
and viz/build_visualiser.py and viz/sky.html draw SVG in a browser.  This is the
missing piece - a plate you can put in a document, a paper or a print.

    python tools/render_png.py --year 1168 --month 4 --day 22 --hour 21 --minute 55 \
        --out /tmp/horoscope_1168.png

    python tools/render_png.py --jd 2461300.2 --kernel de440s.bsp --out today.png

    # with the long-zodiac constraints shown against what the sky actually did
    python tools/render_png.py --year 1168 --month 4 --day 22 --hour 21 --minute 55 \
        --spec examples/dendera_long_dl2_iau.json --out dl2.png

Everything drawn comes from the engine: the positions from the kernel, the sector
boundaries from the boundary set named in the footer, and the specification test
from the same rule objects the search uses.  Nothing is recomputed here, so a
plate cannot disagree with a report.  The date is printed in the calendar that
was in use at that instant - Julian before 15 October 1582, Gregorian after.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                       # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from zodiac_dating.boundaries import get_set                          # noqa: E402
from zodiac_dating.engine import _jd_to_calendar                       # noqa: E402
from zodiac_dating.ephemeris import BODIES, load_ephemeris             # noqa: E402
from zodiac_dating.specifications import Spec                          # noqa: E402

GROUND = "#0a0c10"
INK = "#e9e6e0"
DIM = "#9aa0a8"
FAINT = "#464b53"
ACCENT = "#c9a45c"
COLOUR = {
    "Sun": "#e8b44a", "Moon": "#d3dae1", "Mercury": "#9fb3c8", "Venus": "#d9a7a0",
    "Mars": "#c9705f", "Jupiter": "#c9a45c", "Saturn": "#a9a2b8",
}
GLYPH = {"Sun": "\u2609", "Moon": "\u263d", "Mercury": "\u263f", "Venus": "\u2640",
         "Mars": "\u2642", "Jupiter": "\u2643", "Saturn": "\u2644"}

R_BAND_IN, R_BAND_OUT = 0.855, 0.955
R_TICK = 1.0
R_MARK = 0.905


def polar(lon_deg: float, r: float):
    a = math.radians(180.0 + lon_deg)
    return r * math.cos(a), r * math.sin(a)


def arc(lon_a: float, lon_b: float, r: float, n: int = 240):
    lons = np.linspace(lon_a, lon_b if lon_b >= lon_a else lon_b + 360.0, n)
    return [polar(float(l), r) for l in lons]


def wrap(s: float) -> float:
    return s % 360.0


def sectors_of(bs):
    """[(name, start, end)] in order, from the boundary set's own crossings."""
    out = []
    cr = bs.crossings
    for i, (start, name) in enumerate(cr):
        end = cr[(i + 1) % len(cr)][0]
        if end <= start:
            end += 360.0
        out.append((name, float(start), float(end)))
    return out


def label_for(jd):
    y, m, d, hh, mm, ss = _jd_to_calendar(jd)
    sign = "-" if y < 0 else "+"
    cal = "Julian" if jd + 0.5 < 2299161.0 else "Gregorian"
    return (f"{sign}{abs(y):04d}-{m:02d}-{d:02d} {hh:02d}:{mm:02d}", cal,
            f"{sign}{abs(y):04d}-{m:02d}-{d:02d} {hh:02d}:{mm:02d}:{ss:02d}")


def main() -> int:
    ap = argparse.ArgumentParser(description="Render a horoscope plate to PNG.")
    ap.add_argument("--jd", type=float, help="Julian date, TT scale")
    ap.add_argument("--year", type=int)
    ap.add_argument("--month", type=int, default=1)
    ap.add_argument("--day", type=int, default=1)
    ap.add_argument("--hour", type=int, default=0)
    ap.add_argument("--minute", type=int, default=0)
    ap.add_argument("--out", required=True)
    ap.add_argument("--spec", help="specification JSON, to show its constraints against the sky")
    ap.add_argument("--table", default="iau_j2000",
                    choices=["iau_j2000", "horos_csn_j2000", "horos_cs_j2000"],
                    help="boundary set to draw (default: the real IAU boundaries)")
    ap.add_argument("--kernel", help="ephemeris .bsp (default: $ZODIAC_DATING_KERNEL or the cache)")
    ap.add_argument("--title", help="override the heading")
    ap.add_argument("--dpi", type=int, default=200)
    ap.add_argument("--width", type=float, default=13.0)
    ap.add_argument("--height", type=float, default=8.6)
    args = ap.parse_args()

    if args.kernel:
        e = load_ephemeris(args.kernel)
    else:
        e = load_ephemeris()

    eph = e
    bs = get_set(args.table)

    # ---- the instant ------------------------------------------------------
    # A calendar date goes through Skyfield's own conversion, so Delta T is
    # applied exactly as it is everywhere else in the project: ts.utc() reads a
    # civil (UT) date in the calendar then in use and returns Terrestrial Time.
    if args.jd is not None:
        jd = float(args.jd)
        t = e.ts.tt_jd(jd)
        delta_t = float(t.delta_t)
    else:
        if args.year is None:
            ap.error("give --jd or --year")
        t = e.ts.utc(args.year, args.month, args.day, args.hour, args.minute)
        jd = float(t.tt)
        delta_t = float(t.delta_t)

    coords = eph.coordinates(np.array([jd]))
    lon = {b: float(np.asarray(coords[b][0]).reshape(-1)[0]) for b in BODIES}
    lat = {b: float(np.asarray(coords[b][1]).reshape(-1)[0]) for b in BODIES}

    spec = Spec.from_dict(json.load(open(args.spec))) if args.spec else None

    # ---- figure -----------------------------------------------------------
    fig = plt.figure(figsize=(args.width, args.height), dpi=args.dpi)
    fig.patch.set_facecolor(GROUND)
    ax = fig.add_axes([0.015, 0.055, 0.70, 0.87])
    ax.set_facecolor(GROUND)
    ax.set_aspect("equal")
    ax.set_axis_off()
    ax.set_xlim(-1.30, 1.30)
    ax.set_ylim(-1.22, 1.22)
    ax.set_xticks([])
    ax.set_yticks([])

    # sector band, unequal widths, alternating weight
    secs = sectors_of(bs)
    for i, (name, a, b) in enumerate(secs):
        outer = arc(a, b, R_BAND_OUT)
        inner = arc(b, a, R_BAND_IN)
        xs = [p[0] for p in outer] + [p[0] for p in inner]
        ys = [p[1] for p in outer] + [p[1] for p in inner]
        ax.fill(xs, ys, color="#ffffff", alpha=0.030 if i % 2 == 0 else 0.055,
                linewidth=0, zorder=1)
        ax.plot(*zip(*outer), color="#6f7681", lw=0.5, zorder=2)
        ax.plot(*zip(*inner), color="#40454c", lw=0.4, zorder=2)
        for edge in (a, b):
            p0, p1 = polar(edge, R_BAND_IN), polar(edge, R_BAND_OUT)
            ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color="#4d525a", lw=0.5, zorder=2)
        # name, radially outside, and the width in degrees
        mid = (a + b) / 2.0
        px, py = polar(mid, R_BAND_OUT + 0.022)
        ang = (180.0 + mid) % 360.0
        ax.text(px, py, name, color=DIM, fontsize=7.2, ha="left", va="center",
                rotation=ang if ang < 180 else ang - 180, rotation_mode="anchor",
                zorder=6)
        px2, py2 = polar(mid, R_BAND_IN - 0.028)
        ax.text(px2, py2, f"{b - a:.1f}\u00b0", color=FAINT, fontsize=6.0,
                ha="center", va="center", zorder=5)

    # degree ring: ticks every 10, labels every 30
    for d in range(0, 360, 10):
        major = d % 30 == 0
        p0 = polar(d, R_TICK - (0.018 if major else 0.010))
        p1 = polar(d, R_TICK + 0.004)
        ax.plot([p0[0], p1[0]], [p0[1], p1[1]],
                color="#8b9199" if major else "#3d424a",
                lw=0.7 if major else 0.5, zorder=3)
        if major:
            lx, ly = polar(d, R_TICK - 0.042)
            ax.text(lx, ly, f"{d}", color="#6b7178", fontsize=5.6, ha="center",
                    va="center", zorder=4)

    # bodies: marker on the ring, name and longitude outside.  Labels are placed
    # into radial lanes by a collision test, not by a fixed rotation, or the
    # bodies that pile up in one constellation - 1168 puts three in Aries - write
    # over each other.
    order = [b for b in BODIES if b in lon]
    lanes = [R_TICK + 0.055, R_TICK + 0.105, R_TICK + 0.155, R_TICK + 0.205,
             R_TICK + 0.252]
    min_sep = 16.0
    taken: list[list[float]] = [[] for _ in lanes]
    lane_of: dict[str, int] = {}

    def fits(lane_lons, L):
        return all(abs(((L - u + 180.0) % 360.0) - 180.0) > min_sep for u in lane_lons)

    for b in sorted(order, key=lambda k: lon[k]):
        for i in range(len(lanes)):
            if fits(taken[i], lon[b]):
                taken[i].append(lon[b])
                lane_of[b] = i
                break
        else:
            taken[-1].append(lon[b])
            lane_of[b] = len(lanes) - 1

    for b in order:
        L = lon[b]
        mx, my = polar(L, R_MARK)
        ax.scatter([mx], [my], s=34, facecolor=GROUND, edgecolor=COLOUR[b],
                   linewidths=1.1, zorder=7)
        ax.text(mx, my, GLYPH[b], color=COLOUR[b], fontsize=8.5, ha="center",
                va="center", zorder=8)
        r_lab = lanes[lane_of[b]]
        lead = [polar(L, R_MARK + 0.016), polar(L, r_lab - 0.008)]
        ax.plot([lead[0][0], lead[1][0]], [lead[0][1], lead[1][1]],
                color=COLOUR[b], lw=0.45, alpha=0.42, zorder=6)
        lx, ly = polar(L, r_lab)
        ang = (180.0 + L) % 360.0
        rot = ang if ang < 180 else ang - 180
        ax.text(lx, ly, f"{b} {L:05.1f}\u00b0", color=COLOUR[b], fontsize=7.2,
                ha="left" if ang < 180 else "right", va="center",
                rotation=rot, rotation_mode="anchor", zorder=9)

    # the interior carries the instant itself, which is what the plate is of
    stamp, cal, stamp_full = label_for(jd)
    ax.text(0.0, 0.075, stamp, color=INK, fontsize=17, family="DejaVu Serif",
            ha="center", va="center", zorder=6)
    ax.text(0.0, 0.008, f"{cal} calendar \u00b7 Terrestrial Time", color=DIM,
            fontsize=7.0, ha="center", va="center", zorder=6)
    ax.text(0.0, -0.048, f"\u0394T {delta_t:.0f} s", color=FAINT, fontsize=7.0,
            family="DejaVu Sans Mono", ha="center", va="center", zorder=6)
    ax.plot(*zip(*arc(0.0, 360.0, 0.62)), color="#2b3038", lw=0.5, zorder=2)

    # ---- headings ---------------------------------------------------------
    # The instant is inside the wheel, so the corner carries only what the plate
    # is measured in: the ephemeris, the boundary table with its hash, and the JD.
    if args.title:
        fig.text(0.015, 0.968, args.title, color=INK, fontsize=16,
                 family="DejaVu Serif", va="top")
        top = 0.916
    else:
        top = 0.968
    fig.text(0.015, top,
             "geocentric apparent positions on the ecliptic of J2000 \u00b7 "
             + ("real IAU constellation boundaries" if args.table == "iau_j2000"
                else f"HOROS boundary table ({args.table})"),
             color=DIM, fontsize=7.6, va="top")
    footer = (f"JPL {Path(eph.filename).name} \u00b7 {bs.name} boundaries {bs.digest()} \u00b7 "
              f"JD(TT) {jd:.5f}")
    fig.text(0.015, 0.035, footer, color=FAINT, fontsize=6.2, family="DejaVu Sans Mono",
             va="bottom")

    # ---- readout column ---------------------------------------------------
    x0 = 0.700
    fig.text(x0, 0.968, "position", color=DIM, fontsize=7.0, family="DejaVu Sans Mono",
             va="top")
    fig.text(x0 + 0.112, 0.968, "constellation", color=DIM, fontsize=7.0,
             family="DejaVu Sans Mono", va="top")
    if spec:
        fig.text(x0 + 0.212, 0.968, "required", color=DIM, fontsize=7.0,
                 family="DejaVu Sans Mono", va="top")
        fig.text(x0 + 0.212, 0.943, "\u2713 satisfied / \u2717 violated", color=FAINT,
                 fontsize=6.2, family="DejaVu Sans Mono", va="top")
    y = 0.922
    for b in order:
        sect = bs.sector_name(lon[b])
        fig.text(x0, y, GLYPH[b], color=COLOUR[b], fontsize=9.5, va="top")
        fig.text(x0 + 0.018, y - 0.001, b, color=INK, fontsize=8.0, va="top")
        fig.text(x0, y - 0.024, f"{lon[b]:07.3f}\u00b0", color=DIM, fontsize=7.4,
                 family="DejaVu Sans Mono", va="top")
        fig.text(x0 + 0.112, y - 0.024, f"{lat[b]:+05.2f}\u00b0", color=FAINT,
                 fontsize=6.6, family="DejaVu Sans Mono", va="top")
        fig.text(x0 + 0.152, y - 0.001, sect, color=INK, fontsize=7.6, va="top")
        if spec and b in spec.rules:
            hits = [nm for (nm, a, bb) in sectors_of(bs)
                    if bool(spec.rules[b].contains(np.array([wrap((a + bb) / 2.0)]))[0])]
            ok = bool(spec.rules[b].contains(np.array([lon[b]]))[0])
            txt = ", ".join(hits) if hits else "on the boundary"
            fig.text(x0 + 0.212, y - 0.001, txt, color=INK if ok else DIM, fontsize=7.0,
                     va="top")
            fig.text(x0 + 0.212, y - 0.024, "\u2713 satisfied" if ok else "\u2717 violated",
                     color="#7fae7a" if ok else "#c9705f", fontsize=6.8,
                     family="DejaVu Sans Mono", va="top")
        y -= 0.058

    if spec:
        fig.text(x0, y - 0.012, f"specification: {Path(args.spec).name}", color=FAINT,
                 fontsize=6.4, family="DejaVu Sans Mono", va="top")

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, facecolor=GROUND)
    plt.close(fig)
    print(f"wrote {args.out}  ({Path(args.out).stat().st_size / 1024:.0f} KB)")
    print(f"  {stamp_full} ({cal} calendar, TT)")
    for b in order:
        print(f"  {b:<8} {lon[b]:7.3f}\u00b0  {lat[b]:+6.2f}\u00b0  {bs.sector_name(lon[b])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
