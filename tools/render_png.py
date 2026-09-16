#!/usr/bin/env python3
"""Render a horoscope to PNG, from this engine's own positions.

The engine computes numbers and the visualiser draws SVG in a browser; this is the
raster output, for a document, a paper or a print.

    python tools/render_png.py --year 1168 --month 4 --day 22 --hour 21 --minute 55 \
        --out /tmp/horoscope_1168.png

    python tools/render_png.py --jd 2461300.2 --kernel de440s.bsp --out today.png

    # with the long-zodiac constraints shown against what the sky actually did
    python tools/render_png.py --year 1168 --month 4 --day 22 --hour 21 --minute 55 \
        --spec examples/dendera_long_dl2_iau.json --out dl2.png

Everything drawn comes from the engine: the positions from the kernel, the sector
boundaries from the boundary set named in the footer, and the specification test
from the same rule objects the search uses.  Nothing is recomputed here, so a
plate cannot disagree with a report.  The date is printed in the calendar that was
in use at that instant - Julian before 15 October 1582, Gregorian after.

Text is never rotated.  Labels sit upright in radial lanes, assigned by a
collision test, because rotated labels on the lower half of a circle read upside
down and a chart nobody can read is not a chart.
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

from zodiac_dating.boundaries import get_set                           # noqa: E402
from zodiac_dating.engine import _jd_to_calendar                        # noqa: E402
from zodiac_dating.ephemeris import BODIES, load_ephemeris              # noqa: E402
from zodiac_dating.specifications import Spec                           # noqa: E402

GROUND = "#0a0c10"
INK = "#f1eee8"
TEXT = "#d7dce2"
DIM = "#98a0a9"
FAINT = "#5c636b"
FAINT2 = "#3b4149"
ACCENT = "#c9a45c"
GOOD = "#8dc187"
BAD = "#e0796a"
COLOUR = {
    "Sun": "#f0bd57", "Moon": "#e2e8ef", "Mercury": "#a8bdd2", "Venus": "#e2b0a8",
    "Mars": "#d97f6c", "Jupiter": "#d3ad63", "Saturn": "#b5aec4",
}
GLYPH = {"Sun": "\u2609", "Moon": "\u263d", "Mercury": "\u263f", "Venus": "\u2640",
         "Mars": "\u2642", "Jupiter": "\u2643", "Saturn": "\u2644"}

R_BAND_IN, R_BAND_OUT = 0.74, 0.86
R_MARK = 0.800
# Names and body labels are in separate radial bands with a gap wider than the
# text is tall (0.135 units ~ 24 pt at this figure size), so a name and a planet
# label can never touch however close their angles are - which is what happened
# when the two bands sat 0.07 apart.
NAME_LANES = (0.935, 1.010, 1.085)
BODY_LANES = (1.200, 1.300, 1.400, 1.500)
LIM = 1.58
NAME_MIN_SEP = 24.0
BODY_MIN_SEP = 27.0

# Standard IAU abbreviations, for the columns where the full Latin name does not
# fit: "Aries, Taurus, Gemini" runs off the page at any readable size.
ABBREV = {"Aries": "Ari", "Taurus": "Tau", "Gemini": "Gem", "Cancer": "Cnc",
          "Leo": "Leo", "Virgo": "Vir", "Libra": "Lib", "Scorpio": "Sco",
          "Ophiuchus": "Oph", "Sagittarius": "Sgr", "Capricornus": "Cap",
          "Aquarius": "Aqr", "Pisces": "Psc"}


def wrap_names(names, per_line=3):
    """Short names, at most `per_line` to a line, so the column cannot overflow."""
    short = [ABBREV.get(n, n[:3]) for n in names]
    return [" ".join(short[i:i + per_line]) for i in range(0, len(short), per_line)] or ["-"]


def polar(lon_deg: float, r: float):
    a = math.radians(180.0 + lon_deg)
    return r * math.cos(a), r * math.sin(a)


def arc(lon_a: float, lon_b: float, r: float, n: int = 260):
    lons = np.linspace(lon_a, lon_b if lon_b >= lon_a else lon_b + 360.0, n)
    return [polar(float(l), r) for l in lons]


def sep(a: float, b: float) -> float:
    return abs(((a - b + 180.0) % 360.0) - 180.0)


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


def place(items, lanes, min_sep):
    """Assign each (key, angle) to the innermost lane that is not already taken.

    Collision-tested rather than rotated or fudged: labels here are horizontal, so
    two of them crowd each other when their angles are close.
    """
    taken: list[list[float]] = [[] for _ in lanes]
    lane_of: dict[str, int] = {}
    for key, ang in sorted(items, key=lambda t: t[1]):
        for i in range(len(lanes)):
            if all(sep(ang, u) > min_sep for u in taken[i]):
                taken[i].append(ang)
                lane_of[key] = i
                break
        else:
            taken[-1].append(ang)
            lane_of[key] = len(lanes) - 1
    return lane_of


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

    eph = load_ephemeris(args.kernel) if args.kernel else load_ephemeris()
    bs = get_set(args.table)

    # ---- the instant ------------------------------------------------------
    # Calendar dates go through Skyfield's own conversion, so Delta T is applied
    # exactly as it is elsewhere in the project.
    if args.jd is not None:
        jd = float(args.jd)
        delta_t = float(eph.ts.tt_jd(jd).delta_t)
    else:
        if args.year is None:
            ap.error("give --jd or --year")
        t = eph.ts.utc(args.year, args.month, args.day, args.hour, args.minute)
        jd = float(t.tt)
        delta_t = float(t.delta_t)

    coords = eph.coordinates(np.array([jd]))
    lon = {b: float(np.asarray(coords[b][0]).reshape(-1)[0]) for b in BODIES}
    lat = {b: float(np.asarray(coords[b][1]).reshape(-1)[0]) for b in BODIES}
    spec = Spec.from_dict(json.load(open(args.spec))) if args.spec else None
    order = [b for b in BODIES if b in lon]

    # ---- figure -----------------------------------------------------------
    fig = plt.figure(figsize=(args.width, args.height), dpi=args.dpi)
    fig.patch.set_facecolor(GROUND)
    ax = fig.add_axes([0.005, 0.045, 0.700, 0.905])
    ax.set_facecolor(GROUND)
    ax.set_aspect("equal")
    ax.set_axis_off()
    ax.set_xlim(-LIM, LIM)
    ax.set_ylim(-LIM, LIM)
    ax.set_xticks([])
    ax.set_yticks([])

    secs = sectors_of(bs)

    # the band: real, unequal sectors
    for i, (name, a, b) in enumerate(secs):
        outer = arc(a, b, R_BAND_OUT)
        inner = arc(b, a, R_BAND_IN)
        ax.fill([p[0] for p in outer] + [p[0] for p in inner],
                [p[1] for p in outer] + [p[1] for p in inner],
                color="#ffffff", alpha=0.035 if i % 2 == 0 else 0.065, linewidth=0, zorder=1)
        ax.plot(*zip(*outer), color="#9aa2ab", lw=0.6, zorder=2)
        ax.plot(*zip(*inner), color="#6a717a", lw=0.5, zorder=2)
        for edge in (a, b):
            p0, p1 = polar(edge, R_BAND_IN), polar(edge, R_BAND_OUT)
            ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color="#7c848d", lw=0.55, zorder=2)

    # degree ticks: every 10, longer every 30.  No numbers - the readout has the
    # longitudes, and a ring of small numerals is the first thing to become
    # unreadable when the image is scaled down.
    for d in range(0, 360, 10):
        major = d % 30 == 0
        p0 = polar(d, 0.876)
        p1 = polar(d, 0.906 if major else 0.894)
        ax.plot([p0[0], p1[0]], [p0[1], p1[1]],
                color="#aab2bb" if major else "#5f666e",
                lw=0.8 if major else 0.5, zorder=3)

    # constellation names: upright, in lanes, with a tick to the sector
    names = [(name, (a + b) / 2.0 % 360.0) for (name, a, b) in secs]
    name_lane = place(names, NAME_LANES, NAME_MIN_SEP)
    for name, mid in names:
        r = NAME_LANES[name_lane[name]]
        p0, p1 = polar(mid, R_BAND_OUT + 0.012), polar(mid, r - 0.030)
        ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color=FAINT2, lw=0.5, zorder=3)
        px, py = polar(mid, r)
        ax.text(px, py, name, color=TEXT, fontsize=11.0, ha="center", va="center",
                zorder=6)

    # bodies: glyph on the band, upright label in a lane, leader between them
    body_lane = place([(b, lon[b]) for b in order], BODY_LANES, BODY_MIN_SEP)
    for b in order:
        L = lon[b]
        mx, my = polar(L, R_MARK)
        ax.scatter([mx], [my], s=52, facecolor=GROUND, edgecolor=COLOUR[b],
                   linewidths=1.4, zorder=7)
        ax.text(mx, my, GLYPH[b], color=COLOUR[b], fontsize=12.5, ha="center",
                va="center", zorder=8)
        r = BODY_LANES[body_lane[b]]
        lx, ly = polar(L, r - 0.028)
        ax.plot([mx, lx], [my, ly], color=COLOUR[b], lw=0.5, alpha=0.5, zorder=6)
        tx, ty = polar(L, r)
        ax.text(tx, ty, f"{b} {L:05.1f}\u00b0", color=COLOUR[b], fontsize=12.0,
                ha="center", va="center", zorder=9)

    # the instant, inside the wheel
    stamp, cal, _ = label_for(jd)
    ax.text(0.0, 0.085, stamp, color=INK, fontsize=25, family="DejaVu Serif",
            ha="center", va="center", zorder=6)
    ax.text(0.0, 0.010, f"{cal} calendar", color=DIM, fontsize=10.5, ha="center",
            va="center", zorder=6)
    ax.text(0.0, -0.052, f"Terrestrial Time \u00b7 \u0394T {delta_t:.0f} s", color="#8d949c",
            fontsize=10.5, ha="center", va="center", zorder=6)
    ax.plot(*zip(*arc(0.0, 360.0, 0.535)), color="#252a31", lw=0.6, zorder=2)

    # ---- headings ---------------------------------------------------------
    if args.title:
        fig.text(0.012, 0.972, args.title, color=INK, fontsize=15,
                 family="DejaVu Serif", va="top")
        top = 0.925
    else:
        top = 0.972
    fig.text(0.012, top,
             "geocentric apparent positions on the ecliptic of J2000 \u00b7 "
             + ("real IAU constellation boundaries" if args.table == "iau_j2000"
                else f"HOROS boundary table ({args.table})"),
             color=DIM, fontsize=9.0, va="top")
    fig.text(0.012, 0.030,
             f"JPL {Path(eph.filename).name} \u00b7 {bs.name} boundaries "
             f"{bs.digest()} \u00b7 JD(TT) {jd:.5f}",
             color="#7b838c", fontsize=8.6, family="DejaVu Sans Mono", va="bottom")

    # ---- readout column ---------------------------------------------------
    x0 = 0.700
    fig.text(x0, 0.972, "position", color=DIM, fontsize=9.0, va="top")
    fig.text(x0 + 0.098, 0.972, "constellation", color=DIM, fontsize=9.0, va="top")
    if spec:
        fig.text(x0 + 0.196, 0.972, "required", color=DIM, fontsize=9.0, va="top")
        fig.text(x0 + 0.196, 0.949, "\u2713 satisfied   \u2717 violated", color=FAINT,
                 fontsize=7.8, va="top")
    y = 0.918
    for b in order:
        sect = bs.sector_name(lon[b])
        fig.text(x0 - 0.014, y + 0.004, GLYPH[b], color=COLOUR[b], fontsize=13, va="top")
        fig.text(x0 + 0.014, y, b, color=INK, fontsize=11.5, va="top")
        fig.text(x0, y - 0.030, f"{lon[b]:07.2f}\u00b0", color=TEXT, fontsize=11.0,
                 family="DejaVu Sans Mono", va="top")
        fig.text(x0 + 0.098, y - 0.030, f"{lat[b]:+05.2f}\u00b0", color=DIM,
                 fontsize=9.0, family="DejaVu Sans Mono", va="top")
        fig.text(x0 + 0.140, y, sect, color=INK, fontsize=11.0, va="top")
        if spec and b in spec.rules:
            hits = [nm for (nm, a, bb) in secs
                    if bool(spec.rules[b].contains(np.array([((a + bb) / 2.0) % 360.0]))[0])]
            ok = bool(spec.rules[b].contains(np.array([lon[b]]))[0])
            fig.text(x0 + 0.194, y - 0.001, "\u2713" if ok else "\u2717",
                     color=GOOD if ok else BAD, fontsize=12.5, va="top")
            for i, line in enumerate(wrap_names(hits)):
                fig.text(x0 + 0.212, y + i * 0.022, line, color=TEXT if ok else DIM,
                         fontsize=9.6, va="top")
        y -= 0.080

    if spec:
        fig.text(x0, y + 0.014, f"specification: {Path(args.spec).name}", color=FAINT,
                 fontsize=8.4, family="DejaVu Sans Mono", va="top")

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, facecolor=GROUND)
    plt.close(fig)
    print(f"wrote {args.out}  ({Path(args.out).stat().st_size / 1024:.0f} KB)")
    print(f"  {stamp} ({cal} calendar, TT)")
    for b in order:
        print(f"  {b:<8} {lon[b]:7.3f}\u00b0  {lat[b]:+6.2f}\u00b0  {bs.sector_name(lon[b])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
