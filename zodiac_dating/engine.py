"""
The scanning engine.

Given a :class:`Spec`, walk the requested instant range at the requested step,
evaluate every planet's constraint at every instant, and report the *windows of
time* during which the whole configuration satisfies the specification.  Window
edges are then refined by bisection down to about a second, so the reported
times are limited by the specification and the ephemeris, not by the step.

Three deliberate design choices, all of them the opposite of what the existing
programs do:

* Windows, not instants.  A horoscope constraint is an inequality, so the
  answer is an interval.  Reporting a single date and time hides how wide that
  interval is.
* No hidden slack.  Nothing expands an interval except an explicit
  ``tolerance_deg``, which is printed with the results.
* No statistical claim.  The output states what was searched, what matched, and
  how far each planet sat from its nominated point.  It does not compute a
  p-value, a likelihood or a posterior; see README for why that is left out.
"""
from dataclasses import dataclass
import datetime as _dt

import numpy as np

from .ephemeris import BODIES, load_ephemeris
from .specifications import Spec

# Instants evaluated per block inside scan(); see the comment there.
CHUNK = 200_000


@dataclass
class Window:
    start_jd_tt: float
    end_jd_tt: float
    midpoint_jd_tt: float
    longitudes: dict
    deviations: dict
    mean_abs_deviation: float
    order_ok: bool
    order_found: list

    @property
    def width_days(self):
        return self.end_jd_tt - self.start_jd_tt


def _jd_to_calendar(jd):
    """Calendar date-time for a JD on the TT scale, in the calendar that was in
    use: Julian before 15 October 1582, Gregorian from that day onward.

    Meeus, Astronomical Algorithms, ch. 7, with the switch at JD 2299161. That is
    the same rule Skyfield applies when it parses a calendar date, so a date this
    function prints can be read back in and lands on the same instant - an earlier
    version applied the Gregorian correction to every date while labelling the
    result Julian, which shifted twelfth-century dates by seven days and made
    printing and parsing disagree for anything after 1582.
    """
    jd = float(jd) + 0.5
    z = int(jd)
    f = jd - z
    if z >= 2299161:                        # Gregorian
        alpha = int((z - 1867216.25) / 36524.25)
        a = z + 1 + alpha - alpha // 4
    else:                                   # Julian
        a = z
    b = a + 1524
    c = int((b - 122.1) / 365.25)
    d = int(365.25 * c)
    e = int((b - d) / 30.6001)
    day = b - d - int(30.6001 * e) + f
    month = e - 1 if e < 14 else e - 13
    year = c - 4716 if month > 2 else c - 4715
    di = int(day)
    frac = day - di
    hours = frac * 24.0
    hh = int(hours)
    mm = int((hours - hh) * 60)
    ss = int((((hours - hh) * 60) - mm) * 60)
    return year, month, di, hh, mm, ss


def _calendar_name(jd):
    return "Gregorian calendar" if float(jd) + 0.5 >= 2299161 else "Julian calendar"


def _fmt_jd(jd):
    y, m, d, hh, mm, ss = _jd_to_calendar(jd)
    return f"{y:+05d}-{m:02d}-{d:02d} {hh:02d}:{mm:02d}:{ss:02d} ({_calendar_name(jd)}, TT)"


def _cyclic_order_match(found_order, required):
    """True if ``found_order`` (list of body names, in longitude order) contains
    ``required`` as a cyclic subsequence, where bodies in the same group may
    appear in any order."""
    n = len(found_order)
    if not required:
        return True
    first_group = set(required[0])
    starts = [i for i, body in enumerate(found_order) if body in first_group]
    for s in starts:
        seq = found_order[s:] + found_order[:s]
        pos = 0
        ok = True
        for group in required:
            g = set(group)
            end = pos
            while end < len(seq) and seq[end] in g:
                end += 1
            if end == pos or not any(b in g for b in seq[pos:end + 1]):
                ok = False
                break
            # every member of the group must be inside the block
            block = seq[pos:end]
            if not g.issubset(set(block)):
                ok = False
                break
            pos = end
        if ok:
            return True
    return False


def _evaluate(spec, jds, eph):
    """Return (matches, longitudes) for an array of JD(TT)."""
    lons = eph.longitudes(jds)
    mask = np.ones(len(jds), dtype=bool)
    for body, rule in spec.rules.items():
        mask &= rule.contains(lons[body], spec.tolerance_deg)
    return mask, lons


def scan(spec: Spec, progress=None):
    """Return (windows, stats dict)."""
    eph = load_ephemeris(spec.ephemeris)
    jd0, jd1 = spec.resolve_range()
    _check_range(eph, jd0, jd1, spec.ephemeris)

    step = spec.step_hours / 24.0
    n = int(np.ceil((jd1 - jd0) / step)) + 1
    jds = jd0 + step * np.arange(n)
    jds = jds[jds <= jd1 + 1e-9]

    # Evaluate in blocks.  A one-hour sweep of the whole DE441 span is tens of
    # millions of instants; the longitudes for all of them at once would not fit
    # in memory, and they are not needed - only the per-instant match matters.
    mask = np.zeros(len(jds), dtype=bool)
    done = 0
    for start in range(0, len(jds), CHUNK):
        block = jds[start:start + CHUNK]
        m, _ = _evaluate(spec, block, eph)
        mask[start:start + len(block)] = m
        done += len(block)
        if progress:
            progress(done, int(mask[:done].sum()))

    # The order constraint is reported for every window; it is only matched on
    # when the specification asks for it, because the published tables state an
    # order that not every candidate window satisfies (see README).
    if spec.enforce_order and spec.order:
        for i in np.nonzero(mask)[0]:
            lons_i = eph.longitudes(np.array([jds[i]]))
            present = sorted((b for b in BODIES if b in spec.rules),
                             key=lambda b: float(lons_i[b][0]))
            if not _cyclic_order_match(present, spec.order):
                mask[i] = False

    windows = []
    i = 0
    while i < len(mask):
        if not mask[i]:
            i += 1
            continue
        j = i
        while j + 1 < len(mask) and mask[j + 1]:
            j += 1
        start = _refine_edge(spec, eph, jds[i], -1) if i > 0 else jds[i]
        end = _refine_edge(spec, eph, jds[j], +1) if j < len(mask) - 1 else jds[j]
        windows.append(_make_window(spec, eph, start, end))
        i = j + 1

    stats = {"sampled_instants": int(len(jds)), "matching_instants": int(mask.sum()),
             "windows": len(windows), "step_hours": spec.step_hours}
    return windows, stats


def _check_range(eph, jd0, jd1, filename):
    from .ephemeris import EphemerisRangeError
    try:
        seg_start, seg_end = eph.range_years()
    except Exception:
        return
    if jd1 < seg_start or jd0 > seg_end:
        raise EphemerisRangeError(
            f"requested JD range {jd0:.1f}..{jd1:.1f} lies outside {filename} "
            f"({seg_start:.1f}..{seg_end:.1f})")


def _matches_at(spec, eph, jd):
    mask, _ = _evaluate(spec, np.array([jd]), eph)
    return bool(mask[0])


def _refine_edge(spec, eph, jd_inside, direction, tol=1e-5):
    """Push a window edge outwards until the configuration stops matching."""
    step = 1.0 / 24.0
    lo = jd_inside
    hi = jd_inside + direction * step
    for _ in range(60):
        if _matches_at(spec, eph, hi):
            lo = hi
            hi = lo + direction * step
            step *= 2.0
            if step > 400:
                break
        else:
            break
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if _matches_at(spec, eph, mid):
            lo = mid
        else:
            hi = mid
        if abs(hi - lo) < tol:
            break
    return lo


def _make_window(spec, eph, start, end):
    mid = 0.5 * (start + end)
    lons = {b: float(eph.longitudes(np.array([mid]))[b][0]) for b in BODIES}
    devs, absdevs = {}, []
    for body, rule in spec.rules.items():
        if rule.best_deg is None or rule.free:
            continue
        d = (lons[body] - rule.best_deg + 180.0) % 360.0 - 180.0
        devs[body] = d
        absdevs.append(abs(d))
    order_found = sorted((b for b in BODIES if b in spec.rules),
                         key=lambda b: lons[b])
    return Window(start_jd_tt=start, end_jd_tt=end, midpoint_jd_tt=mid,
                  longitudes=lons, deviations=devs,
                  mean_abs_deviation=float(np.mean(absdevs)) if absdevs else float("nan"),
                  order_ok=_cyclic_order_match(order_found, spec.order),
                  order_found=order_found)


def format_report(spec, windows, stats):
    out = []
    out.append("=" * 78)
    out.append("ZODIAC DATING ENGINE - results")
    out.append("=" * 78)
    out.append(spec.describe())
    out.append("-" * 78)
    out.append(f"instants sampled  : {stats['sampled_instants']}")
    out.append(f"instants matching : {stats['matching_instants']}")
    out.append(f"windows found     : {stats['windows']}")
    out.append("-" * 78)
    if not windows:
        out.append("No instant in the searched range satisfies the specification.")
        return "\n".join(out)

    out.append(f"{'#':>3} {'start (JD TT)':>15} {'end (JD TT)':>15} {'days':>9} "
               f"{'mean|dev|':>9}  order")
    for k, w in enumerate(windows, 1):
        out.append(f"{k:>3} {w.start_jd_tt:>15.5f} {w.end_jd_tt:>15.5f} "
                   f"{w.width_days:>9.4f} {w.mean_abs_deviation:>9.3f}  "
                   f"{'ok' if w.order_ok else 'VIOLATED'}")
    out.append("-" * 78)
    for k, w in enumerate(windows, 1):
        out.append(f"[{k}] {_fmt_jd(w.start_jd_tt)}")
        out.append(f"    to {_fmt_jd(w.end_jd_tt)}")
        out.append(f"    duration {w.width_days * 24.0:.3f} hours "
                   f"({w.width_days:.4f} days)")
        out.append("    geocentric apparent J2000 ecliptic longitudes at midpoint:")
        for body in BODIES:
            if body in w.longitudes:
                dev = w.deviations.get(body)
                sector = spec.boundaries.sector_name(w.longitudes[body])
                line = f"      {body:<8s} {w.longitudes[body]:8.3f} deg  ({sector})"
                if dev is not None:
                    line += f"   deviation from nominated point {dev:+7.3f} deg"
                out.append(line)
        out.append(f"    cyclic order found: {' < '.join(w.order_found)}")
        out.append(f"    order requirement : "
                   f"{'satisfied' if w.order_ok else 'NOT satisfied'}")
    return "\n".join(out)
