"""
Input specifications.

A specification says, in plain JSON:

    {
      "boundaries": "iau_j2000",          # or "horos_csn_j2000", "horos_cs_j2000",
                                          # or {"name": ..., "provenance": ...,
                                          #     "crossings": [{"longitude_deg":..,"sector":..}, ...]}
      "tolerance_deg": 0.0,               # 0.0 means zero. Nothing is added silently.
      "ephemeris": "de441_part-1.bsp",
      "search": {"from_jd": 990000, "to_jd": 2500000, "step_hours": 1.0},
      "planets": {
        "Sun":     {"from": "Psc",        "to": "Ari",      "best": "Ari:0.5"},
        "Moon":    {"from": "Lib",        "to": "Lib",      "best": "Lib:0.5"},
        "Mercury": {"from": "Aqr",        "to": "Psc",      "best": "Psc:0.5"},
        "Venus":   {"from": "Psc",        "to": "Ari",      "best": "Ari:0.25"},
        "Mars":    {"from": "Cap",        "to": "Cap",      "best": "Cap:0.5"},
        "Jupiter": {"from": "Gem",        "to": "Cnc",      "best": "Gem:0.5"},
        "Saturn":  {"from": "Vir",        "to": "Lib",      "best": "Lib:0.5"}
      },
      "order": [["Venus"], ["Jupiter"], ["Saturn"], ["Moon"], ["Mars"],
                ["Mercury"], ["Sun"]],
      "note": "Dendra round zodiac, variant DR9, as described in Fomenko &
               Nosovsky, New Chronology of Egypt, Appendix 2."
    }

Semantics, deliberately strict:

* An interval runs from ``from`` to ``to`` inclusive, counter-clockwise in
  increasing longitude, wrapping through 360 if necessary.
* A sector named once, in both ``from`` and ``to``, means that whole sector -
  which is how the published tables phrase it ("the Sun is in Pisces"), and not
  a zero-width point.  Write ``"Pisces:1.0"`` if you ever mean the far edge.
* ``"Sector:f"`` addresses the point a fraction ``f`` of the way across the
  sector, so ``"Taurus:0.5"`` is its middle.
* A planet may also be declared free: ``{"free": true}``.
* ``best`` is used *only* to report how far the computed position sits from the
  point the decoder nominated on the drawing.  It never influences matching.
  If it is omitted, no deviation is reported for that planet.
* ``tolerance_deg`` widens every interval symmetrically and is quoted in the
  output.  At 0.0 a planet one thousandth of a degree outside its interval does
  not match.
"""
from dataclasses import dataclass, field
import json

from .boundaries import BoundarySet, get_set
from .ephemeris import BODIES


class SpecError(ValueError):
    pass


@dataclass
class PlanetRule:
    body: str
    free: bool = False
    start_deg: float = 0.0
    end_deg: float = 0.0
    best_deg: float = None

    def contains(self, lon, tolerance=0.0):
        import numpy as np
        if self.free:
            return np.ones_like(np.asarray(lon, dtype=float), dtype=bool)
        a = (self.start_deg - tolerance) % 360.0
        b = (self.end_deg + tolerance) % 360.0
        x = np.asarray(lon, dtype=float) % 360.0
        width = (b - a) % 360.0
        if width == 0.0 and tolerance == 0.0:
            # degenerate interval: only an exact hit counts, and floating point
            # makes an exact hit meaningless, so it can never match.
            return np.zeros_like(x, dtype=bool)
        if (self.end_deg - self.start_deg) % 360.0 + 2 * tolerance >= 360.0:
            return np.ones_like(x, dtype=bool)
        return ((x - a) % 360.0) <= width


@dataclass
class Spec:
    boundaries: BoundarySet
    rules: dict                     # body -> PlanetRule
    tolerance_deg: float = 0.0
    ephemeris: str = "de441_part-1.bsp"
    from_jd: float = None
    to_jd: float = None
    step_hours: float = 1.0
    order: list = field(default_factory=list)
    note: str = ""

    # ------------------------------------------------------------------
    @classmethod
    def from_dict(cls, data) -> "Spec":
        b = data.get("boundaries", "iau_j2000")
        if isinstance(b, str):
            bset = get_set(b)
        else:
            bset = BoundarySet(name=b.get("name", "custom"),
                               crossings=[(c["longitude_deg"], c["sector"])
                                          for c in b["crossings"]],
                               provenance=b.get("provenance", ""))

        rules = {}
        for body, rule in data["planets"].items():
            if body not in BODIES:
                raise SpecError(f"unknown body {body!r}; expected one of {BODIES}")
            if rule.get("free"):
                rules[body] = PlanetRule(body=body, free=True)
                continue
            start = bset.degrees(rule["from"])
            end = bset.degrees(rule["to"])
            if rule["from"] == rule["to"] and ":" not in str(rule["from"]):
                # A bare sector named once means the whole sector, not a point.
                start = bset.start_of(rule["from"])
                end = bset.end_of(rule["from"])
            best = bset.degrees(rule["best"]) if rule.get("best") else None
            rules[body] = PlanetRule(body=body, start_deg=start, end_deg=end, best_deg=best)

        search = data.get("search", {})
        return cls(boundaries=bset,
                   rules=rules,
                   tolerance_deg=float(data.get("tolerance_deg", 0.0)),
                   ephemeris=data.get("ephemeris", "de441_part-1.bsp"),
                   from_jd=search.get("from_jd"),
                   to_jd=search.get("to_jd"),
                   step_hours=float(search.get("step_hours", 1.0)),
                   order=data.get("order", []),
                   note=data.get("note", ""))

    @classmethod
    def from_file(cls, path) -> "Spec":
        with open(path) as f:
            return cls.from_dict(json.load(f))

    # ------------------------------------------------------------------
    def resolve_range(self):
        if self.from_jd is None or self.to_jd is None:
            raise SpecError("search.from_jd and search.to_jd are required")
        if self.to_jd <= self.from_jd:
            raise SpecError("search.to_jd must be greater than search.from_jd")
        return float(self.from_jd), float(self.to_jd)

    def describe(self):
        lines = [f"boundary set      : {self.boundaries.name} "
                 f"(sha256 {self.boundaries.digest()})",
                 f"boundary values   : {self.boundaries.describe()}",
                 f"ephemeris         : {self.ephemeris}",
                 f"tolerance         : {self.tolerance_deg:g} degrees "
                 f"({'exact, nothing added' if self.tolerance_deg == 0 else 'user-supplied'})",
                 f"search            : JD {self.from_jd:.1f} to {self.to_jd:.1f}, "
                 f"step {self.step_hours:g} h"]
        for body in BODIES:
            rule = self.rules.get(body)
            if rule is None:
                lines.append(f"  {body:<8s} not constrained")
            elif rule.free:
                lines.append(f"  {body:<8s} free (any longitude)")
            else:
                best = "" if rule.best_deg is None else f", best point {rule.best_deg:.2f}"
                lines.append(f"  {body:<8s} {rule.start_deg:.2f} to {rule.end_deg:.2f} "
                             f"(width {(rule.end_deg - rule.start_deg) % 360:.2f} deg){best}")
        if self.order:
            lines.append("order constraint  : " +
                         " < ".join("=".join(g) for g in self.order))
        if self.note:
            lines.append(f"note              : {self.note}")
        return "\n".join(lines)
