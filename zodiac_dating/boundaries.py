"""
Constellation boundary sets used to turn a longitude into a sector name.

Three kinds of set are supported, and the one in force is always named in the
output:

* ``iau_j2000``      the real constellation boundaries on the ecliptic,
                     derived by ``tools/derive_boundaries.py`` from Skyfield's
                     bundled IAU boundary map (Delporte 1930 boundaries as
                     tabulated by Roman 1987, CDS VI/42).  Thirteen sectors:
                     Ophiuchus really is on the ecliptic and is not merged away.
* user sets          any explicit list of crossings, e.g. from a published
                     source, carried in the input file with its provenance.
* ``horos_cs_j2000`` / ``horos_csn_j2000``
                     the two boundary tables distributed with Fomenko and
                     Nosovsky's HOROS (old table used 2002-2007, new table from
                     17.11.2007), so that dates can be compared with that
                     program on its own terms.  Both omit Ophiuchus; that is
                     recorded in the set, not hidden.

A boundary set is a list of (longitude, sector) pairs, ordered by increasing
longitude, where each longitude marks the *start* of that sector.  Sector
membership is therefore exact: a longitude belongs to the sector whose start
boundary it has most recently passed.  Nothing is expanded or nudged.
"""
from dataclasses import dataclass, field
import json
import hashlib
from importlib import resources

MODERN_IAU = "iau_j2000"

HOROS_CS = [
    (26.0, "Aries"), (51.0, "Taurus"), (89.0, "Gemini"), (118.0, "Cancer"),
    (143.0, "Leo"), (174.0, "Virgo"), (215.0, "Libra"), (236.0, "Scorpio"),
    (266.0, "Sagittarius"), (301.0, "Capricorn"), (329.0, "Aquarius"),
    (346.0, "Pisces"),
]
HOROS_CSN = [
    (31.0, "Aries"), (56.0, "Taurus"), (92.0, "Gemini"), (118.0, "Cancer"),
    (137.0, "Leo"), (172.0, "Virgo"), (215.0, "Libra"), (239.0, "Scorpio"),
    (266.0, "Sagittarius"), (296.0, "Capricorn"), (326.0, "Aquarius"),
    (349.0, "Pisces"),
]

HOROS_PROVENANCE = (
    "Boundary tables distributed with HOROS (Fomenko & Nosovsky). The old "
    "table (CS) was used from 2002 to 17.11.2007; the new table (CSN) from "
    "17.11.2007. Both are given in J2000 ecliptic longitude in the program's "
    "own documentation files (readme and the CSN scale listing shipped in "
    "zodiak02.zip). Both omit Ophiuchus, which the IAU boundaries place on the "
    "ecliptic between 247.64 and 266.24 degrees."
)


@dataclass
class BoundarySet:
    name: str
    crossings: list                    # [(longitude_deg, sector_name), ...] increasing
    provenance: str = ""
    notes: str = ""

    def __post_init__(self):
        self.crossings = sorted(((float(l) % 360.0, s) for l, s in self.crossings),
                                key=lambda x: x[0])
        if len(self.crossings) < 2:
            raise ValueError("a boundary set needs at least two crossings")

    # -- membership -------------------------------------------------------
    def sector_index(self, lon_deg):
        import numpy as np
        starts = np.array([c[0] for c in self.crossings])
        idx = np.searchsorted(starts, np.asarray(lon_deg) % 360.0, side="right") - 1
        idx = np.where(np.asarray(lon_deg) % 360.0 < starts[0], len(starts) - 1, idx)
        return idx

    def sector_name(self, lon_deg):
        idx = self.sector_index(lon_deg)
        names = [c[1] for c in self.crossings]
        if getattr(idx, "shape", ()) == ():
            return names[int(idx)]
        import numpy as np
        return np.array(names, dtype=object)[idx]

    def start_of(self, sector_name) -> float:
        for lon, name in self.crossings:
            if name == sector_name:
                return lon
        raise KeyError(f"sector {sector_name!r} is not in boundary set {self.name!r}")

    def end_of(self, sector_name) -> float:
        """Longitude where this sector ends, i.e. the next sector's start."""
        for i, (lon, name) in enumerate(self.crossings):
            if name == sector_name:
                return self.crossings[(i + 1) % len(self.crossings)][0]
        raise KeyError(f"sector {sector_name!r} is not in boundary set {self.name!r}")

    def degrees(self, point) -> float:
        """Resolve 'Sector' or 'Sector:0.5' (fraction across the sector) to degrees."""
        if isinstance(point, (int, float)):
            return float(point) % 360.0
        if ":" in point:
            sector, frac = point.split(":", 1)
            frac = float(frac)
        else:
            sector, frac = point, 0.0
        a = self.start_of(sector)
        b = self.end_of(sector)
        width = (b - a) % 360.0
        if not 0.0 <= frac <= 1.0:
            raise ValueError(f"fraction in {point!r} must lie between 0 and 1")
        return (a + frac * width) % 360.0

    def digest(self) -> str:
        payload = json.dumps([[round(l, 9), s] for l, s in self.crossings],
                             sort_keys=True).encode()
        return hashlib.sha256(payload).hexdigest()[:16]

    def describe(self) -> str:
        return ", ".join(f"{l:.4f} {s}" for l, s in self.crossings)


def iau_set() -> BoundarySet:
    """The real boundaries on the ecliptic, derived from the IAU boundary data."""
    path = resources.files("zodiac_dating.data").joinpath(
        "ecliptic_constellations_iau_j2000.json")
    raw = json.loads(path.read_text())
    crossings = [(c["longitude_deg"], c["starts_sector"]) for c in raw["crossings"]]
    return BoundarySet(
        name=MODERN_IAU,
        crossings=crossings,
        provenance=raw["source"],
        notes="Thirteen sectors including Ophiuchus (real constellation on the ecliptic).",
    )


def horos_set(which: str) -> BoundarySet:
    table = {"horos_cs_j2000": HOROS_CS, "horos_csn_j2000": HOROS_CSN}[which]
    return BoundarySet(name=which, crossings=table, provenance=HOROS_PROVENANCE,
                       notes="Twelve sectors; Ophiuchus omitted as in the HOROS tables.")


def get_set(name: str) -> BoundarySet:
    if name == MODERN_IAU:
        return iau_set()
    if name in ("horos_cs_j2000", "horos_csn_j2000"):
        return horos_set(name)
    raise KeyError(f"unknown boundary set {name!r}; use '{MODERN_IAU}', "
                   "'horos_cs_j2000', 'horos_csn_j2000', or supply one in the input file")
