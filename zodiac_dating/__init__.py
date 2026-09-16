"""Zodiac dating engine: geocentric horoscope dating with real boundaries."""

from .boundaries import BoundarySet, MODERN_IAU, get_set, iau_set
from .ephemeris import BODIES, Ephemeris, load_ephemeris
from .specifications import PlanetRule, Spec, SpecError
from .engine import Window, scan, format_report

__version__ = "1.0.0"

__all__ = [
    "BoundarySet", "MODERN_IAU", "get_set", "iau_set", "BODIES", "Ephemeris",
    "load_ephemeris", "PlanetRule", "Spec", "SpecError", "Window", "scan",
    "format_report", "__version__",
]
