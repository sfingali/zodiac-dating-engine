# References and attributions

Everything this engine depends on, with its licence and where it came from.

## Ephemeris

**Park, R. S., Folkner, W. M., Williams, J. G., Boggs, D. H. (2021).** *The JPL
Planetary and Lunar Ephemerides DE440 and DE441.* The Astronomical Journal
**161**(3), 105. doi:10.3847/1538-3881/abd414

The underlying ephemeris files (`de441_part-1.bsp`, `de441_part-2.bsp`,
`de440s.bsp`, `de421.bsp`) are produced by NASA/Jet Propulsion Laboratory,
California Institute of Technology, and are in the public domain. They are
served from <https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/>
and <https://ssd.jpl.nasa.gov/ftp/eph/planets/bsp/>. No ephemeris file is
redistributed here; the engine downloads one into its cache directory on first
use, or you point it at a copy you already have.

## Astronomy library

**Skyfield** — Brandon Rhodes. *Skyfield: Elegant Astronomy for Python.*
MIT licence. <https://github.com/skyfielders/python-skyfield>

Skyfield provides: the SPK kernel reader, the geocentric observation chain
(light-time and aberration corrections), the J2000 ecliptic frame, the
constellation boundary map, and the Delta T model. The engine is thin on top of
it by design — reimplementing any of that would only add a place for a mistake
to hide.

## Constellation boundaries

**Delporte, E. (1930).** *Délimitation définitive des constellations.*
International Astronomical Union, Cambridge.

**Roman, N. G. (1987).** *Identification of a Constellation from a Position.*
Publications of the Astronomical Society of the Pacific **99**, 695.
CDS catalogue VI/42.

These are the official constellation boundaries. The boundary map bundled with
Skyfield encodes them, and `tools/derive_boundaries.py` walks the ecliptic
through that map to find where the boundaries actually fall on the zodiac —
that derivation, not a hand-copied list, is what `iau_j2000` uses.

## Calendar and time

**Meeus, J. (1998).** *Astronomical Algorithms*, 2nd ed., Willmann-Bell.
Chapter 7 (Julian Day) is the source of the calendar conversion used when a
window edge is printed as a date.

**Morrison, L. V. & Stephenson, F. R. (2004).** *Historical values of the Earth's
clock error ΔT and the calculation of eclipses.* Journal for the History of
Astronomy **35**, 327.

**Stephenson, F. R., Morrison, L. V. & Hohenkerk, C. Y. (2016).** *Measurement of
the Earth's rotation: 720 BC to AD 2015.* Proceedings of the Royal Society A
**472**, 20160404.

Delta T (TT − UT1) converts civil dates into the Terrestrial Time arguments the
ephemerides require. The model applied is Skyfield's; `tools/delta_t.py` prints
it across the covered span so it is never an unstated assumption. This is the
dominant uncertainty in the engine for any date before 1600, and it is an
uncertainty about the Earth's rotation, not about the planets.

## Published zodiac decipherment tables

**Fomenko, A. T. & Nosovsky, G. V.** *New Chronology of Egypt*, Appendix 2 —
the published decipherment variants for the Egyptian zodiacs (the round and long
Dendera zodiacs, the Esna zodiacs, the Athribis zodiacs), giving each planet's
admissible constellations and, for some variants, a nominated point.
<https://chronologia.org/nx_egypt/pril2.html>

The example specifications in `examples/` are transcriptions of two of those
published variants — round Dendera **DR9** and long Dendera **DL2** — quoted with
their variant codes so any reader can check them against the source. They are
included to give the engine something real to run on, and to let a HOROS result
be compared on HOROS's own terms.

**HOROS** (Fomenko & Nosovsky's dating program, distributed at chronologia.org)
ships two boundary tables: the old **CS** table used from 2002 to 17 November
2007, and the current **CSN** table used from that date. Both are reproduced in
`zodiac_dating/boundaries.py` *as data*, with their provenance and switchover
date stated, so that this engine can reproduce their behaviour and measure it.
No code from HOROS is used or copied; it is a Windows binary built on 1990s
Fortran (PLANETAP, ELP2000-85) and carries no licence granting reuse.

## The paper this engine was built to check

**Baiget Orts, C. (2025).** *Astronomical Refutation of the New Chronology by
Fomenko and Nosovsky: The 1151-Year Planetary Cycle...* arXiv:2504.12962.

The engine exists because of the claim in that paper and the discussion around
it. It is independent of it: it implements neither the paper's method nor its
metric, it computes real geocentric positions from JPL ephemerides under real
constellation boundaries, and it reports windows rather than a single date or a
score. His code, which is MIT-licensed
(<https://github.com/carbaior/1151cycle>, <https://github.com/carbaior/sescc>),
was read and run separately during the review that prompted this; nothing from it
is included here.

## This code

MIT. See `LICENCE`. Written for Stephen Fingleton, 2026.
