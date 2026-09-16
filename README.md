# zodiac-dating-engine

A sky-dating engine. You describe a planetary configuration — each of the seven
classical bodies in a constellation or a range of them, plus optionally the order
they must stand in around the ecliptic — and it returns **every window of time in
which the real sky matched that description**, computed geocentrically from JPL
ephemerides under the real IAU constellation boundaries.

It was written to date ancient sky depictions (the Egyptian zodiacs and the like)
without the two things that make dating software hard to check: hidden slack and
hidden statistics.

## What it deliberately does not do

* **No hidden tolerance.** A body is in its interval or it is not. A tolerance can
  be set in the specification, and it is then printed with the results. The
  default is zero and zero is exact.
* **No score, no p-value, no "most probable date".** It reports the windows that
  match, how wide each is, and how far each body sat from any point the source
  drawing nominated. Whether a coincidence means anything is a question for the
  reader — an engine that answers it for you has smuggled its author's
  assumptions into the answer.
* **No substitute ephemeris.** Real JPL kernels, real constellation boundaries
  derived from the IAU map, geocentric apparent positions including light-time
  and aberration.
* **No silent substitution.** The kernel, its coverage, the boundary table and
  its hash are printed with every result.

## Install

```bash
pip install -r requirements.txt          # skyfield, numpy
# or
pip install -e .                         # gives you the `zodiac-dating` command
```

The engine needs a JPL ephemeris kernel. On first use it downloads DE441 part 1
(about 1.5 GB, public domain) into `~/.cache/zodiac_dating`. To use a copy you
already have:

```bash
export ZODIAC_DATING_KERNEL=/path/to/de441_part-1.bsp
export ZODIAC_DATING_KERNEL_DIR=/path/to/directory   # searched before downloading
export ZODIAC_DATING_CACHE=/path/to/cache            # where downloads go
```

`ZODIAC_DATING_KERNEL` is honoured everywhere, including by the example
specifications, which name no kernel and therefore follow the environment. A
kernel named *inside* a specification overrides it. Nothing is fetched while a
usable kernel is already on disk: the path actually read is printed with every
run (`selftest` shows it explicitly).

DE441 part 1 covers −13200 to 1969, which is what historical dating needs;
`de440s.bsp` or `de421.bsp` serve modern dates.

## Quick start

```bash
python -m zodiac_dating boundaries                       # print every boundary table in use
python -m zodiac_dating selftest --jd 2437000.5          # seven longitudes for one instant
python -m zodiac_dating run examples/dendera_round_dr9_iau.json
python tools/smoke.py                                    # offline checks, exit 0 on pass
python tools/crosscheck_horizons.py                      # compare against JPL Horizons (network)
python tools/delta_t.py                                  # Delta T across the covered span
python tools/render_png.py --year 1168 --month 4 --day 22 --hour 21 --minute 55 \
    --out plate.png                                       # a horoscope as a PNG plate
```

`run` accepts overrides without touching the file:
`--from-jd`, `--to-jd`, `--step-hours`, `--tolerance`, `--ephemeris`,
`--json-out`.

## How a specification is written

```json
{
  "boundaries": "iau_j2000",
  "tolerance_deg": 0.0,
  "ephemeris": "de441_part-1.bsp",
  "search": {"from_jd": 990000, "to_jd": 2440400, "step_hours": 6.0},
  "planets": {
    "Sun":     {"from": "Pisces",    "to": "Pisces"},
    "Moon":    {"from": "Libra",     "to": "Libra"},
    "Mercury": {"from": "Aquarius",  "to": "Pisces"},
    "Venus":   {"from": "Pisces",    "to": "Aries"},
    "Mars":    {"from": "Capricorn", "to": "Capricorn"},
    "Jupiter": {"from": "Gemini",    "to": "Cancer"},
    "Saturn":  {"from": "Virgo",     "to": "Libra"}
  },
  "order": [["Venus"], ["Jupiter"], ["Saturn"], ["Moon"], ["Mars"], ["Mercury", "Sun"]],
  "note": "where this description comes from"
}
```

* **Intervals** run from `from` to `to` inclusive, counter-clockwise in
  increasing longitude, wrapping through 360 when they have to.
* **A sector named once** (`"from": "Pisces", "to": "Pisces"`) means that whole
  sector, because that is how the published tables phrase it. `"Pisces:1.0"` is
  its far edge; `"Taurus:0.5"` is its middle.
* **`free`** releases a body from any constraint: `{"free": true}`; a body not
  mentioned at all is also unconstrained.
* **`best`** is the point the drawing's decoder nominated. It is used only to
  report how far the computed sky sat from it. It never influences matching.
* **`order`** is the sequence of the bodies around the ecliptic by increasing
  longitude, checked cyclically. Bodies in one group (`["Mercury", "Sun"]`) may
  appear in either order, which expresses "these two are drawn interchangeably".
* **`enforce_order`** (default `false`) decides what the order constraint *does*.
  By default the engine reports, for every window it finds, whether the published
  order held — it does not throw the window away for it, and the report says
  `ok` or `VIOLATED`. Set `enforce_order` (or pass `--enforce-order`) and an
  instant only matches if the order holds too. The distinction is not academic:
  the round-zodiac DR9 window at 1836 satisfies the planet intervals but
  **violates** the published order, so it survives one reading and not the other.
* **`step_hours`** is how finely the range is sampled. Window edges are then
  refined by bisection to about a second, so the step does not limit the reported
  times — but a step longer than the shortest window can step over it entirely.
  With constellation-scale constraints, half a day is safe for the outer bodies;
  the Moon moves 13.2° per day, so keep the step well under a day when the Moon
  constrains the answer.

## Boundary sets

The table in force is named and hashed in the output.

| Name | Sectors | What it is |
|---|---|---|
| `iau_j2000` | 13 | The real constellation boundaries on the ecliptic, including **Ophiuchus**, derived by `tools/derive_boundaries.py` from the IAU map (Delporte 1930, as tabulated by Roman 1987, CDS VI/42). |
| `horos_cs_j2000` | 12 | The boundary table distributed with Fomenko & Nosovsky's program HOROS, old **CS** table, used 2002–17.11.2007. |
| `horos_csn_j2000` | 12 | The same program's current **CSN** table, in use from 17.11.2007. |
| custom | any | Any set of crossings you supply in the specification, with its provenance, e.g. to test a published table against the real sky. |

The IAU set as derived (J2000 ecliptic longitude, degrees, each value being where
a sector begins):

```
Aries 28.6853   Taurus 53.4157   Gemini 90.1378   Cancer 117.9850
Leo 138.0353    Virgo 173.8483   Libra 217.8086   Scorpio 241.0300
Ophiuchus 247.6361   Sagittarius 266.2355   Capricorn 299.6531
Aquarius 327.4845    Pisces 351.6493
```

Ophiuchus holds 18.6° of the ecliptic. The two HOROS tables do not contain it:
they assign those degrees to Scorpio and Sagittarius. This engine keeps it,
because it is there in the sky — and both tables are available so that a result
published under HOROS can be compared on HOROS's own terms instead of argued
with. `examples/dendera_round_dr9_iau.json` and
`examples/dendera_round_dr9_horos_csn.json` are the same drawing and the same
published decipherment run under each, which isolates what the boundary table
alone does to the answer.

## What an answer means

```bash
python -m zodiac_dating run examples/dendera_round_dr9_iau.json
```

prints, for each window: the Julian-calendar start and end (TT), the duration in
hours, the seven geocentric apparent J2000 ecliptic longitudes at the midpoint
with the sector each falls in, the deviation from any nominated `best` point, and
whether the order constraint held. A window is a real interval: if the
configuration stood for three days, three days is the answer, and a program that
reports one date and time has thrown that width away.

## Accuracy, stated plainly

* **Ephemeris and frame handling: better than 0.01°.** `tools/crosscheck_horizons.py`
  asks JPL Horizons itself for geocentric vectors and converts them to J2000
  ecliptic longitude independently of Skyfield's frame code; the residual against
  this engine's apparent longitudes runs to about 0.008° at worst, which is the
  known difference between a geometric and an apparent position (light-time and
  aberration), not an error in the ephemeris. `tools/smoke.py` locks that in as a
  regression test with the Horizons values hard-coded.
* **ΔT is the real limit, and only for old dates.** Dates are converted to
  Terrestrial Time using Skyfield's ΔT model, because the ephemerides are
  functions of TT while documents are in civil time. The Moon moves 0.549° per
  hour, so a disagreement of 1,000 s between ΔT models — the order of the spread
  across published models in the first millennium BC — moves the Moon by about
  0.15°, while the planets move by amounts too small to matter. `tools/delta_t.py`
  prints the model's values across the covered span (25,310 s at −1000, 10,430 s
  at +1, 1,650 s at 1000, 109 s at 1600) so that the assumption is visible rather
  than buried. A lunar longitude in antiquity should be read as good to a few
  tenths of a degree, not to 0.01°.
* **What the engine cannot tell you.** Whether a configuration is remarkable. It
  reports windows; the number of days per century on which such a configuration
  occurs is a property of the sky, and anyone claiming a dating is interesting
  owes you that number.

## Worked examples

Three specifications are bundled, all transcribed from published decipherment
tables with their variant codes (see `REFERENCES.md`):

| File | Drawing and variant | Boundary set |
|---|---|---|
| `examples/dendera_round_dr9_iau.json` | Round Dendera, DR9 | real IAU |
| `examples/dendera_round_dr9_horos_csn.json` | Round Dendera, DR9 | HOROS current (CSN) |
| `examples/dendera_long_dl2_iau.json` | Long Dendera, DL2 | real IAU |

Every one of them was run over the whole span of DE441 part 1 — JD 990000 to
2440400, about −13000 to 1969 — sampled every 6 hours, 5,801,601 instants. The
full reports are in `results/`, and `results/run_examples.sh` regenerates them.
No tolerance is applied unless a run says so.

| Specification | Boundary set | Tolerance | Windows found |
|---|---|---|---|
| Long Dendera DL2 | real IAU | 0° | **1**: 22–24 Apr 1168, 45.5 h, order satisfied |
| Round Dendera DR9 | real IAU | 0° | **0** |
| Round Dendera DR9 | HOROS CSN | 0° | **0** |
| Round Dendera DR9 | real IAU | 5° | **2**: 19–21 Feb −169, 6–9 Mar 1836 — both violate the order |
| Round Dendera DR9 | HOROS CSN | 5° | **2**: 16–18 Feb 271 (order satisfied), 6–7 Mar 1836 (violated) |

Dates are printed in the calendar that was in use: Julian before 15 October 1582,
Gregorian from that day on (see *Time conventions*).

Three things are worth reading out of that table.

**The long zodiac reproduces, and it is not a near miss.** Under real boundaries,
with nothing widened, the DL2 decipherment admits exactly one window in thirteen
thousand years: 22–24 April 1168 (Julian). Fomenko and Nosovsky publish the long
zodiac's exhaustive solution as 22–26 April 1168. The window this engine finds sits
*inside* theirs, from a different implementation, a different ephemeris and
independently derived boundary tables.

**The round zodiac's DR9 variant produces nothing at all on an exact reading**,
under the real boundaries or under HOROS's own current table — so the variant as
published is not merely rare, it never happens. It only yields dates when HOROS's
documented preliminary tolerance of ±5° is applied, at which point two appear, and
then the published *order* of the bodies eliminates one of them under their own
table and both of them under the real boundaries.

**The boundary table moves ancient dates by centuries.** With ±5° allowed, the
same drawing and the same published variant date to −169 or to 271 depending
solely on whether the real IAU boundaries or the HOROS table is used; the
nineteenth-century window is unaffected. That is the size of the effect the choice
of boundary table has on a dating, and it is why this engine prints which one it
used, with its hash.


## Visualiser

`viz/build_visualiser.py` writes a self-contained page — no build step, no CDN,
no server-side compute at request time — showing the wheels, the constraint
bands, the real (uneven) constellation sectors including Ophiuchus, a scrubber
over each window, and a timeline of every window found in the 13,000-year span.
It reads the reports in `results/` and recomputes the match state with the
engine's own rule objects, so the page cannot disagree with the reports.

```bash
python viz/build_visualiser.py --out index.html      # needs the kernel, as usual
```

## Sky view, and the positions API

`viz/sky.html` is a depiction of the sky rather than of any decipherment: the
seven bodies at their true geocentric positions, a date and time you can step by
an hour or a century, trails showing where each body has been, and the real
constellations. It runs from a small FastAPI service in `serve/app.py`:

```
GET /api/health
GET /api/range
GET /api/sky?year=&month=&day=&hour=&minute=      # or ?jd=<TT Julian date>
GET /api/track?jd=&span_days=&n=                  # for the trails
```

The API exists because Skyfield is Python and the kernel is 1.5 GB — a browser
cannot compute this, and precomputing a grid would fix the step and the range.
The service loads the ephemeris once at startup and answers in milliseconds
after that. It is supervised by s6 (`/opt/data/services/zodiac/run`, watchdog in
`/opt/data/scripts/`, boot install in `/opt/data/s6-install-zodiac.sh`) and
reached through Caddy at `/zodiac/api/*`.

A copy runs at <https://stephenfingleton.com/zodiac/>.

## Time conventions

The engine is strict about calendar and time scale, and this is the part most
easily got wrong:

* A date you type is read in the calendar that was in use: **Julian before
  15 October 1582, Gregorian from that day on** (Skyfield's switch is set
  explicitly to 2299161.0, because Skyfield's default is proleptic Gregorian,
  which would put a twelfth-century date seven days out).
* Positions are functions of **Terrestrial Time**, so a civil date is converted
  with Skyfield's ΔT model, and ΔT is printed with the result — 25,310 s in
  1000 BC, 1,002 s in 1168, about 69 s today.
* The same switch applies when parsing and when printing, so a date the engine
  prints can be typed back in and land on the same instant. `tools/delta_t.py`
  shows the model across the covered span.

## Plates

`tools/render_png.py` draws a horoscope as a PNG, for a document or a print:
real sector widths, degree ring, the seven bodies with their ecliptic latitudes,
and the instant in the centre with the calendar and Delta T. It reads nothing but
the engine - positions from the kernel, sectors from the boundary set, constraint
tests from the same rule objects the search uses - so a plate cannot disagree with
a report. `--spec` adds a *required* column and marks each body satisfied or
violated, which turns the plate into a picture of a dating.

```sh
python tools/render_png.py --year 1168 --month 4 --day 22 --hour 21 --minute 55 --out plate.png
python tools/render_png.py --year 1168 --month 4 --day 22 --hour 21 --minute 55 \
    --spec examples/dendera_long_dl2_iau.json --out dl2.png
python tools/render_png.py --year 2026 --month 9 --day 16 --kernel de440s.bsp --out today.png
```

Three are in `results/plates/`: the long zodiac as its own dating requires it
(1168), the same constraints at the instant 1,151 years earlier (17 CE, where
three of the seven bodies violate them), and the plain sky.

## Corrections made while building this

Kept here because both changed published numbers, and a reader comparing an
earlier version deserves to know why they moved.

1. **Dates were being printed in the wrong calendar.** The conversion carried the
   Gregorian correction term for every date while labelling the result Julian, and
   Skyfield's own switch was being set to the opposite of what was intended
   (its default is proleptic Gregorian, not Julian-before-1582). Twelfth-century
   windows therefore appeared seven days late: the long zodiac's unambiguous
   window moved from "29 April – 1 May 1168" to **22–24 April 1168**, which is
   inside the 22–26 April that Fomenko and Nosovsky publish. What had looked like
   a seven-day disagreement with the published result was this bug.
2. **The order constraint was reported but never matched on.** It is now explicit
   (`enforce_order`, `--enforce-order`, off by default). Turning it on removes one
   of the two round-zodiac candidates under HOROS's own table and both of them
   under the real boundaries.

## Files

```
zodiac_dating/
  ephemeris.py       kernels, geocentric apparent longitudes, kernel discovery
  boundaries.py      the boundary tables, and the provenance of each
  specifications.py  the input format, and its semantics
  engine.py          the scan, window edges, and the report
  cli.py             command line
  data/              the derived IAU boundary table
tools/
  derive_boundaries.py    regenerate the IAU table from the IAU map
  smoke.py                offline self-test, exit 0 on pass
  crosscheck_horizons.py  compare against JPL Horizons
  delta_t.py              print the Delta T model in use
viz/
  build_visualiser.py     generate the self-contained visualiser page
  sky.html                the sky view (served at /zodiac/ alongside the API)
serve/
  app.py                  positions API behind the sky view
results/plates/           rendered horoscope plates (PNG)
```

## Licence and attribution

MIT — see `LICENCE`. Built on Skyfield (MIT), JPL DE440/DE441 (NASA, public
domain), and the IAU constellation boundaries (Delporte 1930; Roman 1987, CDS
VI/42). The HOROS boundary tables are reproduced as data with their provenance
and switchover date; no code from HOROS, or from any other dating program, is
used here. Full citations, including the published Dendera tables the examples
come from, are in `REFERENCES.md`.
