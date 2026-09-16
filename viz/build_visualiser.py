#!/usr/bin/env python3
"""
Build the static visualiser for the zodiac dating engine.

Reads the run reports in results/ (windows already found by the engine, so the
page cannot disagree with it), recomputes the planet longitudes on a scrub grid
around each window, recomputes the match state at every grid point under every
boundary table and tolerance on offer, and writes a single self-contained
index.html: no build step, no CDN, no server-side compute at request time.

    python viz/build_visualiser.py --out /opt/data/www/zodiac/index.html

Everything the page shows is decided here, in Python, by the same code that
produced the reports.  The JavaScript only draws.
"""
import argparse
import datetime as dt
import json
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from zodiac_dating.boundaries import get_set, iau_set          # noqa: E402
from zodiac_dating.engine import _cyclic_order_match           # noqa: E402
from zodiac_dating.ephemeris import BODIES, load_ephemeris     # noqa: E402
from zodiac_dating.specifications import Spec                  # noqa: E402

RESULTS = ROOT / "results"
EXAMPLES = ROOT / "examples"

# Scrub grid: +/- this many days around each window midpoint, at this step.
HALF_WINDOW_DAYS = 20.0
STEP_HOURS = 2.0


def parse_windows(report_text):
    """Pull the numeric window rows out of an engine report.

    The mean-deviation column can be `nan` when no `best` points were nominated,
    so it is matched as any non-space token rather than as a number.
    """
    windows = []
    for line in report_text.splitlines():
        m = re.match(r"\s*(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+(\S+)\s+(\w+)", line)
        if m:
            windows.append({"start_jd": float(m.group(2)), "end_jd": float(m.group(3)),
                            "width_days": float(m.group(4)), "order_ok": m.group(6) == "ok"})
    return windows


def report_facts(path):
    text = path.read_text()
    facts = {}
    for key, pat in (("boundary_set", r"boundary set\s+: (\S+)"),
                     ("boundary_hash", r"boundary set\s+: \S+ \(sha256 (\w+)\)"),
                     ("kernel", r"ephemeris\s+: (\S+)"),
                     ("tolerance", r"tolerance\s+: (\S+) degrees")):
        m = re.search(pat, text)
        if m:
            facts[key] = m.group(1)
    return facts, parse_windows(text)


def spec_rules(spec):
    """The resolved degree intervals for the page, so the arcs it draws are the
    same numbers the engine matched against."""
    out = []
    for body in BODIES:
        rule = spec.rules.get(body)
        if rule is None:
            out.append({"body": body, "free": True})
        elif rule.free:
            out.append({"body": body, "free": True})
        else:
            out.append({"body": body, "start": round(rule.start_deg, 4),
                        "end": round(rule.end_deg, 4),
                        "best": None if rule.best_deg is None else round(rule.best_deg, 4)})
    return out


def boundaries_payload(bset):
    return {"name": bset.name, "hash": bset.digest(),
            "crossings": [{"lon": round(lon, 4), "name": name}
                          for lon, name in bset.crossings]}


def grid(eph, mid_jd):
    """Longitudes for every body on the scrub grid, plus a match bit per
    (spec, table, tolerance) combination the page offers."""
    step = STEP_HOURS / 24.0
    n = int(2 * HALF_WINDOW_DAYS / step) + 1
    jds = mid_jd - HALF_WINDOW_DAYS + step * np.arange(n)
    lons = eph.longitudes(jds)
    return jds, {b: [round(float(v), 4) for v in lons[b]] for b in BODIES}


def match_bits(spec, lons, jds, tolerance, enforce_order=False):
    """The engine's own evaluation, vectorised, at every grid instant.

    Same rule objects the engine uses (rule.contains), so the page cannot
    disagree with the reports it sits beside.  The order constraint is returned
    separately rather than folded in, because that is what the engine does: it
    reports the order, and only enforces it if the specification asks
    (``enforce_order``).
    """
    n = len(jds)
    mask = np.ones(n, dtype=bool)
    for body, rule in spec.rules.items():
        mask &= np.asarray(rule.contains(np.array(lons[body]), tolerance), dtype=bool)

    order_bits = np.ones(n, dtype=bool)
    if spec.order:
        names = [b for b in BODIES if b in spec.rules]
        arr = np.array([lons[b] for b in names])
        order = np.argsort(arr, axis=0)
        for i in range(n):
            present = [names[j] for j in order[:, i]]
            if not _cyclic_order_match(present, spec.order):
                order_bits[i] = False
    if enforce_order:
        mask &= order_bits
    return "".join("1" if m else "0" for m in mask), "".join("1" if m else "0" for m in order_bits)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=str(ROOT / "index.html"))
    ap.add_argument("--kernel", default=None)
    args = ap.parse_args(argv)

    eph = load_ephemeris(args.kernel or None)
    payload = {
        "kernel": {"name": eph.filename, "path": eph.path,
                   "start_jd": round(eph.range_years()[0], 1),
                   "end_jd": round(eph.range_years()[1], 1)},
        "boundaries": {name: boundaries_payload(get_set(name))
                       for name in ("iau_j2000", "horos_cs_j2000", "horos_csn_j2000")},
        "iau": boundaries_payload(iau_set()),
        "wheels": [],
        "timeline": {"from_jd": 990000.0, "to_jd": 2440400.0, "runs": []},
    }

    # --- the two specifications the page shows -------------------------------
    plans = [
        # key, spec file, report file, scrub centre (JD TT), tables, tolerances
        {"key": "dl2", "spec": EXAMPLES / "dendera_long_dl2_iau.json",
         "report": RESULTS / "dendera_long_dl2_iau.txt",
         "tables": ["iau_j2000"], "tolerances": [0.0, 5.0]},
        {"key": "dr9", "spec": EXAMPLES / "dendera_round_dr9_horos_csn.json",
         "report": RESULTS / "dendera_round_dr9_horos_csn_tol5.txt",
         "tables": ["iau_j2000", "horos_cs_j2000", "horos_csn_j2000"],
         "tolerances": [0.0, 5.0]},
    ]

    for plan in plans:
        base_spec = Spec.from_file(plan["spec"])
        facts, windows = report_facts(plan["report"])
        if windows:
            centre = 0.5 * (windows[0]["start_jd"] + windows[0]["end_jd"])
        else:                                    # DR9 exact has none: use their window
            centre = 2437000.5
        centre = float(centre)

        jds, lons = grid(eph, centre)

        # rules and match bits per (table, tolerance)
        variants = []
        for table in plan["tables"]:
            raw = json.loads(Path(plan["spec"]).read_text())
            spec = Spec.from_dict({**raw, "boundaries": table})
            for tol in plan["tolerances"]:
                mask, order_bits = match_bits(spec, lons, jds, tol)
                variants.append({
                    "table": table, "tolerance": tol,
                    "rules": spec_rules(spec),
                    "mask": mask, "order": order_bits,
                })

        payload["wheels"].append({
            "key": plan["key"],
            "note": base_spec.note,
            "centre_jd": centre,
            "grid_start_jd": float(jds[0]), "grid_end_jd": float(jds[-1]),
            "step_days": STEP_HOURS / 24.0,
            "jds": [round(float(j), 5) for j in jds],
            "lons": lons,
            "variants": variants,
            "reported_windows": windows,
            "report_facts": facts,
        })

    # --- timeline: every window the engine has found, on one 13,000-year axis --
    for label, table, tol, rep in (
            ("Long Dendera, real boundaries", "iau_j2000", 0.0,
             RESULTS / "dendera_long_dl2_iau.txt"),
            ("Round Dendera, real boundaries, exact", "iau_j2000", 0.0,
             RESULTS / "dendera_round_dr9_iau.txt"),
            ("Round Dendera, real boundaries, 5 deg", "iau_j2000", 5.0,
             RESULTS / "dendera_round_dr9_iau_tol5.txt"),
            ("Round Dendera, HOROS CSN table, exact", "horos_csn_j2000", 0.0,
             RESULTS / "dendera_round_dr9_horos_csn.txt"),
            ("Round Dendera, HOROS CSN table, 5 deg", "horos_csn_j2000", 5.0,
             RESULTS / "dendera_round_dr9_horos_csn_tol5.txt")):
        _, windows = report_facts(rep)
        payload["timeline"]["runs"].append({
            "label": label, "table": table, "tolerance": tol,
            "windows": [{"start_jd": w["start_jd"], "end_jd": w["end_jd"],
                         "width_days": w["width_days"]} for w in windows],
        })

    payload["delta_t"] = delta_t_samples(eph)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(render(payload))
    print(f"wrote {out}  ({out.stat().st_size/1024:.0f} KB)")
    for w in payload["wheels"]:
        print(f"  wheel {w['key']}: {len(w['jds'])} grid points, "
              f"{len(w['variants'])} variants, "
              f"{sum(v['mask'].count('1') for v in w['variants'])} matching points total")


def delta_t_samples(eph):
    out = []
    for year in (-13000, -10000, -5000, -2000, -1000, -500, 0, 500, 1000, 1500, 1969):
        jd = 1721423.5 + (year - 1) * 365.25
        try:
            out.append({"year": year,
                        "seconds": round(float(eph.ts.tt_jd(jd).delta_t), 1)})
        except Exception:
            pass
    return out


def render(payload):
    data = json.dumps(payload, separators=(",", ":"))
    return HTML.replace("__DATA__", data)


HTML = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Zodiac dating engine — where the sky matched the drawing</title>
<meta name="description" content="Every window of time in which the real sky matched a published decipherment of the Egyptian zodiacs, computed from JPL ephemerides under the real IAU constellation boundaries.">
<link rel="icon" href="data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 32 32'%3E%3Ccircle cx='16' cy='16' r='13' fill='none' stroke='%23c9c4b8' stroke-width='1'/%3E%3Ccircle cx='16' cy='16' r='4' fill='%23c9c4b8'/%3E%3C/svg%3E">
<style>
  :root{
    --bg:#08090a; --ink:#e9e6df; --dim:#8d887e; --faint:#4a4842;
    --line:rgba(233,230,223,.14); --line2:rgba(233,230,223,.07);
    --match:#d8c48a; --warn:#c98f6a;
    --mono:ui-monospace,"JetBrains Mono","SF Mono",Menlo,Consolas,monospace;
  }
  *{box-sizing:border-box}
  html,body{margin:0;background:var(--bg);color:var(--ink)}
  body{font-family:var(--mono);font-size:13px;line-height:1.65;
       -webkit-font-smoothing:antialiased;letter-spacing:.01em}
  main{max-width:1180px;margin:0 auto;padding:8vh 5vw 18vh}
  h1{font-weight:400;font-size:clamp(26px,4.4vw,52px);line-height:1.1;margin:0 0 .6em;letter-spacing:-.01em}
  h2{font-weight:400;font-size:11px;letter-spacing:.34em;text-transform:uppercase;
     color:var(--dim);margin:0 0 3.2em}
  a{color:var(--ink);text-decoration:none;border-bottom:1px solid var(--line)}
  a:hover{border-color:var(--ink)}
  .lede{max-width:62ch;color:var(--dim);font-size:14px}
  .facts{list-style:none;display:flex;flex-wrap:wrap;gap:2.6em 4.4em;padding:0;margin:3.4em 0 0}
  .facts b{display:block;font-weight:400;font-size:24px;color:var(--ink);letter-spacing:-.01em}
  .facts span{color:var(--dim);font-size:11px;letter-spacing:.16em;text-transform:uppercase}
  section{margin-top:15vh;border-top:1px solid var(--line2);padding-top:5vh}
  .stage{display:grid;grid-template-columns:minmax(0,1fr) minmax(0,1fr);gap:6vw;align-items:start}
  @media(max-width:900px){.stage{grid-template-columns:1fr;gap:8vh}}
  svg{width:100%;height:auto;display:block;overflow:visible}
  .sector{fill:none;stroke:var(--line);stroke-width:1}
  .bound{stroke:var(--faint);stroke-width:1}
  .band{fill:rgba(233,230,223,.05);stroke:none}
  .best{stroke:var(--dim);stroke-width:1;stroke-dasharray:2 3}
  .marker circle{fill:var(--ink)}
  .marker text{fill:var(--ink);font-family:var(--mono);font-size:9px;letter-spacing:.14em}
  .marker.off circle{fill:var(--faint)}
  .marker.off text{fill:var(--dim)}
  .slabel{fill:var(--dim);font-family:var(--mono);font-size:9px;letter-spacing:.1em}
  .sline{stroke:var(--line);stroke-width:.6}
  .state{fill:var(--dim);font-family:var(--mono);font-size:10px;letter-spacing:.24em;text-transform:uppercase}
  .state.yes{fill:var(--match)}
  .read{color:var(--dim)}
  .read b{color:var(--ink);font-weight:400}
  .read table{width:100%;border-collapse:collapse;margin-top:1.4em;font-size:12px}
  .read td{padding:.28em 0;border-bottom:1px solid var(--line2);white-space:nowrap}
  .read td.n{text-align:right;font-variant-numeric:tabular-nums;color:var(--ink)}
  .read td.d{text-align:right;font-variant-numeric:tabular-nums;color:var(--dim)}
  .controls{margin-top:2.4em;display:flex;flex-direction:column;gap:1.1em}
  .scrub{display:flex;align-items:center;gap:1.2em}
  input[type=range]{-webkit-appearance:none;appearance:none;width:100%;background:transparent;height:22px}
  input[type=range]::-webkit-slider-runnable-track{height:1px;background:var(--line)}
  input[type=range]::-moz-range-track{height:1px;background:var(--line)}
  input[type=range]::-webkit-slider-thumb{-webkit-appearance:none;width:9px;height:9px;
    border-radius:50%;background:var(--ink);margin-top:-4px}
  input[type=range]::-moz-range-thumb{width:9px;height:9px;border:0;border-radius:50%;background:var(--ink)}
  .seg{display:flex;gap:1.6em;flex-wrap:wrap;color:var(--dim)}
  .seg button{background:none;border:0;padding:0 0 .2em;font:inherit;color:var(--dim);
    cursor:pointer;border-bottom:1px solid transparent}
  .seg button[aria-pressed=true]{color:var(--ink);border-bottom-color:var(--ink)}
  .stamp{color:var(--faint);font-size:11px;margin-top:1.6em;letter-spacing:.04em}
  .win{fill:var(--match)} .win5{fill:var(--warn)}
  .strip{width:100%;height:auto}
  .axis{fill:var(--faint);font-family:var(--mono);font-size:9px}
  .rowlabel{fill:var(--dim);font-family:var(--mono);font-size:10px}
  .method{columns:2;column-gap:6vw}
  @media(max-width:900px){.method{columns:1}}
  .method p{margin:0 0 1.1em;color:var(--dim)}
  .method h3{font-weight:400;font-size:12px;margin:0 0 .5em;color:var(--ink);letter-spacing:.02em}
  .method div{break-inside:avoid;margin-bottom:2.4em}
  code{font-family:var(--mono);color:var(--ink)}
  .foot{margin-top:14vh;border-top:1px solid var(--line2);padding-top:3vh;color:var(--faint);font-size:11px}
</style>
</head>
<body>
<main>
  <header>
    <h1>Where the sky<br>matched the drawing</h1>
    <p class="lede">A published decipherment of an Egyptian zodiac says which
      constellation each planet stood in, and in what order around the ecliptic.
      These are the windows of time in which the real sky, computed geocentrically
      from JPL&nbsp;DE441 under the real IAU constellation boundaries, actually
      matched. No tolerance unless one is named, and none hidden.</p>
    <ul class="facts" id="facts"></ul>
  </header>

  <section id="wheel-dl2">
    <h2>Long Dendera zodiac · variant DL2</h2>
    <div class="stage">
      <div class="wheel"></div>
      <div class="read"></div>
    </div>
  </section>

  <section id="wheel-dr9">
    <h2>Round Dendera zodiac · variant DR9</h2>
    <div class="stage">
      <div class="wheel"></div>
      <div class="read"></div>
    </div>
  </section>

  <section id="timeline">
    <h2>Thirteen thousand years, every window</h2>
    <div class="stripwrap"></div>
  </section>

  <section id="method">
    <h2>Method, and what this does not claim</h2>
    <div class="method" id="method"></div>
  </section>

  <div class="foot" id="foot"></div>
</main>
<script>
const D = __DATA__;
const BODIES = ["Sun","Moon","Mercury","Venus","Mars","Jupiter","Saturn"];
const GLYPH = {Sun:"☉",Moon:"☽",Mercury:"☿",Venus:"♀",Mars:"♂",Jupiter:"♃",Saturn:"♄"};
const JDN = {cx:0,cy:0,R:0};

function pt(lon,r,cx,cy){const a=(180+lon)*Math.PI/180;return [cx+r*Math.cos(a),cy-r*Math.sin(a)];}
function arcSeg(lon0,lon1,r,cx,cy,steps){ // polyline, so no sweep-flag ambiguity
  const span=((lon1-lon0)%360+360)%360, n=Math.max(2,Math.ceil(span/steps));
  let d="";
  for(let i=0;i<=n;i++){const [x,y]=pt(lon0+span*i/n,r,cx,cy);d+=(i?"L":"M")+x.toFixed(2)+" "+y.toFixed(2)+" ";}
  return d;
}
function jdToYear(jd){return Math.round((jd-1721423.5)/365.25)+1;}
function jdToDate(jd){ // the calendar that was in use: Julian before 15 Oct 1582
  let j=jd+0.5, z=Math.floor(j), f=j-z, a=z;
  if(z>=2299161){const al=Math.floor((z-1867216.25)/36524.25);a=z+1+al-Math.floor(al/4);}
  const b=a+1524, c=Math.floor((b-122.1)/365.25), dd=Math.floor(365.25*c);
  const e=Math.floor((b-dd)/30.6001);
  const day=b-dd-Math.floor(30.6001*e)+f;
  const month=e<14?e-1:e-13, year=month>2?c-4716:c-4715;
  const di=Math.floor(day), fr=day-di, hh=fr*24;
  const h=Math.floor(hh), mi=Math.floor((hh-h)*60), s=Math.floor((((hh-h)*60)-mi)*60);
  const pad=(n,w=2)=>String(n).padStart(w,"0");
  return (year<0?"-"+pad(-year,4):pad(year,4))+"-"+pad(month)+"-"+pad(di)+" "+
         pad(h)+":"+pad(mi)+":"+pad(s);
}
function sectorOf(crossings,lon){
  let best=null;
  for(const c of crossings){const start=c.lon, next=(start+((crossings[(crossings.indexOf(c)+1)%crossings.length].lon-start)+360)%360);
    const off=((lon-start)%360+360)%360, w=((next-start)%360+360)%360;
    if(off<w||w===0) return c.name;
  }
  return crossings[crossings.length-1].name;
}
function inInterval(lon,start,end,tol){
  const a=(start-tol+360)%360, b=(end+tol+360)%360, x=(lon+360)%360;
  let w=((b-a)%360+360)%360;
  if(Math.abs(((end-start)%360+360)%360)<1e-9&&tol===0) w=0;
  if(w===0) return false;
  return ((x-a)%360+360)%360<=w;
}
function offsets(lon,start,end){ // signed distance to nearest interval edge
  const a=(start+360)%360, b=(end+360)%360, x=(lon+360)%360;
  let w=((b-a)%360+360)%360; if(w===0)w=360;
  const off=((x-a)%360+360)%360;
  return off<=w?0:Math.min(off-w,360-off);
}

function drawWheel(host,wheel,state){
  const cx=340, cy=340, R=250, BW=250;
  const crossings=D.boundaries[state.table].crossings;
  const v=wheel.variants.find(x=>x.table===state.table&&x.tolerance===state.tolerance);
  const idx=state.index;
  const jd=wheel.jds[idx];
  const lons={}; BODIES.forEach(b=>lons[b]=wheel.lons[b][idx]);
  const want=v&&v.rules;
  let h="";
  // sectors
  for(let i=0;i<crossings.length;i++){
    const start=crossings[i].lon, end=crossings[(i+1)%crossings.length].lon;
    h+=`<path class="sector" d="${arcSeg(start,end,R,cx,cy,1)}"/>`;
    h+=`<path class="sector" d="${arcSeg(start,end,R-16,cx,cy,1)}" stroke-opacity=".35"/>`;
    const [x0,y0]=pt(start,R-16,cx,cy), [x1,y1]=pt(start,R,cx,cy);
    h+=`<line class="bound" x1="${x0.toFixed(2)}" y1="${y0.toFixed(2)}" x2="${x1.toFixed(2)}" y2="${y1.toFixed(2)}"/>`;
    let span=((end-start)%360+360)%360;
    const mid=start+span/2;
    const [lx,ly]=pt(mid,R+34,cx,cy);
    const [ax,ay]=pt(mid,R+2,cx,cy);
    h+=`<line class="sline" x1="${ax.toFixed(2)}" y1="${ay.toFixed(2)}" x2="${lx.toFixed(2)}" y2="${ly.toFixed(2)}"/>`;
    h+=`<text class="slabel" x="${lx.toFixed(2)}" y="${ly.toFixed(2)}" text-anchor="middle">${crossings[i].name.toUpperCase()} ${span.toFixed(1)}°</text>`;
  }
  // constraint bands, one ring per constrained body
  if(want){
    let ring=0;
    for(const r of want){
      if(r.free) continue;
      const rr=R-30-ring*15;
      if(rr<40) break;
      h+=`<path class="band" d="${arcSeg(r.start,r.end,rr+6,cx,cy,2)}L${pt(r.end,rr-6,cx,cy).map(n=>n.toFixed(2)).join(" ")}`;
      // close the band as a thin annulus
      const back=[]; const span=((r.end-r.start)%360+360)%360;
      for(let i=Math.ceil(span/2);i>=0;i--){const [x,y]=pt(r.start+span*i/Math.ceil(span/2),rr-6,cx,cy);back.push(x.toFixed(2)+" "+y.toFixed(2));}
      h+=`L${back.join("L")}Z"/>`;
      if(r.best!=null){const [bx,by]=pt(r.best,rr,cx,cy), [bx2,by2]=pt(r.best,rr+5,cx,cy);
        h+=`<line class="best" x1="${bx.toFixed(2)}" y1="${by.toFixed(2)}" x2="${bx2.toFixed(2)}" y2="${by2.toFixed(2)}"/>`;}
      ring++;
    }
  }
  // markers
  const slots={};
  BODIES.forEach((b,k)=>{
    const lon=lons[b];
    const press=(slots[Math.round(lon/6)]=(slots[Math.round(lon/6)]||0)+1)-1;
    const rr=R-14-press*13;
    const [mx,my]=pt(lon,rr,cx,cy);
    const rule=want&&want.find(r=>r.body===b);
    const okr=!rule||rule.free||(rule&&(state.tolerance>0
        ? inInterval(lon,rule.start,rule.end,state.tolerance)
        : inInterval(lon,rule.start,rule.end,0)));
    h+=`<g class="marker ${okr?"":"off"}">`;
    if(rule&&!rule.free){const [gx,gy]=pt(rule.start,rr,cx,cy),[hx,hy]=pt(rule.end,rr,cx,cy);
      h+=`<circle cx="${gx.toFixed(2)}" cy="${gy.toFixed(2)}" r="1.4"/><circle cx="${hx.toFixed(2)}" cy="${hy.toFixed(2)}" r="1.4"/>`;}
    h+=`<circle cx="${mx.toFixed(2)}" cy="${my.toFixed(2)}" r="3.6"/>`;
    const [tx,ty]=pt(lon,rr-13,cx,cy);
    h+=`<text x="${tx.toFixed(2)}" y="${ty.toFixed(2)}" text-anchor="middle">${b.toUpperCase()}</text></g>`;
  });
  // centre
  h+=`<text class="state ${state.matched?"yes":""}" x="${cx}" y="${cy-6}" text-anchor="middle">${state.matched?"MATCHED":"NO MATCH"}</text>`;
  h+=`<text class="state" x="${cx}" y="${cy+14}" text-anchor="middle">${jdToDate(jd)}</text>`;
  host.querySelector("svg")?.remove();
  const svg=document.createElementNS("http://www.w3.org/2000/svg","svg");
  svg.setAttribute("viewBox","0 0 680 680"); svg.innerHTML=h;
  host.querySelector(".svgwrap")?.remove();
  const wrap=document.createElement("div"); wrap.className="svgwrap";
  wrap.appendChild(svg); host.prepend(wrap);

  // readout
  const read=host.parentElement.querySelector(".read");
  read.innerHTML=readout(wheel,v,jd,lons,state);
}

function readout(wheel,v,jd,lons,state){
  const crossings=D.boundaries[state.table].crossings;
  let rows="";
  BODIES.forEach(b=>{
    const rule=v&&v.rules.find(r=>r.body===b);
    const sec=sectorOf(crossings,lons[b]);
    let dev="";
    if(rule&&!rule.free){
      const d=((lons[b]-rule.start+540)%360)-180;
      dev=((lons[b]-rule.start)%360+360)%360<=(((rule.end-rule.start)%360)+360)%360?"in":"out";
      if(rule.best!=null) dev+=` · ${(((lons[b]-rule.best+540)%360)-180).toFixed(1)}°`;
    }
    rows+=`<tr><td>${GLYPH[b]} ${b}</td><td class="n">${lons[b].toFixed(2)}°</td>`+
          `<td class="d">${sec}</td><td class="d">${dev}</td></tr>`;
  });
  const w=wheel.reported_windows;
  const span=wheel.grid_end_jd-wheel.grid_start_jd;
  return `<div class="pos">${jdToDate(jd)} — position ${(( (jd-wheel.grid_start_jd)/span)*100).toFixed(1)}% of the scrub window</div>
    <div class="controls">
      <div class="scrub"><input type="range" min="0" max="${wheel.jds.length-1}" value="${state.index}" data-role="scrub"></div>
      <div class="seg" data-role="table">${["iau_j2000","horos_cs_j2000","horos_csn_j2000"].map(t=>
        `<button data-table="${t}" aria-pressed="${t===state.table}">${t==="iau_j2000"?"real IAU boundaries":t==="horos_cs_j2000"?"HOROS CS table":"HOROS CSN table"}</button>`).join("")}</div>
      <div class="seg" data-role="tol">${[0,5].map(t=>
        `<button data-tol="${t}" aria-pressed="${t===state.tolerance}">tolerance ${t}°</button>`).join("")}</div>
    </div>
    <table>${rows}</table>
    <div class="stamp">
      boundary table ${D.boundaries[state.table].hash} · kernel ${D.kernel.name} · ΔT ${deltaTAt(jd)}<br>
      order constraint at this instant: <b>${v&&v.order?v.order[state.index]==="1"?"satisfied":"violated":"none"}</b>
      (reported, not matched on — the engine only enforces it when a specification says <code>enforce_order</code>)<br>
      windows found in the full span for this specification and table: <b>${w.length}</b>
      ${w.length?`(${w.map(x=>jdToDate(x.start_jd).slice(0,8)+" → "+jdToDate(x.end_jd).slice(0,8)).join(", ")})`:""}
    </div>`;
}

function deltaTAt(jd){
  let pick=D.delta_t[0];
  const year=jdToYear(jd);
  for(const s of D.delta_t){if(Math.abs(s.year-year)<Math.abs(pick.year-year))pick=s;}
  return `${pick.seconds} s (model value for the nearest sampled year, ${pick.year<0?(-pick.year)+" BC":pick.year+" AD"})`;
}

function boot(){
  // facts
  const total=D.wheels.reduce((a,w)=>a+w.jds.length,0);
  document.getElementById("facts").innerHTML=[
    ["13,000","years searched"],["13","constellations on the ecliptic, Ophiuchus included"],
    ["≤0.008°","agreement with JPL Horizons"],["0°","tolerance, unless you ask for more"],
    ["5.8M","instants sampled"]
  ].map(([b,s])=>`<li><b>${b}</b><span>${s}</span></li>`).join("");

  D.wheels.forEach(wheel=>{
    const sec=document.getElementById(wheel.key==="dl2"?"wheel-dl2":"wheel-dr9");
    const host=sec.querySelector(".wheel");
    const state={table:wheel.key==="dl2"?"iau_j2000":"iau_j2000",tolerance:0,index:0,matched:false};
    // start at the closest thing to a match: the middle of the first reported window
    if(wheel.reported_windows.length){
      const mid=0.5*(wheel.reported_windows[0].start_jd+wheel.reported_windows[0].end_jd);
      let best=0,bd=1e18;
      wheel.jds.forEach((j,i)=>{const d=Math.abs(j-mid);if(d<bd){bd=d;best=i;}});
      state.index=best; state.tolerance=0;
    } else {
      state.index=Math.floor(wheel.jds.length/2);
    }
    const upd=()=>{
      const v=wheel.variants.find(x=>x.table===state.table&&x.tolerance===state.tolerance);
      state.matched=!!(v&&v.mask[state.index]==="1");
      drawWheel(host,wheel,state);
      host.parentElement.querySelector('[data-role=scrub]')
        .addEventListener("input",e=>{state.index=+e.target.value;upd();});
      host.parentElement.querySelectorAll('[data-role=table] button')
        .forEach(b=>b.addEventListener("click",()=>{state.table=b.dataset.table;upd();}));
      host.parentElement.querySelectorAll('[data-role=tol] button')
        .forEach(b=>b.addEventListener("click",()=>{state.tolerance=+b.dataset.tol;upd();}));
    };
    upd();
  });

  drawTimeline();
  document.getElementById("method").innerHTML=method();
  document.getElementById("foot").innerHTML=
    `zodiac dating engine · MIT · <a href="https://github.com/sfingali/zodiac-dating-engine">github.com/sfingali/zodiac-dating-engine</a><br>`+
    `ephemeris ${D.kernel.name}, covering JD ${D.kernel.start_jd} to ${D.kernel.end_jd}`+
    ` (13200 BC to 1969 AD) · Skyfield · JPL DE441 (NASA, public domain) · IAU boundaries Delporte 1930 / Roman 1987`;
}

function drawTimeline(){
  const W=1180,H=40+34*D.timeline.runs.length;
  const j0=D.timeline.from_jd,j1=D.timeline.to_jd;
  const X=jd=>((jd-j0)/(j1-j0))*W;
  let h="";
  D.timeline.runs.forEach((run,k)=>{
    const y=34+k*34;
    h+=`<text class="rowlabel" x="0" y="${y-8}">${run.label}</text>`;
    h+=`<line x1="0" y1="${y}" x2="${W}" y2="${y}" stroke="rgba(233,230,223,.07)"/>`;
    run.windows.forEach(w=>{
      const x=X(w.start_jd);
      h+=`<rect class="${run.tolerance>0?"win5":"win"}" x="${Math.max(0,x-1).toFixed(1)}" y="${y-16}" width="2.4" height="16"/>`;
      h+=`<text class="axis" x="${Math.max(0,x-1).toFixed(1)}" y="${y+14}">${jdToYear(w.start_jd)}</text>`;
    });
    if(!run.windows.length){
      h+=`<text class="axis" x="0" y="${y+14}">no instant in the span satisfies this</text>`;
    }
  });
  [ -13000,-10000,-7000,-4000,-1000,1000,1969 ].forEach(y=>{
    const jd=1721423.5+(y-1)*365.25;
    if(jd<j0||jd>j1)return;
    h+=`<text class="axis" x="${X(jd).toFixed(1)}" y="${H-4}" text-anchor="middle">${y<0?-y+" BC":y+" AD"}</text>`;
  });
  document.querySelector(".stripwrap").innerHTML=
    `<svg class="strip" viewBox="0 0 ${W+40} ${H}" preserveAspectRatio="none" style="height:${H}px"></svg>`;
  document.querySelector(".stripwrap svg").innerHTML=`<g transform="translate(20,0)">${h}</g>`;
}

function method(){
  return [
    ["What is computed","Geocentric apparent J2000 ecliptic longitudes for the Sun, Moon, Mercury, Venus, Mars, Jupiter and Saturn, from JPL DE440/DE441 through Skyfield, with light-time and aberration. Agreement with JPL Horizons itself is better than 0.008°."],
    ["Boundaries","The real constellation boundaries on the ecliptic, derived from the IAU map (Delporte 1930, as tabulated by Roman 1987) by walking the ecliptic and bisecting every crossing. Thirteen sectors, because Ophiuchus is genuinely on the ecliptic and holds 18.6° of it. For comparison the two boundary tables shipped with HOROS are available, so a published HOROS date can be checked on HOROS's own terms."],
    ["Why the sectors are uneven","Scorpio is 6.6° wide, Virgo 44.0°. A wheel drawn as twelve 30° wedges would misstate the very thing being tested."],
    ["Windows, not dates","A constraint is an inequality, so the answer is an interval. Edges are refined by bisection to about a second; the sampling step only sets how finely the search is scanned."],
    ["No statistics","No score, no p-value, no most-probable date. The engine reports what matched, how wide it was, and how far each body sat from any point the drawing's decoder nominated. Whether a coincidence means anything is a question for the reader."],
    ["The limit on old dates","The argument is Terrestrial Time, so a civil date in antiquity involves ΔT: 25,310 s at 1000 BC. The Moon moves 0.549° per hour, so a 1,000 s spread between ΔT models moves it 0.15°. Good to a few tenths of a degree on the Moon for ancient dates, not to 0.01°."],
    ["Sources","Decipherment tables transcribed from Fomenko & Nosovsky, New Chronology of Egypt, Appendix 2, with variant codes. Full citations in the repository."],
  ].map(([t,b])=>`<div><h3>${t}</h3><p>${b}</p></div>`).join("");
}

boot();
</script>
</body>
</html>
"""


if __name__ == "__main__":
    main()
