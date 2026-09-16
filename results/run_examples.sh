#!/bin/bash
# Run every bundled example over the full DE441 part-1 span and keep the reports.
cd /opt/data/repos/zodiac-dating-engine || exit 1
export ZODIAC_DATING_KERNEL=/opt/data/fomenko_check/sescc/de441_part-1.bsp
mkdir -p /opt/data/tmp/examples
for f in examples/*.json; do
  b=$(basename "$f" .json)
  echo "=== $b  $(date -u +%H:%M:%S) ==="
  start=$(date +%s)
  env -u PYTHONPATH /tmp/eph3/bin/python -m zodiac_dating run "$f" \
      --json-out "/opt/data/tmp/examples/$b.json" > "/opt/data/tmp/examples/$b.txt" 2>&1
  rc=$?
  dur=$(( $(date +%s) - start ))
  echo "  exit $rc in ${dur}s -> $(grep -E 'windows found|instants matching' "/opt/data/tmp/examples/$b.txt" | tr '\n' ' ')"
done
echo "ALL EXAMPLES DONE"
