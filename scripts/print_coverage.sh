#!/bin/bash

set -euo pipefail

echo -e "strain\tpct_coverage"

for file in $@
do
	NAME=$(basename $file .tsv)
	PCT=$(awk -F "\t" 'NR > 1 {t += $3} {s += $5} END {print s/t}' $file)
	echo -e "${NAME}\t${PCT}"
done
