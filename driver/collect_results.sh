#!/bin/bash

datestr=$(date +"%Y-%m-%dT%H-%M")

CMD="python3 ../scripts/merge_to_sql.py --db ../../out_db/$datestr-run.db /scratch/alaric/out/percolate --cleanup"
echo $CMD
$CMD

