#!/bin/bash

datestr=$(date +"%Y-%m-%dT%H-%M")

CMD="python3 ../scripts/merge_to_sql.py --db ../../out_db/$datestr-run.db ../../out --cleanup"
echo $CMD
$CMD

