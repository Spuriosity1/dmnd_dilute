
tmp="../../tmp"

outfile="$(date -I)_$(hostname)_benchmark.csv"

mkdir -p $tmp

for L in `seq 1 20`; do
    start=`date +%s.%N`
    
    ../build/dmnd_dilute $L 0 0 0 $L 0 0 0 $L -p 0.01 -o $tmp --seed 2bd1dde03c3db836 --neighbours 2 4

    end=`date +%s.%N`
    runtime=$( echo "$end - $start" | bc -l )

    echo $L $runtime >> $outfile
done
