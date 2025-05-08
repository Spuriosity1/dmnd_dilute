
tmp="../../tmp"

outfile="$(date -I)_$(hostname)_benchmark.csv"

mkdir -p $tmp
rm $tmp/*

for L in `seq 1 30`; do
    start=`perl -MTime::HiRes=time -e 'printf "%.9f\n", time'`
    
    ../build/dmnd_dilute $L 0 0 0 $L 0 0 0 $L -p 0.1 -o $tmp --seed 2bd1dde03c3db836 -n 2 4 -f >/dev/null

    end=`perl -MTime::HiRes=time -e 'printf "%.9f\n", time'`
    runtime=$( echo "$end - $start" | bc -l )
    echo $runtime
    echo $L,$(( L*L*L )),$runtime >> $outfile
done
