#!/bin/bash
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )


# Check if the user provided an input file
if [ "$#" -lt 4 ]; then
    #                1                 2              3       4   5
    echo "Usage: $0 <seeds_per_point> <inital_seedid> <L> <db_repo> <auxiliary args>"
    exit 1
fi


seeds_per_point="$1"
index="$2"
L="$3"
DB_REPO="${4}"

seedfile="${SCRIPT_DIR}/seeds.txt"

if [ ! -d "${DB_REPO}" ]; then
    echo "Error: specified db repo does not exist. Create?"
    read -p "Continue (y/n)?" choice
    case "$choice" in 
      y|Y ) echo "yes"; mkdir -p "${DB_REPO}";;
      n|N ) echo "no"; exit 1;;
      * ) echo "invalid"; exit 1;;
    esac
fi

max_seed=`wc -l < "${seedfile}"`

for p in `seq 0 0.01 0.1`; do
    for i in `seq 0 "${seeds_per_point}"`; do

        if [ $index -ge "${max_seed}" ]; then
            echo "Reached end of seed file!"
            exit 1
        fi
        seed=`sed "${index}q;d" "${seedfile}" | sed 's/ //g'` 
        echo ../build/dmnd_dilute $L 0 0 0 $L 0 0 0 $L -p $p -o "${DB_REPO}" --seed "${seed:16}" $5
        let index+=1
    done
done
