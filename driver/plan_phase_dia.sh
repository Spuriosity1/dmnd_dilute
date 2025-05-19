#!/bin/bash
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )


# Check if the user provided an input file
if [ "$#" -lt 3 ]; then
    #                1                 2              3       4   5
    echo "Usage: $0 <seeds_per_point>  <L> <db_repo> <auxiliary args>"
    exit 1
fi


seeds_per_point="$1"
L="$2"
DB_REPO="${3}"
auxargs="${4}"


if [ ! -d "${DB_REPO}" ]; then
    echo "Error: specified db repo does not exist. Create?"
    read -p "Continue (y/n)?" choice
    case "$choice" in 
      y|Y ) echo "yes"; mkdir -p "${DB_REPO}";;
      n|N ) echo "no"; exit 1;;
      * ) echo "invalid"; exit 1;;
    esac
fi


for p in `seq 0 0.01 0.1`; do
    for i in `seq 0 "${seeds_per_point}"`; do

	seed=`openssl rand -hex 4`
        echo ../build/dmnd_dilute $L 0 0 0 $L 0 0 0 $L -p $p -o "${DB_REPO}" --seed "${seed}" $auxargs
        let index+=1
    done
done
