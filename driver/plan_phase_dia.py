#!/usr/bin/env python3
import os
import sys
import numpy as np
import secrets
import argparse

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Run dmnd_dilute with varying parameters')
    parser.add_argument('seeds_per_point', type=int, help='Number of seeds per probability point')
    parser.add_argument('db_repo', type=str, help='Database repository directory')
    parser.add_argument('--aux', type=str, default="", help='Auxiliary arguments for dmnd_dilute')
    parser.add_argument('-n', '--delete_nn', nargs='+', type=str, metavar='ARG',
                       help='Arguments to be appended as "-n ARG1 [ARG2...]"')
    parser.add_argument('-P', '--prob_range', nargs=3, type=float, default=[0.0, 0.01, 0.1],
                       metavar=('START', 'STEP', 'STOP'),
                       help='Probability range as START STEP STOP (inclusive of both terminals)')
    parser.add_argument('-L', '--L_range', nargs=3, type=int, default=[0.0, 0.01, 0.1],
                       metavar=('START', 'STEP', 'STOP'),
                       help='System linear size range as START STEP STOP (inclusive of both terminals)')
    
    # Parse arguments
    args = parser.parse_args()
    
    # Get the directory of the script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Check if the db_repo directory exists, create if needed
    if not os.path.isdir(args.db_repo):
        print(f"Error: specified db repo '{args.db_repo}' does not exist. Create?")
        while True:
            choice = input("Continue (y/n)? ")
            if choice.lower() in ['y', 'yes']:
                print("yes")
                os.makedirs(args.db_repo, exist_ok=True)
                break
            elif choice.lower() in ['n', 'no']:
                print("no")
                sys.exit(1)
            else:
                print("invalid")
     

    p_start, p_step, p_stop = args.prob_range
    p_values = np.arange(p_start, p_stop + 0.5 * p_step, p_step)

    L_start, L_step, L_stop = args.L_range
    L_values = np.arange(L_start, L_stop + L_step, L_step, dtype=int)

    # Prepare delete_nn arguments if provided
    delete_nn_args = ""
    if args.delete_nn:
        delete_nn_args = "-n " + " ".join(args.delete_nn)

    # Loop over probability values from min_p to max_p with steps of p_step
    for L in L_values:
        for _p in p_values:
            p = "%.6f" % _p
            
            # For each probability value, run the command seeds_per_point+1 times
            for i in range(args.seeds_per_point + 1):
                # Generate a random seed (equivalent to openssl rand -hex 4)
                seed = secrets.token_hex(4)
                
                # Construct the command
                cmd = f"../build/dmnd_dilute {L} 0 0 0 {L} 0 0 0 {L} -p {p} -o {args.db_repo} --seed {seed}"
                if delete_nn_args:
                    cmd += f" {delete_nn_args}"

                if args.aux:
                    cmd += f" {args.aux}"
                
                # Print the command
                print(cmd)
                

if __name__ == "__main__":
    main()
