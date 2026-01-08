#!/usr/bin/env python3
"""
run_ranger_no_logs.py - Ranger-DTL runner without log files, ensuring 100 replicates per combination
"""

import argparse
import os
import subprocess
import sys
from multiprocessing import Pool
import time

def generate_l_values(l_center, l_range, l_step):
    """Generate L values list"""
    l_min = l_center - l_range
    l_max = l_center + l_range
    l_list = list(range(l_min, l_max + 1, l_step))
    l_list = [l for l in l_list if l >= 1]
    l_list.sort()
    
    if not l_list:
        raise ValueError("No valid L values (>=1) generated!")
    
    return l_list

def run_bootstrap_task(args):
    """Run a bootstrap task without log files"""
    d_val, l_val, rep, seed, idx, total, init_file, t_val, ranger_bin, output_prefix = args
    
    output_file = f"{output_prefix}.D{d_val}.T{t_val}.L{l_val}.{rep:03d}"
    
    print(f"[{idx}/{total}] D={d_val} L={l_val} Rep={rep:03d}")
    
    cmd = [
        ranger_bin,
        "-i", init_file,
        "-o", output_file,
        "-D", str(d_val),
        "-T", str(t_val),
        "-L", str(l_val),
        "--seed", str(seed)
    ]
    
    try:
        # Run without capturing output (no log files)
        result = subprocess.run(cmd, capture_output=False, check=False)
        
        if result.returncode == 0:
            return True
        else:
            print(f"  WARNING: Task failed with exit code {result.returncode}")
            return False
            
    except Exception as e:
        print(f"  ERROR: Exception running task: {e}")
        return False

def ensure_replicates(args):
    """Ensure each parameter combination has exactly the required number of replicates"""
    d_values = [int(d.strip()) for d in args.d_values.split(',')]
    l_list = generate_l_values(args.l_center, args.l_range, args.l_step)
    
    # Check existing files to see which replicates are missing
    tasks = []
    total_needed = 0
    total_completed = 0
    
    for d_val in d_values:
        for l_val in l_list:
            # Check how many replicates already exist for this combination
            existing_reps = set()
            pattern = f"{args.output_prefix}.D{d_val}.T{args.t_val}.L{l_val}."
            
            # Look for existing output files
            for filename in os.listdir('.'):
                if filename.startswith(pattern):
                    try:
                        # Extract replicate number from filename
                        rep_num = int(filename.split('.')[-1])
                        existing_reps.add(rep_num)
                    except:
                        pass
            
            # Add tasks for missing replicates
            needed_reps = set(range(1, args.replicates + 1)) - existing_reps
            total_needed += len(needed_reps)
            total_completed += len(existing_reps)
            
            for rep in needed_reps:
                seed = d_val * 1000000 + l_val * 1000 + rep
                tasks.append((d_val, l_val, rep, seed))
    
    print(f"Replicates status: {total_completed} completed, {total_needed} needed")
    return tasks, total_needed

def main():
    parser = argparse.ArgumentParser(
        description="Ranger-DTL runner without log files - ensures 100 replicates per combination"
    )
    
    # Required arguments
    parser.add_argument('-i', '--input', required=True, help='Initial input file')
    parser.add_argument('-l', '--l_center', type=int, required=True, help='L center value')
    
    # Optional arguments with defaults
    parser.add_argument('-j', '--jobs', type=int, default=8, help='Number of parallel jobs (default: 8)')
    parser.add_argument('-r', '--replicates', type=int, default=100, 
                       help='Number of bootstrap replicates per combination (default: 100)')
    parser.add_argument('--l_range', type=int, default=75, 
                       help='L value range around center (default: 75)')
    parser.add_argument('--l_step', type=int, default=25, 
                       help='L value step size (default: 25)')
    parser.add_argument('--d_values', type=str, default='400,600,800,1000',
                       help='Comma-separated D values (default: 400,600,800,1000)')
    parser.add_argument('-t', '--t_val', type=int, default=600,
                       help='T parameter value (default: 600)')
    parser.add_argument('--ranger_bin', default='Ranger-DTL.linux',
                       help='Path to Ranger-DTL executable (default: Ranger-DTL.linux)')
    parser.add_argument('--output_prefix', default='ranger',
                       help='Output file prefix (default: ranger)')
    parser.add_argument('--resume', action='store_true',
                       help='Resume from existing files (skip completed replicates)')
    
    args = parser.parse_args()
    
    # Validate input file
    if not os.path.isfile(args.input):
        print(f"ERROR: Input file '{args.input}' not found!")
        sys.exit(1)
    
    # Validate Ranger executable
    if not os.path.isfile(args.ranger_bin) or not os.access(args.ranger_bin, os.X_OK):
        which_result = subprocess.run(['which', args.ranger_bin], capture_output=True, text=True)
        if which_result.returncode != 0:
            print(f"ERROR: Ranger executable '{args.ranger_bin}' not found or not executable!")
            print("Please specify full path using --ranger_bin")
            sys.exit(1)
        args.ranger_bin = which_result.stdout.strip()
    
    # Parse D values
    try:
        d_values = [int(d.strip()) for d in args.d_values.split(',')]
    except ValueError:
        print("ERROR: Invalid D values format. Use comma-separated integers.")
        sys.exit(1)
    
    # Generate L values
    try:
        l_list = generate_l_values(args.l_center, args.l_range, args.l_step)
    except ValueError as e:
        print(f"ERROR: {e}")
        sys.exit(1)
    
    # Display configuration
    print("=========================================")
    print(f"Input file   : {args.input}")
    print(f"L center     : {args.l_center} (±{args.l_range}, step {args.l_step})")
    print(f"L values     : {l_list}")
    print(f"D values     : {d_values}")
    print(f"T value      : {args.t_val}")
    print(f"Replicates   : {args.replicates} per combination")
    print(f"Total combos : {len(d_values) * len(l_list)}")
    print(f"Total tasks  : {len(d_values) * len(l_list) * args.replicates}")
    print(f"Parallel     : {args.jobs} jobs")
    print(f"Ranger bin   : {args.ranger_bin}")
    print(f"Resume mode  : {args.resume}")
    print("=========================================")
    print()
    
    # Prepare task list
    if args.resume:
        print("Resume mode: checking existing files...")
        tasks, total_needed = ensure_replicates(args)
        if total_needed == 0:
            print("All replicates already completed! Nothing to do.")
            return
    else:
        tasks = []
        for d_val in d_values:
            for l_val in l_list:
                for rep in range(1, args.replicates + 1):
                    seed = d_val * 1000000 + l_val * 1000 + rep
                    tasks.append((d_val, l_val, rep, seed))
        total_needed = len(tasks)
    
    print(f"Running {total_needed} tasks...")
    
    # Prepare arguments for parallel processing
    task_args = []
    for idx, (d_val, l_val, rep, seed) in enumerate(tasks, 1):
        task_args.append((
            d_val, l_val, rep, seed, idx, total_needed,
            args.input, args.t_val, args.ranger_bin, args.output_prefix
        ))
    
    # Run tasks in parallel
    success_count = 0
    with Pool(processes=args.jobs) as pool:
        results = pool.map(run_bootstrap_task, task_args)
        success_count = sum(results)
    
    # Final status
    failed_count = total_needed - success_count
    print(f"\nCompleted: {success_count}/{total_needed} successful, {failed_count} failed")
    
    if args.resume:
        # Calculate overall completion rate
        d_values = [int(d.strip()) for d in args.d_values.split(',')]
        l_list = generate_l_values(args.l_center, args.l_range, args.l_step)
        total_expected = len(d_values) * len(l_list) * args.replicates
        
        # Count actually completed files
        actual_completed = 0
        for d_val in d_values:
            for l_val in l_list:
                pattern = f"{args.output_prefix}.D{d_val}.T{args.t_val}.L{l_val}."
                for filename in os.listdir('.'):
                    if filename.startswith(pattern):
                        actual_completed += 1
        
        print(f"Overall completion: {actual_completed}/{total_expected} ({actual_completed/total_expected*100:.1f}%)")
    
    print("=========================================")

if __name__ == "__main__":
    main()