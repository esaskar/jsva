#!/usr/bin/env python2
"""
Run a command in the cluster using SLURM. 
"""

import sys, os, subprocess, argparse, datetime

SLURM_TEMPLATE = """#!/bin/bash
#SBATCH --job-name=%(job_name)s
#SBATCH --cpus-per-task=%(threads)s
#SBATCH --mem=%(vmem)s
#SBATCH --partition=%(queue)s
#SBATCH --workdir=%(work_dir)s
#SBATCH --output=%(output_file)s
#SBATCH --error=%(error_file)s
#SBATCH --open-mode=truncate
#SBATCH --requeue
%(extra)s
%(job)s
"""

DEFAULT_MEM = 10        # maximum memory in MBs
DEFAULT_QUEUE = "all"
DEFAULT_THREADS = 1

def run_cmd(cmd, job_name = None, threads = DEFAULT_THREADS, vmem = DEFAULT_MEM, queue = DEFAULT_QUEUE,
            out_fn = None, err_fn = None, job_fn = None):
    """
vmem in MBs
"""
    now = datetime.datetime.now()
    tag = now.strftime("%Y%m%d_%H%M%S_%f")

    if not job_name:
        job_name = tag

    if not out_fn:
        out_fn = "%s.out" % (tag)
    if not err_fn:
        err_fn = "%s.err" % (tag)

    work_dir = os.getcwd()

    params = {"job_name" : job_name,
              "threads" : threads,
              "vmem" : vmem,
              "queue" : queue,
              "work_dir" : work_dir,
              "output_file" : out_fn,
              "error_file" : err_fn,
              "extra" : "#",
              "job" : cmd}

    if not job_fn:
        job_fn = "%s.job" % (tag)
    
    jobf = open(job_fn, "w")
    jobf.write(SLURM_TEMPLATE % (params))
    jobf.flush()

    rc = subprocess.call(["sbatch", jobf.name])
    if rc != 0:
        sys.stderr.write("sbatch failed: %d\n" % (rc))
        sys.exit(rc)

def __main(args):
    run_cmd(" ".join(args.cmd))

if __name__ == "__main__":
     p = argparse.ArgumentParser()
     p.add_argument("cmd", nargs = argparse.REMAINDER, help = "Command to run in the cluster")
     args = p.parse_args()
     if len(args.cmd) == 0:
         p.print_help(sys.stderr)
         sys.exit(2)
     __main(args)