#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import pp_eqb_1 as eqb
import os, subprocess

pdb_id = "mA18"
new_job = "test_"+pdb_id
jobid, new_job = eqb.equillibration(pdb_id, new_job)    # Start Equlibration
os.chdir('/home/jji110/parch_platform/')
cmd = "sbatch --depend=afterok:{0} annealing_py.sh {1} {2}".format(jobid, new_job, jobid)  # Start annealing job with dependance on eqb
status, ann_jobid = subprocess.getstatusoutput(cmd)

print("Eqb: ", jobid)
print("Ann: ", ann_jobid)

