#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 4 14:12:19 2022

@author: jji110
"""


import os, subprocess
import numpy as np
import math
import shutil
import json

import MDAnalysis as mda
import MDAnalysis.analysis.distances as Mdistance


def equillibration(pdb_id, new_job):
    user_account = "/home/jji110/"
    #path         = user_account + "parch_platform/"
    path_pp_bank = user_account + "parch_platform/backup_pp_bank/"
    path_src     = user_account + "parch_platform/backup_eqb/"    # path_src: source directory
    path_dst     = user_account + "parch_platform/"+new_job+"/"   # path_dst: destination directory.    The directory "new_job" is being created when copying the content (files and directories) of backup_dir 
    pdb_fname = pdb_id+'.pdb'    
    shutil.copytree(path_src, path_dst)   # At the same time, create the new_job directory. Copy the files and directories of path_src into path_dst.  shutil.copy() & shutil.copyfile()
    shutil.copy(path_pp_bank+pdb_fname, path_dst)   # file => directory
    os.chdir(path_dst)  
    
    os.system("echo q|gmx make_ndx -f "+pdb_fname+" -o c1.ndx")   # to get rid of ions (CL, NA) and H2O (HOH)
    os.system("echo 1|gmx editconf -f "+pdb_fname+" -o pp_clean.pdb -n c1.ndx")  
    
    ############# check the number of protein moleucles ##############
    
    fin = open('pp_clean.pdb', "r")
    content = fin.readlines()
    fin.close()
    
    ids = []       # chain_id info of each atom
    resids = []    # resid info of each atom
    for line in range(0,len(content)):
        if (content[line][:4] == "ATOM"):
            row = content[line].strip().split()
            if (len(row[4])==1):       
                ids.append(row[4])           # the 4th column is chain_ID (character) info
                resids.append(int(row[5]))   # the 5th column is resid (integer) info
            else:
                ids.append(row[4][0])
                resids.append(int(row[4][1:]))
    res,ind = np.unique(ids, return_index=True)
    chain_id = res[np.argsort(ind)]
    
    resid_dict = {}     
    for ii in range(0,len(ids)):
        cha = ids[ii]
        if cha in resid_dict.keys():
            resid_dict[cha] = resid_dict[cha] + [resids[ii]]   # [resids[ii]]: one integer in list.
        else:
            resid_dict[cha] = [resids[ii]]     # Initially, "key(chain_id, character)" and "value(resid, integer)" of "atom_based"
    
    len_dict = len(list(resid_dict.keys()))
    for jj in range(0,len_dict):
        cha = list(resid_dict.keys())[jj]
        resid_dict[cha] = np.unique(np.array(resid_dict[cha]))   # resid_dict is one array
        resid_dict[cha] = list(resid_dict[cha])                  # resid_dict is one list    
    print("resid_dict = ",resid_dict)          # Then, "key(chain_id, character)" and "value(resid, integer)" of "residue_based"     
    
    ###### save the dict data (must be in character or string) ######
    #################################################################
    
    min_aa = resid_dict[list(resid_dict.keys())[0]][0]   # min_aa: the resid of "the 1st residue" in the 1st Protein molecule
    max_bb = resid_dict[list(resid_dict.keys())[0]][-1]  # max_bb: the resid of "the last residue" in the 1st Protein molecule
    
    for kk in range(0,len_dict):
        cha = list(resid_dict.keys())[kk]
        min_aa  = min(min_aa,resid_dict[cha][0])
        max_bb  = max(max_bb,resid_dict[cha][-1])
        for mn in range(0,len(resid_dict[cha])):
            resid_dict[cha][mn] = str(resid_dict[cha][mn])       # this commandline only works when "resid_dict[cha]" is one "list", instead of one "array"
    print("resid_dict = ",resid_dict)          # Finally, "key(chain_id, character)" and "value(resid, character)" of residue_based      
        
    fout = open("resid_dict.txt", "w")
    fout.write(json.dumps(resid_dict))         # json.dumps(): convert "python dict (character or string, integer does not work)" into "json data"
    fout.close()
    
    ################################################################
    
    os.system("echo 2|gmx pdb2gmx -f pp_clean.pdb -o pp_clean.gro -water tip3p -ignh")   
    os.system("gmx editconf -f pp_clean.gro -o pp_newbox.gro -c -d 1.5 -bt cubic")
    os.system("gmx solvate -cp pp_newbox.gro -cs spc216.gro -o pp_solv.gro -p topol.top")
    
    os.system("gmx grompp -f ions.mdp -c pp_solv.gro -p topol.top -o ions.tpr")  # neutralize the system
    os.system("echo 13|gmx genion -s ions.tpr -o pp_solv_ions.gro -p topol.top -pname NA -nname CL -neutral")   # 13: SOL, 12: water
    os.system("gmx grompp -f ions.mdp -c pp_solv_ions.gro -p topol.top -o ions.tpr")  # 0.15 mol NaCl
    os.system("echo 13|gmx genion -s ions.tpr -o final_solv_ions.gro -p topol.top -conc 0.15")   # 13: SOL, 12: water
    os.system("gmx make_ndx -f final_solv_ions.gro -o index.ndx < INP_ion")   # based on the ions.tpr info, add one loop
    os.system("mv topol.top old_topol.top")     
     
       
    if (len(chain_id)==1):        
        os.system("sed '1,7d' posre.itp > p_1.itp")
        os.system("awk '{printf(\"%6s%6s\\n\",$1,$2)}' p_1.itp > p_2.itp")   # formated output the first two columns
        os.system("awk '{printf(\"%8s%8s%8s\\n\",$3,$4,$5)}' p_1.itp > p_3.itp")    # formated output the 3rd, 4th, 5th columns
        os.system("sed 's/1000/POSRES_FC_BB/g' p_3.itp > p_4.itp")  # s: substitute or replace 1000 using the variable POSRES_FC_BB
        os.system("paste -d '' p_2.itp p_4.itp > p_5.itp")  
        os.system("sed '1i [ position_restraints ]' p_5.itp > p_6.itp")   # insert the content in the first row
        os.system("sed '1i #ifdef POSRES' p_6.itp > p_7.itp")
        os.system("sed '$a #endif' p_7.itp > p_8.itp")   # add the content at the end of file
    
        os.system("sed '/; Include Position restraint file/,$d' old_topol.top > t_1.top") 
        os.system("sed '1,21d' t_1.top > t_2.top")                     
        os.system("cat t_2.top p_8.itp > Protein_chain_A.itp")    
        os.system("echo '#include \"Protein_chain_A.itp\"' >> temp_topol.top")
                  
        os.system("echo '  ' >> temp_topol.top")                  
        os.system("sed -n '/system/,$ p' old_topol.top >> temp_topol.top")   # append the number of molecules in temp_topol.top. p: print
        os.system("cp temp_topol.top topol.top")   
                  
        status, jobid = subprocess.getstatusoutput("sbatch run.sh")
        
    elif (len(chain_id)>1):    
        for i in range(0,len(chain_id)):
            cx = chain_id[i]
            os.system("sed '1,7d' posre_Protein_chain_"+cx+".itp > p_"+cx+"_1.itp")
            os.system("awk '{printf(\"%6s%6s\\n\",$1,$2)}' p_"+cx+"_1.itp > p_"+cx+"_2.itp")   # formated output the first two columns
            os.system("awk '{printf(\"%8s%8s%8s\\n\",$3,$4,$5)}' p_"+cx+"_1.itp > p_"+cx+"_3.itp")    # formated output the 3rd, 4th, 5th columns
            os.system("sed 's/1000/POSRES_FC_BB/g' p_"+cx+"_3.itp > p_"+cx+"_4.itp")  # s: substitute or replace 1000 using the variable POSRES_FC_BB
            os.system("paste -d '' p_"+cx+"_2.itp p_"+cx+"_4.itp > p_"+cx+"_5.itp")  
            os.system("sed '1i [ position_restraints ]' p_"+cx+"_5.itp > p_"+cx+"_6.itp")   # insert the content in the first row
            os.system("sed '1i #ifdef POSRES' p_"+cx+"_6.itp > p_"+cx+"_7.itp")
            os.system("sed '$a #endif' p_"+cx+"_7.itp > p_"+cx+"_8.itp")   # add the content at the end of file
            
            os.system("sed '/; Include Position restraint file/,$d' topol_Protein_chain_"+cx+".itp > t_"+cx+"_1.top")
            os.system("sed '1,18d' t_"+cx+"_1.top > t_"+cx+"_2.top")
            os.system("cat t_"+cx+"_2.top p_"+cx+"_8.itp > Protein_chain_"+cx+".itp")   # the final protein.itp file 
            os.system("echo '#include \"Protein_chain_"+cx+".itp\"' >> temp_topol.top")
            
        os.system("echo '  ' >> temp_topol.top")   # insert blank lines    
        os.system("sed -n '/system/,$ p' old_topol.top >> temp_topol.top")   # append the number of molecules in temp_topol.top. p: print
        os.system("cp temp_topol.top topol.top")
        
        status, jobid = subprocess.getstatusoutput("sbatch run.sh")
    jobid = jobid.split()[-1]
    return jobid, new_job

