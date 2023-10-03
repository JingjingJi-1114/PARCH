#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os, sys
import numpy as np
import MDAnalysis as mda
import shutil, json


user_account = "/home/jji110/"
new_job      = sys.argv[1]

pp_shell = "shell_2"

path_src  = user_account + "parch_platform/backup_analysis/" 
path      = user_account + "parch_platform/"+new_job+"/"+pp_shell+"/"

###############################################################################
print(path)
pv = []   
 
for iii in range(1,4):   # 3 jobs: mid_1, mid_2, mid_3   

    path_dst = path+"mid_"+str(iii)+"/"
    path_dst_sub = path_dst+"3-15A_test"+"/"

    os.chdir(path_dst)     

    ############### compute the parch value (pv) ################

    shutil.copy(path_src+"4_pv_correlation.py", path_dst)    
    os.system("./4_pv_correlation.py")

    #################### write pv into .pdb #####################

    TBF = np.loadtxt(path_dst+'hp_type_for_color_LYS_ref_3-15A.txt')   
    pv.append(TBF)  
    
    ref_pdb = 'pp_for_bfactor_0.pdb'    # refer to the specified frame, since protein's positions change a little during the annealing process.
    ref_ind = 'correlation_pv.pdb'
                
    pdb1 = path_dst_sub+ref_pdb   #!! reference    
    u    = mda.Universe(pdb1)
    ref  = mda.Universe(pdb1)  
        
    pdbtrj = path_dst_sub+ref_ind
    
    num = []    
    with mda.Writer(pdbtrj, multiframe=True, bonds=None, n_atoms=u.residues.n_atoms) as PDB:  # multiframe=False:  single frame            
        ref.trajectory[0]   # reference coordinates: set to first frame
        for jj in range(0,len(u.residues)):        
            for kk in range(0,len(u.residues.tempfactors[jj])):  # u.residues.tempfactors: is one m*n array. m: the # of residues. n: the # of atoms of each residue            
                num.append(TBF[jj])
        u.atoms.tempfactors = num   # displacement for each atom    
        PDB.write(u.atoms)                             

######################### compute the average_pv ##############################

os.chdir(path+"mid_1/")
 
ave_pv = np.average(np.array(pv), axis=0)     # axis=0, along the rows; axis=1, along the columns
std_pv = np.std(np.array(pv), axis=0)

#print ('ave_pv = ', ave_pv)
#print ('std_pv = ', std_pv)

np.savetxt(path+'mid_1/ave_pv_3-15A.txt', ave_pv)
np.savetxt(path+'mid_1/std_pv_3-15A.txt', std_pv)

######################### write average_pv into .pdb ##########################      
        
ave_TBF = np.loadtxt(path+'mid_1/ave_pv_3-15A.txt')   

ref_pdb = 'pp_for_bfactor_0.pdb'    # refer to the specified frame, since protein's positions change a little during the annealing process.
ref_ind = 'ave_correlation_pv.pdb'
            
pdb1 = path+'mid_1/3-15A_test/'+ref_pdb   #!! reference    
u    = mda.Universe(pdb1)
ref  = mda.Universe(pdb1)  
    
pdbtrj = path+'mid_1/3-15A_test/'+ref_ind
num = []    
with mda.Writer(pdbtrj, multiframe=True, bonds=None, n_atoms=u.residues.n_atoms) as PDB:             
    ref.trajectory[0]   
    for jj in range(0,len(u.residues)):        
        for kk in range(0,len(u.residues.tempfactors[jj])):              
            num.append(ave_TBF[jj])
    u.atoms.tempfactors = num       
    PDB.write(u.atoms)                             

################################# Correcting pdb ##############################
    
def correct_pvfile(resid, pdb):
    replacements = []
    for key, value in resid.items():
        for rid in value:
            replacements.append(key + rid.rjust(4, " "))   # .rjust: right adjusted
            print(rid)

    replacements_ptr = 0
    prev_rid = None

    for i in range(0,len(pdb)):         
        if pdb[i][:4] == "ATOM":
            row = pdb[i].split()
            
            if (len(row[4])==1):                        
                now_rid = int(row[5])   
            else:                            
                now_rid = int(row[4][1:])     
                
            if (prev_rid==None):        
                prev_rid = now_rid
                pdb[i] = pdb[i][:21] + replacements[replacements_ptr] + pdb[i][26:]
            else:
                if (prev_rid == now_rid):    
                    pdb[i] = pdb[i][:21] + replacements[replacements_ptr] + pdb[i][26:]
                else:                       
                    prev_rid = now_rid       
                    replacements_ptr += 1
                    pdb[i] = pdb[i][:21] + replacements[replacements_ptr] + pdb[i][26:]
                    if (replacements_ptr != 0 and replacements[replacements_ptr][0] != replacements[replacements_ptr-1][0]):
                        pdb[i-1] = pdb[i-1] + 'TER\n'   # pdb is a list now.  "sss\n"+"TER\n" = "sss\nTER\n", does not change the length of pdb list. But, when pdb is output to a file, "TER\n" will be in next line.    
    return pdb


os.chdir(path+"mid_1/3-15A_test/")    
os.system("cp ave_correlation_pv.pdb ave_correlation_pv.txt")    


fin_1 = open(user_account+"parch_platform/"+new_job+"/resid_dict.txt","r")
data = fin_1.read()
fin_1.close()
resid = json.loads(data)

fin_2 = open("ave_correlation_pv.txt","r")
pdb = fin_2.readlines()
fin_2.close()
    
corrected_pdb = correct_pvfile(resid, pdb)    
    
fout = open('ave_correlation_pv_correct.pdb', 'w')
for line in corrected_pdb:
    fout.write(line)
fout.close()

os.chdir(user_account+"parch_platform/")        

