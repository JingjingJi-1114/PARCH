#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os, sys
import numpy as np
import MDAnalysis as mda


user_account = "/home/jji110/"
new_job      = sys.argv[1]

pp_shell = "shell_2"

path  = user_account + "parch_platform/"+new_job+"/"+pp_shell+"/"   # path_dst: destination directory.    The directory "new_job" is being created when copying the content (files and directories) of backup_dir 

os.chdir(path)   

dd = 3.15   # the cutoff value for water calculations.

print(os.getcwd())

for iii in range(1,4):   # 3 jobs: mid_1, mid_2, mid_3   
    path_dst = path+"mid_"+str(iii)+"/"

    os.chdir(path_dst)      

    os.system("mkdir 3-15A_test")
    path_dst_sub = path_dst+"3-15A_test"+"/"
    
    ######## compute the initial H2O number attaching to each residue #########
    
    U0 = mda.Universe(path_dst+'em.tpr', path_dst+'em.gro')  # resid is 1-based.  dimer/multi-mer is fine. 
    
    aa = np.unique(U0.select_atoms("all").resids)[0]    # all: Protein+SOL
    bb = np.unique(U0.select_atoms("resname SOL").resids)[0] - 1 
    
    
    ww_num = []
    for k in range(aa,bb+1):    
        nn = str(k)    
        WI = U0.select_atoms('(around '+ str(dd) +' (resid '+ nn +')) and (resname SOL)')
        RI = np.unique(WI.resids)
        ww_num.append(len(RI))    
    NUM_I = np.array(ww_num)
    hp0 = len(np.argwhere(NUM_I==0))
    print ('hp0=',hp0)    
    np.savetxt(path_dst+'init_ww_3-15A.txt', NUM_I)    

    ################################## NUM_w.txt ##############################

    NUM = []
    U1 = mda.Universe(path_dst+'w_h.tpr', path_dst+"w_h.xtc")  # resid is 1-based.  Dimer is fine.
    for ts in U1.trajectory:    
        w_num = []    
        for k in range(aa,bb+1):        
            nn = str(k)        
            WW = U1.select_atoms('(around '+ str(dd) +' (resid '+ nn +')) and (resname SOL)')
            RR = np.unique(WW.resids)
            w_num.append(len(RR))
        NUM.append(np.array(w_num))        
    np.savetxt(path_dst_sub+'NUM_w.txt', NUM)    
    
    ###########################################################################
    ###### remove the fluctuations of H2O number attaching to each residue #####
        
    wb0 = np.loadtxt(path_dst_sub+'NUM_w.txt')   # len(wb0)=1001

    num_ww_tot = np.transpose(wb0)         # wb0 ~ row: number of time-frames;    num_ww_tot ~ row: number of residues

    temp_point_spec = []
    for i in range(0,len(num_ww_tot[0]),100):          # selecting the specified temperature data [0, 100, 200, 300,..., 1000]   
        temp_point_spec.append(num_ww_tot[:,i])        # temp_point ~ row: number of frames (11)
    ori_num_ww_temp = np.transpose(temp_point_spec)    # ori_num_ww_temp ~ row: number of residues  
    num_ww_temp     = np.transpose(temp_point_spec)    # row : residues,  column : specified temperature points
    
    for i in range(0,len(num_ww_temp)):    
        for j in range(1,len(num_ww_temp[0])):    # remove the fluctuations and keep the decreasing trend        
            cp_a = num_ww_temp[i][j-1]
            cp_b = num_ww_temp[i][j]
            if (cp_a<cp_b):            
                num_ww_temp[i][j] = cp_a
        #end of for(j)
    #end of for(i)    
    np.savetxt(path_dst+'num_ww_temp_3-15A.txt', num_ww_temp)
    #num_ww_temp = np.loadtxt(path+'num_ww_temp.txt')
  
    ############################## PP_ww.gro ##################################      

    ct = 0
    temp_points = 11   # the number of selected temperature points
    step_size = int((len(U1.trajectory)-1)/(temp_points-1))    # (1001-1)/(11-1)
    for ts in U1.trajectory:   # ts: iterate the number of frames    
        if (ct%step_size==0):        
            sol_shell_prot = U1.select_atoms('(resid '+str(aa)+':'+str(bb)+') or ((around '+ str(dd) +' (resid '+str(aa)+':'+str(bb)+')) and (resname SOL))')
            SP = sol_shell_prot        
            name = str(int(ct/step_size))        
            SP.residues.atoms.write(path_dst_sub+'pp_ww_'+name+'.gro', reindex=True)   #  reindex=True: bynum is 1-based (pp_ww_0.gro). reindex=False: bynum is 0-based (pp_ww_0.gro).             

            os.chdir(path_dst_sub)    # navigation

            os.system("gmx editconf -f pp_ww_"+name+".gro -o pp_ww_"+name+".pdb")             # pp_ww_0.pdb 
            os.system("echo q|gmx make_ndx -f pp_ww_"+name+".pdb -o pp_ww_"+name+".ndx")      # pp_ww_0.ndx
            os.system("echo 1|gmx editconf -f pp_ww_"+name+".pdb -o pp_for_bfactor_"+name+".pdb -n pp_ww_"+name+".ndx")    # pp_for_bfactor_0.pdb        
        ct = ct + 1  
    
    ###########################################################################
    ################### Coloring protein based on SOL number ##################

    evp_st_time   = 0   
    evp_ed_time   = 5000    # 5000 ps = 5 ns. Please check .mdp file
    evp_step_time = 500
    
    TT = []
    for i in range(evp_st_time,evp_ed_time+1,evp_step_time):    
        evp_fre = int(i/5)       # frame: 0, 100, 200, ..., 1000.   Temp: 300K, 350K, 400K, ..., 800K    
        TT.append(evp_fre)
            
    for t in range(0,len(TT)):   # TT=[0, 100, 200, 300, ..., 1000]
        TBF = wb0[TT[t]]
        ref_pdb = 'pp_for_bfactor_'+str(t)+'.pdb'    # refer to the specified frame, since protein's positions change a little during the annealing process.
        ref_ind = 'ww_factor_'
                
        pdb1 = path_dst_sub+ref_pdb   #!! reference    
        u = mda.Universe(pdb1)
        ref = mda.Universe(pdb1)  
        
        pdbtrj = path_dst_sub+ref_ind+str(t)+".pdb"
    
        num = []       
        with mda.Writer(pdbtrj, multiframe=True, bonds=None, n_atoms=u.residues.n_atoms) as PDB:  # multiframe=False:  single frame                
            ref.trajectory[0]          
            # iterate through our trajectory, if multifrane=True. (of course, multiframe=True works for single frame)
            # project displacement on structure via bfactor ("tempfactor") field       
            for jj in range(0,len(u.residues)):            
                for kk in range(0,len(u.residues.tempfactors[jj])):  # u.residues.tempfactors: is one m*n array. m: the # of residues. n: the # of atoms of each residue                
                    num.append(TBF[jj])
            u.atoms.tempfactors = num   # displacement for each atom
            
            PDB.write(u.atoms)                             

    ###########################################################################
    ######################### bfactor_ww + SOL ################################

    bynum_a = 1   # the 1st atom number of protein
    bynum_b = U1.select_atoms("resname SOL").indices[0]   # the last atom number of protein (index or indices is 0-based).
    
    sa = str(bynum_a+4)
    sb = str(bynum_b+4)
    
    na = str(bynum_b+6)
    nb = str(bynum_b+11)
    
    os.chdir(path_dst_sub)    
    
    #for t in range(0,len(TT)):
    for t in range(0,11):    
        name = str(t)
        os.system("sed ' "+sa+","+sb+" d' pp_ww_"+name+".pdb > pep_"+name+".pdb")    # pp_ww_0.pdb   
        os.system("cat ww_factor_"+name+".pdb pep_"+name+".pdb > media_ww_"+name+".pdb")    
        os.system("sed ' "+na+","+nb+" d' media_ww_"+name+".pdb > bfactor_ww_"+name+".pdb")
    
os.chdir(path)  


