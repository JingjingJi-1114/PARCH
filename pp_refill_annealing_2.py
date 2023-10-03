#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import matplotlib.pyplot as plt
import matplotlib.pylab as py
from mpl_toolkits.mplot3d import Axes3D     
import sys, subprocess
import numpy as np
import math
import os
import shutil


import MDAnalysis as mda
import MDAnalysis.analysis.distances as Mdistance


def annealing(new_job, jobid):
    user_account = "/home/jji110/"
    #new_job      = "test_1ycr"
    
    pp_shell = "shell_2"
    
    path_src  = user_account + "parch_platform/backup_annealing/"    # path_src: source directory
    path_dst  = user_account + "parch_platform/"+new_job+"/"   # path_dst: destination directory.    The directory "new_job" is being created when copying the content (files and directories) of backup_dir 
    
    
    os.chdir(path_dst)   
     
    
    U_11 = mda.Universe(path_dst+"final_solv_ions.gro")   # check if the protein is charged or not
    num_na = len(np.unique(U_11.select_atoms("resname NA").resids))
    num_cl = len(np.unique(U_11.select_atoms("resname CL").resids))
    pp_net_charge = num_cl - num_na
    print('pp_net_charge = ', pp_net_charge)
    np.savetxt(path_dst+'pp_net_charge.txt',np.array([[pp_net_charge]]))
    
    os.system("echo 19 19|gmx trjconv -f md.xtc -s md.tpr -b 10 -e 10 -tu ns -o pp_10ns.gro -n index.ndx -pbc mol -center")
    
    
    os.system("mkdir "+pp_shell)
    os.chdir(path_dst+pp_shell)    
    
    os.system("cp ../pp_10ns.gro .")
    os.system("cp ../Protein_chain_* .")
    os.system("cp -r ../charmm36-jul2021.ff .")
    
    os.system("cp ../topol.top .")
    os.system("mv topol.top old_topol.top")
    os.system("sed '/NA/ d' old_topol.top > old_topol_1.top") 
    os.system("sed '/CL/ d' old_topol_1.top > old_topol_2.top")
    os.system("sed '/SOL/ d' old_topol_2.top > old_topol_3.top")
    os.system("cp old_topol_3.top topol.top")
    os.system("echo 'SOL                 1' >> topol.top")  
    
    os.system("cp old_topol_3.top old_template.top")
    os.system("sed '/ions.itp/ d' old_template.top > old_template_1.top")
    os.system("sed '/system/i #include \"SOD.itp\"' old_template_1.top > old_template_2.top")  # insert (#include "SOD.itp") above [ system ]
    os.system("sed '/system/i #include \"CLA.itp\"' old_template_2.top > old_template_3.top")
    os.system("sed '/system/{x;p;x;}' old_template_3.top > template.top")  # insert one blank line above [ system ]
              
    shutil.copy(path_src+"SOD.itp", path_dst+pp_shell)     # file => directory
    shutil.copy(path_src+"CLA.itp", path_dst+pp_shell)     # CL ions are fixed
    shutil.copy(path_src+"1_SOL.gro", path_dst+pp_shell)
    
    ### gmx insert-molecules : re-arrange the resids of several proteins following the order
    os.system("gmx insert-molecules -f pp_10ns.gro -ci 1_SOL.gro -o PP_1_SOL.gro -nmol 1")  # gmx insert can help "re-arrange resids" of "several proteins" following the order
    os.system("gmx solvate -cp PP_1_SOL.gro -cs spc216.gro -o box_solv.gro -p topol.top")   # refill

    #######################  refilled protein-water layer ##############################
    
    dd = 4.15     # angstrom, 2nd water shell thickness
    
    input_file = 'box_solv.gro'   # Protein+SOL
    U_2 = mda.Universe(path_dst+pp_shell+"/"+input_file)
    
    aa = np.unique(U_2.select_atoms("all").resids)[0]    # all: Protein+SOL
    bb = np.unique(U_2.select_atoms("resname SOL").resids)[0] - 1   
    
    sol_shell = U_2.select_atoms('(around '+ str(dd) +' (resid '+str(aa)+':'+str(bb)+')) and (resname SOL)')
    sol_shell_prot = U_2.select_atoms('(resid '+str(aa)+':'+str(bb)+') or ((around '+ str(dd) +' (resid '+str(aa)+':'+str(bb)+')) and (resname SOL))')
    
    SP = sol_shell_prot
    SP.residues.atoms.write(path_dst+pp_shell+'/W_shell.gro', reindex=False)    #!!  https://www.mdanalysis.org/docs/documentation_pages/coordinates/GRO.html 
    
    SS = sol_shell
    os.system("echo 'SOL         "+str(len(SS.residues))+"' >> template.top")    # >> : append
    os.system("cp template.top W_topol.top")
        
    #############################################################################################
    ################ adding ions (or not) to prepare the systems for annealing ##################
    #############################################################################################
    
    ############# (define functions) ##############
    
    def random_vector(shell_radius):    #  the function that is used for generating spherical surface points randomly
        rand_i, rand_j = np.random.rand(2)                # np.random.rand(2) = array([0.35709415, 0.85351321]).    np.random.rand(3,2) = gererate a (3 rows * 2 columns) array, the elements inside are from a uniform distribution in the range (0, 1).   
        theta          = np.pi * rand_i                   #  theta = [0,180]
        phi            = (2.0 * np.pi) * rand_j           #  phi = [0,360]
        x              = shell_radius * np.sin(theta) * np.cos(phi)        # Cartesian coordinate x
        y              = shell_radius * np.sin(theta) * np.sin(phi)        # Cartesian coordinate y
        z              = shell_radius * np.cos(theta)                      # Cartesian coordinate z
        vector         = np.array([x, y, z])
        return vector
        
    def cal_distance(p1, p2):   #  the function that is used for computing the distane between any two points.
        dd_ion = math.sqrt(math.pow((p2[0]-p1[0]),2) + math.pow((p2[1]-p1[1]),2) + math.pow((p2[2]-p1[2]),2))    #  math.pow( x, y ) : x^y
        return dd_ion
    
    ############### (preparation) ################
        
    path = path_dst+pp_shell+"/"
    
    os.chdir(path)   # navigation
    
    box_buffer  = 3.0    # nm, to avoid the mirror issue
    dis_ion_ps  = 3.0    # nm, the distance of ion from the protein surface
    dis_ion_ion = 3.0    # nm, the distance between any two points
    
    os.system("gmx editconf -f W_shell.gro -center 0 0 0 -o W_init_c.gro > inform.txt")    # the geometrical center of protein is at (0,0,0)
    
    U0 = mda.Universe(path+'W_init_c.gro')
    PP_COG = U0.select_atoms('resid '+str(aa)+':'+str(bb)).center_of_geometry()
    dist = []
    for i in range(aa,bb+1):
        name = str(i)
        res  = U0.select_atoms('resid '+name).center_of_geometry()
        DD   = Mdistance.distance_array(PP_COG,res)   # unit: angstrom
        dist.append(DD[0][0])    # dist: angstrom 
    dist_max = np.max(np.array(dist))/10.0   # dist_max : nm
    
    print(os.system("pwd"))   #  pwd : the work directory
    
    #############################################
    ####### (check if protein is charged) #######
    
    if (pp_net_charge!=0):    # protein has net charge(s)
    
        shutil.copyfile(path_src+"INP_ion", path+"INP_ion")
        shutil.copyfile(path_src+"hydrated_na.gro", path+"hydrated_na.gro")
        shutil.copyfile(path_src+"hydrated_cl.gro", path+"hydrated_cl.gro")
    
        if (pp_net_charge<0):
            ion_name = "hydrated_na.gro"
            ion_resname = "NA"
        elif (pp_net_charge>0):
            ion_name = "hydrated_cl.gro"
            ion_resname = "CL"
    
        shell_radius =  dist_max + dis_ion_ps   #   # unit: nm
        
        position = []
        for i in range(0,abs(pp_net_charge)*10000):   #  to avoid that two points are close
            point = random_vector(shell_radius)   # the function that is used for generating spherical surface points randomly
            count = 0
            if(i==0):        
                position.append(point)   # generate the 1st point
            elif(i>0):
                for j in range(0,len(position)):    #  compare with the points generatd
                    temp_dis = cal_distance(position[j],point)   #  the function that is used for computing the distane between any two points.
                    if (temp_dis<dis_ion_ion):
                        break
                    elif (temp_dis>=dis_ion_ion):
                        count = count + 1     # make sure that the distances between the current point and previous points are all >= 5.0 nm
                #end of for
                if (count==len(position)):
                    position.append(point)
            #end of if
            if (len(position)==abs(pp_net_charge)):
                break
        #end of for
        
        np.savetxt(path+'position.txt', np.array(position))
        os.system("cp position.txt position.dat")
        os.system("gmx insert-molecules -f W_init_c.gro  -ci "+ion_name+" -o pp_ion.gro -ip position.dat")
        
        print('shell_radius(nm) = ', shell_radius)
        np.savetxt(path+'ion_shell_radius.txt',np.array([[shell_radius]]))
        
        ################## Some edits about files ###################
    
        os.system("cp pp_ion.gro ori_pp_ion.gro")
        
        U1 = mda.Universe(path+'pp_ion.gro')
        
        pp_sol = U1.select_atoms("(resid "+str(aa)+":"+str(bb)+") or (resname SOL)")
        pp_sol.residues.atoms.write(path+'PP_SOL.gro', reindex=True)
        
        ions = U1.select_atoms("resname "+ion_resname)
        ions.residues.atoms.write(path+'only_'+ion_resname+'.gro', reindex=True)
        
        na = str(len(pp_sol.indices)+3)   # the number of atoms of PP_SOL + additional 3 lines
        nb = str(len(pp_sol.indices)+5)   # delete the middle 3 lines
        
        sa = str(len(pp_sol.indices))
        sb = str(len(pp_sol.indices)+abs(pp_net_charge))
        
        os.system("cat PP_SOL.gro only_"+ion_resname+".gro > tot_1.gro")
        os.system("sed ' "+na+","+nb+" d' tot_1.gro > tot_2.gro")    # delete the middle 3 lines
        os.system("sed '2s/"+sa+"/"+sb+"/' tot_2.gro > tot_3.gro")   # replace the old number of atoms in the 2nd line
        
        box_length = str((shell_radius+box_buffer)*2.0)   # reset the size of box
        os.system("gmx make_ndx -f tot_3.gro -o W_ind.ndx < INP_ion")   # update the W_ind.ndx
        os.system("echo 19 19|gmx editconf -f tot_3.gro -o W_init.gro -box "+box_length+" "+box_length+" "+box_length+" -n W_ind.ndx -c")   # 19: PP + Solv (SOL+ions)
        
        
        sol = U1.select_atoms("resname SOL")
        os.system("cp W_topol.top old_W_topol.top")     # update the W_topol.top
        os.system("sed -i '/SOL/ d' W_topol.top")
        os.system("echo 'SOL              "+str(len(sol.residues))+"' >> W_topol.top")
        os.system("echo '"+ion_resname+"                "+str(abs(pp_net_charge))+"' >> W_topol.top")
        
        ##################### files copy&paste #####################
    
        os.chdir(path)     # navigation
        os.system("mkdir mid_1")
        
        os.system("cp -r charmm36-jul2021.ff mid_1")
        os.system("cp Protein_chain_* mid_1")
        os.system("cp SOD.itp mid_1")
        os.system("cp CLA.itp mid_1")
        os.system("cp W_topol.top mid_1")
        os.system("cp W_init.gro mid_1")
        os.system("cp W_ind.ndx mid_1")
        shutil.copy(path_src+"em.mdp", path+"mid_1")   # file ==> directory
        
        os.system("cp -r mid_1 mid_2")
        os.system("cp -r mid_1 mid_3")
        
    elif (pp_net_charge==0):  # protein is neutral
        
        shutil.copyfile(path_src+"INP_no_ion", path+"INP_no_ion")
        
        shell_radius =  dist_max + dis_ion_ps  # though no ions are inserted into box.
        box_length = str((shell_radius+box_buffer)*2.0)
        
        os.system("gmx make_ndx -f W_init_c.gro -o W_ind.ndx < INP_no_ion")   # update the W_ind.ndx
        os.system("echo 17 17|gmx editconf -f W_init_c.gro -o W_init.gro -box "+box_length+" "+box_length+" "+box_length+" -n W_ind.ndx -c")   # 17: PP + Solv (SOL)
    
        ##################### files copy&paste #####################
    
        os.chdir(path)     # navigation
        os.system("mkdir mid_1")
        
        os.system("cp -r charmm36-jul2021.ff mid_1")
        os.system("cp Protein_chain_* mid_1")
        os.system("cp SOD.itp mid_1")
        os.system("cp CLA.itp mid_1")
        os.system("cp W_topol.top mid_1")
        os.system("cp W_init.gro mid_1")
        os.system("cp W_ind.ndx mid_1")
        shutil.copy(path_src+"em.mdp", path+"mid_1")   # file ==> directory
        
        os.system("cp -r mid_1 mid_2")
        os.system("cp -r mid_1 mid_3")    
    #end of if
           
    for i in range(1,4):   # generate mid_1, mid_2, mid_3 
        shutil.copy(path_src+"heat_nvt_"+str(i)+".mdp", path+"mid_"+str(i))
        shutil.copy(path_src+"evp_"+str(i)+".sh", path+"mid_"+str(i))
    
    ann_jobs = []
    for i in range(1,4):
        os.chdir(path+"mid_"+str(i))
        cmd = "sbatch --depend=afterok:{0} evp_{1}.sh".format(jobid, i)
        status, jid = subprocess.getstatusoutput(cmd)
        ann_jobs.append(jid.split()[-1])
    
    print(os.system("pwd"))    
    os.chdir(path)    
    return ann_jobs, new_job

def start_water_analysis(jobids, new_job):
    j1, j2, j3 = jobids
    print("********************Starting Water Analysis*******************************")
    cmd = "sbatch --depend=afterok:{0}:{1}:{2} wt_analysis.sh {3}".format(j1, j2, j3, new_job)
    print(cmd)
    status, wa_jobid = subprocess.getstatusoutput(cmd)
    print(wa_jobid)
    print("************************Submitted Water Analysis Job*********************************")
    return wa_jobid.split()[-1]
    
def calculate_ave_pv(jobid, new_job):
    print("*************************Calculating PARCH Values*******************************")
    cmd = "sbatch --depend=afterok:{0} parch_values.sh {1}".format(jobid, new_job)
    print(cmd)
    status, pv_jobid = subprocess.getstatusoutput(cmd)
    print(pv_jobid)
    return pv_jobid.split()[-1]


ann_jobs, new_job = annealing(sys.argv[1], sys.argv[2])

os.chdir("/home/jji110/parch_platform/")

wa_jobid = start_water_analysis(ann_jobs, new_job)

pv_jobid = calculate_ave_pv(wa_jobid, new_job)

