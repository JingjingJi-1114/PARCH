#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 4 14:12:19 2022

@author: jji110
"""


import numpy as np


const_kb_si_units = 1.3806490e-23
const_h_si_units  = 6.62607015e-34
const_kb_by_h_si_units  = const_kb_si_units/const_h_si_units  
const_ps_to_sec = 1.0e-12


def get_monotonic_series(full_wt_vec):  ###### obtain a "monotonically" decreasing series
    ndata = np.size(full_wt_vec)   
    Lwt   = []
    Lwt.append(full_wt_vec[0])
    for i in range(1,ndata):
        if full_wt_vec[i] < full_wt_vec[i-1]:
            Lwt.append(full_wt_vec[i])
        #end of if
    #end of for
    mono_wt_vec   = np.asarray(Lwt)  
    return mono_wt_vec


def get_no_recross_series(full_wt_vec):    ###### allow a "platform" series
    ndata = np.size(full_wt_vec)    
    norecross_wt_vec = np.zeros_like(full_wt_vec)    
    norecross_wt_vec[0] = full_wt_vec[0]
    for i in range(1,ndata):
        norecross_wt_vec[i] = min(full_wt_vec[i],full_wt_vec[i-1])   
    #end of for
    return norecross_wt_vec


def get_hp_type(x):       ######  defining the strength of hydrophobicity:  1 (hydrophobic) ~ 10 (hydrophilic)
    assert x >= 0.0e0
    assert x <= 1.0e0
    x1 = 0.0e0
    hp_type = -9

    if (x==0.0e0):
        hp_type = 0
    else:
        for i in range(1,12):
            print('i=',i)
            x2 = x1 + 0.10e0
            if (x > x1 and x <= x2) :
                hp_type = i
                print('x1=',x1)
                print('x2=',x2)
                break        
            else:
                x1 = x2
        #end of for
    if hp_type > 10: hp_type = 10
    assert hp_type != -9
    return hp_type


def calc_time_corr_fn(skip_step,fvec):   ###### defining correlation function, return "one" value
    ndata   = np.size(fvec)              
    npt     = 0
    sum_fij = 0.0e0
    for i in range(ndata):
        j  = i + skip_step
        fj = 0.0e0
        if j < ndata: fj = fvec[j]    
        npt += 1    
        fij  = fvec[i]*fj           
        sum_fij += fij                      
        #print("npt fij",npt,fij)       
    #end of for
    if npt > 0:
        avg_fij = sum_fij/npt
    else:
        avg_fij = 0.0e0
    #end of if
    corr_fn = avg_fij
    return corr_fn         


def calc_time_corr_series(ndata,temp_vec,time_vec,wt_vec):   
    #--check that the data is consistent--#
    assert ndata == np.size(temp_vec)
    assert ndata == np.size(time_vec)
    assert ndata == np.size(wt_vec)
    
    ctau_vec = np.zeros(ndata)
   
    norecross_wt_vec = get_no_recross_series(wt_vec) 
    for itau in range(ndata):
        ctau_vec[itau] = calc_time_corr_fn(itau,norecross_wt_vec)                   
    return ctau_vec      

    
def calc_time_integral_corr_fn (ndata,temp_vec,time_vec,ctau_vec,norm_factor):
    intval   = 0.0e0
    x1       = time_vec[0]
    fx1      = ctau_vec[0] 
    sum_area = 0.0e0
    for i in range(1,ndata):
        x2   = time_vec[i]
        fx2  = ctau_vec[i]
        area = (x2-x1)*(fx1+fx2)*0.50e0
        sum_area = sum_area + area
        x1       = x2
        fx1      = fx2
    #end of for
    intval = sum_area/norm_factor
    return intval


def reader_wtmat_all_res(infile):    
    Ldata   = []
    fin = open(infile,"r")
    for line in fin:        
        buff  = line.strip().split()  
        v = []                          
        for x in buff:
            v.append(float(x))    
        #end of for
        Ldata.append(v)
    fin.close()
    mat_wt_series = np.asarray(Ldata)
    num_res,ndata   = np.shape(mat_wt_series)   
    print("num_res_from_file= ",num_res)
    print("num_data_from_file= ",ndata)
    return num_res,ndata,mat_wt_series


def writer_ctau_mat_all_res(outfile,num_res,ndata,ctau_mat):  
    fout = open(outfile,"w")
    for ires in range(num_res):
        v = ctau_mat[ires,:]    
        s = ''
        for x in v:
            s = s + str(x) + "  "    
        fout.write(s+" \n")
    fout.close()
    return None


def find_first_zero_idx(vec):
    n = np.size(vec)
    ians = 0
    for i in range(n):
        if abs(vec[i]) <= 1.0e-100:
            ians = i
            break
    return ians


def main_program():
    #-------user_input----------#
    infile  = "num_ww_temp_3-15A.txt"
    outfile = "ctau_all_res_3-15A.txt"
    
    print("---reading file--- ",infile)
    
    num_res,ndata,mat_wt_series = reader_wtmat_all_res(infile)
    time_vec = np.array([0, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000])
    temp_vec = np.array([300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800])
    assert ndata == np.size(time_vec)
    assert ndata == np.size(temp_vec)

    norm_factor = 1.0e0
    hard_ref    = 743340.909090909

    ctau_mat    = np.zeros((num_res,ndata))
    ctau_intval = np.zeros(num_res)
    max_ctau_intval = -1
    for ires in range(num_res):
        wt_vec = mat_wt_series[ires,:]
        ctau_mat[ires,:]  = calc_time_corr_series(ndata,temp_vec,time_vec,wt_vec)   
        ctau_intval[ires] = calc_time_integral_corr_fn(ndata,temp_vec,time_vec,ctau_mat[ires,:],norm_factor)   
        max_ctau_intval   = max(max_ctau_intval,ctau_intval[ires])    
    writer_ctau_mat_all_res(outfile,num_res,ndata,ctau_mat)
    #-------print out ----------#
    fout = open("hp_type_for_color_LYS_ref_3-15A.txt","w")
    for ires in range(num_res):
        ctau_intval[ires] = (ctau_intval[ires]/hard_ref)*10.0e0
        hp_type = round(ctau_intval[ires],2)     
        print(ctau_intval[ires],ires+1,hp_type)
        fout.write(str(hp_type)+"\n")        
    fout.close()
    print('max_ctau_intval = ',max_ctau_intval)
    return None


######################################### test ##############################################

def my_fake_data_test1():
    ndata    = 11
    time_vec = np.array([0, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000])
    temp_vec = np.array([300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800])
    wt_vec   = np.array([9,    9,   8,   8,   8,   7,   5,   4,   4,   3,   3])

    calc_time_corr_series(ndata,temp_vec,time_vec,wt_vec)

    return None


if __name__ == "__main__" :
    print("--starting driver--")
    #my_fake_data_test1()
    main_program()

