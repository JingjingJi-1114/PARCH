#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 09:10:00 2023

@author: jji110
"""


import numpy as np
import pandas as pd
import csv
import ast
from itertools import permutations, combinations
import json
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns


class Triad():     
    def __init__(self, pdbId, chainIds, residIds, residues):
        self.pdbId = pdbId
        self.chainIds = chainIds
        self.residIds = residIds
        self.residues = residues 
        
    def __eq__(self, triad2):    # def equal(self,triad2):
        if self.pdbId != triad2.pdbId:
            return False        
        unique_residues_1 = list(zip(self.chainIds, self.residues, self.residIds))  
        unique_residues_2 = list(zip(triad2.chainIds, triad2.residues, triad2.residIds))         
        for i in unique_residues_1:
            if i not in unique_residues_2:
                return False            
        return True 
    
    def __str__(self):   # str(triad) by default.
        return str(self.pdbId)+"-"+str(list(zip(self.chainIds, self.residIds, self.residues)))
        
        
class Triads():   # all triads    
    def __init__(self):
        self.explored_triads = {}  
    
    def add_triad(self, newTriad):
        if newTriad.pdbId in self.explored_triads.keys():
            self.explored_triads[newTriad.pdbId] = np.append(self.explored_triads[newTriad.pdbId], newTriad)
        else:
            self.explored_triads[newTriad.pdbId] = np.array([newTriad], dtype=object)            
            
    def exists(self, newTriad):
        if newTriad.pdbId not in self.explored_triads.keys():
            return False
        else:
            for triad in self.explored_triads[newTriad.pdbId]:
                if newTriad == triad:     # newTriad.equal(triad). same as newTriad.__eq__(triad) => 
                    return True
            return False
        
    def find_all_triads(self, residues):
        all_triads = []
        for pdb, triads in self.explored_triads.items():
            for triad in triads:
                if sorted(triad.residues) == sorted(residues):
                    all_triads.append(triad)
        return all_triads
        

def pair_final():
    pair = []
    for i in range(0,20):
        sub_pair = []
        for j in range(0,20):
            ss = ref_LL[i]+'_'+ref_LL[j]
            sub_pair.append(ss)
        pair.append(sub_pair)
    pair_arr = np.array(pair)
    return pair_arr


def get_neighbors(sort_nbr_info):    
    neighbors = []
    ori_resid = []
    chain_id  = []    
    for ele in sort_nbr_info:
        ss = ele.split('_')[2]  # 3-letter
        neighbors.append(ref_LL[np.argwhere(ref_AA==ss)[0][0]])  
        ori_resid.append(int(ele.split('_')[1]))
        chain_id.append(ele.split('_')[0]) 
        
    neighbors_arr = np.array(neighbors)
    ori_resid_arr = np.array(ori_resid)
    chain_id_arr  = np.array(chain_id)
    
    neighbors_new = np.sort(neighbors_arr)  # single-letter
    ori_resid_new = ori_resid_arr[np.argsort(neighbors_arr)]
    chain_id_new  = chain_id_arr[np.argsort(neighbors_arr)]        
    return list(neighbors_new), list(ori_resid_new), list(chain_id_new)


def remove_nan(matrix):
    return matrix[~np.isnan(matrix)]


def make_symmetric(matrix):
    return matrix + matrix.T - np.diag(matrix.diagonal())


def get_permutations(neighbors, nbr_pv):
    return list(permutations(neighbors, 2)), list(permutations(nbr_pv, 2))


def get_combinations(neighbors, nbr_pv, ori_resid, chain_id):
    return list(combinations(neighbors, 2)), list(combinations(nbr_pv, 2)), list(combinations(ori_resid, 2)), list(combinations(chain_id, 2))


def find_nbr_pv(pdb_id, chain_id, ori_resid, neighbors, df_arr):
    pvs = []
    for i in range(len(neighbors)):
        pv = df_arr[np.argwhere((df_arr[:,2] == pdb_id) & (df_arr[:,3] == chain_id[i]) & (df_arr[:,4] == ori_resid[i]) & (df_arr[:,1] == neighbors[i])), 6]
        if len(pv) != 0:
            pvs.append(pv[0][0])        
    return pvs


def save_data_bank(data_bank, filename):
    with open(filename, 'w') as fout:
        csv_data = csv.writer(fout, dialect='excel')
        csv_data.writerows(data_bank)   
        fout.close()


def unique_counts(all_cmb):
    unique_counts = {}
    for pdbId, groups in all_cmb.explored_triads.items():
        for gg in groups:
            cmb_name = ''.join(sorted(gg.residues))
            if cmb_name in unique_counts.keys():
                unique_counts[cmb_name] = unique_counts[cmb_name] + 1
            else:
                unique_counts[cmb_name] = 1
    return unique_counts


path = "/home/jji110/parch_platform/Analysis_script/pp_1000_t2/AA_nbr/new_triad_file/"
ref_LL = np.array(['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'])
ref_AA = np.array(['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR'])
fname = 'LYS'  # input: which amino acid
fn = fname+'_surf'
clt = 2   # label

df0 = pd.read_csv(path+'v2_1000_pp_cluster.csv', index_col=False, header=0)
df_arr = df0.to_numpy()


def main_program(df_arr):
    triple_nbr = np.zeros((20, 20), dtype=np.int32)
    triple_nbr_pv = np.zeros((20, 20), dtype=float)
    triple_nbr_list_pv = np.empty((20,20),dtype=object)
    allTriads = Triads()
    for row in df_arr:
        if (row[5]==fname and row[10]==clt):
            neighbors, ori_resid, chain_id = get_neighbors(ast.literal_eval(row[-1]))  
            pdb_id = row[2]
            nbr_pv = find_nbr_pv(pdb_id, chain_id, ori_resid, neighbors, df_arr)
            nbr_perms, nbr_pv_perms, resid_perms, chain_perms = get_combinations(neighbors, nbr_pv, ori_resid, chain_id)    
            common_resid_pv = row[6]

            assert(len(nbr_pv_perms) == len(nbr_perms))
            
            for idx in range(0, len(nbr_perms)):
                triad = Triad(pdb_id, [row[3]] + list(chain_perms[idx]), [row[4]] + list(resid_perms[idx]), [row[1]] + list(nbr_perms[idx]))
                if not allTriads.exists(triad):                      
                    row_idx = np.where(ref_LL == nbr_perms[idx][0])[0][0]
                    col_idx = np.where(ref_LL == nbr_perms[idx][1])[0][0]                    
                    triple_nbr[row_idx, col_idx] += 1
                    ave_pv_triple = (nbr_pv_perms[idx][0]+nbr_pv_perms[idx][1]+common_resid_pv)/3.0
                    triple_nbr_pv[row_idx, col_idx] += ave_pv_triple
                    triple_nbr_list_pv[row_idx, col_idx] = np.append(triple_nbr_list_pv[row_idx, col_idx], ave_pv_triple) 
                    allTriads.add_triad(triad)                                                                    
    return triple_nbr, triple_nbr_pv, triple_nbr_list_pv, allTriads

triple_nbr_final, triple_nbr_pv_final, triple_nbr_list_pv_final, allTriads = main_program(df_arr)

################################## data output ##############################
#############################################################################

unique_counts = unique_counts(allTriads)
df = pd.DataFrame(unique_counts, index=[0]).T
df.to_excel(path+'triad_all_cmb.xlsx')

for i in range(0,20):
    for j in range(0,20):
        if triple_nbr_list_pv_final[i,j] is not None:
            triple_nbr_list_pv_final[i,j] = triple_nbr_list_pv_final[i,j][np.where(triple_nbr_list_pv_final[i,j] != None)]
        else:
            triple_nbr_list_pv_final[i,j] = triple_nbr_list_pv_final[j,i]
            
data_list = []
for i in range(0,20):
    for j in range(0,20):
        data_list.append(list(triple_nbr_list_pv_final[i,j]))
        
with open(path+fn+"_pv_list.json", "w") as json_file:
    json.dump(data_list, json_file)
                      
with open(path+fn+"_triads.txt", "w") as file:
    for pdbId, triads in allTriads.explored_triads.items():
        for triad in triads:
            file.write(str(triad)+"\n")

triple_nbr_symm    = make_symmetric(triple_nbr_final)
triple_nbr_pv_symm = make_symmetric(triple_nbr_pv_final)

df = pd.DataFrame(triple_nbr_symm)
df.to_excel(path+"triad_nbr_"+fn+".xlsx", index=False)

df_pv = pd.DataFrame(triple_nbr_pv_symm)
df_pv.to_excel(path+"triad_nbr_pv_"+fn+".xlsx", index=False)

################################### plot ######################################
###############################################################################

data_arr = triple_nbr_symm
plot_arr = np.diag(np.diag(data_arr))+np.tril(data_arr,k=-1)
#plot_arr = triple_nbr_symm

plt.style.use('default')
font_0 = {'size' : 26, 'weight' : 'bold'}  
font   = {'size' : 28, 'weight' : 'bold'}  
matplotlib.rc('font', **font_0)
matplotlib.rc('axes', linewidth=1.5)

matplotlib.rcParams['font.family'] = "sans-serif"   
matplotlib.rcParams['font.sans-serif'] = "Arial"

fig = plt.figure(figsize=(10,10))
ax  = fig.add_subplot(111)

sns.heatmap(plot_arr, annot=False, fmt='.0f', mask=(plot_arr==0), vmin=0, vmax=1500, cmap='PiYG_r', cbar=False,
            cbar_kws={"shrink": 1, "orientation": "vertical", "pad": 0.05, "format": "%.0f", "ticks": [0, 375, 750, 1125, 1500]})

x_label = ref_LL
y_label = ref_LL

ax.set_xticklabels(x_label, fontdict=font_0, rotation=0)
ax.set_yticklabels(y_label, fontdict=font_0, rotation=0)

ax.set_ylabel('Residue', fontdict=font)
ax.set_xlabel('Residue', fontdict=font)

#for spine in ax.spines.values():
#    spine.set(visible=True, lw=1.0, edgecolor="black")

#ax.spines['top'].set_visible(False)

plt.savefig(path+'pic_nbr_'+fn+'.png', transparent=True, bbox_inches='tight', dpi=600)
#plt.show()
