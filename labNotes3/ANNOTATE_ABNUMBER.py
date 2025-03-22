#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 20 11:32:31 2025

@author: toma
"""

import numpy as np
import pandas as pd
from abnumber import Chain
import glob
from pandarallel import pandarallel
import time
import os


dfs=pd.read_csv('FILE_WITH_SEQUENCES.csv')

def ABN(row):    
    GF_DATA=[]
    hc=row['[HC]']
    lc=row['[LC]']
    try:
            varH = Chain(hc, scheme='imgt', allowed_species='human', assign_germline=True)
            varL = Chain(lc, scheme='imgt', allowed_species='human', assign_germline=True)
            GF_DATA.append([varH.seq,varH.fr1_seq,varH.cdr1_seq,varH.fr2_seq,varH.cdr2_seq,varH.fr3_seq,varH.cdr3_seq,varH.fr4_seq,varH.v_gene,varH.j_gene,varL.seq,varL.fr1_seq ,varL.cdr1_seq,varL.fr2_seq ,varL.cdr2_seq,varL.fr3_seq,varL.cdr3_seq,varL.fr4_seq,varL.v_gene,varL.j_gene])
    except:
              print("Gen pair not recognized ")
              print(hc)
              print(lc)
              
              
    return GF_DATA           



res_par=dfs.parallel_apply(ABN,axis=1)

DATA_N=[]
for ik in range(len(res_par)):
    ti=res_par[ik]
    if len(ti)>0:
       DATA_N.append(ti[0])
print(f'Recognized  HC/LC are {len(DATA_N)}')
c_names=['HC_V','hfr1','hcdr1','hfr2','hcdr2','hfr3','hcdr3','hfr4','heavy_V','heavy_J','LC_V','lfr1','lcdr1','lfr2','lcdr2','lfr3','lcdr3','lfr4','light_V','light_J',]
dfF=pd.DataFrame(DATA_N,columns=c_names)
dfF.to_csv('ANNOTATED_FILE.csv',index=False)