# -*- coding: utf-8 -*-
"""
Chi square gender test in all targets from HTS screen
"""

import pyodbc
import pandas as pd
from math import floor, log10
from statsmodels.sandbox.stats.multicomp import multipletests

Connstr_icd = '' # constructor for SQL connection
Cursor = '' # constructor for SQL connection
TABLE = '' # SQL table
TARGETS = '' # targets list

def gender_perc (target):
    # get table for target
    table = TABLE
    # Screened females
    df_female = pd.read_sql("SELECT * FROM " + table + " WHERE  Gender = 'W'", pyodbc.connect(Connstr_icd))
    # Hit females
    df_hits_female = pd.read_sql("SELECT * FROM " + table + " WHERE hits = 'positive' AND Gender = 'W'", pyodbc.connect(Connstr_icd))   
    # Screened males
    df_male = pd.read_sql("SELECT * FROM " + table + " WHERE  Gender = 'M'", pyodbc.connect(Connstr_icd)) 
    # Hit males
    df_hits_male = pd.read_sql("SELECT * FROM " + table + " WHERE hits = 'positive' AND Gender = 'M'", pyodbc.connect(Connstr_icd))   
    df_hits = pd.read_sql("SELECT * FROM " + table + " WHERE hits = 'positive'", pyodbc.connect(Connstr_icd))   
    df_nonhits = pd.read_sql("SELECT * FROM " + table + " WHERE hits = 'negative'", pyodbc.connect(Connstr_icd))   
 
    # Calculate percentage    
    perc_male = 100*len(df_male)/len(df_nonhits)
    perc_male_hits = 100*len(df_hits_male)/len(df_hits)
    perc_female = 100*len(df_female)/len(df_nonhits)
    perc_female_hits = 100*len(df_hits_female)/len(df_hits)
    
    p_chi = None
    if len(df_hits_female)>=5 and len(df_hits_male)>=5:
        obs=[[len(df_hits_female), len(df_hits_male)],[len(df_female)-len(df_hits_female), len(df_male)-len(df_hits_male)]]
        chi2, p_chi, dof, expected = stats.chi2_contingency (obs)
        
    return perc_male, perc_male_hits, perc_female, perc_female_hits, p_chi, target

p_values = []
perc_male_hits_list = []
perc_female_hits_list = []
perc_male_list = []
perc_female_list = []
target_list = []
for target in TARGETS:
    perc_male, perc_male_hits, perc_female, perc_female_hits, p_chi, target_p = gender_perc (target)
    if p_chi is not None:
        perc_male_hits_list.append (round(perc_male_hits,2))
        perc_female_hits_list.append (round(perc_female_hits,2))
        perc_male_list.append (round(perc_male,2))
        perc_female_list.append (round(perc_female,2))
        p_values.append (p_chi)
        target_list.append (target_p)
      
df_gender = pd.DataFrame ({'Target':target_list, 'P_values':p_values})

#Bonferroni correction
alpha = 0.05
lista = []
p_adjusted = multipletests(df_gender['P_values'].tolist(), alpha = alpha,  method='bonferroni')
exponent = [int(floor(log10(abs(i)))) for i in p_adjusted[1]]
coeff = [round(i/float(10**j),2) for i,j in zip(p_adjusted[1],exponent)]
p_val_format = [r"{}e{}".format(i, j)  if int(j)!=0 else r"{}".format(i) for i,j in zip(coeff,exponent)]

df_gender['Bonferroni'] =  p_val_format
df_gender['Signif'] = p_adjusted[0]
df_gender['Perc male hits'] = perc_male_hits_list
df_gender['Perc female hits'] = perc_female_hits_list
df_gender['Perc male nonhits'] = perc_male_list
df_gender['Perc female nonhits'] = perc_female_list

significance = []
for pvalue in df_gender['Bonferroni'].tolist():
    if float(pvalue) > 0.05:
        significance.append ('ns')
    else:
        if float(pvalue) <= 0.0001:
            significance.append ('****')   
        elif float(pvalue) <= 0.001:
            significance.append ('***')  
        elif float(pvalue) <= 0.01:
            significance.append ('**')    
        elif float(pvalue) <= 0.05:
            significance.append ('*')   
df_gender['significance'] = significance

print (df_gender)
 