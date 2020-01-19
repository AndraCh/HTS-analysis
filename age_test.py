# -*- coding: utf-8 -*-
"""
Median age test in all targets from HTS screen
"""
import pyodbc
import pandas as pd
from math import floor, log10
import scipy
import numpy as np
from statsmodels.sandbox.stats.multicomp import multipletests

Connstr_icd = '' # constructor for SQL connection
TABLE = '' # SQL table
TARGETS = '' # targets list

def get_age_target (table, type_hits):
    db_query = "SELECT Age FROM "+ table + " WHERE hits = '" + type_hits + "'"
    data_frame = pd.read_sql(db_query, pyodbc.connect(Connstr_icd)) 
    age_list = data_frame['Age'].tolist()
    return age_list

def median_test_age (targets):
    median_age_hits = []
    median_age_non_hits = []
    p_val = []
    p_value_format = []
    no_hits = []
    no_non_hits = []
    for target in targets:
        table = TABLE
        age_hits = get_age_target (table, 'positive')
        no_hits.append (len(age_hits))
        age_non_hits = get_age_target (table, 'negative')
        no_non_hits.append (len(age_non_hits))
        stat, p_value, med, tbl = scipy.stats.median_test(age_hits, age_non_hits)
        median_age_hits.append (np.median(age_hits))
        median_age_non_hits.append (np.median(age_non_hits))
        p_val.append (p_value)
        exponent = int(floor(log10(abs(p_value))))
        coeff = round(p_value / float(10**exponent), 2)
        p_value_f = r"{}e{}".format(coeff, exponent) #r"{}.10^{}".format(coeff, 2), exponent)
        p_value_format.append (p_value_f)
        
    df = pd.DataFrame ({'Target':targets, 'Median hits':median_age_hits, 'Median non hits':median_age_non_hits, 'P value':p_val,
                       '# hits':no_hits, '# non hits':no_non_hits})
    alpha = 0.05
    p_adjusted = multipletests(df['P value'].tolist(), alpha = alpha,  method='bonferroni')
    exponent = [int(floor(log10(abs(i)))) for i in p_adjusted[1]]
    coeff = [round(i/float(10**j),2) for i,j in zip(p_adjusted[1],exponent)]
    p_val_format = [r"{}e{}".format(i, j)  if int(j)!=0 else r"{}".format(i) for i,j in zip(coeff,exponent)]
    df['Bonferroni'] = p_val_format
    df['Signif'] = p_adjusted[0]
    df['P value'] = p_value_format 
    return df

df = median_test_age (TARGETS)

significance = []
for pvalue in df['Bonferroni'].tolist():
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
df['significance'] = significance

print (df)
