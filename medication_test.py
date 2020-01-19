# -*- coding: utf-8 -*-
"""
Medication test

Calculate the occurence of all medications in negative and positive patient 
samples. 
"""

import pyodbc
import pandas as pd
from collections import Counter
from scipy import stats
from statsmodels.sandbox.stats.multicomp import multipletests
from math import floor, log10

Connstr = '' # constructor for SQL connection
TABLE = ''
MEDICATION_TABLE = ''

def get_medication (type_hits):
    db_query =  '' # SQL query to extract data
    df_med =  pd.read_sql(db_query, pyodbc.connect(Connstr))
    return df_med

def count_medication (df):
    all_medication = []
    for index, row in df.iterrows():
        all_medication = all_medication + [row['Medication']]
    return Counter(all_medication)

def process_medication (type_hits):
    df = get_medication (type_hits)
    counter = count_medication (df) 

    medication = []
    counts = []
    for i in counter:
        medication.append (i)
        counts.append (counter[i])
    df = pd.DataFrame ({'medication':medication, 'counts':counts})
    return df

def calc_no_patients (type_hits):
    db_query = (" SELECT COUNT(*) as TotalNo FROM " + TABLE 
                + " WHERE hits =  '" + type_hits + "' ")
    df =  pd.read_sql(db_query, pyodbc.connect(Connstr))
    return df

df_pos = process_medication ('positive')
df_pos.head()

df_neg = process_medication ('negative')
df_neg.head()

no_pos = calc_no_patients ('positive')
no_neg = calc_no_patients ('negative')

df_pos_count_old = df_pos.copy()
df_neg_count_old = df_neg.copy()

df_pos['perc'] = (df_pos['counts']*100)/no_pos['TotalNo'].iloc[0]
df_pos.sort_values(by='counts', ascending=False)[:5] 

df_neg['perc'] = (df_neg['counts']*100)/no_neg['TotalNo'].iloc[0]
df_neg.sort_values(by='counts', ascending=False) [:5]

name = []
pos_perc = []
neg_perc = []
pos_no = []
neg_no = []
for med in df_pos['medication'].tolist():
    name.append (med)
    pos_perc.append (df_pos[df_pos['medication']==med]['perc'].iloc[0]) 
    pos_no.append (df_pos[df_pos['medication']==med]['counts'].iloc[0]) 
    if len(df_neg[df_neg['medication']==med]) > 0:
        neg_perc.append (df_neg[df_neg['medication']==med]['perc'].iloc[0])
        neg_no.append (df_neg[df_neg['medication']==med]['counts'].iloc[0])
    else:
        neg_perc.append (0)
        neg_no.append (0)
    
df_comp = pd.DataFrame ({'medication':name, 'pos_perc':pos_perc, 'pos_no':pos_no, 'neg_perc':neg_perc, 'neg_no':neg_no})
df_comp.sort_values(by='pos_perc', ascending=False) [:5]

# Calculate significance hits
def process_pvalue(p_value):
    exponent = [int(floor(log10(abs(i)))) for i in p_value]
    coeff = [round(i/float(10**j),2) for i,j in zip(p_value,exponent)]
    p_val_format = [r"{}e{}".format(i, j)  if int(j)!=0 else r"{}".format(i) for i,j in zip(coeff,exponent)]
    return p_val_format

df_comp_val = df_comp[(df_comp['pos_no'] >=5) & (df_comp['neg_no'] >= 5)]
p_value = []
for index, row in df_comp_val.iterrows ():
    obs = [[row['pos_no'], row['neg_no']],[no_pos['TotalNo'].iloc[0]-row['pos_no'],no_neg['TotalNo'].iloc[0]-row['neg_no']]]
    chi2, p_chi, dof, expected = stats.chi2_contingency (obs)
    p_value.append (p_chi)
df_comp_val['p value'] = p_value
df_comp_val.head()

# Adjust for multiple comparisons
df_comp_val = df_comp_val.sort_values (['p value'], ascending=True)
methods = ['bonferroni'] 
alpha = 0.05
for method in methods:
    p_adjusted = multipletests(df_comp_val['p value'].tolist(), alpha = alpha, method=method)
    df_comp_val[method]=p_adjusted[1]
    df_comp_val['significant ' +method]=p_adjusted[0]
    
df_comp_val['pos_perc'] = [round(i,2) for i in df_comp_val['pos_perc']] 
df_comp_val['neg_perc'] = [round(i,2) for i in df_comp_val['neg_perc']]  
df_comp_val['p value'] = process_pvalue(df_comp_val['p value'].tolist())
df_comp_val[method] = process_pvalue(df_comp_val[method].tolist())
    
significance = []
for pvalue in df_comp_val['bonferroni'].tolist():
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
df_comp_val['significant bonferroni'] = significance

print (df_comp_val[df_comp_val['significant bonferroni']==True])

df_comp_val.to_excel ('Medication.xlsx')
