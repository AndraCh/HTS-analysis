# -*- coding: utf-8 -*-
"""
ICDs test

Calculate the occurence of all ICD codes in negative and positive patient 
samples. 
Select the first description fro ICD description table as the ICD code 
description. 
This script is using high evel diagnoses codes (A40, A41, ..).
"""

import pyodbc
import pandas as pd
from math import floor, log10
import scipy
import numpy as np
from statsmodels.sandbox.stats.multicomp import multipletests

Connstr_icd = '' # constructor for SQL connection
Cursor = '' # constructor for SQL connection
TABLE = '' # SQL table
TARGETS = '' # targets list

def common_codes_positive_negative (type_patients, table):
    # Create query
    db_query = '' # SQL query to extract data
    # Execute query and get all codes
    list_codes = pd.read_sql(db_query, pyodbc.connect(Connstr_icd))['ICDCode'].tolist() 
    # Select unique counts
    unique_list_codes = set (list_codes)
    count_appearance_list = []
    unique_id_codes = []
    for unique_code in unique_list_codes:
        count_appearance_list.append (list_codes.count(str(unique_code)))
        unique_id_codes.append (unique_code)
    df_occ = pd.DataFrame ({'code':unique_id_codes, 
                            'occurrence':count_appearance_list}).sort_values(['occurrence'], ascending=False)
    # Get ICD code description 
    code = []
    description = []
    group = []
    indexes_list = []
    for index, row in df_occ.iterrows():
        indexes_list.append (index)
        db_query = '' # SQL query to extract data
        Cursor_icds.execute (db_query)
        for row_sql in Cursor_icds:
            code.append (row_sql[0])
            description.append (row_sql[1])
            group.append (row_sql[2])            
    df_occ['code_name'] = code
    df_occ['description'] = description
    df_occ['group'] = group
    # Return data frame
    return df_occ.sort_values(["occurrence"],ascending=False)  

def get_patients_code (type_patients, code, no_patients):
    """
    Get the total number of positive/negative sample for a specific ICD code.
    """
    # Create query
    db_query = '' # SQL query to extract data
    Cursor.execute (db_query)
    # Get no codes
    for row in Cursor:
        no_code =  row[0]
    # Create array codes
    list_pos = [1 for i in range(0,int(no_code))]
    list_neg = [0 for i in range(0,no_patients - int(no_code))]
    list_values_code = list_pos + list_neg
    return list_values_code , no_code

def calc_codes_significance (df_occ_pos, df_occ_neg, type_hits, table):
    """
    Calculate the significance of each ICD code in positive versus negative 
    samples based on p value computed using a chi square test.  
      
    Chiq-square test Python  
    *obs = [[no of positive patient samples with the ICD code (N1), 
             no of positive patient samples with the ICD code (N2)],    
           [total number of positive patients - N1,
           total number of negative patients - N2]]*  
     Test:  
     *chi2, p_chi, dof, expected = stats.chi2_contingency (obs)*
     """
    # Get screened patients
    db_query = (" SELECT * FROM " + table)   
    data_frame = pd.read_sql(db_query, pyodbc.connect(Connstr_icd)) 
    
    df_codes = df_occ_pos
    #type_hits = 'positive'  
    hits_codes =  df_codes["code"].tolist()
    hits_codes_names = df_codes["code_name"].tolist()
    hits_description = df_codes["description"].tolist() 
    ind = 0
    p_value_list = []
    code_list = []
    code_list_desc = []
    n_code_1_list = []
    n_code_2_list = []
    perc_pos = []
    perc_neg = []
    for code in  hits_codes:
        if ind < len(hits_codes):
           
            pos_patients, no_code1 = get_patients_code ("positive",code,len(data_frame[data_frame['hits']=='positive']))
            neg_patients, no_code2 = get_patients_code ("negative",code,len(data_frame[data_frame['hits']=='negative']))
            obs = [[sum(pos_patients), sum(neg_patients)],[len(pos_patients)-sum(pos_patients),len(neg_patients)-sum(neg_patients)]]
            #print ( hits_codes_names[ind], obs)
            chi2, p_chi, dof, expected = stats.chi2_contingency (obs)
            table = [[sum(pos_patients),sum(neg_patients)],[len(pos_patients)-sum(pos_patients),len(neg_patients)-sum(neg_patients)]]
            oddsratio, p_chi = scipy.stats.fisher_exact(table, alternative='two-sided') 
            p_value_list.append (p_chi) #(float("%0.4f" % (p_chi)))  
            code_list.append ( hits_codes_names[ind])
            code_list_desc.append (hits_description[ind])
            n_code_1_list.append (no_code1)
            n_code_2_list.append (no_code2)
            perc_pos.append (round(sum(pos_patients)*100/len(pos_patients),2))
            perc_neg.append (round(sum(neg_patients)*100/len(neg_patients),2))
            ind = ind + 1
    # Create data frame
    df_rank_codes = pd.DataFrame ({'code':code_list, 'description':code_list_desc, 
                                   'p value':p_value_list, 'no hits':n_code_1_list, 
                                   'no non hits':n_code_2_list,
                                   'perc hits': perc_pos,
                                   'perc non hits': perc_neg})
    
    # Adjust for multiple comparisons
    df_rank_codes = df_rank_codes.sort_values (['p value'], ascending=True)
    methods = ['bonferroni']  
    alpha = 0.05 
    for method in methods:
        p_adjusted = multipletests(df_rank_codes['p value'].tolist(), alpha = alpha, method=method)
        df_rank_codes[method]=p_adjusted[1]
        df_rank_codes['significant ' +method]=p_adjusted[0]
    return df_rank_codes

def get_ranked_codes (table):
    """
    Select all ICD codes for negative and positive patients. 
    Use a threshold for a minimum no of observations for each ICD code. 
    Many statistical tests require at least 5 observations per sample.
    """
    # Extract code for positive and negative 
    df_occ_pos = common_codes_positive_negative ("positive", table)
    df_occ_neg = common_codes_positive_negative ("negative", table)
    # Minimum number of codes 
    min_codes = 5
    # Select ICD codes which appear at least e.g. min_codes = 5 times
    df_occ_pos = df_occ_pos[df_occ_pos['occurrence'] >=min_codes]
    df_occ_neg = df_occ_neg[df_occ_neg['occurrence'] >= min_codes]
    return df_occ_pos, df_occ_neg

# Calculate significance hits
def process_pvalue(p_value):
    exponent = [int(floor(log10(abs(i)))) for i in p_value]
    coeff = [round(i/float(10**j),2) for i,j in zip(p_value,exponent)]
    p_val_format = [r"{}e{}".format(i, j)  if int(j)!=0 else r"{}".format(i) for i,j in zip(coeff,exponent)]
    return p_val_format

df_occ_pos, df_occ_neg = get_ranked_codes (TABLE)
df_rank_codes = calc_codes_significance (df_occ_pos, df_occ_neg, 'positive', TABLE)
df_rank_codes['p value'] = process_pvalue(df_rank_codes['p value'].tolist())
df_rank_codes['bonferroni'] = process_pvalue(df_rank_codes['bonferroni'].tolist())

significance = []
for pvalue in df_rank_codes['bonferroni'].tolist():
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
df_rank_codes['significance'] = significance

df_rank_codes.head()