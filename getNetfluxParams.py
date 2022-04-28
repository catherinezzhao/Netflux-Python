# getNetfluxParams.py
# Based on getNetfluxParams.m
# CZZ on 12/13/2021

import pandas as pd
import numpy as np
import regex

# PARAMETER MISMATCH FUNCTIONS
def mismatchParams(rxnRules, rxnID, w, n, EC50):
    # Returns the reactionID and reaction rule strings that don't have
    # defined parameters, and changes parameters for those reactions to the
    # default values
    noparams = []
    for i in range(len(rxnRules)):
        # if value in n, w, EC50 is NaN (empty) or a string
        if pd.isna(n.iloc[i]) or pd.isna(w.iloc[i]) or pd.isna(EC50.iloc[i]) or isinstance(n.iloc[i], str) or isinstance(w.iloc[i], str) or isinstance(EC50.iloc[i], str):
            noparams.append(rxnID.iloc[i] + ": " + rxnRules.iloc[i])
            w[i] = 1
            n[i] = 1.4
            EC50[i] = 5
    return noparams, w, n, EC50


def mismatchSpecParams(specID, y0, ymax, tau):
    # % Returns the species IDs that don't have defined parameters and sets the
    # % parameters to the default values for those species. 
    noSpecParams = []
    for i in range(len(specID)):
        if pd.isna(y0.iloc[i]) or pd.isna(ymax.iloc[i]) or pd.isna(tau.iloc[i]) or isinstance(y0.iloc[i], str) or isinstance(ymax.iloc[i], str) or isinstance(tau.iloc[i], str):
            noSpecParams.append(specID.iloc[i])
            y0[i] = 0
            ymax[i] = 1
            tau[i] = 1
    return noSpecParams, y0, ymax, tau


def getNetfluxParams(xlsfilename):
    # Extracts the parameters from Netflux Excel file
    
    # READ IN SPECIES SHEET
    spec_df = pd.read_excel(xlsfilename, sheet_name='species')    
    spec_df.columns = spec_df.iloc[0] # use row 0 as headers
    y0 = spec_df['Yinit'][1:] # exclude 1st row
    ymax = spec_df['Ymax'][1:]
    tau = spec_df['tau'][1:]
    
    # READ IN REACTIONS SHEET
    rxn_df = pd.read_excel(xlsfilename, sheet_name='reactions')    
    rxn_df.columns = rxn_df.iloc[0] # use row 0 as headers
    w = rxn_df['Weight'][1:]
    n = rxn_df['n'][1:]
    EC50 = rxn_df['EC50'][1:]
    
    specID = spec_df['ID'][1:]
    rxnID = rxn_df['ID'][1:]
    rxnRules = rxn_df['Rule'][1:]
    
    # make sure all reactions have parameters. mismatchParams returns the ID
    # and reaction rule strings which do not have corresponding reaction
    # parameters, and sets default parameters to those reactions.

    noparams = '';
    (noparams, w, n, EC50) = mismatchParams(rxnRules, rxnID, w, n, EC50)
    
    # make sure all species have defined parameters, and set default
    # parameters for those that don't. 
    noSpecParams = '';
    (noSpecParams, y0, ymax, tau) = mismatchSpecParams(specID, y0, ymax, tau)
    
    # 1) remove reaction and ID from respective lists if reaction rules aren't
    # present, if rules present but no ID, nothing is changed.
    # 2) make sure that all reactions have defined paramaters
    
    temp_ind_lst = []
    for i in range(len(rxnRules)):
        if pd.isna(rxnRules.iloc[i]):
            temp_ind_lst.append(i)
    rxnRules = rxnRules.drop(rxnRules.index[temp_ind_lst])
    rxnID = rxnID.drop(rxnID.index[temp_ind_lst])
    
    # delete species with no specID
    specID = specID.dropna()
    
    # Create CNA error list
    CNAerror = []
    if len(noparams)>0:
        CNAerror.append('Missing parameter(s) for following reactions')
        CNAerror.append('(default parameters assigned):')
        for i in range(len(noparams)):
            CNAerror.append(noparams[i])
        CNAerror.append('')
        raise ValueError('ParamErr:MissingParams','Reaction or species parameters may be missing');

    if len(noSpecParams)>0:
        CNAerror.append('Missing parameter(s) for following species')
        CNAerror.append('(default parameters assigned):')
        for i in range(len(noSpecParams)):
            CNAerror.append(noSpecParams[i])
        CNAerror.append('')
        raise ValueError('ParamErr:MissingParams','Reaction or species parameters may be missing');

    if specID.duplicated().any():
        CNAerror.append('Warning: Duplicate species detected')
        raise ValueError('DuplicateError:DuplicateSpecies', 'Duplicate species detected')

    paramList = [w,n,EC50,tau,ymax,y0] # parameter list for Netflux2pythonODE.py
    return paramList, CNAerror

getNetfluxParams("exampleNet.xlsx")

