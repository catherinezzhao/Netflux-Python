# xls2Netfluxpython.py
# Based on xls2Netflux.m
# CZZ on 12/13/2021

def xls2Netfluxpython(networkID, xlsfilename):
    """
    Reads the Netflux Excel spreadsheet and converts them to differential 
    equations
    
    Converts an Excel spreadsheet of specified format into a Normalized-Hill 
    differential equation system
    
    Notes:
        1) xls2Netflux makes a number of assumptions about the xlsfile
           structure and will throw errors if there are problem: 'species' 
           sheet: ID's in column B, starting in row 3, names in column C 
           'reactions' sheet: ID's in column B starting in row 3, reaction 
           rule in column C
    
    Dependencies: cna2Netflux load XLS network reconstruction and read
    species, reaction rules
    """
    import pandas as pd
    import numpy as np
    import regex
    from Netflux2pythonODE import Netflux2pythonODE
    
    
    # READ IN SPECIES SHEET
    spec_df = pd.read_excel(xlsfilename, sheet_name='species')    
    spec_df.columns = spec_df.iloc[0] # use row 0 as headers
    specID = spec_df['ID'][1:] # In 'ID' col, get everything after 1st value (index 0)
    specNames = spec_df['name'][1:]
    numSpecs = len(specID)
    
    # READ IN REACTIONS SHEET
    rxn_df = pd.read_excel(xlsfilename, sheet_name='reactions')    
    rxn_df.columns = rxn_df.iloc[0] # use row 0 as headers
    rxnID = rxn_df['ID'][1:] # exclude 1st row
    rxnRules = rxn_df['Rule'][1:]
    numRxns = len(rxnID)

    # CREATE interMat AND notMat
    interMat = pd.DataFrame(0, index=range(numSpecs), columns=range(numRxns)) # creates df of 0s that's num of specs by num of rxns
    notMat = pd.DataFrame(0, index=range(numSpecs), columns=range(numRxns))
    rMat = pd.DataFrame(0, index=range(numSpecs), columns=range(numRxns))
    pMat = pd.DataFrame(0, index=range(numSpecs), columns=range(numRxns))
    
    # POPULATING interMat
    for i in range(numRxns):
        # Calculate interMat column for reaction 'i'
        # specIDReg is a list of indices or NaN
        # For each species, if species is in the rxn string, the index will be stored
        # e.g., If '=> A' is the rxn string, then specIDReg contains [3, NaN, NaN, ...]
        specIDReg = pd.Series(0, index=range(numSpecs))
        for j in range(numSpecs):
            # If the species is in the reaction
            if (specID.iloc[j] in rxnRules.iloc[i]):
                # Store the index of where the species was found in the rxn string
                specIDReg[j] = rxnRules.iloc[i].find(specID.iloc[j])
            else:
                specIDReg[j] = None
        
        # Position of the reaction arrow
        eqPos = rxnRules.iloc[i].find('=>')
        
        # Parse out the reactants and products from the string
        reactants = pd.Series(0, index=range(numSpecs))
        products = pd.Series(0, index=range(numSpecs))
        for j in range(numSpecs):
            array = specIDReg.iloc[j]
            # MATLAB needed another for loop because it's matrix based, but Python does need the following line:
            # for k in range(array.length)):
            # If index of species in rxn is before rxn arrow, it's a reactant
            if array < eqPos:
                reactants.iloc[j] = -1
            # If index of species in rxn is after rxn arrow, it's a product
            elif array > eqPos:
                products.iloc[j] = 1
        
        for j in range(numSpecs):
            # .at() allows access to individual cells
            # It's j,i because j is each spec (row) and i is each rxn (col)
            interMat.at[j,i] = reactants.iloc[j] + products.iloc[j]
            
            rMat.at[j,i] = reactants.iloc[j]
            pMat.at[j,i] = products.iloc[j]
        
        for j in range(numSpecs):
            # If the species is in the reaction
            # "!"+specID.iloc[j] checks if ! is immediately proceeded by
            # the specific species (row), and that's in the rxn rule we're
            # looking at
            # e.g., "!D" in "C & !D => E"
            if "!"+specID.iloc[j] in rxnRules.iloc[i]:
                notMat.at[j,i] = 0
            else:
                notMat.at[j,i] = 1
        
        
    # Create model dict CNAmodel and then create ODE files via cna2Netflux
    # CNAmodel is a dict
    CNAmodel = {}
    CNAmodel["specID"] = specID # keys are strings; in .m, these are fields
    CNAmodel["interMat"] = interMat
    CNAmodel["reactantMat"] = rMat
    CNAmodel["productMat"] = pMat
    CNAmodel["notMat"] = notMat
    CNAmodel["net_var_name"] = networkID # taken in as argument
    # CNAmodel["type"] = spec_df['type'][1:]
    # CNAmodel["location"] = spec_df['location'][1:]
    # CNAmodel["module"] = spec_df['module'][1:]
    # CNAmodel["ref"] = spec_df['ref'][1:]
    
    
    (paramList, ODElist, CNAerror) = Netflux2pythonODE(CNAmodel,xlsfilename)

    ReadError = CNAerror
    
    return specID, rxnID, rxnRules, paramList, ODElist, CNAmodel, ReadError

