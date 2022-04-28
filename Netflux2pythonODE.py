# Netflux2pythonODE.py
# Based on Netflux2ODE.m
# CZZ on 12/13/2021

def Netflux2pythonODE(CNAmodel,xlsfilename):
    """
    Generates the system of differential equations
    
    Netflux2pythonODE takes the CNAmodel dict generated by xls2Netflux and
    generates the system of normalized Hill differential equations. Outputs
    ODElist and paramList are sent to ODE to solve. 
    10/22/21 by CZZ
    """
    from getNetfluxParams import getNetfluxParams
    import pandas as pd
    
    # %% Pull in info from CNAmodel

    interMat = CNAmodel["interMat"]
    reactantMat = CNAmodel["reactantMat"]
    productMat = CNAmodel["productMat"]
    notMat = CNAmodel["notMat"]
    speciesNames = CNAmodel["specID"]
    modelname = CNAmodel["net_var_name"]
    
    numspecies = len(interMat.index) # number of rows
    numrxns = len(interMat.columns)
    
    # %% read parameters from Excel file
    [paramList, CNAerror] = getNetfluxParams(xlsfilename)
    
    # %% Create ODElist, a cell array of strings containing the ODEs
    reactions = calcReactions(reactantMat,productMat,notMat,speciesNames)
    mapVar = [0] * numspecies
    for i in range(numspecies):
        mapVar[i] = '\t' + speciesNames.iloc[i] + " = " + str(i); # downshifted for Python
    pythonODElist = []
    pythonODElist.append('\tdydt = np.zeros(%i)' %numspecies) #Create matrix of Strings with Differential Equations
    
    # %Create matrix of Strings with Differential Equations
    for i in range(numspecies):
        s = speciesNames.iloc[i]
        if pd.isna(paramList[4].iloc[i]):
            yMaxStr = paramList[4].iloc[i]
        else:
            yMaxStr = 'ymax[' + s + ']'
        pythonODElist.append(f'\tdydt[{s}] = ({reactions[i]}*{yMaxStr} - y[{s}])/tau[{s}]')
    pythonODElist.append('\treturn dydt')
    for each in pythonODElist:
        mapVar.append(each)
    pythonODElist = mapVar
    
    return paramList, pythonODElist, CNAerror

    
# %% function calcReactions
def calcReactions(reactantMat, productMat,notMat,speciesNames):
    numspecies = len(reactantMat.index)
    numrxns = len(reactantMat.columns)
    rxnString = []
    
    # % generate reaction strings for each reaction of the form: act(y(A),rpar(:,1))
    for i in range(numrxns):
        # If index of species in rxn is before rxn arrow, it's a reactant (-1). If not, 0.
        # print(reactantMat)
        reactants = [i for i, x in enumerate(reactantMat.iloc[:,i]) if x == -1]
        if len(reactants)==0:                  # input reaction, of the form: w[3]
            my_str = 'w[' + str(i) + ']'
        elif len(reactants)==1:            # single reactant, of the form: act(y[A],w[3],n[3],EC50[3])
            if notMat.at[reactants[0],i] == 0:
                my_str = 'inhib(y[' + speciesNames.iloc[reactants[0]] + '],w[' + str(i) + '],n['+ str(i) +'],EC50[' + str(i) + '])'
            else:
                my_str = 'act(y[' + speciesNames.iloc[reactants[0]] + '],w[' + str(i) + '],n[' + str(i) + '],EC50[' + str(i) + '])'
        else: # multiple reactants, of the form: AND(w[3],[act(y[A],w[3],n[3],EC50[3]),inhib(y[B],w[3],n[3],EC50[3])])
            my_str = 'AND(w[' + str(i) + '],['
            for j in range(len(reactants)):         # loop over each reactant
                if notMat.at[reactants[j],i] == 0:
                    my_str = my_str + 'inhib(y[' + speciesNames.iloc[reactants[j]] + '],w[' + str(i) + '],n[' + str(i) + '],EC50[' + str(i) + '])'
                else:
                    my_str = my_str + 'act(y[' + speciesNames.iloc[reactants[j]] + '],w[' + str(i) + '],n[' + str(i) + '],EC50[' + str(i) + '])'
                if j < len(reactants): # more reactants to come
                    my_str = my_str + ','
            my_str = my_str + '])' # cap with close parentheses
        rxnString.append(0)
        rxnString[i] = my_str.replace(",)",")")
    
    
    # combine relevant reactions for each species, using nested ORs if needed
    reactions = []
    for i in range(numspecies):
        rxnList = [i for i, x in enumerate(productMat.iloc[i,:]) if x == 1] # find reactions that affect specie 'i'
        if len(rxnList)==0:                  # 0 reactions for that specie
            reactions.append(0)    
            reactions[i] = '0'
        elif len(rxnList)==1:                # 1 reaction for that specie
            reactions.append(0)
            reactions[i] = rxnString[rxnList[0]]
        else: # combine reactions with nested 'OR's
            my_str = ''
            for j in range(len(rxnList)-1):
                my_str = my_str + 'OR(' + rxnString[rxnList[j]] + ','
            my_str = my_str + rxnString[rxnList[-1]] + ')' #+ repmat(')',1,j)];
            reactions.append(0)
            if my_str.count("(") != my_str.count(")"):
                my_str += ")"
            reactions[i] = my_str
            # print(reactions[i])
                
    return reactions

