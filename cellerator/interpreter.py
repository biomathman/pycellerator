#!/usr/bin/env python
#
# Cellerator Interpreter Module
# B.E.Shapiro 21 Feb 2012 
# bug fixes 11/11/14 did not parse X + X -> X correctly 
# rev 8/30/2015
#
#
#****************************************************************************
#
#    This file is part of Cellerator
#
#    Cellerator converts reactions, expressed in a text-formatted arrow-based 
#    notation into differential equations and performs numerical simulations.
#    Copyright (C) 2012 Bruce E Shapiro.
#    
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>
#    
#****************************************************************************

#
#Linux Command Line Syntax:
#
#  python interpreter.py
#       -in inputfilename
#       [-out outputfile]
#       [-dump]    # to dump the parser
#

# from pyparsing import *
import os.path
import sys
from parser import Reaction,StripComments
from sympy import *
from expander import expand
import datetime
from utils import *
from reader import readmodel
from time import clock




#*******************************************************************************
        
def invokeParser(lines, dump=False):
    """
    input: list of lines of text data in the format of Cellerator reactions
    output: list of parsed reactions (parser data structures)
    """
    results=[]
    for line in lines:
        line = StripComments(line,"#")
        if line=="": continue # ignore blanks and/or comments
        r=Reaction(line)
        if dump: r.Dump()
        results.append(r)
    return (results)

#*******************************************************************************

def runparser(datafile, trace=False):
    if ("-in" in sys.argv):
        i = (sys.argv).index("-in")+1
        if i<len(sys.argv):
            infile=sys.argv[i]
        else:
            raise SystemExit("Error: expecting filename after -in")
    elif os.path.isfile(datafile):
        infile=datafile
    else:
        raise SystemExit("Error: interpreter (runparser): No Input File Given; use -in option.")
    if not os.path.isfile(infile):
        raise SystemExit("Error: "+os.path.abspath(infile)+" File not found.")
    
    rates = {}
    if os.path.splitext(infile)[1]==".model":
        (lines, ics, rates, frozenvars, functions,assignments, fullfile)=readmodel(INFILE=infile)
        
    else:
        f=open(infile)         
        lines = f.readlines()
        f.close()
    results=invokeParser(lines,dump=trace)
    # print "runparser: rates=", rates
    return (results, rates, functions, frozenvars, assignments, ics)

#*******************************************************************************

def makeSymbolDictionary(reactions, inputrates, IC):
    """
    generates a dictionary of Sympy symbols from a list of reactions
        {"variable":symbol, "variable":symbol, ...}
    "variable" - string value of variable in input (text file)
    symbol - symbol used by Sympy to represent the variable
    --------------------------------------------------------------
    In fact, the name of the symbol generated will be identical to
    the name of the string used in the input file. The purposed of
    generating the dictionary is to generate all the initial Sympy
    unique Symbol assignments so that they can subsequently be
    referenced by their string name. 
    
    """
    
    # print "makeSymbolDictionary: inputrates:", inputrates
    
    symbolDictionary = {"Nil":1, "t":0}

    def ismath (x):
        if type(x) != type("Hello World"): return False
        if len(x)<3:      return False
        if x[0]!= "\"":   return False
        if x[-1] != "\"": return False
        expr = x[1:-1]
        try:
            math=sympify(expr)
        except:
            print "Error: Interpreter: Unable to sympify the expression "+x
            symbolDictionary["Error"]=1
            return False
        return True
        
        
    def _addToDict(lis):
     

       for x in lis:
            # print "_addToDict: x: ", x
            if type(x) == type([]):
                map(_addToDict, x)
                continue
            if x not in symbolDictionary:
                if not is_number(x):
                    if ismath(x):
			            symbolDictionary[x] = sympify(x[1:-1])
                    else:
                        symbolDictionary[x]=Symbol(x)    
                    
    
    for reaction in reactions:
        _addToDict(flatten(list(reaction.LHS())))
        _addToDict(list(reaction.RHS()))
        _addToDict(flatten(list(reaction.MHS())))        
        rates = flatten(list(reaction.Rates()))
        _addToDict(rates)
    
    # print "makeSymbolDictionary: rates = ", inputrates
    for x in inputrates:
        if x in symbolDictionary:
            # print x, " is in the dictionary already"
            continue 
        else: 
            # print x, " is not in the dictionary, so add it"
            symbolDictionary[x] = Symbol(x)
       
    # print "makeSymbolDictionary: IC = ", IC  
    for x in IC: 
        if x in symbolDictionary:
            continue
        else:
            symbolDictionary[x] = Symbol(x)
            
    if "-DICTIONARY" in map(lambda x:x.upper(), sys.argv):
        print "***************** Symbol Dictionary *****************\n",symbolDictionary,\
            "\n****************************************************"
    return symbolDictionary


#*******************************************************************************


def make_mass_action_ODE_Terms(r, ODETerms, symbolDictionary, frozen):

        DBG = False
        if DBG: print "make_mass_action_ODE_Terms"

        def rint(x,digs):
            if abs(round(x,digs)-int(x))<10**-digs:
                return int(x)
            return round(x,digs)

        LS = zip(listify(r.LHS()), \
                 map(lambda x: rint(float(x),12), list(r.LH_STOIC())))
        RS = zip(listify(r.RHS()), \
                 map(lambda x: rint(float(x),12), list(r.RH_STOIC())))
                 
        if DBG: print "LS = ", LS      
             
        LD = list2Dict(LS)
        RD = list2Dict(RS)
        # print "LD = ", LD, " RD = ", RD

        species = list ( set([x for (x,s) in LS]) |  \
                         set([x for (x,s) in RS]) |  \
                         set( list(r.MHS())))
#        species = [x for (x, s) in LS] + [x for (x, s) in RS] + list(r.MHS())

        if DBG: print "species = ", species
        if DBG: print "list(r.LHS()):", list(r.LHS())
        NS = {}
        baseterm = 1
        for s in set(r.LHS()):
            sym = symbolDictionary[s]
            exp = LD[s]
            if exp == 1.0:
                baseterm *= sym
            else:
                baseterm *= (sym**exp)
        if DBG: print "LD: ", LD, " RD: ", RD
        # stoich = 0
        for s in species:
            if DBG: print "loop species = ", s
            stoich=0
            if s in LD: stoich -= LD[s]
            if s in RD: stoich += RD[s]
            NS[s] = stoich
            if s in LD: 
                if DBG: print "LD[s]=",LD[s]
            if s in RD: 
                if DBG: print "RD[s]=", RD[s]
            if DBG: print "loop: s=", s, " NS[s]= ", NS[s]
         # print "NS = ", NS

        if DBG: print "make_mass_action_ODE_Terms: baseterm: ", baseterm

        if DBG: print "make_mass_action_ODE_Terms: NS: ", NS

        k = list(r.Rates())
        if r.Arrow() in ["->","-->"]:
            if len(k)!= 1:
                print "Warning: Only expecting one rate constant "+\
                  "in reaction "+r.Input()+" but found "+str(len(k))+\
                  " - this may cause unexpected results."
        rate = k[0]
        if is_number(rate):
            rate = float(rate)
        else:
            rate = symbolDictionary[rate]
        
        if DBG: print "make_mass_action_ODE_Terms: rate: ", rate
  
        for s in NS:
            if DBG: print "make_mass_action_ODE_Terms: checking s:", s
            if s not in frozen:
                STUFF = baseterm
                sto = NS[s]
                if is_minusone(sto): 
                    STUFF = - STUFF
                elif not is_one(sto):
                    STUFF = sto*STUFF
                if is_number(rate):
                     if float(rate) != 1: STUFF = STUFF*rate
                else:
                    STUFF = STUFF * rate               
                if s in ODETerms:
                    ODETerms[s] += STUFF
                else:
                    ODETerms[s] = STUFF
                    
        return ODETerms
#*******************************************************************************

def make_hill_ODE_Terms(r, ODETerms, symbolDictionary, frozen):
    def setval(val):
        if is_number(val): return(float(val))
        if not(val in symbolDictionary):
            print "Error: undefined symbol: ", val
        return(symbolDictionary[val])
    
    # print "make_hill_ODE_Terms(r, ODETerms, symbolDictionary, frozen)"
    X = listify(r.LHS())
    substrates = X
    X = map(setval,X)
    P = listify(r.RHS())
    P = P[0]
    if P in frozen: return ODETerms
    E = list(r.MHS())
    nonregulatory=False
    if len(E)==1:
        E = E[0]
        nonregulatory=True
    else:   
        E = 1
    rates = list(r.Rates())
    v,n,K,alpha,T=1,1,1,0,1
    nrates=len(rates)
    if nrates>0: v=setval(rates[0])
    if nrates>1: n=setval(rates[1])
    if nrates>2: K=setval(rates[2])
    if nrates>3: alpha=setval(rates[3])
    if nrates>4: T=rates[4]   
    nx=len(X)
    if not(type(T)==type([])): T=[T]
    while len(T) < len(X): T.append(1)
    T = map(setval, T)
    # print "T=",T
    # T = map(lambda x: round(float(x),7), T)
    X=[t*x for (t,x) in zip(T,X)]
    # print "X=", X
    sum_of_x=alpha
    for u in X: 
        sum_of_x = sum_of_x + u    
    # print "sum_of_x=",sum_of_x
    dPdt = (v*setval(E)*(sum_of_x**n))/(K**n + sum_of_x**n)
    # print dPdt   
  
    if P in ODETerms:
        ODETerms[P] += dPdt
    else:
        ODETerms[P] = dPdt
        
    if nonregulatory:
        for S in substrates:
            if S in ODETerms:
                ODETerms[S] -= dPdt
            else:
                ODETerms[S] = -dPdt
    
    return ODETerms

#*******************************************************************************

def make_rational_ODE_Terms(r, ODETerms, symbolDictionary, frozen):
    
    DBG = False
    def setval(val):
        # if type(val)==type([]):
        if isinstance(val, list):
            return(prodlist(map(setval, val)))
        if is_number(val): return(float(val))
        if val in symbolDictionary:
            q = symbolDictionary[val]
            return(q)
        else:
            print "Error: undefined symbol: ", val

    # print symbolDictionary
    

    if DBG: print "make_rational_ODE_Terms: r.Input: ", r.Input()
    if DBG: print "LHs: ", r.LHS(), "len(r.LHS()):", len(r.LHS())
    N,D=r.LHS()
    if DBG: print "N:", N,"len(N):", len(N)
    #N,D = listify(r.LHS())
    N = eval(N)
    D = eval(D)
    if DBG: print "make_rational_ODE_Terms: N:", N, "len(N):", len(N)
    
    N = map(setval, N); ''
    if DBG: print "make_rational_ODE_Terms: setval:N:", N
    
    if DBG: print "make_rational_ODE_Terms: D:", D

    D = map(setval, D)
    if DBG: print "make_rartional_ODE_Terms: D:", D

    P = listify(r.RHS())[0]
    if P in frozen: return ODETerms
    numcoef, denomcoef, numexps, denomexps = list(r.Rates())
    numcoef = map(setval, numcoef)
    denomcoef = map(setval, denomcoef)
    numexps = map(setval, numexps)
    denomexps = map(setval, denomexps)
    # print "Numerator=",N, " Denominator=",D, " Product=",P
    # print "numcoef = ", numcoef
    # print "denomcoef = ", denomcoef
    # print "numexps = ", numexps
    # print "denomexps = ", denomexps
    
    N = [a**b for (a,b) in zip(N, numexps[1:])] ; 
    D = [a**b for (a,b) in zip(D, denomexps[1:])]
    
    
    
    nbase = numcoef[0]
    nexp = numexps[0]
    if not(nexp == 1): nbase = nbase**nexp
    dbase = denomcoef[0]
    dexp = denomexps[0]
    if not(dexp == 1): dbase = dbase**dexp    
    
    numerator = nbase + sumlist([a*b for (a,b) in zip(numcoef[1:],N)])
    denominator = dbase + sumlist([a*b for (a,b) in zip(denomcoef[1:], D)])
    # print "numerator:", numerator, " denominator: ", denominator
    
    dPdt = numerator/denominator
    # print "dPdt = ", dPdt
    
    if P in ODETerms:
        ODETerms[P] += dPdt
    else:
        ODETerms[P] = dPdt
    
    return ODETerms

#*******************************************************************************

def make_MMH_ODE_Terms(r, ODETerms, symbolDictionary, frozen):
    def setval(val):
        if is_number(val): return(float(val))
        if val in symbolDictionary:
            q = symbolDictionary[val]
            return(q)
        else:
            print "Error: undefined symbol: ", val
    
    S = listify(r.LHS())[0]
    P = listify(r.RHS())[0]
    if P in frozen: return ODETerms
    E = listify(r.MHS())
    if len(E)>0:
        E = setval(E[0])
    else:
        E = 1.0
        
    
    # print S, ":-->", P, "mod[",E,"]"

    rates = list(r.Rates())
    if len(rates)==2:
        K,v = map(setval, list(rates))
    elif len(rates)==3:
        k1,k2,k3 = map(setval, list(rates))        
        K = (k2+k3)/k1
        v = k3
    else:
        print len(rates), " rates found in reaction ",r.Input(),\
            " precisely 2 or 3 rates must be specified in an MMH reaction; ",\
            " unable to parse rate constants. Reaction ignored."
        return ODETerms
    # print "rates: ", K, v
    
    dPdt = v*E*setval(S)/(K+setval(S))
    if not (P in frozen):
        if P in ODETerms:
            ODETerms[P] += dPdt
        else:
            ODETerms[P] = dPdt
    if not (S in frozen):
        if S in ODETerms:
            ODETerms[S] -= dPdt
        else:
            ODETerms[S] = -dPdt


    return ODETerms

#*******************************************************************************
    
def make_GRN_ODE_Terms(r, ODETerms, symbolDictionary, frozen):
    def setval(val):
        if is_number(val): return(float(val))
        if val in symbolDictionary:
            q = symbolDictionary[val]
            return(q)
        else:
            print "Error: undefined symbol: ", val
    
    S = listify(r.LHS())
    P = listify(r.RHS())[0]
    E = listify(r.MHS())
    if len(E)>0:
        E = setval(E[0])
    else:
        E = 1.0
        
    v,T,n,h=1,1,1,0
    rates = list(r.Rates())
    nrates=len(rates)
    if nrates>0: v=setval(rates[0])
    if nrates>1: T=rates[1]
    if nrates>2: n=setval(rates[2])
    if nrates>3: h=setval(rates[3])
    nS=len(S)
    if not(type(T)==type([])): T=[T]
    while len(T) < nS: T.append(1)
    T = map(setval, T)
    # print "v,T,n,h=",[v,T,n,h]
    X=-h-sumlist([t*x for (t,x) in zip(T,map(setval,S))])
    # print "X=",X
    dPdt = v/(1+exp(X))
    if P in ODETerms:
        ODETerms[P] += dPdt
    else:
        ODETerms[P] = dPdt

    
    return ODETerms

#*******************************************************************************

def make_NHCA_ODE_Terms(r, ODETerms, symbolDictionary, frozen):
    def setval(val):
        if is_number(val): return(float(val))
        if val in symbolDictionary:
            q = symbolDictionary[val]
            return(q)
        else:
            print "Error: undefined symbol: ", val
    
    S = listify(r.LHS())
    P = listify(r.RHS())[0]
    E = listify(r.MHS())
    
    if len(E)>0:
        E = setval(E[0])
    else:
        E = 1.0
        
    v,TP,TM,n,m,k=1,1,1,1,1,1
    rates = list(r.Rates())
    nrates=len(rates)
    if nrates>0: v=setval(rates[0])
    if nrates>1: TP=rates[1]
    if nrates>2: TM=rates[2]
    if nrates>3: n=rates[3]
    if nrates>4: m=setval(rates[4])
    if nrates>5: k=setval(rates[5])
    
    nS=len(S)
    if not(type(TP)==type([])): TP=[TP]
    TP1=TP[-1]
    while len(TP) < nS: TP.append(TP1)
    TP = map(setval, TP)

    if not(type(TM)==type([])): TM=[TM]
    TM1=TM[-1]
    while len(TM) < nS: TM.append(TM1)
    TM = map(setval, TM)
    
    if not(type(n)==type([])): n=[n]
    n1=n[-1]
    while len(n) < nS: n.append(n1)
    
    
    # print "v, VP, TM, n, m, k=", v, TP, TM, n, m, k
    U = [x**y for (x,y) in zip(map(setval,S), map(setval,n))]
    # print "U=",U
       
    xplus=  prodlist([(1+t*x)**m for (t,x) in zip(TP, U)])
    xminus= prodlist([(1+t*x)**m for (t,x) in zip(TM, U)])
    # print "XP,XM=",xplus, xminus
    
    dPdt = v*xplus/(k*xminus + xplus)
    
    if P in ODETerms:
        ODETerms[P] += dPdt
    else:
        ODETerms[P] = dPdt

    
    return ODETerms

#*******************************************************************************

def make_MWC_ODE_Terms(r, ODETerms, symbolDictionary, frozen):
             
    def setval(val):
        if is_number(val): return(float(val))
        if val in symbolDictionary:
            q = symbolDictionary[val]
            return(q)
        else:
            print "Error:undefined symbol: ", val
    
    
    S = list(r.LHS())
    P = list(r.RHS())[0]
    modifiers = map(listify, list(r.MHS()))
    
    
    E,A,I = modifiers[:3]
    # print "make_MWC_ODE_Terms: E=", E, " A=" ,A, " I=", I    
    if len(E)>0:
        E = setval(E[0])
    else:
        E = 1.0
        
    # print "make_MWC_ODE_Terms:  E=", E, " A=" ,A, " I=", I    
    rates = list(r.Rates())
    k,n,c,L,K,KA,KI=1,1,1,1,1,1,1
    if len(rates)>0: k=setval(rates[0])
    if len(rates)>1: n=setval(rates[1])
    if len(rates)>2: c=setval(rates[2])
    if len(rates)>3: L=setval(rates[3])
    if len(rates)>4: K=rates[4]
    if len(rates)>5: KA=rates[5]
    if len(rates)>6: KI=rates[6]

    if type(K)==type("K"): K=[K]
    if type(KA)==type("KA"): KA=[KA]
    if type(KI)==type("KI"): KI=[KI]
    
    # print "make_MWC_ODE_Terms:  rates=",rates
    nS = len(S)
    K1 = getlast(K, 1)
    # print "make_MWC_ODE_Terms:  K=",K, " K1=",K1
    while len(K) < nS: K.append(K1)
    nA=len(A)
    nI=len(I)
    KA1 = getlast(KA, 1)
    KI1 = getlast(KI, 1)
    while len(KA) < nA: KA.append(KA1)
    KA = map(setval, KA)
    while len(KI) < nI: KI.append(KI1)
    KI = map(setval, KI)
    A = map(setval, A)
    I = map(setval, I)
    
    inhibitors = modifiers[3:]
    # print "make_MWC_ODE_Terms: inhibitors=", inhibitors
    rates = rates[7:]
    # print "make_MWC_ODE_Terms: rates=",rates
    
    
    SBAR=[]
    ABAR=[]
    for i in range(nS):
        # print "make_MWC_ODE_Terms: Substrate=",S[i]
        CS=map(setval, inhibitors[i])
        CSK=map(setval, rates[i])
        sqbar=sumlist([c*x/y for (x,y) in zip(CS,CSK)])
        # print "make_MWC_ODE_Terms: Competitors=",CS, " rates=",CSK," sqbar=", sqbar
        SBAR.append(sqbar)
    for i in range(nI):
        # print "make_MWC_ODE_Terms: Inhibitors=",I[i]
        CI=map(setval, inhibitors[i+nS])
        CIK=map(setval, rates[i+nS])
        aqbar=sumlist([c*x/y for (x,y) in zip(CI,CIK)])
        # print "make_MWC_ODE_Terms: Competitors=",CI, " rates=", CIK, " aqbar=", aqbar
        ABAR.append(aqbar)

    # print "S=", S, " P=",P," E=",E, " A=",A, " I=",I," K=",K, " KA=",KA, " KI=", KI
        
    s = [s/k for (s,k) in zip(map(setval,S),map(setval,K))]
    a = [a/k for (a,k) in zip(A, KA)]
    i = [i/k for (i,k) in zip(I, KI)]
    
    # print "s=",s , " a=", a, " i=",i
    
    S_plus_SBAR = [x+y for (x,y) in zip(s,SBAR)]
    A_plus_ABAR = [x+y for (x,y) in zip(a,ABAR)]
    
    afactor = prodlist([(1+aval) for aval in A_plus_ABAR])
    ifactor = prodlist([(1+ival) for ival in i])
    sfactor = prodlist([(1+sval) for sval in s])
    ssbarfactor = prodlist([(1+sval) for sval in S_plus_SBAR])
    csfactor = prodlist([(1+c*sval) for sval in S_plus_SBAR])
    sprod = prodlist(s)
    csprod = prodlist([c*sval for sval in s])
    
    numerator = (afactor**n)*(ssbarfactor**(n-1))*(sprod) + \
        L * csprod * (csfactor ** (n-1)) * (ifactor**n)
        
    denominator = (afactor**n)*(sfactor**n) + L * (csfactor**(n-1)) * (ifactor**n)
    
    dPdt=E * numerator/denominator
    
    # print "afactor=",afactor
    # print "ifactor=",ifactor
    # print "sfactor=",sfactor
    # print "csfactor=", csfactor
    # print "sprod=",sprod
    # print "csprod=", csprod
    # print "numerator=", numerator
    # print "denominator=",denominator
    # print "dPdt=",dPdt    
    
    if not (P in frozen):
        if P in ODETerms:
            ODETerms[P] += dPdt
        else:
            ODETerms[P] = dPdt
        
    for U in S:
        if U in frozen: 
            continue
        elif U in ODETerms:
            ODETerms[U] -= dPdt
        else:
            ODETerms[U] = -dPdt
            

        
    return ODETerms
 

#*******************************************************************************

def make_SSYSTEM_ODE_Terms(r, ODETerms, symbolDictionary, frozen,eps=10**(-37)):
    

    if "-epsilon" in sys.argv:
        i = (sys.argv).index("-epsilon")+1
        if i<len(sys.argv):
            try:
                eps = sys.argv[i]
            except:
                print "Error: SSystem: expecting a numerical value after -epsilon option."
                raise SystemExit("Invalid command syntax.")
         
    def setval(val):
        if is_number(val): return(float(val))
        if val in symbolDictionary:
            q = symbolDictionary[val]
            return(q)
        else:
            print "Error:undefined symbol: ", val
            
    # print "SSystem: input: ",r.Input()
    
    
    S = map(setval, listify(r.LHS()))   
    P = listify(r.RHS())[0]
    if P in frozen: return ODETerms
    
    rates = listify(r.Rates())
       
    # print "SSystem: S=",S, " P=",P
    # print "SSystem: Rates=",rates
    ns = len(S)
    
    tau, alpha, beta,nplus,nminus = 1, 1, 1,[],[]
    if len(rates)>0: tau = setval(rates[0])
    if tau == 0:
        print "Error: Interpreter: SSystem time constant cannot be zero."
        print "Input reaction: ",r.Input()
        raise SystemExit("Unable to interpret reaction: Please correct reaction and resubmit job.")
    if len(rates)>1: alpha = setval(rates[1])
    if len(rates)>2: beta = setval(rates[2])
    if len(rates)>3: nplus = map(setval, listify(rates[3]))
    if len(rates)>4: nminus = map(setval, listify(rates[4]))
    
    
    if len(nplus)>ns: 
        nplus=nplus[:ns]
    else:
        np1=0
        if len(nplus)>0: np1=nplus[0]
        while len(nplus)<ns: nplus.append(np1)
        
    if len(nminus)>ns: 
        nminus=nminus[:ns]
    else:
        nm1=0
        if len(nminus)>0: nm1=nminus[0]
        while len(nminus)<ns: nminus.append(nm1)
    
    
    # print "SSystem: tau=",tau," alpha=", alpha, " beta = ", beta
    # print "SSystem: nplus = ", nplus, " nminus = ",nminus
    
    
     
    SP = prodlist([max(x,eps)**y for (x,y) in zip(S, nplus)])
    SM = prodlist([max(x,eps)**y for (x,y) in zip(S, nminus)])
 

    # print "SSystem: SP = ", SP
    # print "SSystem: SM = ", SM
    dPdt = (alpha * SP - beta * SM)/tau
    # print "SSystem: dPdt = ", dPdt
    
    if P in ODETerms:
        ODETerms[P] += dPdt
    else:
        ODETerms[P] = dPdt

    
    return ODETerms
    
    
#*******************************************************************************
def make_using_ODE_Terms(r, ODETerms, symbolDictionary, frozen):
    S = listify(r.LHS())
    SS = listify(r.LH_STOIC())
    P = listify(r.RHS())
    PS = listify(r.RH_STOIC())
    
    rate = r.rates[0]
    species = list(set(S+P))
    
    for X in species: 
        # 
        # Note that "frozen" variables may be changed by a "using" expression
        #
        SL = 0
        if X in S: SL = int(SS[S.index(X)])
        SR = 0
        if X in P: SR = int(PS[P.index(X)])
        dXdt = str(SR-SL)+"*"+rate[1:-1]
        dXdt = sympify(dXdt)
        
        if X in ODETerms:
            ODETerms[X] += dXdt
        else:
            ODETerms[X] = dXdt

    return ODETerms    
#*******************************************************************************
    
def make_USER_ODE_Terms(r, ODETerms, symbolDictionary, frozen):
    def setval(val):
        if is_number(val): return(float(val))
        if val in symbolDictionary:
            q = symbolDictionary[val]
            return(q)
        else:
            print "Error: undefined symbol: ", val
    
    S = listify(r.LHS())
    P = listify(r.RHS())[0]
    E = listify(r.MHS())
    if len(E)>0:
        if E[0]=="Nil":
            E = 1
        else:        
            E = setval(E[0])
    else:
        E = 1.0
        
    # print "make_USER_ODE_Terms: S = ", S, " P = ", P, " E = ", E
   
        
    rates = list(r.Rates())
    nrates=len(rates)
    if not (nrates == 5):
        print "Error: USER Reaction must have precisely five rates representing "\
        + "v, T, n, h and a lambda-expression. The reaction `"+r.Input()+"' has "\
        + str(nrates)+ " rates."
        raise SystemExit("Correct reaction and re-submit job.")
    v, T, n, h, f = rates
    # print "make_USER_ODE_Terms: v=",v," T=",T," n=", n," h=",h, " f=", f
    
    #
    # this section blocked off because functions should be referenced through 
    # the function section now 7-29-16
    #
    if false:
		f = f[1:-1].strip() # drop surrounding quotes
		fparts = f.split(":")
		if len(fparts) != 2:
			print "Error: Incorrectly formated lambda expression in User Reaction '"+r.Input()+"'"
			raise SystemExit("Correct reaction and re-submit job.")
		fargs,fdef=fparts
		fargs = fargs.strip()
		fdef = fdef.strip()
		# print "make_User_ODE_Terms: fargs='"+fargs+"' fdef='"+fdef+"'"
		if not fargs.startswith("lambda"):
			print "Error: Expecting a lambda expression as the final rate in User " \
			+ "Reaction '"+r.Input()+"'"
			raise SystemExit("Correct reaction and re-submit job.")
		fargs = flatten(map(lambda x: x.split(","), fargs.split(" ")))
		fargs = filter(lambda x: len(x)>0, map(lambda x:x.strip(), fargs[1:]))         
		if len(fargs)!= 1:
				print "Error: The lambda expression in the User reaction '"+r.Input()+\
				"' should have precisely one argument; this one has "+str(len(fargs)-1)
				raise SystemExit("Correct reaction and re-submit job.")
    

    # print "returning from make_USER_ODE_Terms"

    v = setval(v)
    h = setval(h)
        
    nS=len(S)
    if not(type(T)==type([])): T=[T]
    T0 = 1
    if len(T)>0: T0=T[-1]
    while len(T) < nS: T.append(T0)
    T = map(setval, T)
    
    if not(type(n)==type([])): n=[n]
    n0=0
    if len(n)>0: n0=n[-1]
    while len(n) < nS: n.append(n0)
    n = map(setval, n)
 
    # print "make_USER_ODE_Terms: v=",v," T=",T," n=", n," h=",h, " f=", f

    s_to_the_ns = [s**q for (s,q) in zip(map(setval,S), n)]
    X = -h - sumlist([t*x for (t,x) in zip(T,s_to_the_ns)])
 
    # X=-h-sumlist([t*x for (t,x) in zip(T,map(setval,S))])
    
    # print "make_USER_ODE_Terms: X = ", X
    
    f = sympify(f)
    # print "make_USER_ODE_Terms: f = ", f
    
    dPdt = v*E*f(X)
    # print "make_user_ODE_Terms: dPdt=",dPdt

 

    if P in ODETerms:
        ODETerms[P] += dPdt
    else:
        ODETerms[P] = dPdt

    
    return ODETerms    

#*******************************************************************************

ODEHandler={
            "Mass Action":make_mass_action_ODE_Terms,
            "Hill":make_hill_ODE_Terms,
            "rational":make_rational_ODE_Terms,
            "MMH": make_MMH_ODE_Terms, 
            "GRN": make_GRN_ODE_Terms,
            "NHCA": make_NHCA_ODE_Terms,
            "MWC":  make_MWC_ODE_Terms,
            "SSystem": make_SSYSTEM_ODE_Terms, 
            "USER": make_USER_ODE_Terms,
            "using": make_using_ODE_Terms
            }
#*******************************************************************************


def makeODETerms(reactions, symbolDictionary, frozen=[]):
    """
    Convert list of reactions to a dictionary of terms in Differential
    Equations.
    
    The expander should be called first, to convert higher order mass-action
    reactions to basic mass-action reactions

    So far, this only works on basic, low-level, mass-action reactions

    input:  reactions - list of reactions (output of ParseArrow)
            s - a symbolDictionary {"symbol":symbol, ...} where each
                symbol is the Sympy symbol used for the string "symbol"
            frozen - optional list of frozen variables (same as in Cellerator)
                for which the ODEs will be d/dt = 0.

    Output: dictionary of ode terms {"symbol": term
        where term is the symbolic representation of the term in Sympy
        notation. Similar to Basic Input Form in Mathematica. This gives the RHS
        of the ODE for each variable. 
        
    """
    
    DBG = False
    
    ODETerms = {}

    for s in frozen:
        ODETerms[s] = 0
        

    for r in reactions:
        if DBG: print "makeODETerms: reaction is ", r.Input()
        at = r.ArrowType()
        if at in ODEHandler:
            odemaker=ODEHandler[at]
            ODETerms = odemaker(r, ODETerms, symbolDictionary, frozen)
        else:
            print "interpreter not implemented for reaction of type `"+ str(r.ArrowType())\
            +"' input reaction: ", r.Input()
            
    if DBG: print "makeODETerms: removing NIL"
    if len(ODETerms)>0:
        if "Nil" in ODETerms:
            del ODETerms["Nil"]
            # print "removing Nil"
              
    return ODETerms

#*******************************************************************************

def jacobian(odeterms):
    J=[]
    for y in odeterms:
        JY=[]
        for x in odeterms:
            XI=diff(odeterms[y],Symbol(x))
            JY.append(XI)
        J.append(JY)
    output = Matrix(J)
    return output

#*******************************************************************************

def steadystate(odeterms, rates, symbolic):
    J = jacobian(odeterms)
    variables = list(odeterms)
    symbols = map(Symbol,variables)
    equations = [odeterms[v] for v in variables]
    

    # print "rates = ", rates
    # print "equations: ", equations
    
    J1 = J
    for v in rates:
        symbol=Symbol(v)
        value = rates[v]
        J1=J1.subs(symbol,value)
        equations = map(lambda x: x.subs(symbol,value), equations)
        
    
    
    # print "equations: ", equations
    # print "symbols: ", symbols
    
    ss = solve(equations, symbols)
    
    # print "ss = ", ss   
    # DJ = J1.det()
    
    # print "J1 = ",J1, " determinant: ", DJ
    ev = J1.eigenvals()
    # print "eigenvalues: ", ev
    
    eigenvalues = [u for u in ev]
    multiplicity = []
    for m in ev: multiplicity.append(ev[m])
    
    # print "eigenvalues = ", eigenvalues
    
    evalsss=[]
    for SS in ss:
        # print "examining steady state at ", SS
        evss = eigenvalues
        for i in range(len(symbols)):
            name = symbols[i]
            ssvalue = SS[i]
            # print name, " = ", ssvalue
            evss = map(lambda x:x.subs(name,ssvalue),evss)
        # print "evss = ",evss
        evalsss.append(evss)
            

        
    # print "eigenvalues: ", eigenvalues
    
    if symbolic:
        return (symbols, ss, evalsss)
    else:
        output=[]
        output.append("There are "+str(len(ss))+" steady states")
        for i in range(len(ss)):
            output.append(str(tuple(symbols))+" = "+str(ss[i])+" has eigenvalues "+str(evalsss[i]))
        # output.append("The Jacobian is \n"+str(J))
        # output.append("The Jacobian with constants is \n"+str(J1))
        # output.append("The general formula for the eigenvalues is:\n"+str(eigenvalues))
        output = "\n".join(output)
        
        return (output)
        
    
#*******************************************************************************
    
def interpret(filename, ofile="", dumpparser=False, format="ODE", frozen=[], \
    symbolic=False):
    
    if "-time" in sys.argv: 
        start = clock()
    
    DBG = False
    if DBG: print "interpret: Starting"
    
    output=[]
    dump = dumpparser or ("-dump" in sys.argv)
    SYMBOLIC = symbolic or "-symbolic" in sys.argv
        
    freeze=[]
    if ("-frozen" in sys.argv):
         i=(sys.argv).index("-frozen")+1
         while i<len(sys.argv):
            next = sys.argv[i]
            if str(next)[0]=="-": break
            freeze.append(next)
            i+=1
    elif type(frozen)==type([]):
        freeze = frozen
    if DBG: print "interpret: calling Parser"
       
    (results, rates, functions, frozen, assignments, ICS) = runparser(filename, trace=dump)
    # 
    # combine frozen variables on command line with those in model file
    #
    freeze = list (set (freeze + frozen))
    # print rates
    if DBG: print "interpret: expanding"   
    if "-time" in sys.argv: 
        parsetime = clock()-start
        print "(parse time): %20.10f" %(parsetime)
        start = clock()
    results = expand(results)
    if DBG: print "interpret: making Symbol Dicitionary"
    s=makeSymbolDictionary(results, rates, ICS)
    if DBG: print "interpret: making ODE Terms"    
    odeterms = makeODETerms(results, s, frozen=freeze)
    if DBG: print "interpret: printing terms"
    
    # print odeterms
    if "-format" in sys.argv:
        i = (sys.argv).index("-format")+1
        if i < len(sys.argv):
            format = sys.argv[i]
        else:
            print "Error: interpret: expecting format string after -format keyword."
            format = "ODE"
    fmt = str(format).upper()

    if not fmt in ["ODE","DICT","CODE", "PYTHON", "LATEX", "JACOBIAN", "STEADYSTATE"]:
        print "Error: interpret: invalid format='"+str(fmt)+"' requested."
        raise SystemExit("Resubmit job with correct value for keyword -format")
    k=0
    if fmt in ["CODE","PYTHON"]:
        output.append("from math import *")
        if len(functions)>0:
            for f in functions:
                output.append(f)
        output.append("def f(y,t):")
        for x in rates:
            output.append("    "+str(x)+" = "+str(rates[x]))
        for x in assignments:
            output.append("    "+str(x))
        for x in odeterms:
            line = "    y["+str(k)+"] = " + str(x)
            output.append(line)
            k += 1
        k = 0
    elif fmt == "DICT":
        output=odeterms
        return (output)
        # if symbolic: return output
    elif fmt == "JACOBIAN":
        output = jacobian(odeterms)
        return (output)
        # if symbolic: return output
    elif fmt == "STEADYSTATE":
        output = steadystate(odeterms, rates, SYMBOLIC)

      
    elif fmt == "LATEX":
        output.append("\\documentclass[12pt,letterpaper]{article}")
        output.append("\\usepackage[latin1]{inputenc}")
        output.append("\\usepackage{amsmath, amsfonts, amssymb}")
        output.append("\\author{py[cellerator]}")
        output.append("\\date{\\today}")
        output.append("\\title{Automatically Generated Equations}")
        output.append("\\begin{document}")
        output.append("\\maketitle")
        output.append("\\begin{align*}")
        for f in functions:
            output.append(latex(f))
    if fmt == "ODE":
        for f in functions: 
            output.append(f)
        for f in assignments:
            output.append(f)
    for x in odeterms:
        if fmt == "ODE":
            line =  str(x)+"' = "+str(odeterms[x])
        elif fmt == "DICT":
            continue
        elif fmt == "JACOBIAN":
            continue
        elif fmt == "STEADYSTATE":
            continue
        elif fmt in ["CODE", "PYTHON"]:
            line = "    yp["+str(k)+"] = " + str(odeterms[x])
            k+=1
        elif fmt == "LATEX":
            line =x+"' &= "+latex(odeterms[x])[1:-1]+"\\\\"
        else:
            print "Error: interpret: invalid format requested: " + str(format)
            raise SystemExit("Resubmit job with correct value for keyword.")
        output.append(line)

    if fmt =="ODE":
        output = "\n".join(output)
    elif fmt == "LATEX":
        output.append("\\end{align*}")
        output.append("\\end{document}")
        output = "\n".join(output)    
    elif fmt in ["CODE", "PYTHON"]: 
        output.append("    return (yp)")
        output = "\n".join(output)
    
    
    outfile=ofile    
    if ("-out" in sys.argv):
        i = (sys.argv).index("-out")+1
        if i<len(sys.argv):
            outfile=sys.argv[i]
        else:
            print("Error: interpret: expecting filename after -out option.")
            outfile="interpreted-equations.out"

    if len(outfile)>0:
        outputfile = outfile
        (left,right)=outfile.split(".")
        j=1
        while os.path.isfile(outputfile):
            outputfile=left+str(j)+"."+right
            j+=1
        f=open(outputfile,"w")
        if fmt == "DICT":
            output = str(output)
            f.write(output)
        else:
            f.write(output)
            # output = map(lambda x:x+"\n",output)
            # f.writelines(output)
        f.close()
        if "-time" in sys.argv: 
            print "(interpret time): %20.10f" %(clock() - start)

        
        return os.path.abspath(outputfile)
    
    if "-time" in sys.argv: 
        print "(interpret time): %20.10f"  %(clock() - start)
       
    return(output)

#*******************************************************************************    

def interpreter():
    if not ("-quiet" in sys.argv):
        d=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print "py[cellerator]: interpreter ("+d+")"
    rv=interpret("", ofile="")
    if "-quiet" in sys.argv: return
    print rv
#*******************************************************************************    


if __name__ == "__main__":
    interpreter()



