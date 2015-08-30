#!/usr/bin/env python
#
# Cellerator Shorthand Parsing Module
# B.E.Shapiro 3 April 2011 
#
# 11/8/14 add diffusion arrow
# rev 8/30/15
#
#
#
#****************************************************************************
#    Cellerator converts reactions, expressed in a text-formatted arrow-based 
#    notation into differential equations and performs numerical simulations.
#
#****************************************************************************
#
#    Copyright (C) 2012-14 Bruce E Shapiro
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

from pyparsing import *
import os.path
import sys, datetime
from utils import *
from reader import readmodel
from sympy import *


dummy=Word(alphas)
parseType=type(dummy.parseString("teststring"))
 
def isnumber(s):
    """
    Checks a string to determine if it encodes a number. 
    """
    try:
        float(s)
        return True
    except ValueError:
        return False


#**********************************************************************

class Reaction(object):
    """
    Object to store a cellerator arrow.

    To invoke: 
    x=ParseArrow(arrow)
    arrow = string representation
    
    methods:
    x.RHS() - gives the RHS of the reaction (products)
    X.RH_STOIC() - gives the stoichiometry of RHS (for Mass Action only; 
        may give unexpected values for other reactions)
    x.LHS() - gives the LHS of the reaction (reactants)
    X.LH_STOIC() - gives the stoichiometry of the LHS (for Mass Action only; 
        may give unexpected values for other eactions)
    x.Arrow() - returns the arrow, e.g., "->", "-->", etc.
    x.Rates() - returns a list of rate constants
    x.MHS() - returns the list of modifiers
    x.Number() - a unique reaction counter is incremented for each objecte
    x.Input() - returns the original input string
    x.Parsed() - returns the original parsed object
    x.Dump() - prints a dump of the object
    """
    counter = 0
    
 
    def __init__(self, raw):

        def _needsModifier(x):
            needsmod={"->":False,  "-->":True, "<->":False,
                      "=>":False,  "<=>":True, "|->":False,
                      "|-->":True, "==>":True, ":=>":True,
                      ":->":False, ":-->":True }
            if x in needsmod:
                return needsmod[x]
            else:
                return False

        def _CheckRates(rates,arrow):
            """
            Checks to make sure that the right number of rates are provided.

            If too few rates are given, defaults are filled in.
            If too many are given, the extras are chopped off.
            """
            # print "DBG: CheckRates: ", rates, arrow
            arrow2check=arrow
            if len(rates)>0:
                if rates[0]=='using':
                    r = rates[1:]
                    r = ["using", "".join(r[0])]
                    return(r)
                if rates[0]=='Flux':
                    arrow2check='Flux'
                if rates[0]=='Enz':
                    return(rates)
                if rates[0]=='Diffusion':
                    return(rates)
 
            # print "DBG: CheckRates: ", rates, arrow2check
          
          
            def fluxrate():
                number = self.num
                v="v"+str(self.num)
                return (v)
                               
            DEFAULT_RATES={
                "->":  [1],
                "-->": [1],
                "<->": [1,1],
                "=>":  [1,1,1,0],
                "<=>": [1,1,1,0,1,1,1,0],
                ":=>": [1,1,1,0,1,0],
                "Flux": ["Flux", -float("inf"), fluxrate(), float("inf"),1, 0]
                }
            OTHER_DEFAULT_RATES={ 
                ":->": [1,1],
                ":-->":[1,1]
                }
            N_MAX_RATES={
                ":->":5,
                ":-->":5 }
            Descriptions={
                "->":  "Simple Mass Action",
                "-->": "Simple Mass Mass Action with a Modifier",
                "<->": "Bi-directional Simple Mass Action",
                "=>":  "Catalyzed Mass Action with Intermediate Complex Formation",
                "<=>": "Reversible Catalyzed Mass Actiona with Intermediate Complex Formation",
                ":=>": "Catalyzed Mass Action with Sequential Complex Formation",
                ":->": "Michaelis-Menten-Henri Rate Law",
                ":-->":"Michaelis-Menten-Henri Rate Law with explicit Enzyme dependence"
                }
            r=rates[:]
            nrates=len(r)

            def _copy_rates(default):
                nexpected=len(default)
                while len(r) < nexpected:
                    nextrate = default[len(r)]
                    r.append(str(default[len(r)]))
                return nexpected

            def _chop_rates(r, nexpected):
                if len(r)>nexpected:
                   r = r[:nexpected]
                   print "Warning: too many rates specified in the reaction ",\
                         raw
                return r       
            
            # print "DBG: CheckRates: rates: ", rates
            # print "DBG: CheckRates: arrow: ", arrow2check
          
            if arrow2check in DEFAULT_RATES:
                default = DEFAULT_RATES[arrow2check]
                nexpected=_copy_rates(default)
                r = _chop_rates(r, nexpected)
                return r[:]
            elif arrow2check in OTHER_DEFAULT_RATES:
                default=OTHER_DEFAULT_RATES[arrow2check]
                _copy_rates(default)
                if arrow2check in N_MAX_RATES:
                    nexpected = N_MAX_RATES[arrow2check]
                    r=_chop_rates(r, nexpected)
                return r[:]
            else:
                # make sure lists get turned into lists
                def f(x):
                    if type(x)==type("hello world"): return x
                    return list(x)
                r = map(f,r)
                # pass
                #
                # NEED to put code here to handle non-mass action rates
                #
            return r

       # END OF CHECKRATES
       ############################################################################

        NUM_MODS_NEEDED={
            "->":  0,
            "-->": 1,
            "<->": 0,
            "=>" : 1,
            "<=>": 2,
            ":->": 0,
            ":-->":1,
            ":=>": 1
            }

        def _CheckMods(mods, arrow, arowtype):
            if arrow == "<=>":
                # print (arrow,arrowtype)
                if arrowtype == "kMech":
                    return mods
            if arrow in NUM_MODS_NEEDED:
                need = NUM_MODS_NEEDED[arrow]
                have = len(mods)
                m = mods[:]
                if have != need:
                    print "Error: Expecting ", need," modifiers; found ", \
                          have, " in reaction ", raw
                    if have > need:
                        m = m[:need]
                return m
            else: 
                return mods
        #*************************************************************************
                        
        def _DetermineType(rates, arrow):
            Types = {
                "->":  "Mass Action",   # mass action, Flux, or Diffusion using arrowtypes
                "-->": "Mass Action",
                "<->": "Mass Action",
                "=>" : "Mass Action",
                "<=>": "Enzymatic",     # mass action or kMech arrowtype
                ":->": "MMH",
                ":-->":"MMH",
                ":=>": "Mass Action",
                "|->": "GRN",           # multiple subtypes given by arrowtype
                "|-->":"GRN",           # multpile subtypes given by arrowtype
                "==>": "Multiuse-Arrow" #MWC, rational
                }
            

            # print "entering determineType: rates: ", rates
            
            # arrowtypes that are **NOT** given by "Types"

            if len(rates)>0:
                if rates[0]=='using':
                    arrowtype="using"
                    r=rates[1:]
                    return (r, arrowtype)
                if rates[0]=='Flux':
                    arrowtype="Flux"
                    r=rates[1:]
                    return(r, arrowtype)
                if rates[0]=="Diffusion":
                    arrowtype="Diffusion"
                    
                    r=rates[1:]
          
                    return(r,arrowtype)
            
            # arrowtypes that **ARE** given by "Types"
            
            if arrow in Types:
                arrowtype = Types[arrow]                
            else:
                arrowtype = "Unknown"
            
            # special modifications to types that **ARE** in "TYPES"
            
            if arrowtype == "GRN":
                arrowtype = rates[0]
                r = rates[1:]
            elif arrowtype == "Enzymatic":
                arrowtype="Mass Action"
                r=rates
                if len(rates)>0:
                    if rates[0]=="Enz":
                        arrowtype="kMech"
                        r=rates[1:]
            elif arrowtype == "Multiuse-Arrow":
                if len(rates)>0 and rates[0]=="MWC":
                    arrowtype = "MWC"
                    r = rates[1:]
                else:
                    r = rates              
            else:
                r = rates
            return ( r, arrowtype)
        #*************************************************************************
            
        def _errorCheckUSER():
            if self.arrowtype=="USER":
                if len(self.rates)<5:
                    print "Error: User Arrow "+str(self.raw).rstrip()+" has only "+\
                        str(len(self.rates))+" User fields, 5 expected. Expected "+\
                        "format is User[r, T, n, h, f]"
            return
        #*************************************************************************
            
        def _errorCheckRATIONAL():
        #
        # check and if possible fix rational reactions by adding default
        # exponents and coefficients
        #
        # rational species are always going to be recognized as products
        # even if they are not products
        #
            if self.arrowtype=="rational":
                A,D,N,M=self.rates
                [nA,nD,nN,nM]=map(len,self.rates)
                num,den=self.lhs
                num = list(num) 
                den = list(den)
                numerator=[]; denominator=[]
                for X in num: numerator.append(X[:])
                NUM = [(u if len(u)>1 else u[0]) for u in numerator]
                for X in den: denominator.append(X[:])
                DEN = [(u if len(u)>1 else u[0]) for u in denominator]
                # print "numerator=",numerator, " denominator:",denominator
                # print "NUM=",NUM, " DEN = ", DEN


                S = [NUM,DEN]
                # print "S=",S
                self.lhs = S
                nn = len(num); nd=len(den)
                # check number of coefficients - default is zero
                if nA>nn+1:
                    print "Error: Too many Numerator coefficients in reaction "\
                        + str(self.raw).rstrip()\
                        +"; expecting at most "+str(nn+1)+" but "+str(nA)+" found"
                    A=A[:nn+1]
                else:
                    while len(A)<nn+1:
                        A.append(0)
                if nD>nd+1:
                    print "Error: Too many Denominator coefficients in reaction "\
                        + str(self.raw).rstrip()\
                        +"; expecting at most "+str(nd+1)+" but "+str(nD)+" found"
                    D=D[:nd+1]
                else:
                    while len(D)<nd+1:
                        D.append(0)
                
                # check exponents - default is 1
                if nN>nn+1:
                    print "Error: Too many Numerator exponents in reaction "\
                        + str(self.raw).rstrip() \
                        +"; expecting at most "+str(nn+1)+" but "+str(nN) +" found"
                    N = N[:nn+1]
                else:
                    while len(N)<nn+1:
                        N.append(1)                      
                if nM>nd+1:
                    print "Error: Too many Denominator exponents in reaction "\
                        + str(self.raw).rstrip() \
                        +"; expecting at most "+str(nd+1)+" but "+str(nM) +" found"
                    M = M[:nd]
                else:
                    while len(M)<nd+1:
                        M.append(1)
                self.rates = [A,D,N,M]
            return
   
        #*************************************************************************
             
        def _errorCheckMWC():
        #
        # Expected format of rates 
        # MWC[kcat, n, c, L, KS] or
        # MWC[kcat, n, c, L, KS, KA, KI] 
        # where KA and KI are lists, e.g., 
        # MWC[kcat, n, c, L, KS, [KA1, KA2, KA3,...], [KI1, KI2,...]]
        # too many KA's or KI's is considered and error and a message is printed
        # too few KA's or KI's is considered lazy input and 1-filled
        # if n, c, L, KS, KA's, or KA's are missing, they are replaced with 1
        # the number of KA's and KA's are expanded to equal the number of 
        # Activators and Inhibitors in the reaction
        #
        # If there is competitive inhibition, the modifiers look like:
        # mod[E, [A1,..], [I1,..], [CS11,..],..,[CA11,...]]
        # mod[enzyme, activators, inhibitors, CI for S1, CI for S2, ...,
        #       CI for SN, CI for A1, CI for A2, ..., CI for last A 
        #
            if self.arrowtype == "MWC":
                # parse modifiers into [Enzyme, [Activators], [Inhibitors]] for GMWC
                # even or [Enzyme, [], []] for the standard form of MWC
                modifiers = self.mods
                enzyme = modifiers[0]
                modifiers = modifiers[1:]
                activators = []
                inhibitors = []
                if len(modifiers)>0: activators = listify(modifiers[0])
                if len(modifiers)>1: inhibitors = listify(modifiers[1])
                    
                self.mods = [enzyme, activators, inhibitors]
                substrates = listify(self.lhs)
                
                nS = len(substrates)
                nA = len(activators)
                nI=len(inhibitors)
                
                # each activator and substrate is allowed to have
                # its own separate list of competitive inhibitors
                # there may be empty lists along the way.
                
                CIS = map(listify, modifiers[2:])
                # print "CIS(1)=",CIS
                
                CS=[]
                for i in range(nS):
                    if len(CIS)==0: break
                    CI = CIS[0]
                    CIS = CIS[1:]
                    CS = CS + [CI]
                # print "CIS(2)=", CIS
                while len(CS)<nS: CS.append([])
                 
                CA=[]
                for i in range(nA):
                    if len(CIS)==0: break
                    CI = CIS[0]
                    CIS = CIS[1:]
                    CA = CA + [CI]
                # print "CIS(3)=", CIS
                while len(CA)<nA: CA.append([])
                # print "CS=",CS," CA=",CA
                competitors = CS + CA
                # print "competitors = ", competitors
                self.mods = self.mods + competitors
                                
                
                R = self.rates
                    
                k,n,c,L,K,KA,KI,KCS,KCA=1,1,1,1,1,1,1,1,1

                if len(R)>0: k=(R[0])
                if len(R)>1: n=(R[1])
                if len(R)>2: c=(R[2])
                if len(R)>3: L=(R[3])
                if len(R)>4: K=rates[4]
                if len(R)>5: KA=rates[5]
                if len(R)>6: KI=rates[6]
                
                # make sure the right number of rates are aquired for each species

                K=deatomize(K)
                KA=deatomize(KA)
                KI=deatomize(KI)
      
                # print "K=",K, " KA=", KA," KI=", KI
             
                K1 = getlast(K, 1)
                while len(K) < nS: K.append(K1)
                while len(K) > nS: K = K[:-1]
                
                KA1 = getlast(KA, 1)
                while len(KA) < nA: KA.append(KA1)
                while len(KA) > nA: KA = KA[:-1]
 
                KI1 = getlast(KI, 1)
                while len(KI) < nI: KI.append(KI1)
                while len(KI) > nI: KI = KI[:-1]
                
                # print "K=",K, " KA=", KA," KI=", KI
                
                
                self.rates = R[:4] + [K] + [KA] + [KI]
                
                # acquire the rates for the competitive inhibitors
                
                KC=[]
                if len(rates)>7: KC=map(deatomize, rates[7:])
                # print "KC=", KC
                
                nC = map(len, competitors)
                nKC = map(len, KC)
                # print "nC= ",nC
                # print "nKC=",nKC
                
                # SA=substrates+activators
                # print "SA=",SA
                
                # each substrate and activator is allowed to have 
                # a list of competitive inhibitors of undetermined length
                # the number of rates for each inhibitor should be the same 
                # as the number of species
                
                KCS=[]
                for i in range(len(competitors)):
                    kvalues=[]
                    if i<len(KC): kvalues=KC[i]
                    nkstart=len(kvalues)
                    # print i," ", kvalues," ",nC[i]
                    if len(kvalues)<nC[i]:
                        knew=getlast(kvalues,1)
                        while len(kvalues)<nC[i]:
                            kvalues.append(knew)
                    while len(kvalues)>nC[i]:
                        kvalues=kvalues[:-1]
                    # if len(kvalues)!=nkstart:
                    #    print "modified kvalues:", kvalues
                    KCS = KCS + [kvalues]
                # print "KCS=",KCS
                self.rates = self.rates + KCS

            return
 
        #*************************************************************************
        
        

        def _errorCheckMathRates():
        #
        # make sure that math-strings retain their leading and trailing quotation marks
        # these will be used to indicate the "math-ness" of the expression
        #
        
        #   The parser returns the string as a list with the quotes as separate
        #   entities, joinmath reattaches them
        # 
            def joinmath(l):
                # print "joinmath: l=",l
                if type(l)==type([]):
                    if len(l)>0:
                        # print "joinmath: l=",l
                        if l[0]=="\"" and l[-1]=="\"": 
                            return "".join(l)
                return l
                
            def mathcheck(s):
                if type(s)==type("hello world"):                 
                    if len(s)>2:
                        if(s[0]=="\"") and (s[-1]=="\""):                        
                            try:
                                # print "s[1:-1]: ",s[1:-1]
                                # expression = sympify(s[1:-1],locals=locals())
                                expression = sympify(s, locals=locals())
                            except:
                                print "Error: Parser: Unable to sympify the math "\
                                    "expression " + s
                                
                return
                
            r = self.rates
            r = r if type(r)==type("abc") else list(r)
            # print "errorCheckRates: rates: ", r
            
            r = map((lambda x: x if type(x)==type("abc") else list(x)), r)
            # print "errorCheckRates: rates: ", r

                    
            r = joinmath(r)
            if type(r)==type([]): r = map(joinmath, r)
            
            # print "errorCheckRates: rates: ", r
            # or k in r:   print "k = ", k, " type: ", type(k)
            
            
            mathcheck(r)
            if type(r)==type([]): map(mathcheck,r)
                
            self.rates = r
            return

        #*************************************************************************
 
        def _DetermineStoichiometries(oldvars):
        
            newvars = []; 
            sto = []
            n = len(oldvars)
            i = 0
            
            while i < n:
                this = (oldvars)[i]
                if i<n-1:
                    that = (oldvars)[i+1]
                else:
                    that = "DONE"
                if isnumber(this):
                    sto.append(this)
                    newvars.append(that)
                    i+=2
                else:
                    sto.append("1")
                    newvars.append(this)
                    i+=1                   
            return (newvars, sto)
            
        #*************************************************************************
        
        # TOLIST added 11/18/14 because under certain conditions self.lhs or 
        # self.rhs could be a single species (str) which would be converted to 
        # a list of charcters by list(X)
        
        def TOLIST(X):
            if isinstance(X, str):
                return ([X])
            return (list(X))
        
        
        def _SetStoichiometries():
        
            if self.arrowtype == "rational":
                self.slhs=[1]
                self.srhs=[1]
                return
            
            (newlhs, sleft)  = _DetermineStoichiometries(TOLIST(self.lhs))
            (newrhs, sright) = _DetermineStoichiometries(TOLIST(self.rhs))
            
            self.lhs = newlhs
            self.rhs = newrhs
            self.slhs = sleft
            self.srhs = sright
            
            return
        #***************************************
        
        def _checkcascade():
            r=self.parsed
            fields=len(r)
            if len(r)<6:
                # print "checkcascade: ", r , " is not a cascade"
                return False
            # print "checkcascade:",r
            k = list(r[-1])
            r = r[:-1]
            if evenq(fields):
                mods=[] 
            else:
                mods=list(r[-1])
                r=r[:-1]

            if len(k)>0:
                if k[0]=="rates": k=k[1:]
            if len(mods)>0:
                    if mods[0]=="mod": mods=mods[1:]
         
            # print "checkcascade:  k:",k," mods: ", mods
 
            n=len(r)        
            arrows =  [r[i] for i in range(1,n,2)]
            stages = len(arrows)
            species = map(list,[r[i] for i in range(0,n,2)])            
            
            iscascade = len(set(arrows))==1
            arrow = str(stages)+"*"+arrows[0] 
            if arrows[0] in ["|->", "|-->"]:
                clas = k[0]+"-";
                k=k[1:] 
            elif arrows[0] in [":->",":-->"]:
                clas = "MMH-"
                k=k[1:]
            else:
                clas = ""
            arrowtype = clas+"Cascade"
            if iscascade:
 
             
                self.arrowtype = arrowtype
                self.arrow = arrow
                self.lhs = species[0]
                self.slhs = []
                self.rhs = flatten(list(species[1:]))
                self.srhs = []
                self.mods = flatten(mods)
                self.rates=k
                
                
            else:
                print "Error: CheckCascade: The reaction '" + self.raw +\
                    "' appears to be an improperly formatted cascade. Check to "+\
                    " make sure all the arrows are the same." 
            # print "checkcascade: ", arrow, " iscascade: ", iscascade
                
            return iscascade
        #***************************************    


        def _reformatAllSpeciesIndices():
       
            self.lhs = map(deindex, self.lhs)
            self.rhs = map(deindex, self.rhs)
            self.mods =  map(deindex, self.mods)
        
            return
                
        #***************************************    
               
        # parse the input arrow        
        # print "DBG: arrow:", arrow
        # print "DBG: raw:", raw
        parsed = parse_cellerator_arrow(raw)
        # print "DBG: parsed: ", parsed

        # save the basic object
        inputstring=raw[:]
        self.raw = inputstring.strip()
        parsedobject=parsed[:]
        self.parsed = parsedobject
        Reaction.counter += 1
        self.num = Reaction.counter

        #decode the fields
        l = len(parsed)
        
 
        if l==5:
            lhs,arrow,rhs,mods,rates=parsed
        elif l==4:
            lhs,arrow,rhs,rates=parsed
            mods=[]
        elif l==3:
            first,second,third=parsed
            lhs,arrow,rhs=first,second,third
            rates=[]
            mods=[]
        else:
            data=str(raw.strip())
            matched,lb,rb=bracketsmatch(data)
            if not matched: 
                print "ARROW ", Reaction.counter," ERROR: the reaction " + data + " does not "\
                + "have matching brackets; there are "+str(lb)+" [  and "+str(rb)+" ]."\
                + "\nCorrect data file and re-submit job."
                raise SystemExit(1)
            if _checkcascade(): 
                _reformatAllSpeciesIndices()
                return
           
            
            
            print "ARROW ", Reaction.counter, "ERROR (Parser): Unable to identify reaction: '"+data+"'"
            # itemnumber=0
            # for item in parsed:
            #    print itemnumber,":",item
            #    itemnumber+=1
            # print "Counter = ", Reaction.counter
            lhs=rhs=rates=mods=[]
            arrow='Unknown'
            
        
        # print ">>>>> PARSED ARROW: ",parsed
        # print "HELLO", len(parsed)
        # print "lhs=",lhs
        # print "rhs=",rhs
        # print "arrow=", arrow
        # print "rates=",rates
        # print "mods=",mods
        # itemnumber=0
        # for item in parsed:
        #        print itemnumber,":",item
        #        itemnumber+=1
         

         
        needs_a_modifier = _needsModifier(arrow)
        if needs_a_modifier and l==4:
            (rates, mods) = (mods, rates)
        rates = rates[:]
        mods = mods[:]
 
        # print "**** rates=", rates, " mods = ", mods 

       
        if type(rates) is str:
            rates=[ rates ]
        if type(mods) is str:
            mods=[ mods ]

        if len(rates)>1 and rates[0]=="mod":
            (rates,mods)=(mods,rates)
        if len(mods)>1 and mods[0]=="rates":
            (rates,mods)=(mods,rates)

        if len(rates)>0 and rates[0]=="rates":
            rates=rates[1:]
        if len(mods)>0 and mods[0]=="mod":
            mods=mods[1:]

        if len(rates)>0 and rates[0]=="MMH":
            rates=rates[1:]
        
        
        #if len(rates)>0 and rates[0]=="Diffusion":
        #   rates=rates[1:]
        
        # print "**** rates=", rates, " mods = ", mods 
        rates = _CheckRates(rates,arrow)
        (rates,arrowtype) = _DetermineType(rates, arrow)
       
        # print "************** ARROWTYPE: ", arrowtype
        # print "arrow: ", arrow
        mods = _CheckMods(mods,arrow,arrowtype)
        
        if len(mods)>0 and mods[0]=="rational":
            rates = map(list, mods[1:])
            arrowtype="rational"
            mods=[]
           
           
        
        self.arrow=arrow
        self.lhs=lhs
        self.rhs=rhs
        self.rates=rates
        self.mods=mods
        self.arrowtype = arrowtype
        
        #
        #  Put Additional Error Checking Here
        #
        
        
        _errorCheckUSER()
        _errorCheckRATIONAL()
        _errorCheckMWC()
        _errorCheckMathRates()
        
        #
        # Set the Stoichiometries
        #
        _SetStoichiometries()
        _reformatAllSpeciesIndices()

        
    def Input(self):
        "Returns the Input string originally provided to a Reaction object."
        return self.raw

    def Parsed(self):
        "Returns the original output of pyparsing provided to a Reaction object."
        return self.parsed

    def Number(self):
        """
        Returns the reaction number assigned to a Reaction object.

        Reaction numbers are incremented by one each time a Reaction
        object is instantiated.
        """
        return self.num

    def Length(self):
        "Returns the length of the object returns by pyparsing to the Reaction object."
        return len(self.parsed)

    def ArrowType(self):
        "Returns the cellerator arrow type as a string."
        return self.arrowtype
    
    def Arrow(self):
        "Returns the type of arrow in the reaction as a string."
        return self.arrow

    def LHS(self):
        "Returns a list of reactants as strings."
        return self.lhs
        
    def LH_STOIC(self):
        "Retursn a list of the stoichiometries of all the reactants as strings"
        return self.slhs
       
    def RH_STOIC(self):
        "Returns a list of the stoichiometries of all of the products as strings"
        return self.srhs

    def RHS(self):
        "Returns a list of products as strings."
        return self.rhs

    def MHS(self):
        "Returns a list of modifiers as strings."
        return self.mods

    def Rates(self):
        "Returns a list of rate constants as strings."
        return self.rates

    def Dump(self):
        "Prints a dump of the arrow object to standard output."
        print "======================= Reaction Number ",\
              self.Number()," ======================="
        print "Input:    ", self.Input()
        print "Parsed:   ", self.Parsed()
        print "Arrow:    ", self.Arrow()
        print "ArrowType:", self.ArrowType()
        print "LHS:      ", self.LHS()
        print "LH_STOIC: ", self.LH_STOIC()
        print "RHS:      ", self.RHS()
        print "RH_STOIC: ", self.RH_STOIC()
        print "modifiers:", self.MHS()
        print "rates:    ", self.Rates()
        




def parse_cellerator_arrow(text_reaction):
    """
    Basic parser for Cellerator Shorthand
    """

#    newline = Literal("\n")
#    anything=Regex(".*")
#    USKeys="!@#$%^&*()_+-={}[]|\\\"\?/<>,.'"+ alphanums

#    junkword =Word(USKeys)
#    junkline = junkword+Optional(ZeroOrMore(junkword))


#    result = "hello world"

    arrow_one_dash =    Literal("->")
    arrow_two_dash =    Literal("-->")
    arrow_two_way =     Literal("<->")
    arrow_one_equal =   Literal("=>")
    arrow_colon_equal = Literal(":=>")
    arrow_two_equal =   Literal("==>")
    arrow_equal_2way =  Literal("<=>")
    arrow_tee =         Literal("|->")
    arrow_tee_dash =    Literal("|-->")
    arrow_colon_dash =  Literal(":->")
    arrow_colon_2dash=  Literal(":-->")
    
    an_arrow = arrow_one_dash | arrow_two_dash | arrow_two_way | arrow_one_equal |\
        arrow_colon_equal | arrow_two_equal | arrow_equal_2way | arrow_tee | \
        arrow_tee_dash | arrow_colon_dash | arrow_colon_2dash

    LB = Suppress("[")
    RB = Suppress("]")
    lb = Literal("[")
    rb = Literal("]")

    #identifier  = Word(alphas+'$', alphanums+'$')
    #
    # don't include $ because some Python pattern matching libraries 
    # use it to emulate linux pattern matching
    #
    identifier  = Word(alphas+"$_", alphanums+"_")
    mathstring = Word(alphanums+"_+-/* %^&~<>()|=:[],.")
    
    MATH = Keyword("Math")
    RATES = Keyword("rates")|Keyword("Rates")
    MOD = Keyword("mod")|Keyword("Mod")
    MMH = Keyword("MMH")
    MWC = Keyword("MWC")
    HILL = Keyword("Hill")
    GRN = Keyword("GRN")
    SSYSTEM=Keyword("SSystem")
    NHCA = Keyword("NHCA")
    USER = Keyword("USER")
    RATIONAL = Keyword("rational")
    USING = Keyword("using")
    DIFFUSION = Keyword("Diffusion")
    
    INFINITY = Keyword("inf") | Keyword("Inf") | Keyword ("INF") | Keyword("infinity") | \
        Keyword("Infinity") | Keyword("INFINITY")
    NEGINFINITY =Combine("-"+INFINITY)
    FLUX = Keyword("Flux")
    
   
    
    keywords=[ RATES, MOD, MMH, MWC, HILL, GRN, SSYSTEM, NHCA, USER, RATIONAL, \
        USING, DIFFUSION ]
    grykeywords=[HILL, GRN, SSYSTEM, NHCA, USER]
    multikeywords=[MWC,RATIONAL]
    
    

    for k in keywords:
        identifier.ignore(k)

    plus        = Suppress("+")
    comma       = Suppress(",")
    star        = Suppress("*")
    LT          = Suppress("<")
    GT          = Suppress(">")
    quote       = "\""
    unsignedinteger   = Word(nums)
    negativeinteger   = Combine("-"+unsignedinteger)
    integer     = unsignedinteger | negativeinteger
    number      = Combine(integer+"."+ unsignedinteger) | \
                  Combine(integer+".") | \
                  Combine("."+unsignedinteger) | \
                  integer
    asterisk    = Suppress("*")
    
    index = identifier|unsignedinteger
    indexes = index + ZeroOrMore(Literal(",")+index)  
    indexedIdentifier = Combine(identifier+"("+indexes+")")
    species = indexedIdentifier | identifier

    sto_species = species | number + asterisk + species | number + species
    product_of_species = Group(species + ZeroOrMore(star + species))
    
    
    sspecies = Group(species)
    sum_of_species=Group(species + ZeroOrMore(plus + species))
    sum_of_sto_species = Group(sto_species + ZeroOrMore(plus + sto_species))
    
    seq_of_species=species + ZeroOrMore(comma + species)
    seq_of_products=product_of_species + ZeroOrMore(comma + product_of_species)
    
    list_of_species = Group(LB+seq_of_species+RB)|Group(LB+RB)
    list_of_products_of_species = Group(LB+seq_of_products+RB)|Group(LB+RB)
    
    list_of_modifiers= Group(MOD + LB + seq_of_species + RB)
    modifier_group= Group(species) | list_of_modifiers 
    nested_species = (species | list_of_species) + ZeroOrMore(comma + (species | list_of_species))
    nested_modifier_group = Group(MOD+LB+nested_species+RB)
        
    math_rate = Group(quote + mathstring + quote) 
    rate_constant = ( indexedIdentifier | identifier | number | math_rate )
    
    sequence_of_rates = rate_constant + (ZeroOrMore(comma + rate_constant))
    rate_list =  (LB + sequence_of_rates + RB) | LB+RB
    list_of_rates = Group(rate_list)
    multi_rate = rate_constant | Group(rate_list)
    rates = RATES + LB + sequence_of_rates + RB
    MMHrates = Group(MMH + LB + sequence_of_rates + RB)
    rate_group=Group( rate_constant | rates )
    
    single_rate = Group(rate_constant | (RATES+  LB + rate_constant + RB))
    
     
    usingrates = Group(USING+LB+math_rate+RB)
    USING_arrow = sum_of_sto_species + Literal("->") + sum_of_sto_species + comma + usingrates
    
    diffusionrate = Group(DIFFUSION+LB+rate_constant + RB)
    DIFFUSION_arrow = species + Literal("->") + species + comma + diffusionrate


    lowlimit = number | NEGINFINITY
    uplimit = number | INFINITY
    fluxrange = lowlimit + LT + identifier + LT + uplimit
    fluxrates = Group(FLUX + LB + fluxrange + Optional(comma + number + Optional(comma + number)) + RB)
    flux_arrow = sum_of_sto_species + arrow_one_dash + sum_of_sto_species \
        + comma + fluxrates 


    simple_mass_action_arrow = (sum_of_species + arrow_one_dash + sum_of_species) \
        + Optional(comma + rate_group)
        

    simple_sto_mass_action_arrow = (sum_of_sto_species + arrow_one_dash + sum_of_sto_species) \
        + Optional(comma + rate_group)


    modified_sto_mass_action_arrow = sum_of_sto_species + arrow_two_dash + sum_of_sto_species \
        + comma +  modifier_group + comma + rate_group 


    modified_mass_action_arrow = sum_of_species + arrow_two_dash + sum_of_species \
        + comma +  modifier_group + comma + rate_group 


    mass_action_bidirectional = sum_of_species + arrow_two_way + sum_of_species\
        + comma + rate_group
        
    mass_action_sto_bidirectional = sum_of_sto_species + arrow_two_way + sum_of_sto_species\
        + comma + rate_group
        

    mass_action_catalytic = sspecies + arrow_one_equal + sspecies  \
        + comma + modifier_group + comma + rate_group

    mass_action_catalytic_bidirectional = \
        sspecies + arrow_equal_2way + sspecies  \
        + comma + modifier_group + comma + rate_group

    mass_action_catalytic_sequential = \
        sspecies + arrow_colon_equal + sspecies \
        + comma + modifier_group + comma + rate_group
                       
    MMH_simple = sspecies + arrow_colon_dash + sspecies + comma + MMHrates
    MMH_explicit = sspecies + arrow_colon_2dash + sspecies + comma + \
                   modifier_group + comma + MMHrates

    GRNKEY = HILL | GRN | SSYSTEM | NHCA | USER
    GRNrates = Group(GRNKEY + LB + sequence_of_rates + RB)
    GRNMultirates = Group(GRNKEY + LB + multi_rate + ZeroOrMore(comma + multi_rate) + RB)
    GRN_Basic = sspecies + arrow_tee + sspecies + comma + GRNrates
    GRN_Multi = list_of_species + arrow_tee + sspecies + comma + GRNMultirates 
    GRN_Catalytic = sspecies + arrow_tee_dash + sspecies + comma + \
                    modifier_group + comma + GRNrates
    GRN_Multi_Catalytic = list_of_species + arrow_tee_dash + sspecies + comma + \
                    modifier_group + comma + GRNMultirates               
                   
    RATIONAL_arrow = Group( LB+list_of_products_of_species+comma+\
        list_of_products_of_species+RB )  + \
        arrow_two_equal + sspecies + comma +   Group( RATIONAL + LB + list_of_rates + comma + \
        list_of_rates + comma + list_of_rates + comma + list_of_rates +  RB)
    


    MWCRates = Group(MWC + LB + multi_rate + ZeroOrMore(comma + multi_rate) + RB)
    
    MWC_arrow = (sspecies | list_of_species) + arrow_two_equal + sspecies + comma + \
      (modifier_group | nested_modifier_group) + comma + MWCRates 
       

    modified_cascade =  sspecies + an_arrow + sspecies + OneOrMore(an_arrow+ sspecies) + \
        (comma + modifier_group) + \
        (comma + rate_group)   
    simple_cascade =    sspecies + an_arrow + sspecies + OneOrMore(an_arrow+ sspecies) + \
        (comma + rate_group) 
 
    GRN_modified_cascade =    sspecies + Literal("|-->") + sspecies + OneOrMore(Literal("|-->")+ sspecies) + \
        (comma + modifier_group) + \
        (comma + GRNMultirates) 

    GRN_cascade =    sspecies + Literal("|->") + sspecies + OneOrMore(Literal("|->")+ sspecies) + \
        (comma + GRNMultirates) 

    MMH_modified_cascade =    sspecies + Literal(":-->") + sspecies + \
        OneOrMore(Literal(":-->") + sspecies) + \
        (comma + modifier_group) + \
        (comma + MMHrates) 

    MMH_cascade =    sspecies + Literal(":->") + sspecies + \
        OneOrMore(Literal(":->")+ sspecies) + \
        (comma + MMHrates) 

    cascade = simple_cascade | modified_cascade | GRN_cascade | GRN_modified_cascade  \
    | MMH_cascade  | MMH_modified_cascade
  
    plain_arrow = USING_arrow \
                | DIFFUSION_arrow \
                | flux_arrow \
                | simple_mass_action_arrow \
                | simple_sto_mass_action_arrow \
                | modified_mass_action_arrow \
                | modified_sto_mass_action_arrow \
                | mass_action_bidirectional \
                | mass_action_sto_bidirectional \
                | MMH_simple \
                | MMH_explicit \
                | mass_action_catalytic \
                | mass_action_catalytic_bidirectional \
                | mass_action_catalytic_sequential \
                | GRN_Basic \
                | GRN_Catalytic \
                | GRN_Multi \
                | GRN_Multi_Catalytic \
                | RATIONAL_arrow \
                | MWC_arrow 
                
    # kMech Arrows
    
    # this particular choice is actually idential to
    #  arrow_equal_2way which is also used for enzymatic arrows
    arrow_kMech=  Literal("<=>")
   
    Enz = Keyword("Enz")
    BiBi=Keyword("BiBi")
    BiTer=Keyword("BiTer")
    BiUni=Keyword("BiUni")
    CI=Keyword("CI")
    CI1=Keyword("CI1")
    CI2=Keyword("CI2")
    En=Keyword("En")
    Enx=Keyword("Enx")
    MulS=Keyword("MulS")
    NC=Keyword("NC")
    NC2=Keyword("NC2")
    NCI=Keyword("NCI")
    NCI1=Keyword("NCI1")
    NCI2=Keyword("NCI2")
    NCmI=Keyword("NCmI")
    OrderedBiBi=Keyword("OrderedBiBi")
    OrderedBiUni=Keyword("OrderedBiUni")
    PingPong=Keyword("PingPong")
    PingPongTerTerF=Keyword("PingPongTerTerF")
    PingPongTerTerR=Keyword("PingPongTerTerR")
    RandomBiBi=Keyword("RandomBiBi")
    TerBi=Keyword("TerBi")
    TerTer=Keyword("TerTer")
    UCI=Keyword("UCI")
    UniBi=Keyword("UniBi")
    UniUni=Keyword("UniUni")
    kMechKeyword = BiBi|BiTer|BiUni|CI|CI1|CI2|En|Enx|MulS|NC|NC2|NCI|NCI1|NCI2|NCmI|\
        OrderedBiBi|OrderedBiUni|PingPong|PingPongTerTerF|PingPongTerTerR|\
        RandomBiBi|TerBi|TerTer|UCI|UniBi|UniUni
        
    kMechRate=Group(kMechKeyword+LB+multi_rate+ZeroOrMore(comma+multi_rate)+RB)
      
    kMechRateList=kMechRate+OneOrMore(comma + kMechRate)
    
    
    kMechSpecies = Group(species | LB + species + ZeroOrMore(comma + species) + RB)
    kMechRates = Group(Enz + LB + kMechRateList+RB)
    kMechArrow = kMechSpecies + arrow_kMech + kMechSpecies  \
        + comma + modifier_group + comma + kMechRates
    bracketedkMechArrow = LB + kMechArrow + RB

    # or some other type of arrow

    bracketed_arrow = (LB + plain_arrow + RB)
    

    arrow = plain_arrow | bracketed_arrow   \
            | kMechArrow | bracketedkMechArrow \
            | (LB + cascade + RB) \
            | restOfLine
    
    result = arrow.parseString(text_reaction)

    # print result

    return result
#*****************************************************************************

def ParseArrow(text):
    return Reaction(text)

#*****************************************************************************
    
def StripComments(lineOfText, commentchar):
    #
    # input: line of text, comment delimiter
    # return value: line of text, with all characters starting from the
    #   comment delimiter to the end of the line removed, and all 
    #   additional white space either at the beginning or end of the 
    #   line (after the comment has been removed) also stripped off
    #
    # example:  StripComments(" hello world! # this is a test", "#")
    # should return "hello world!"
    #
    i = lineOfText.find(commentchar)
    if i>= 0:
        #print "stripping comment on ",lineOfText
        newlineoftext=lineOfText[:i]
    else:
        newlineoftext=lineOfText
    newlineoftext=newlineoftext.strip()
    return newlineoftext
    
def parser(inputfile="", trace=False):
    if not (inputfile==""):
        infile = str(inputfile)
    elif ("-in" in sys.argv):
        i = (sys.argv).index("-in")+1
        if i<len(sys.argv):
            infile=sys.argv[i]
        else:
            raise SystemExit("Error: expecting filename after -in")
    else:
        raise SystemExit("No Input File Given; use -in option.")
    if not os.path.isfile(infile):
        raise SystemExit("Error: "+os.path.abspath(infile)+" File not found.")
        
        
    if os.path.splitext(infile)[1]==".model":
        (lines, ics, rates, frozenvars, functions, assignments,fullfile)=readmodel(INFILE=infile)
    else:
        f=open(infile)         
        lines = f.readlines()
        f.close()    
    # f=open(infile)
    # lines = f.readlines()
    # f.close()
    results=[]
    rs=[]

    for line in lines:
        line = StripComments(line,"#")
        if line=="": continue

        r=Reaction(line)
        if "-dump" in sys.argv: r.Dump()
        if not r.Arrow() == "Unknown":
            results.append(r.Input().strip())
            rs.append(r)
            if trace or ("-trace" in sys.argv):
                print "input: " + str(r.Input()).strip()
        elif trace or ("-trace" in sys.argv):
                print "input ignored: " + line.strip()
    # print results
    return rs

def runpyx():
    if not ("-quiet" in sys.argv):
        d=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print "py[cellerator]: parser ("+d+")"
    rrr = parser()
    if ("-quiet" in sys.argv): return
    print rrr

if __name__ == "__main__":
    runpyx()

