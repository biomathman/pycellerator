#!/usr/bin/env python
#
#  Cellerator Text Input --> Matheamtica List Conversion
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


from parser import Reaction
import os.path
import re
import sys
import utils
from reader import readmodel
import datetime
from expander import expand

#******************************************************************************

def _cleanup(stuff):
    "Removes quotes and apostrophes from a string"
    cleanstuff = re.sub('[\"\']', '', stuff)
    return cleanstuff

#******************************************************************************

def _MList(head, somelist):
    """
    returns a string representing a Mathematica expression
        head[arguments]
    whose head is given by "head" and whose arguments are
    given by the elements of "somelist". If somelist is not
    a list or a tuple it is converted into a single element
    list, e.g.,
        MList("foo",("a","b")) --> "foo[a,b]"
        MList("foo",17) --> 'foo[17]'
    """
    l=somelist[:]
    if type (l)==list:
        pass
    elif type(l) == tuple:
        l = list(l)
    else:
        l = list([l])
    return _cleanup(str(head)+str(l))
    
#******************************************************************************
    
def _MListRates(arrow, TYPE):

    rates = []
    for r in arrow.Rates():
        if type(r) == type([]):
            rates.append(_MList('List',r))
        else:
            rates.append(r)   
    rates = _MList(TYPE, rates)
    return rates

#******************************************************************************

def _Simple_Mass_Action(arrow):
    """
    [ A + B + C -> P + Q, k ] or
    [ A + B + C -> P + Q, rates[k] ] 
    """
    # print arrow.Arrow(), " = ", arrow.ArrowType()
    lhs = _MList('Plus',arrow.LHS())
    rhs = _MList('Plus',arrow.RHS())
    r = _MList('ShortRightArrow', [lhs, rhs])
    rates = arrow.Rates()
    if arrow.ArrowType()=="Flux": rates = [_MList('Flux',rates)]
    r = _MList('List', [r] + rates)
    return r

#******************************************************************************

def _Modified_Mass_Action(arrow):
    """
    [ A + B  --> P, E, k ] or
    [ A + B  --> P, mod[E], rates[k] ]
    equivalent to:
    [ A + B + E -> P + E, k], etc.
    """
    # print arrow.Arrow(), " = ", arrow.ArrowType()
    lhs = _MList('Plus',arrow.LHS())
    rhs = _MList('Plus',arrow.RHS())
    r = _MList('Rule', [lhs, rhs])
    modifiers = arrow.MHS()
    r = _MList('Overscript', [r] + modifiers)
    r = _MList('List',[r] + arrow.Rates())
    return r

#******************************************************************************

def _Mass_Action_Bidirectional(arrow):
    """
    [ A + B + C <-> P + Q, rates[kf, kr]]
    """
    # print arrow.Arrow(), " = ", arrow.ArrowType()
    lhs = _MList('Plus',arrow.LHS())
    rhs = _MList('Plus',arrow.RHS())
    r = _MList('RightArrowLeftArrow', [lhs, rhs])
    rates = arrow.Rates()
    r = _MList('List', [r] + rates)
    return r

#******************************************************************************

def _Mass_Action_Catalytic(arrow):
    """
    [ X => Y, mod[E], rates[k1,k2,...] ]
    """
    # print arrow.Arrow(), " = ", arrow.ArrowType()
    lhs = arrow.LHS()
    rhs = arrow.RHS()
    mhs = arrow.MHS()
    rates=arrow.Rates()
    r = _MList('RightArrowLeftArrow', lhs + rhs)
    r = _MList('Overscript', [r]+mhs)
    r = _MList('List', [r]+rates)
    return r

#******************************************************************************

def _Mass_Action_Catalytic_Bidirectional(arrow):
    """
    [ X<=>Y, mod[F,R], rates[k1,k2,..] ]
    """
    # print arrow.Arrow(), " = ", arrow.ArrowType()
    lhs = arrow.LHS()
    rhs = arrow.RHS()
    mhs = arrow.MHS()
    rates=arrow.Rates()
    mhs.reverse()
    r = _MList('RightArrowLeftArrow', lhs + rhs)
    r = _MList('Underoverscript',[r]+mhs)
    r = _MList('List', [r]+rates)
    return r    

#******************************************************************************

def _Mass_Action_Catalytic_Sequential(arrow):
    """
    [ X :=> Y, mod[F], rates[k1,k2,...]]
    """
    # print arrow.Arrow(), " = ", arrow.ArrowType()
    lhs = arrow.LHS()
    rhs = arrow.RHS()
    mhs = arrow.MHS()
    rates=arrow.Rates()
    r = _MList('Equilibrium', lhs + rhs)
    r = _MList('Overscript', [r]+mhs)
    r = _MList('List', [r]+rates)
    return r

#******************************************************************************

def _MMH_Basic(arrow):
    """
    [x:->Y, MMH[parameters]]
    """
    # print arrow.Arrow(), " = ", arrow.ArrowType(), " not yet implemented: ", arrow.Input()
    S = arrow.LHS()
    P = arrow.RHS()
    rates = _MList('MM',arrow.Rates())
    reaction=_MList('DoubleLongRightArrow',S+P)
    r = _MList('List', [reaction,rates])
    return r

#******************************************************************************

def _MMH_Explicit(arrow):
    """
    [X:-->Y, mod[cat], MMH[parameters]]
    """
    # print arrow.Arrow(), " = ", arrow.ArrowType(), " not yet implemented: ", arrow.Input()
    S = arrow.LHS()
    P = arrow.RHS()
    E = arrow.MHS()
    r = _MList('DoubleLongRightArrow',S+P)
    r = _MList('Overscript',[r]+E)
    rates = _MList('MM', arrow.Rates())
    r = _MList('List',[r,rates])
    return r

#******************************************************************************
    
def _Format_Basic_GRN_Arrow(arrow,rates):
    S = arrow.LHS()
    P = arrow.RHS()
    if len(S)>1: 
        S=[_MList('List',S)]
        r = _MList('RightTeeArrow',S+[P[0]])
    else:
        r = _MList('RightTeeArrow',S+P)
    r = _MList('List',[r, rates])
    return r    

#******************************************************************************
    
def _Hill_Basic(arrow):
    """
    [X|->Y,Hill[parameters]]
    """
    # print arrow.Arrow(), " = ", arrow.ArrowType(), " not yet implemented: ", arrow.Input()
    rates = _MListRates(arrow, 'Hill')
    r = _Format_Basic_GRN_Arrow(arrow, rates)
    return r

#******************************************************************************

def _Hill_Catalytic (arrow):
    """
    [X|->Y,mod[cat],Hill[parameters]]
    """
    # print arrow.Arrow(), " = ", arrow.ArrowType(), " not yet implemented: ", arrow.Input()

    rates = _MListRates(arrow, 'Hill')
    S = arrow.LHS()
    P = arrow.RHS()
    E = arrow.MHS()
    if len(S)>1:
        S=[_MList('List',S)]
        # print "S=",S
        r = _MList('RightTeeArrow',S+[P[0]])
    else:
        r = _MList('RightTeeArrow',S + P)
    r = _MList('Overscript', [r]+E)
    r = _MList('List', [r, rates])
    return r

#******************************************************************************
    
def _GRN_Catalytic (arrow):
    """
    [X|->Y,mod[cat],GRN[parameters]]

    """
    # print arrow.Arrow(), " = ", arrow.ArrowType(), " not yet implemented: ", arrow.Input()
 
    rates = _MListRates(arrow, 'GRN')
    S = arrow.LHS()
    P = arrow.RHS()
    E = arrow.MHS()
    if len(S)>1:
        S=[_MList('List',S)]
        # print "S=",S
        r = _MList('RightTeeArrow',S+[P[0]])
    else:
        r = _MList('RightTeeArrow',S + P)
    r = _MList('Overscript', [r]+E)
    r = _MList('List', [r, rates])
    return r

#******************************************************************************
    
def _GRN_Basic(arrow):
    """
    [X|->Y,GRN[parameters]]
    """
    # print arrow.Arrow(), " = ", arrow.ArrowType(), " not yet implemented: ", arrow.Input()
    
    rates = _MListRates(arrow, 'GRN')
    r = _Format_Basic_GRN_Arrow(arrow, rates)
    return r

#******************************************************************************
    
def _SSystem_Basic(arrow):
    """
    [X|->Y,SSystem[parameters]]
    """
    # print arrow.Arrow(), " = ", arrow.ArrowType(), " not yet implemented: ", arrow.Input()
    
    rates = _MListRates(arrow, 'SSystem')
    r = _Format_Basic_GRN_Arrow(arrow, rates)
    return r

#******************************************************************************
    
def _NHCA_Basic(arrow):
    """
    [X|->Y,NHCA[parameters]]
    """
    # print arrow.Arrow(), " = ", arrow.ArrowType(), " not yet implemented: ", arrow.Input()
    
    rates = arrow.Rates()
    v,TP,TM,n,m,k=rates
    if type(TP)==type([]): TP = _MList("List",TP)
    if type(TM)==type([]): TM = _MList("List",TM)
    T=[TP,TM]
    T=_MList("List", T)
    if type(n)==type([]): n = _MList("List",n)
    rates = _MList('NHCA', [v, T, n, m, k])
    r = _Format_Basic_GRN_Arrow(arrow, rates)
    return r

#******************************************************************************    

GRNBasicHandler={"Hill": _Hill_Basic,
                 "GRN":  _GRN_Basic,
                 "SSystem": _SSystem_Basic,
                 "NHCA": _NHCA_Basic}

                 
GRNCatalyticHandler={"Hill": _Hill_Catalytic, 
                    "GRN": _GRN_Catalytic}

#******************************************************************************
                 
def _GRN_Arrow_Handler(arrow):
    # print arrow.Arrow(), " = ", arrow.ArrowType(), " not yet implemented: ", arrow.Input()    
    a = arrow.ArrowType()
    # print "_GRN_ARROW_Handler:", arrow.Input()
    if not a in GRNBasicHandler:
        print "A converter for the GRN reaction ", arrow.Input(), " has not been implemented."
        return "Error (_GRN_Arrow_Handler)"
    h=GRNBasicHandler[a]
    r = h(arrow)
    return r

#******************************************************************************

def _GRN_Catalytic_Arrow_Handler(arrow):
    # print arrow.Arrow(), " = ", arrow.ArrowType(), " not yet implemented: ", arrow.Input()
    a = arrow.ArrowType()
    # print "_GRN_Catalytic_Arrow_Handler:", arrow.Input()
    if not a in GRNCatalyticHandler:
        print "A converter for the GRN reaction ", arrow.Input(), " has not been implemented."
        return "Error (_GRN_Catalytic_Arrow_Handler)"
    h = GRNCatalyticHandler[a]
    r = h(arrow)
    return r

#******************************************************************************

def _Rational_Arrow(arrow):
    # print "Rational arrow not implemented: ", arrow.Input()
    num = _MList("List",  list(arrow.LHS())[0])
    den = _MList("List", list(arrow.LHS())[1] )
    r = _MList("List", [num,den])
    prod = list(arrow.RHS())
    r = _MList("DoubleLongRightArrow", [r]+prod)
    rates = arrow.Rates()
    rates = map(lambda x: _MList("List",x),  rates)
    rates = _MList("rational", rates)
    r = _MList("List",[r, rates])
    return r
#******************************************************************************

def _MWC_Arrow(arrow):
    # print "MWC arrow not implemented: ", arrow.Input()
    S = list(arrow.LHS())
    P = list(arrow.RHS())
    mhs = list(arrow.MHS())
    E = mhs[0]
    A = mhs[1]
    I = mhs[2]
    
    # print "*** MWC: S = ", S, " S = ", P, " P = ", P
    # print "*** MWC: E = ", E, " A = ", A, " I = ", I
    C = mhs[3:]
    AREC = len(utils.flatten(C))>0
    
    # print "*** MWC: competitors: ", AREC," = ", C
    C = map(lambda x:_MList("List",x), C)
    # print C
    
    if len(S)>1:
        S = [_MList("List",S)]
    r = _MList("DoubleLongRightArrow", S+P)
    # print "*** MWC: r:", r
    A=_MList("List",A)
    I=_MList("List",I)
    AREAI = len(A)>0 and len(I)>0
    
    
    if AREC:
        AI = _MList("List",[A,I]+C)
    else:
        AI = _MList("List",[A,I])
    # print "*** MWC: AI: ", AI

    r = _MList("Underoverscript",[r, AI, E])
    rates = list(arrow.Rates())
    # if not (len(rates) == 7):
    #    print "MWC Error: Converter Not Implemented: _MWC_Arrow: ", arrow.Input()
    #    # raise sys.exit("Parsing Error: Unable to Handle Input")
    #    return 
    KVALS = rates[4:]
    # print "*** MWC: K = ", KVALS
    K,KA,KI,KC = KVALS[0], KVALS[1],KVALS[2],KVALS[3:]
    # print "*** MWC: K = ", K, " KA=", KA, " KI=", KI, " KC=",KC

    K = _MList("List",K)
    KA = _MList("List",KA)
    KI = _MList("List",KI)
    KC = map(lambda x:_MList("List",x),KC)
    if AREC:
        K = _MList("List", [K, _MList("List",[KA,KI] + KC)])
    else:
        K = _MList("List", [K, _MList("List",[KA,KI])])
    rates = rates[:4]
    rates.append(K)
    rates = _MList("MWC", rates)
    r=_MList("List",[r,rates])
    return r

#******************************************************************************

DoubleEqualHandler={"rational": _Rational_Arrow,
                    "MWC": _MWC_Arrow
                    }

#******************************************************************************
    
def _Double_Equal_Arrow_Handler(arrow):
    # print "_Double_Equal_Arrow_Handler:",arrow.Input()
    a = arrow.ArrowType()
    if not a in DoubleEqualHandler:
        print "A converter for the ==> reaction " + str(arrow.Input()).rstrip()+" has not been implemented."
        return "Error (_Double_Equal_Arrow_Handler)"
    h = DoubleEqualHandler[a]
    r = h(arrow)
    return r

Handler = {"->": _Simple_Mass_Action,
           "-->":_Modified_Mass_Action,
           "<->":_Mass_Action_Bidirectional,
           "=>": _Mass_Action_Catalytic,
           "<=>":_Mass_Action_Catalytic_Bidirectional,
           ":=>":_Mass_Action_Catalytic_Sequential,
           ":->":_MMH_Basic,
           ":-->":_MMH_Explicit,
           "|->": _GRN_Arrow_Handler,
           "|-->":_GRN_Catalytic_Arrow_Handler,
           "==>":_Double_Equal_Arrow_Handler
           }
#******************************************************************************
MathArrowGlyph = {    "->":"ShortRightArrow",
                      "<->":"RightArrowLeftArrow",
                      "=>":"RightArrowLeftArrow",
                      "<=>":"RightArrowLeftArrow",
                      ":=>":"Equilbirium",
                      "|->":"RightTeeArrow",
                      "|-->":"RightTeeArrow",                      
                      ":->":"DoubleLongRightArrow",
                      ":-->":"DoubleLongRightArrow"
                      }


def CascadeHandler(r):
    DBG = False
    casc = r.Arrow()
    stages,arrow = casc.split("*")
    modifiers = r.MHS()
    
    if DBG: print "Cascade Handler: stages,arrow: ", stages, arrow
    if DBG: print "Cascade Handler: modifiers: ", modifiers
    #
    # most of the cascades in py[cellerator] don't actually exist in cellerator
    # so expand to the bottom line
    #
    
    er = expand(r)
    
    er = map(_Handle_Reaction, er)
    
    if DBG: print "Cascade Handler: er: ", er
    
    okr = map(lambda u:u[-1], filter(lambda u:u[0], er))
    notok = map(lambda u:u[-1], filter(lambda u:not(u[0]), er))
    
    if DBG: print "okr: ", okr
    if DBG: print "notok: ", notok
    for r in notok:
        print "Error: unable to convert the reaction: "+str(r)
    
    
    return okr

#******************************************************************************

def _Handle_Reaction(parsed):
    a = parsed.Arrow()
    arrowtype=parsed.ArrowType()
    if a in Handler: 
        h = Handler[a]
        r = h(parsed)
        return (True, str(r))
    else:
        return (False, "")
    

def _Handle_Cellerator_Reaction(raw, parsed):
    
    OK, r = _Handle_Reaction(parsed)
    if OK: return [ str(r) ]
    elif (parsed.ArrowType()).endswith("Cascade"):
        # print "This is a Cascade"
        r = CascadeHandler(parsed)
        return r
    else:
        print "The reaction ", raw, " could not be identified."
        return [ "Error (Handle_Cellerator_Reaction) '"+raw+"'"]

#******************************************************************************
def listspecies(reactions):
    # print "listspecies:"
    species=[]
    for r in reactions:
        lhs = list(r.LHS())
        rhs = list(r.RHS())
        mhs = list(r.MHS())
        s=list(set(utils.flatten([lhs,rhs,mhs])))
        # print "s=",s
        species.append(s)
    species=list(set(utils.flatten(species))) 
    if "Nil" in species: species.remove("Nil")  
    # print "species=",species
    return species
    
def findSplitLocation(arrow):
    level=0
    start=True
    i=-1
    for character in arrow:
        i+=1
        if character=="[": level+=1
        if character=="]": level-=1
        if start and level>1: 
            start = False
        if (not start) and level==1:
            break
    return i+1        
            
            
    
#******************************************************************************

def functionate(arrows, species):
    results=[]
    for a in arrows:
        print "DBG: functionate: a = ", a
        if not(a.startswith("Error")):
            #      "\n     01234567890123456789012345678901234567890123456789012345678901234567890"
            isplit = findSplitLocation(a)
            # print "isplit=",isplit
            # first,rest = a.split(",",1)
            first = a[:isplit]
            rest=a[isplit:]
            # print "split = ", first, rest
            for s in species:
                rest = rest.replace(s, s+"[t]")
            a =first+rest
            a = a.replace("Nil","\\[EmptySet]")
        results.append(a)
        # print "revised = ", a    
    return results 
#******************************************************************************

def GetCelleratorNetwork(filename, dump=False):
    if not os.path.isfile(filename):
        print "Error: readPLY: unable to locate file ", os.path.abspath(filename)
        raise SystemExit("File Not Found")

    rates={}
    ics={}
    if os.path.splitext(filename)[1]==".model":
        (linesOfData, ics, rates, frozenvars, functions,assignments, fullfile)=readmodel(INFILE=filename)
    else:
        f=open(filename)         
        linesOfData = f.readlines()
        f.close()

    lines=[]
    for line in linesOfData:
        line=(line.split("#"))[0].strip()     
        if len(line)>0:  lines.append(line)
    linesOfData = lines
    # f = open(filename)
    # linesOfData = f.readlines()
    # f.close()
    n = len(linesOfData)
    if not ("-quiet" in sys.argv):
        print n," lines read from ", os.path.abspath(filename)
    result=[]
    reactions = []
    for i in range(n):
        raw=linesOfData[i]
        r=Reaction(raw)
        reactions.append(r)  
        if dump: r.Dump()
        next = _Handle_Cellerator_Reaction(raw,r)
        # result.append()
        # print "next = ", next
        result = result + next 
        
    species = listspecies(reactions)
    # result = functionate(result, species)    
    # print "species = ", species
    # print "result = ", result
        
    return (result,rates,ics)

#******************************************************************************

def _test_New_Reaction(raw, dump=True):
    r=Reaction(raw)
    if dump: 
        r.Dump()
    print _Handle_Cellerator_Reaction(raw, r)
#******************************************************************************

def convert(inputfile="", outputfile=""):
    if ("-in" in sys.argv):
        i=(sys.argv).index("-in")+1
        if i < len(sys.argv):
            infile = sys.argv[i]
        else:
            print "Error: converter: expecting file name after -in option."
            raise SystemExit("Correct and resubmit job.")
    elif not inputfile == "":
        infile = inputfile
    else:
        print "Error: converter: expecting `python convert.py -in filename'"
        raise SystemExit("Correct and resubmit job.")
        
        
    dumpdb = "-dump" in sys.argv
    
    (output,rates,ics) = GetCelleratorNetwork(infile, dump=dumpdb)
 
    if ("-out" in sys.argv):
        i=(sys.argv).index("-out")+1
        if i < len(sys.argv):
            outfile = sys.argv[i]
        else:
            print "Warning: converter: expecting file name after -out option."
            outfile = "translated-model.nb"
    elif not outputfile == "":
        outfile = outputfile
    else:
        # print "Warning: converter: no output file specified."
        outfile = "translated-model.nb" 
    
    outfile = utils.uniqueFileName(outfile, type="nb")
    
    
    output = ",\n ".join(output)
    output = "<<xlr8r.m;\nmodel = {"+output+"}/.{Nil->\[EmptySet],inf->\[Infinity]};\n"
    output = output.replace("_","\[UnderBracket]")
        
    def dictpair(d,x):
        return str(x).replace("_","\[UnderBracket]")+"->"+str(d[x])
        
    therates = [dictpair(rates,x) for x in rates]
    therates ="therates = {"+",\n ".join(therates)+"};\n"

    theics = [dictpair(ics,x) for x in ics]
    theics = "theics = {" + ",\n ".join(theics) + "};\n"
    
    of = open(outfile,"w")
    of.write(output)
    of.write(therates)
    of.write(theics)
    of.close()
    
    fout = os.path.abspath(outfile)
    return fout

#******************************************************************************
    
    
def converter():
    if not ("-quiet" in sys.argv):
        d=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print "py[cellerator]: converter ("+d+")"
    u = convert()
    if "-quiet" in sys.argv: return
    print "Output written to ", u
    
if __name__ == "__main__":
    converter()

   
