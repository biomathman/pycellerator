#!/usr/bin/env python
#
# Cellerator Expander Module
# B.E.Shapiro 3 April 2012 
# rev. 30 Aug 2015
#
#================================================
#
#****************************************************************************
#
#    This file is part of Cellerator
#
#    Cellerator converts reactions, expressed in a text-formatted arrow-based 
#    notation into differential equations and performs numerical simulations.
#    Copyright (C) 2012 - 2015 Bruce E Shapiro.
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
#
# Expand a reaction into its base constituent reactions 
#
import datetime
from parser import *
from utils import *
from reader import readmodel
from kMech import expandkMech

#*******************************************************************************

def notZero(s):
    # input: string s
    # output: True if the string is different from a number representing 0, e.g., "0" or "0.0"
    #         False if s=="0"
    try:
        float(s)
        return not(float(s)==0)
    except ValueError:
        return True

#*******************************************************************************
    
def expand_no_expansion_necessary(reaction):
    return [reaction]

#*******************************************************************************
    
def expand_mass_action_bidirectional_reaction(reaction):
    # input: ParsedReaction with <-> arrow
    #  
    # output: List of two component ParsedReaction objects with -> arrows
    #  i.e., splits the A+B+... <-> P+Q+... into 
    #                   A+B+...  -> P+Q+... AND
    #                   P+Q+...  -> A+B 
    #  both as ParsedReaction objects.
    #
    # print "expand_mass_action_bidirectional_reaction"
    #
    lhs = list(reaction.LHS())
    lhss = list(reaction.LH_STOIC())
    for i in range(len(lhss)):
        if lhss[i]=="1": lhss[i]=""
         
    rhs = list(reaction.RHS())
    rhss = list(reaction.RH_STOIC())
    for i in range(len(rhss)):
        if rhss[i]=="1": rhss[i]=""

    rates = list(reaction.Rates())
    l = "+".join(map(lambda x:"".join(x), zip(lhss, lhs)))
    r = "+".join(map(lambda x:"".join(x), zip(rhss, rhs)))
    
    r1 = "["+l+"->"+r+","+rates[0]+"]"
    r2 = "["+r+"->"+l+","+rates[1]+"]"
    
    # print "bidirectional r1=", r1, " r2=", r2
    
    r1 = Reaction(r1)
    r2 = Reaction(r2)
    # print r1,r2
    
    # print "r1=", r1.Input()
    # print "r2=",r2.Input()
    return [r1, r2]

#*******************************************************************************

def expand_mass_action_trivially_catalyzed(reaction):
    # input: ParsedReaction for [ LHS-->RHS, mod[E], rates[k]] arrow
    #  
    # output: List of single componentsimple ParsedReaction objects with -> arrows
    #  [LHS + E -> RHS+E, k]
    #   as a ParsedReaction object.
    #
    lhs = list(reaction.LHS())
    lhss = list(reaction.LH_STOIC())
    for i in range(len(lhss)):
        if lhss[i]=="1": lhss[i]=""
         
    rhs = list(reaction.RHS())
    rhss = list(reaction.RH_STOIC())
    for i in range(len(rhss)):
        if rhss[i]=="1": rhss[i]=""

    rates = list(reaction.Rates())

    X = list(reaction.MHS())[0]

    l = "+".join(map(lambda x:"".join(x), zip(lhss, lhs)))
    r = "+".join(map(lambda x:"".join(x), zip(rhss, rhs)))
    
    r1 = "["+l+"+"+X+"->"+r+"+"+X+","+rates[0]+"]"
    r1 = Reaction(r1)
    # r1.Dump()
    return [r1]

#*******************************************************************************

def joinSpecies(A,B,by="_"):
    # intermediate species name generation "A" + "B" --> "A_B"
    # the "_" can be changed to any other string with keyword by
    #
    
    species = str(A) + str(by) + str(B)
    return species
    
    sA = str(A).split("{") # check for indices
    sB = str(B).split("{")
    
    nA = sA[0]
    nB = sB[0]
    index = ""
    if len(sA)>1:
        index = ((index+sA[1])[:-1])
    if len(sB)>1:
        if (len(index))>0: index = index+","
        index = ((index+sB[1])[:-1])
    if len(index)>0: index = "{"+index+"}"
    
    species = nA + str(by) + nB + index
    
    
    
    return(species)
#*******************************************************************************

def expand_mass_action_catalytic(r):
    # input: ParsedReaction for [ A=>B, mod[E], rates[k1,k2,k3,k4]] arrow
    #  
    # output: List of single componentsimple ParsedReaction objects with -> arrows
    #  [A + E -> A_E, k1]
    #  [A_E -> A + E, k2]
    #  [A + E -> B + E, k3]
    #  [B + E -> A_E, k4]
    # if any of the rate constants are zero, the corresponding reactions are omitted
    #
    S = (r.LHS())[0]
    P = (r.RHS())[0]
    E = (r.MHS())[0]
    XS = joinSpecies(S,E,"_")
    [k1,k2,k3,k4]=r.Rates()
    rnew = []
    if notZero(k1): rnew.append(reac([S,E],[XS],k1))
    if notZero(k2): rnew.append(reac([XS],[S,E],k2))
    if notZero(k3): rnew.append(reac([XS],[P,E],k3))
    if notZero(k4): rnew.append(reac([P,E],[XS], k4))
    rnew = map(Reaction, rnew)

    return rnew

#*******************************************************************************

def expand_mass_action_catalytic_bidirectional(r):
    # input: ParsedReaction for [ A<=>B, mod[E, R], rates[k1,k2,k3,k4,k5,k6,k7,k8]] arrow
    #   note: unlike cellerator, k4 must be specified if k5,k6,k7 are specified
    #  
    # output: List of single componentsimple ParsedReaction objects with -> arrows
    #
    # initially they are converted to [A => B, mod[E], rates[k1,k2,k3,k4] and
    #                                 [B => A, mod[R], rates[k5,k6,k7,k8]
    # and subsequently to:
    #  [A + E -> A_E, k1]
    #  [A_E -> A + E, k2]
    #  [A + E -> B + E, k3]
    #  [B + E -> A_E, k4]
    #  [B + R -> B_R, k5]
    #  [B_R -> B+R, k6]
    #  [B_R -> A + R, k7]
    #  [A + R -> B_R, k8]
    # if any of the rate constants are zero, the corresponding reactions are omitted
    #
    # print "expand_mass_action_catalytic_bidirectional: ", r.Input()
    S = list(r.LHS())[0]
    P = list(r.RHS())[0]
    catalysts = list(r.MHS())
    if len(catalysts) < 2:
        print "Error: reaction: "+r.Input()+" must have precisely two catalysts"\
        + " but "+str(len(catalysts))+" were found."
        raise SystemExit("Please correct input file and resubmit job.")
    [F,R] = catalysts 
    rates = list(r.Rates())
    r1 = enzreac(S,P,F,rates[:4])
    r2 = enzreac(P,S,R,rates[4:])
    r1 = Reaction(r1)
    r2 = Reaction(r2)
    rr = expand([r1, r2])
    
    return rr

#*******************************************************************************

def expand_mass_action_catalytic_sequential(r):
    # print "expand_mass_action_catalytic_sequential: ", r.Input()
    # r.Dump()

    S = list(r.LHS())[0]
    P = list(r.RHS())[0]
    E = list(r.MHS())[0]
    rates=list(r.Rates())
    if (len(rates)>6):
        print "Error: expand_mass_action_catalytic_sequential: there are "\
        + str(len(rates)) + " rates in the reaction " + str(r.Input()) + " but "\
        + "only 6 rates are allowed; the extra will be ignored."
        rates=rates[:6]
    [k1,k2,k3,k4,k5,k6] = rates
    ES = joinSpecies(S,E)
    EP = joinSpecies(P,E)
    r1 = revreacs([S,E],[ES], [k1,k2])
    r2 = revreacs([ES], [EP], [k3,k4])
    r3 = revreacs([EP],[P, E], [k5, k6])
    rr = r1+r2+r3
    
    rr = map(Reaction,rr)
    rr = expand(rr)
  
    return rr


#*******************************************************************************
def handle_two_way_equal_arrow(r):
    if r.ArrowType()=="Mass Action":
        return expand_mass_action_catalytic_bidirectional(r)
    elif r.ArrowType()=="kMech":
        expanded = expandkMech(r)
        
        # print "in expander: expanded: ", expanded
        
        expanded = map(Reaction,expanded)
        
        return expanded
    else:
        print "Error: unrecognized ArrowType in "+r.Input().strip()
        raise SystemExit(1)
    return

#*******************************************************************************


handlers={
      "->" :expand_no_expansion_necessary,
      "<->":expand_mass_action_bidirectional_reaction,
	  "-->":expand_mass_action_trivially_catalyzed,
      "=>" : expand_mass_action_catalytic,
      "<=>":handle_two_way_equal_arrow,
      # "<=>":expand_mass_action_catalytic_bidirectional,
      ":=>":expand_mass_action_catalytic_sequential
        }


#*******************************************************************************

ExpectedNumberOfRates={"<->":2,"=>":4,"<=>":8,":=>":6}
        
def expand_Cascade(reaction):

    DBG = False
    # print "Cascade reation (expansion not implemented): "+reaction.Input()
    
    casc = reaction.Arrow()
    stages,arrow=casc.split("*")
    stages = int(stages)
    arrowtype = reaction.ArrowType()
    if DBG: print "Cascade: arrow: ",stages," stages of ",arrow, " r=",reaction.Input()
    if DBG: print "DBG: expandCascade: arrowType:", arrowtype
    modifiers = reaction.MHS()
    #
    # Ensure that there are either zero modifiers or the same number of 
    # modifiers as there are arrows, by repeating the last modifier for each
    # arrow. Otherwise we have changed arrows in the middle of the cascade and 
    # this violates the definition of what we will allow as mixed arrows are
    # not allowed here
    #
    nmods = len(modifiers)
    # the <=> arrow has two modifiers
    # pair them up into a list ["M1,M2","M1,M2","M1,M2",...] if not enough
    # are given (using the last ones given" 
    if arrow == "<=>":
        firstmod=modifiers[-2:]
        while len(modifiers)<2*stages:
            modifiers = modifiers + firstmod
        modifiers = [modifiers[i:i+2] for i in range(0,len(modifiers),2)] 
        modifiers = map(lambda x:str(x[0]+","+x[1]),modifiers) 
    # all the other cascadable arrows have at most one modifier
    # otherwise just expand the list with the last modifier given
    elif (nmods>0) and (nmods <stages):
        firstmod = modifiers[-1]
        while len(modifiers)<stages:
            modifiers.append(firstmod)
    # print "modifiers = ", modifiers
    rates = reaction.Rates()
    #
    # put code here to expand rates to cascade    
    #
    
    nrates = len(rates)
    if DBG: print "DBG: ExpandCascade: nrates=", nrates, " rates=", rates
    if arrow in ["|->","|-->"]:
        dash=arrowtype.index("-")
        name=arrowtype[:dash]
        if DBG: print "DBG: ExpandCascade: GRN '"+name+"'"
        r1 =name+"["+(",".join(rates))+"]"
        rates=[]
        while len(rates)<stages: rates.append(r1)
        if DBG: print "DBG: ExpandCascade: rates=",rates
        
    elif nrates>0:
        if arrow in ExpectedNumberOfRates:
            nrexpected = ExpectedNumberOfRates[arrow]
        else:
            nrexpected = 1
        if DBG: print "DBG: ExpandCascade: nrexpected:", nrexpected
        while len(rates)<nrexpected: rates.append('0')  # zero-fill
        lastset = rates[-nrexpected:]
        while len(rates)< stages*nrexpected:
            rates = rates + lastset
        lr = len(rates)
        rates = [rates[i:i+nrexpected] for i in range(0,lr,nrexpected)] 
        if DBG: print "DBG: ExpandCascade: rates: ", rates
        rates = map(lambda x: "rates["+(",".join(x))+"]", rates)
        if DBG: print "DBG: ExpandCascade: rates: ", rates
       
    #print "Cascade: modifiers: ", modifiers, " rates=", rates
    
    
    
    S = reaction.LHS()
    P = reaction.RHS()
    SP = flatten([S,P])
    SPpairs = zip(SP, SP[1:])
    reactions = map(lambda x: x[0]+arrow+x[1], SPpairs)
    if nmods>0:       
        reactions = zip(reactions, modifiers)
        reactions = map(lambda x: x[0]+", mod["+x[1]+"]", reactions)
    #classification = (reaction.ArrowType()).split("-")
    # print "Cascade: classification", classification
    #if len(classification)>1:
    #    classification = classification[0]
    #    rates = map(lambda x: classification + "rates["+x+"]", rates)        


        
    rr = zip(reactions,rates)    
        
    reactions = map(lambda x:"["+x[0]+","+x[1]+"]" , rr)
         
    # print "Cascade: rates = ", rates
    
    # print "Cascade: reactions = ", reactions
    
    rr = map(Reaction, reactions)

    rr = expand(rr)
    
    
    return rr     

#*******************************************************************************

def expand(reaction):
    #
    #   if a list of reactions is given, apply to each reaction separately
    #   and then concatenate the final lists together and return the 
    #   concatenated results
    #
    # print "expand: type:", type(reaction)
    if type(reaction)==type([]):
        rlist = map(expand,reaction)
        result = []
        for r in rlist: result += r
        return result
    #
    #   otherwise just process a single reaction
    #
    arrow=reaction.Arrow()
    arrowtype=reaction.ArrowType()
    # print "expand: arrow: ", arrow
    if arrow in handlers:
        expander = handlers[arrow]
        r = expander(reaction)
        # map(lambda x:x.Dump(), r)
        return r
    elif arrowtype.endswith("Cascade"):
        # print "************* CASCADE "+reaction.Input()
        r = expand_Cascade(reaction)          
        return r
    else:
        return expand_no_expansion_necessary(reaction)
    
    
def expandReactions(r, dump=False, text=False):
    inputfile=r
    if ("-in" in sys.argv):
        i = (sys.argv).index("-in")+1
        if i<len(sys.argv):
            inputfile=sys.argv[i]
        else:
            raise SystemExit("Error: expecting filename after -in")
    elif inputfile=="":
        raise SystemExit("No Input File Given; use -in option.")
    if not os.path.isfile(inputfile):
        raise SystemExit("Error: "+os.path.abspath(inputfile)+" File not found.")
        
    if os.path.splitext(inputfile)[1]==".model":
        (lines, ics, rates, frozenvars, functions,assignments, fullfile)=readmodel(INFILE=inputfile)
    else:
        f=open(inputfile)
        lines = f.readlines()
        f.close()
    expanded=[]
    expandedInput=[]
    
    
    for line in lines:
        line = StripComments(line,"#")
        if line=="": continue # ignore blanks
        p = Reaction(line)
        # print "Expander: input: ", line
        p = expand(p)
        # print "Expander: expanded: ", p
        expanded += p
        tp = map(lambda x:x.Input(),p)
        expandedInput += tp

    # if requested, dump the expanded text reactions to the screen
    if ("-dump" in sys.argv): dump=True
    if dump:
        for line in expandedInput:
            print line 
    # if requested, dump the expanded text reactions to a file
    # check to see if file is requested

    if ("-out" in sys.argv):
        i = (sys.argv).index("-out")+1
        if i<len(sys.argv):
            outputfile=sys.argv[i]
        else:
            # print "Error: expecting filename after -out"
            outputfile = "expanded-reactions.out"

        outputfile = uniqueFileName(outputfile, type ="out")
        
        f=open(outputfile,"w")
        for line in expandedInput:
             f.write(line+"\n")
        f.close()
        print "Output written to: "+os.path.abspath(outputfile)         
    
    if text == True:
        return expandedInput
      
    return expanded 
    
def expander():
    if not ("-quiet" in sys.argv):
        d=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print "py[cellerator]: expander ("+d+")"
    istext = "-text" in sys.argv
    rs = expandReactions("", text=istext)
    if "-quiet" in sys.argv:
        return
    for r in rs:
        if istext:
            print r
        else:
            print r.Input() 
       
if __name__ == "__main__":
    expander()
    
 
