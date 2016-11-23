from utils import *
from sys import exit
from reader import *
import interpreter
import expander
from sympy import *
from parser import Reaction
#
#
# SBML Writer -- does not handle complete output
# still to do: 

def determinefiles(args):
    if "-in" in args:
          i=args.index("-in")+1
          if i < len(args):
            inputfile = args[i]
          else:
            exit("Error: an input model must be specified.")
    else:
        exit("Error: an input model must be specified.")
    if "-out" in args:
          i=args.index("-out")+1
          if i < len(args):
            outfile = args[i]
          else:
            outfile = uniqueFileName("SBMLFILE", type="xml")
    else:
        outfile  = uniqueFileName("SBMLFILE", type="xml")

    return(inputfile, outfile)
try:    
	from libsbml import *
except:
	print "Cellerator.sbmlwrite: unable to load libsbml. Will not be able write SBML files."
	
def reduceStoichiometry(species, stoichiometries):
    d={}
    s=list(species)
    for X,Z in zip(species,stoichiometries):
        if X in d:
            d[X] += float(Z)
        else:
            d[X] = float(Z)
    snew = list(d)
    stnew = [str(d[X]) for X in snew]

    return (snew, stnew)
    


def genmodel(infile, outfile, verbose=False):
    (reactions, ic, rates, frozenvars, functions, assignments, filename) = readmodel(infile)
    # Expand reactions
    #
    r = interpreter.invokeParser(reactions, dump=False)
    er = expander.expand(r)
    s = interpreter.makeSymbolDictionary(er, rates, ic)
    
    try:
        SBML=SBMLDocument(3,1)
    except ValueError:
        exit("Could not create the SBMLDocument.")
    m = SBML.createModel()
    m.setTimeUnits("dimensionless")
    m.setSubstanceUnits("item")
    m.setExtentUnits("item")
    c = m.createCompartment()
    c.setId('compartment')
    c.setConstant(True)
    c.setSize(1)
    c.setSpatialDimensions(3)
    c.setUnits("dimensionless")
    
    for function in functions:
        if verbose:
            print "function: ", function
        fnew=m.createFunctionDefinition()
        fname, formula = function.split("=")
        fname=fname.strip()
        formula=formula.strip()
        formula=formula.replace("**","^")
        functionvar, functiondef=formula.split(":")
        # print "functiondef:", functiondef
        functionvar=functionvar.split()[-1]
        # print "functionvar:", functionvar
        lambdadef="lambda("+functionvar+","+functiondef+")"
        # print "lambdadef: ", lambdadef
        math_ast = parseL3Formula(lambdadef)
    
        
        # print "math_ast=",math_ast
        fnew.setId(fname)
        fnew.setMath(math_ast)
    
    species=list(ic)
     
    sp=[]
    for specie in species:
        snew=m.createSpecies()
        snew.setId(specie)
        if specie in frozenvars:
            snew.setConstant(True)
        else:
            snew.setConstant(False)
        snew.setBoundaryCondition(False)
        snew.setHasOnlySubstanceUnits(True)
        amount = ic[specie]
        snew.setInitialAmount(amount)
        snew.setCompartment("compartment")
        sp.append(snew)
    
    # add nil species
    snew=m.createSpecies()
    snew.setId("Nil")
    snew.setConstant(False)
    snew.setBoundaryCondition(False)
    snew.setHasOnlySubstanceUnits(True)
    snew.setInitialAmount(0)
    snew.setCompartment("compartment")
    sp.append(snew)    
        
    RATES=list(rates)
    ks=[]
    for rate in rates:
        k = m.createParameter()
        k.setId(rate)
        k.setConstant(True)
        k.setValue(rates[rate])
        k.setUnits('dimensionless')
        ks.append(k)
    
    #print assignments
    
    i=1
    for assignment in assignments:
        lhs, rhs = assignment.split("=")
        lhs=lhs.strip()
        rhs=rhs.strip()
        if verbose:
            print "assignment: to:", lhs, "from:",rhs
        arule = m.createAssignmentRule()
        patchedlaw=rhs.replace("**","^")    
        math_ast = parseL3Formula(patchedlaw)
        arule.setMath(math_ast)
        arule.setVariable(lhs)
        aid = "arule"+str(i)
        i+=1
        arule.setId(aid)
    
    reacts=[]
    i = 1
    for line in er:
       new_reaction=m.createReaction()
       rid = "r"+str(i)
       new_reaction.setId(rid)
       new_reaction.setReversible(False)
       new_reaction.setFast(False)
       odeterm = interpreter.makeODETerms([line], s, frozen=frozenvars)
       termkeys=list(odeterm)
       # print "termkeys:", termkeys
  
       # print "**** Reaction ", i,"\n", line.LHS(),"\n",line.MHS(),"\n",line.RHS()
  
       S, ST = reduceStoichiometry(line.LHS(), line.LH_STOIC())
       # print "reactant: ", S, ST
       for X,Z in zip(S, ST):
          if X in termkeys:
              sref=new_reaction.createReactant()
              sref.setSpecies(X)
              sref.setStoichiometry(float(Z))
              if X in frozenvars:
                  sref.setConstant(True)
              else:
                  sref.setConstant(False) 
          #
          # some Cellerator "Reactants" are really modifiers
          #
          elif str(X) != "Nil":
              sref=new_reaction.createModifier()
              sref.setSpecies(X)

       S, ST = reduceStoichiometry(line.RHS(), line.RH_STOIC())
       # print "product: ", S, ST
       for X,Z in zip(S, ST):
          sref=new_reaction.createProduct()
          sref.setSpecies(X)
          sref.setStoichiometry(float(Z))
          if X in frozenvars:
            sref.setConstant(True)
          else:
            sref.setConstant(False)
       for X in line.MHS():
          sref=new_reaction.createModifier()
          sref.setSpecies(X)
       
       terms = [str(odeterm[X]) for X in termkeys]
       # print "terms: ", terms
       terms = map(lambda x: x.lstrip("-"), terms)
       terms = set(terms) - set(["0"])
       if len(terms)>0:
          law = list(terms)[0]
       else:
          law = "0"
       if verbose:
           print "using kinetic law ", law, " for "+str(rid)+": ", line.Input()
       
       patchedlaw=law.replace("**","^")
       
       math_ast = parseL3Formula(patchedlaw)
       kinetic_law = new_reaction.createKineticLaw()
       kinetic_law.setMath(math_ast)
       # print "kinetic law ", i, " is ", law, odeterm
       i +=1 
       

    return (writeSBMLToString(SBML))
 
def SBMLWRITE(args):
    # print "SBMLWRITE: ", args
    infile, outfile = determinefiles(args)
    # print "inputfile = ", infile
    # print "outfile =   ", outfile
    s = genmodel(infile,outfile)
    f = open(outfile, "w")
    f.write(s)
    f.close()
    return "model "+infile+" written to "+outfile



