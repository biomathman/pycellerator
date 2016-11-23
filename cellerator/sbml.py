#
#****************************************************************************
#    Cellerator converts reactions, expressed in a text-formatted arrow-based 
#    notation into differential equations and performs numerical simulations.
#
#    This module reads and writes SBML files
#
#****************************************************************************
#
#    Copyright (C) 2015 Bruce E Shapiro
#    last rev. 8/30/2015
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
from sympy import *
from utils import *
from math import *
try:
    from libsbml import *
except:
    print "Warning: Cellerator.sbml: unable to import libsbml"
    

def fstring(rval, digits=14):
    if int(rval)==rval:
        return str(rval)
    x=abs(rval)
    if x<1:
        digs = abs(int(log(x,10)))+digits
    else:
        digs = digits
    s = '{0:.{precision}f}'.format(rval,precision=digs)
    # remove trailing zeroes in decimal part
    if '.' in s: s = s.rstrip('0')
    # remove decimal part if it turns out to be integer (still)
    if s[-1]=='.': s = s.rstrip('.')
    return s
    
#****************************************************************************
def make_formula(rlist):
    #
    # input: reactant dictionary list [ {"stoichiometry", "species"}, {"stoichiometry", "species"}, ...]
    # returns: LHS or RHS formula s1*X1 + s2*X2 + ...
    #
    TEENY=1e-10
    def intok(u):
        if abs(int(u)-u)<TEENY:
            return True
        else:
            return False
    def makeint(u):
        if intok(u):
            return int(u)
        else:
            return u
    stoic = map(lambda x:x["stoichiometry"],rlist)
    stoic = map(makeint, stoic)
    stoic = map(fstring, stoic)
    # species = map(Symbol, map(lambda x:x["species"], rlist))
    species = map(lambda x:x["species"], rlist)
    SSP = zip(stoic,species)
    formula = 'Nil'
    for (st,sp) in SSP:
        if st=='1':
            formula = formula + " + " + sp
        else:
            formula = formula + " + " + st + " " + sp
    # remove Nil but only if it is followed by a + sign, not if it is alone
    formula = formula.replace('Nil +', '')
    return formula
    
def getValue(parameter,parameters, default):
    #
    # given a parameter dictionary looks to see if the value of the parameter
    # is in the dictionary, and returns that value, unless the flag "None" has
    # been placed there, in which case the default value is used
    #
    if parameter in parameters:
        val = str(parameters[parameter])
        if val == "None":
            val = default
    else:
        val = default
    return val

def COBRA_REACTION(reaction):
    #
    #   formats a Cellerator reaction from an SBML file
    #   assumes a COBRA-FLUX mile format file
    #
    id = reaction["id"]
    R = reaction["reactants"]
    P = reaction["products"]
    parameters = reaction["parameters"]
    LHS = make_formula(R)
    RHS = make_formula(P)
    upper = getValue("UPPER_BOUND", parameters, "inf")
    lower = getValue("LOWER_BOUND", parameters, "-inf")
    obj = getValue("OBJECTIVE_COEFFICIENT", parameters, "0")
    val = getValue("FLUX_VALUE", parameters, "0")   
    reaction =" [ "+LHS+" -> "+RHS+ ", Flux["+lower+" < " + id + " < " + upper + ", "+obj+" ," + val+"] ]"

    return reaction
#****************************************************************************

def writemodel(filename, functions, compartments, species, parameters,\
    rules, reactions, COBRA):
    
    DEBUG=False
   

    #***************************************************
    def speciesToStoichs(r):
        # R is LHS speces
        R = r['reactants']
        RS = map(lambda x:x['species'],R)
        RST = map(lambda x:x['stoichiometry'],R)
        # P is RHS species
        P = r['products']
        PS = map(lambda x:x['species'], P)
        PST = map(lambda x:x['stoichiometry'], P)
        # S is all species
        S = list(set(RS+PS))
        ST = []
        
        if DEBUG:
			print "speciesToStoics: R: ",R
			print "speciesToStoics: RS: ",RS
			print "speciesToStoics: PST: ",PST
			
		#
		# add check for non-numeric stoich. 7/1/16
		#	
        for s in S:
			try:
				st = 0; stdelta=-1
				if s in RS: st = st - RST[RS.index(s)] 
				stdelta=1  
				if s in PS: st = st + PST[PS.index(s)]
				
			except: 
				st=stdelta
				print "Warning:  unable to convert non-numeric stoichiometry in reaction",\
					nreactions,"This may produce unexpected results in the model."
			ST.append(st) 
        return (S, ST)
    #****************************************************  
    def getSpeciesCompartment(s):
        c = "UNKNOWN"
        for S in species:
            if S['id']==s:
                c = S['compartment']
                break
        return c
    #****************************************************
    output = []
    rcounter = 0
    acounter = len(reactions)
    if len(rules)>0 or (len(reactions)>0 and not COBRA):
        algcounter=0
        for r in rules:
            if r['type'] == 'Rate':    # these will become "using" reactions
                rcounter += 1
            elif r['type'] == 'Assignment':
                acounter += 1
            elif r['type'] == 'Algebraic':
                algcounter += 1
            else:
                print "Unknown rule type: ", r
        if acounter>0:
            output.append("$ASSIGNMENTS")
            for r in rules:
                if r['type'] == 'Assignment':
                    V = r['variable']
                    maths = sympify(r['maths'])
                    # print parameters
                    for p in parameters:
                        x = Symbol(p)
                        v = parameters[p]
                        maths = maths.subs(x,v)
                    assignment =" " + V + "=" + str(maths)
                    output.append(assignment)
            for r in reactions:
                KL = sympify(r['kineticLaw'])
                for p in parameters:
                    x = Symbol(p)
                    v = parameters[p]
                    KL = KL.subs(x,v)
                rid = r['id']
                KL = " " + rid + "=" + str(KL)
                output.append(KL)
        if algcounter>0:
            print "Warning: The SBML model contains "+str(algcounter)+\
                " algebraic rules, which are not supported by this version of "\
                + " py[cellerator] and will be ignored."
    
    
    if len(reactions)>0 or ((len(rules)>0) and (rcounter>0)):
        output.append("$REACTIONS")
        nreactions=0
        
        for r in reactions:
            nreactions += 1
            
            if DEBUG:
                print "************* SBML reaction:", nreactions, "\n", r
           
            if COBRA:
                reaction=COBRA_REACTION(r)
                output.append(reaction)
            else:
                (S, ST) = speciesToStoichs(r)
                reactionID = r['id']
                M = r['modifiers']  # not actually used !! 
                KL = reactionID
                for (s,st) in zip(S,ST):
                    
                    if DEBUG:
                        print "(s,st)=", (s,st)
                        
                    if s=='EmptySet': continue
                    if s=="Nil":continue
                    if is_number(st):
                        if st>0:
                            reaction = "Nil -> "+str(s)
                        elif st < 0:
                            reaction = str(s) + " -> Nil"
                        else:
                            continue
                        st = abs(st)
                        intst = round(st,0)
                        errst = abs(intst - st)
                        small = 10**-7
                        if errst < small:
                            if intst == 1:
                                st = ""
                            else:
                                st = str(intst)+"*"
                        else:
                            st = str(st)+"*"

                    else:
                        st = str(st)+"*"
                        reaction = "Nil -> " + str(s) 
                    
                    compartment = getSpeciesCompartment(s)
                    kl = st + KL + "/" + compartment
                    reaction=" [" + reaction + ", using[\""+kl+"\"]]"
                     
                    if DEBUG:  
                        print "Converted reaction:\n", reaction
                    
                    output.append(reaction)            
 
            
        if  rcounter > 0:
            for r in rules:
                if r['type']=='Rate':
                    V = r['variable']
                    maths = r['maths']
                    reaction = " [Nil -> "+V+", using[\""+maths+"\"]]"
                    output.append(reaction)  

    if len(parameters)>0 or len(compartments)>0:
        output.append("$RATES")
        for p in parameters:
            output.append(" "+str(p)+"="+str(parameters[p]))
        # need to fix the following to allow for compartments
        # to change size in rules
        for c in compartments:
            id = str(c['id'])
            size = str(c['size'])
            output.append(" "+id+"="+size)

    if len(species)>0:
        output.append("$IC")
        frozen=[]
        for s in species:
            specie = s["id"]
            if specie != "EmptySet":
                ic = s["initialvalue"]
                output.append(" "+str(specie)+"="+str(ic)) 
            if s["BC"] or s["constant"]: frozen.append(specie)
        frozen = list(set(frozen)-set(["EmptySet"]))
        if len(frozen)>0:
            output.append("$FROZEN")
            
            for f in frozen:
                output.append(" "+f)

    if len(functions)>0:
        output.append("$FUNCTIONS")
        for f in functions: 
            args, fdef = (functions[f]).split(":")
            args = filter(lambda x:len(x)>0, (args.strip()).split("lambda"))[0].strip()
            fdef = fdef.strip()
            output.append(" "+str(f)+"("+args+")="+fdef)     
         


    output = '\n'.join(output)
    f = open(filename,"w")
    f.writelines(output)
    f.close()
    return  


#****************************************************************************
   

def readfile(filename):
    reader = SBMLReader()
    SBML = reader.readSBML(filename)
    nerrors = SBML.getNumErrors()
    if nerrors > 0:
        # errorLog = SBML.getErrorLog()
        SBML.printErrors()
        print "************************************************************"
        print " WARNING: Because of the failure to validate the SBML FILE,"
        print " unexpected behavior may occur with this model."
        print "************************************************************"
        

    level = SBML.getLevel()
    version = SBML.getVersion()
    return (SBML, level, version)
    
#****************************************************  
    
def getCompartments(model):
    LOC = model.getListOfCompartments()
    n = len(LOC)
    compartments=[]
    for i in range (n): 
        c = model.getCompartment(i)
        comp = \
        {"id":c.getId(),
         "name":c.getName(),
         "size":c.getSize(),
         "dim":c.getSpatialDimensions(),
         "units":c.getUnits()
        }
        compartments.append(comp)
    return compartments

#****************************************************  
def getRules(model):
    LOR = model.getListOfRules()
    n = len(LOR)
    rules=[]
    for i in range(n):
        r = model.getRule(i)
        v = r.getVariable()
        m = r.getMath()
        maths = formulaToString(m)

        
        if r.isAssignment(): 
            t = "Assignment" 
        elif r.isRate():
            t = "Rate"
        elif r.isAlgebraic():
            t = "Algebraic"
        else:
            t = "Error"
        rule ={"variable":v,
               "maths":maths,
               "type":t
               }
        rules.append(rule)
    
    return rules
    
#****************************************************  
    
def getSpecies(model):
    LOS = model.getListOfSpecies()
    n = len(LOS)
    species = []
    for i in range(n):
        s = model.getSpecies(i)
        hasOnlySubstanceUnits = s.getHasOnlySubstanceUnits()
        # this assumes that the spatial dimension of the compartment is nonzero
        IA = s.getInitialAmount()
        IC = s.getInitialConcentration()
        ID = s.getId()
        if hasOnlySubstanceUnits or IA !=0:
            IV = IA
        else:
            IV = IC
            
        if isnan(IV):  # if there is a bad model put something there 4.25.15
            IV=0.0
       
        specie=\
        {"id":ID,
         "name":s.getName(),
         "initialvalue":IV,
         "constant":s.getConstant(),
         "BC":s.getBoundaryCondition(),
         "compartment":s.getCompartment()       
        }
        species.append(specie)    
    
    return species
#****************************************************  

def getParameters(model):
    LOP = model.getListOfParameters()
    n = len(LOP)
    parameters=[]
    for i in range(n):
        p = model.getParameter(i)
        PV=p.getValue()
        if isnan(PV):
            PV=0.0 # put a zero there if no value
        parameter=\
        {"id":p.getId(),
         "name":p.getName(),
         "value":PV,
         "units":p.getUnits()
         }
        parameters.append(parameter)
        
    plist = {}
    for p in parameters:
        plist[p['id']]=p['value']    
    return plist

#****************************************************  
def getSpeciesReferences(listOfSpeciesReferences):
    #
    # returns a list of species references as a list of dictionaries 
    #
    n = len(listOfSpeciesReferences)
    srefs=[]
    for j in range(n):
        speciesref=listOfSpeciesReferences.get(j)
        stoich = speciesref.getStoichiometry()
        stoichMath=speciesref.getStoichiometryMath()
        if str(stoichMath).strip() != "None":
            math = stoichMath.getMath()
            maths = formulaToString(math)
            stoich=maths
        species = {"species":speciesref.getSpecies(),
            "stoichiometry":stoich
             }
        srefs.append(species)
    return srefs
#**************************************************** 
def getModifierSpeciesReferences(listOfModifierSpeciesReferences):
    #
    # returns the modifier species references as a list of species names
    #
    n = len(listOfModifierSpeciesReferences)
    srefs=[]
    for j in range(n):
        speciesref=listOfModifierSpeciesReferences.get(j)
        species = speciesref.getSpecies()
        srefs.append(species)
    return srefs
    
#**************************************************** 

def replaceTerms(formula, parameterdict):
		
    #print "DBG: formula: ", formula, " parameterdict: ", parameterdict
    #print "DBG: formula is ", type(formula)
    # 
    # replace all terms in string formula with 
    # values given in parameter dictionary
    #
    #print "DBG sympify:", sympify(formula)
    f = sympify(formula)
    #print "DBG: f is", type(f)
    for x in parameterdict:
        val = parameterdict[x]
       
        f = f.subs(x,val)
    formula = str(f)    

    return formula
    
#**************************************************** 
#
#     
def getTheKineticLaw(k):
    Level = k.getLevel()
    if Level<3:
        lop = k.getListOfParameters()
        n = k.getNumParameters()
    else:
        lop = k.getListOfLocalParameters()
        n = k.getNumLocalParameters()
    parameters={}

    for i in range(n):
        if Level<3:
            parameter = k.getParameter(i)
        else:
            parameter = k.getLocalParameter(i)
        parameterValueIsSet = parameter.isSetValue()
        # print "DBG: getTheKineticLaw: ", i, parameter, type(parameter)
        ID = parameter.getId()
        value=parameter.getValue()
        if parameterValueIsSet:
            parameters[ID]=value
        else:
            parameters[ID]="None"
    math = k.getMath()
    maths = formulaToString(math)
    #
    # replace all local parameters with values
    # since py[cellerator] does not have "local" parameters
    #
    
    #print "DBG: math = ", math
    # print "DBG: maths=", maths, " parameters=", parameters
    
    maths = replaceTerms(maths, parameters)
    
    klaw = maths   

    # print "local parameters: ", parameters
    return (klaw, parameters)

#**************************************************** 

def getReactions(model):
   
    LOR = model.getListOfReactions()
    n = len(LOR)
    reactions=[]
    for i in range(n):
        r = model.getReaction(i)
        # print "DBG: ", r
        rev = r.getReversible()
        reactants = getSpeciesReferences(r.getListOfReactants())
        products  = getSpeciesReferences(r.getListOfProducts())
        modifiers = getModifierSpeciesReferences(r.getListOfModifiers())
        reactionid = r.getId()
        tkl = r.getKineticLaw()
        # print "DBG:", tkl
        klaw,localparameters = getTheKineticLaw(tkl)
        
        reactions.append({"reactants":reactants, "products":products, \
                          "modifiers":modifiers, "kineticLaw":klaw, \
                          "reversible":rev, "id":reactionid, "parameters":localparameters})
        
    return reactions
#****************************************************  

def find_nth_comma(haystack, n):
# modified from 
# http://stackoverflow.com/questions/1883980/find-the-nth-occurrence-of-substring-in-a-string
    start = haystack.find(",")
    while start >= 0 and n > 1:
        start = haystack.find(",", start+1)
        n -= 1
    return start

def getFunctions(model):
    n = model.getNumFunctionDefinitions()
    functions = {}
    for i in range(n):        
        f = model.getFunctionDefinition(i)
        fid = f.getId()
        math = f.getMath()
        maths = formulaToString(math)
        #
        # this returns lambda(arg, arg, arg,..., arg, formula) 
        # convert to normal python
        #
        nargs = f.getNumArguments()
        j = find_nth_comma(maths,nargs)
        if j<0: 
            print "Error: getFunctions: libSBML says there are ",nargs, " arguments",\
                " but there are fewer commas in " + maths
        else: 
            
            LP = maths.find("(")
            #print "LP=", LP
            #print "nargs=",nargs
            #print "j=", j
            #print "maths=",maths
            arglist = maths[LP+1:j].strip()
            #print "arglist=",arglist
            if arglist[-1]==",": arglist=arglist[:-1]
            formula = maths[j+1:-1].strip()
            # print "arglist = ", arglist
            # print "formula = ", formula
            maths = "lambda "+arglist+" : "+formula
            
        
        functions[str(fid)]=maths
        # print "function: ", f, fid, " math = ", math, " maths=", maths
        
    return functions
    
#****************************************************  
def SBMLREAD(args, INFILE="", OUTFILE=""):
    DBG = False
    if DBG: print "DBG: SBMLREAD: ", args
    ARGS = map(lambda x: x.upper(), args)
    
        

    if OUTFILE=="":
        if "-MODEL" in ARGS:
            i = (ARGS).index("-MODEL")+1
            if i<len(args):
                outfile = args[i] # use one that is not nec. all uppercase
            else:
                raise SystemExit("Error: Cellerator.sbml.SBMLREAD: expecting filename after -model")
        else:
            outfile="tmp.model"
    else:
        outfile=OUTFILE
    outfile = uniqueFileName(outfile, type="model")
    outfile = os.path.abspath(outfile)
    
    COBRA = True if "-COBRA" in ARGS else False
    
 
 
    if INFILE=="":
        if "-in" in args:
            i = (args).index("-in")+1
            if i<len(args):
                infile=args[i]
            else:
                raise SystemExit("Error: Cellerator.sbml.SBMLREAD: expecting filename after -in")
        else:
            raise SystemExit("Error: Cellerator.sbml.SBMLREAD: No Input File Given; use -in option.")
    else:
        infile=INFILE
        
    if not os.path.isfile(infile):
        raise SystemExit("Error (Cellerator.sbml.SBMLREAD): "+os.path.abspath(infile)+" File not found.")
    fullfile = os.path.abspath(infile)  
    if DBG: print "DBG: SBMLREAD: input file ="+fullfile
    SBML, L, V=readfile(fullfile)
    if DBG: print "DBG: SBMLREAD: L, V = ", L, V
    model = SBML.getModel()
    modelid = model.getId()
    modelname = model.getName()
    if DBG: print "DBG: SBMLREAD: model = ", modelid, modelname
    functions = getFunctions(model)
    compartments = getCompartments(model) 
    species = getSpecies(model)
    parameters=getParameters(model)
    rules = getRules(model)
    reactions = getReactions(model)
    
    if DBG:
        for f in functions:
            print "DBG: SBMLREAD: function ",f," = ", functions[f]
        for c in compartments: 
            print "DBG: SBMLREAD: compartment:", c
        for s in species:
            print "DBG: SBMLREAD: species: ", s 
        # for p in parameters:
        print "DBG: SBMLREAD: parameters: ", parameters
        print "DBG: SBMLREAD: rules: ", rules
        for r in reactions:
            print "DBG: SBMLREAD: reaction:", r

    writemodel(outfile, functions, compartments, species, parameters,\
        rules, reactions, COBRA)

    return outfile

#****************************************************    

from sbmlwrite import *
#****************************************************    


def run(args):
    if type(args)==type(""):
        arguments = filter((lambda x:len(x)>0), args.split(" "))

        sys.argv=[]
        for arg in arguments:
            sys.argv.append(arg)
    else:
        arguments = args 

    ARGS = map(lambda x:x.upper(), arguments)
    
    handlers = {'READ':SBMLREAD,
                'WRITE':SBMLWRITE}

    OK = False
    result = "Error Exit"
    for option in handlers:
        if option in ARGS:
            OK = True
            f = handlers[option]
            result = f(arguments)
            break
    if not OK:       
        print "Error: Cellerator.sbml: no option specified. Options are: ",
        for opt in handlers: 
            print opt+" ",
        print ""
        
    
    return result
    
def runsbml():
    result = run(sys.argv)
    print result
    return
