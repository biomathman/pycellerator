#!/usr/bin/env python
#
# B.E.Shapiro 24 Feb 2012
# solver module for Cellerator
# 11/11/14 add parameter for mxiterations 
# revised 4/1-4/24/15 various B.E.S.
# revised 8/30/15 for Github release and to reflect new direc. structure.
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
#
# ****************************************************************************
# Cellerator-like solver
# ****************************************************************************
# invoking from command line (in Linux):
#
# pythonnsolver.py -in filename -run duration stepsize
#
# Entry Point:
#
# runSimulation(filename, duration=value, step=value)
#
import os.path, sys, datetime
from parser import Reaction
from sympy import *
import interpreter, utils
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from expander import expand
from reader import readmodel
import sys
import re
import time

#*******************************************************************************  

def applyReplacementRules(expression, rules, symboltable):
    #
    # Zero is added to each ode because some odes may be equal to
    # zero and a number does not have an attribute of .subs(...)
    # so we convert the expression to symbolic to force it to be
    # susceptible to .subs(symbol, value)
    #
    # print "expression:", expression
    # print "rules: ", rules
    # print "symboltable: ", symboltable
    ZERO = Symbol("ZERO")
    expr=expression+ZERO
    vlist=list(rules)
    for v in vlist:
        # print "applyReplacementRules: expr:", expr, "v:", v, " rules[v]:", rules[v]
        symbol=symboltable[v]
        value=rules[v]
        if utils.is_number(value):
            if round(value,0)==value:
                value = int(round(value,0))
        expr = expr.subs(symbol,value)
        #expr = expr.replace(symbol, value)
        # print "expr: ", expr, type(symbol), type(value), type(expr)
    expr = expr.subs(ZERO,0)
    #expr = expr.replace("ZERO","0")
    return expr

#*******************************************************************************  

def substituteRateConstants(odeterms, rates, symboltable, ratesToLeaveAsIs):
    """
    input: odeterms: dictionary of odeterms by string symbol such as 
            {'Ce': -1.0*Ce*k5 + 2.0*BrO3*HBrO2*k3,k3*Br*Ce**2, ...}
           rates: dictionary of rate values by string symbol such as
            {'k3': 34.0, 'k2': 2000000.0, 'k1': 1.3, 'k5': 0.02, 'k4': 3000.0}
           symboltable: dictionary of symbols by string {"x":x,"y":y,..}
    output: new dictionary of odeterms
    """
    vlist = list(odeterms)
    newterms = {}

    modifiedrates = {}
    for rate in rates:
        if rate not in ratesToLeaveAsIs:
            modifiedrates[rate]=rates[rate] 
    #print rates
    #print modifiedrates
    #sys.exit()
    
    #print "rates: ", rates
    #print "symboltable: ", symboltable
    for v in vlist:
        expr = odeterms[v]
        #print v," ", expr, "(before)"
        expr = applyReplacementRules(expr, modifiedrates, symboltable)
        #print v,":", expr, "(after)"
        newterms[v] = expr
        
    #print newterms    
    #sys.exit()
    return newterms

#*******************************************************************************  

def makeODEfunc(odeterms, symbolDictionary,ics):
    # print "symbolDictionary in makeODEfunc: ", symbolDictionary
    odelist=[]
    iclist=[]
    varlist=[]
    vlist = list(odeterms)
    for v in vlist:
        # print "v=",v, " type=", type(v)
        oderhs = odeterms[v]
        odelist.append(oderhs)
        variable = symbolDictionary[v]
        varlist.append(variable)
        if v in ics:
            ic = ics[v]
        else:
            ic = '0'  # there may be auto-generated variables
        iclist.append(ic)
        
    # print "iclist in makeODEfunc: ", iclist
    return (varlist, odelist, iclist)
    
#*******************************************************************************  

def makeODEs(odeterms, symbolDictionary):
    # print "symbolDictionary in makeODES: ", symbolDictionary
    odelist=[]
    varlist=[]
    vlist = list(odeterms)
    for v in vlist:
        oderhs = odeterms[v]
        odelist.append(oderhs)
        variable = symbolDictionary[v]
        varlist.append(variable)
    return (varlist, odelist)
    
 #*******************************************************************************  
   
def failIfNegative(x):
    if x<0:
        raise SystemExit("Unexpected negative value "+str(x))
       
#*******************************************************************************  
    
def generatePythonFunction(reactions,  rates, ics, frozenvars, functions, 
    assignments, holdrates=[], substituteRates=False):
    """
    input: reactions - list of reactions in cellerator text form
           inputrates - list of rate constants in dictionary form
    output: writes a function to file tmp.py that is solver-compatible
        return values: (y, y0)
        y: list of string names of variables ["y1", "y2",...]
        y0: list of values of variables at initial time
    """    
  
    
    d=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    codefile=utils.uniqueFileName("tmp.py")
    
    f=open(codefile,"w")
    f.write("from math import *\n")
    if len(functions)>0:
        f.write("\n")
        for funct in functions: f.write(funct+"\n")
        f.write("\n")
    f.write("def ode_function_rhs(y,t):\n")
    f.write("  # \n")
    f.write("  # this odeint(..) compatible function was \n")
    f.write("  # automatically generated by Cellerator "+d+" \n")
    ver = sys.version.replace("\n"," ")
    f.write("  # " + ver + "\n")
    f.write("  # " + sys.platform + "\n")
    f.write("  # \n")
    f.write("  # =============================================================\n")
    f.write("  # Model: \n  #\n")
    #
    # there may be a whole boatload of reactions to print 
    # define a generator
    def printable_reaction(k):
        for j in xrange(k):
            reaction=reactions[j]
            yield "  # " + reaction.strip() + "\n"
    for r in printable_reaction(len(reactions)):
        f.write(r)
    #for reaction in reactions:
    #    f.write("  # " + reaction.strip() + "\n")
    #
    # write a whole shitload of comments with the rate constants, unless
    # the rate constants are themselves written out
    #
    if substituteRates:
        f.write("  #\n  # Parameter values used: \n  #\n")
        for r in rates:
            f.write("  # " + r + "=" + str(rates[r]) + "\n")
    if len(holdrates)>0:
        f.write("  global "+holdrates[0]+"\n")
    if len(frozenvars) >0:
        f.write("  #\n  # Frozen Variables (species with fixed derivatives):\n  #\n") 
        for v in frozenvars: f.write("  # "+str(v)+ "\n")
    f.write("  # =============================================================\n")
    
    
    results = interpreter.invokeParser(reactions, dump=False)
    results = expand(results)
    s=interpreter.makeSymbolDictionary(results, rates, ics)
    odeterms = interpreter.makeODETerms(results, s, frozen=frozenvars)
    
    if substituteRates:
        newterms = substituteRateConstants(odeterms, rates, s, holdrates)
        na = []
       # rate constants also appear in assignment rules
        for a in assignments:
            value,assignment = a.split("=")
            assignment = sympify(assignment)
            newassignment = value+"="+str(applyReplacementRules(assignment, rates, s))
            na.append(newassignment)
        assignments = na
    else:  # better this way because applyReplacementRules grows superquadratically
        newterms = odeterms
        ks=rates.keys()
        ks.sort()
        f.write("  # rate constants\n")
        for k in ks:
            if k in holdrates:
                continue
            v=rates[k]
            nextline="  "+k+" = " + str(v)+"\n"
            f.write(nextline)
        
    (y, yprime, ics) = makeODEfunc(newterms, s, ics)
    
    y = map(str, y)
    yprime = map(str,yprime)
    f.write("# pick up values from previous iteration\n")  
    n = len(y)
    for i in range(n):
        nextline = "  "+y[i]+" = "+"max(0, y["+str(i)+"])\n"
        #
        nextline = nextline.replace("{","[").replace("}","]")
        #
        f.write(nextline)
    if (len(assignments)>0):
        f.write("# apply boundary conditions / assignment rules \n")
    for a in assignments:
        f.write("  "+a+"\n")
    f.write("# calculate derivatives of all variables\n")
    f.write("  yp=[0 for i in range("+str(n)+")]\n")
    for i in range(n):
        nextline = "  yp["+str(i)+"] = " + yprime[i]+"\n"
        #
        #
        nextline = nextline.replace("{","[").replace("}","]")
        #
        #
        f.write(nextline)

    
    
    f.write("  return yp\n")
    f.close()
    
    #sys.exit()
    return (y, ics, codefile)

#*******************************************************************************  
    
def makeplot(x,y,title, fig, i, j, k):
    # input: arrays x, y, string title
    # output: displays a plot using matplotlib
    ax=fig.add_subplot(i,j,k)
    ax.plot(x,y)
    ax.set_title(title)
  
#*******************************************************************************  
    
def samePlot(sol, t, variables, filename, plotvars):
    # plots all variables in sol on a single plot
    # input: sol = output of odeint
    #       t = list of times input to odeint
    # print "samePlot"
    fig = plt.figure()
    d=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    fig.suptitle(filename+" ("+d+")")
    ax = fig.add_subplot(1,1,1)

    plots=[]
    for v in plotvars:
        i = variables.index(v)
        y=[x[i] for x in sol]
        p,= ax.plot(t, y)
        plots.append(p)
    
    ax.legend(plots, plotvars)
    plt.show()
        
    return
  
#*******************************************************************************  
        
def plotVariables(sol, t, variables, filename, plotvars, columns=3):
    # plots all variables in sol
    # input: sol = output of odeint
    #       t = list of times input to odeint
    #       variables = list of string variables naming columns in sol
    if ("-sameplot" in sys.argv):
        samePlot(sol, t, variables, filename, plotvars)
        return
    fig = plt.figure()
    d=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    fig.suptitle(filename+" ("+d+")")
    n = len(plotvars)
    cols = int(columns)
    rows = n/cols
    if cols*rows < n:
        rows = rows + 1

    k = 1
    for v in plotvars:
        i = variables.index(v)
        y = [x[i] for x in sol]
        makeplot(t,y,v,fig, rows,cols,k)
        k+=1
    plt.show()
    return

#*******************************************************************************  
    
def getRunParameters(default_duration, default_step): 
    # input: default duration and step size for output
    # output: if invoked from command line, can override those with the -run 
    # option and those values are returned; otherwise the default values
    # are returned
    if ("-run" in sys.argv):
        i = (sys.argv).index("-run")+1
        j = i+1
        if j<len(sys.argv):
            duration=float(sys.argv[i])
            step=float(sys.argv[j])
            return (duration, step)          
    return(default_duration,default_step)
     
  
    
#*******************************************************************************  

def generateSimulation(inputfile="", step=1, duration=100, plot=False, format="CSV", \
    output="", vars=[], plotcolumns=3, sameplot=True):
    formattype = {"CSV":"CSV","TABLE":"TXT", "TSV":"TSV"}
    # input: model file
    # return value: none
    # output: generates simulation plots from model file 
    # if invoked from command line, "solver -in filename -run duration step"
    
    
    holdrates=[]
    if "-scan" in sys.argv:
        i = (sys.argv).index("-scan")+1
        
        if i+3<len(sys.argv):
            holdrates.append(sys.argv[i])
            scanstart = sys.argv[i+1]
            scanstop = sys.argv[i+2]
            scandelta = sys.argv[i+3]
        else:
            sys.exit("Error: expecting 'K min max delta' after keyword -scan")
        
     
     
    
    (r, ic, rates, frozenvars, functions,assignments, filename)=readmodel(INFILE=inputfile)
    (variables, y0, tmpdotpy) = generatePythonFunction(r, rates, ic, frozenvars, 
        functions, assignments, holdrates)

    
    f=open(tmpdotpy,"r")
    rhs=f.readlines()
    f.close()
    
    try:
        os.remove(tmpdotpy)
    except:
        print "Warning: unable to delete "+os.path.abspath(tmpdotpy)
    
     
    
    
    pyfile = ""
    if ("-pyfile" in sys.argv):
        i = (sys.argv).index("-pyfile")+1
        if i<len(sys.argv):
            pyfile=sys.argv[i]
        else:
            print "Waring: solver: expecting file name after -pyfile option."
            pyfile = "newsolver.py"
    
    if pyfile=="":
        pyfile=utils.timed_file_name("solver-for-"+os.path.basename(filename),"py")

    
    #pyfile = utils.uniqueFileName(pyfile, type="py")   
    f=open(pyfile,"w")
    
    f.write("import numpy as np\n")
    f.write("import cellerator.solver\n")
    f.write("from scipy.integrate import odeint\n")
    f.write("\n")
    for line in rhs: f.write(line)
    f.write("\n\n\n")
    f.write("def thesolver():\n")
    f.write("    filename =\""+filename+"\"\n")
    svars = str(map(utils.deindex, map(str, variables)))

    f.write("    variables="+svars+"\n")
    (runtime, stepsize)=getRunParameters(default_duration=duration,default_step=step)
    f.write("    runtime = "+str(runtime)+" \n")
    f.write("    stepsize = "+str(stepsize)+" \n")
    f.write("    t = np.arange(0,runtime+stepsize,stepsize)\n")
    f.write("    y0 = "+str(y0)+"\n")
    
    if "-mxstep" in sys.argv:
        i=(sys.argv).index("-mxstep")+1
        if i<len(sys.argv):
            mxstep=sys.argv[i]
    else:
        mxstep = "500000"
    
    if len(holdrates)>0:
        f.write("    global " + holdrates[0] + "\n")
        f.write("    results=[]\n")
        f.write("    "+holdrates[0]+"="+str(scanstart)+"\n")
        f.write("    while "+ holdrates[0]+" <= " + str(scanstop) + ":\n")       
        f.write("        sol = odeint(ode_function_rhs, y0, t, mxstep="+mxstep+")\n") 
        f.write("        "+holdrates[0]+"+="+str(scandelta)+"\n")
        f.write("        res=["+ holdrates[0] + "] + list(sol[-1])\n")
        f.write("        results.append(res)\n")
        f.write("        print res\n")
        f.write("    return results\n")
        
    
        f.write("if __name__==\"__main__\":\n")
        f.write("    thesolver()\n\n")
        f.close()
        
           
        return os.path.abspath(pyfile)    
 
            
    f.write("    sol = odeint(ode_function_rhs, y0, t, mxstep="+mxstep+")\n") 
     
    fmt = format
    if "-format" in sys.argv:
        i = (sys.argv).index("-format")+1
        if i<len(sys.argv):
            fmt = sys.argv[i]
    fmt = fmt.upper()
    if not fmt in formattype:
        print "Error: solver: invalid format = ", fmt, " requested."
        fmt = "CSV"
    fmt = formattype[fmt]
    
    outfile=output
    if ("-out" in sys.argv):
        i = (sys.argv).index("-out")+1
        if i<len(sys.argv):
            outfile=sys.argv[i]
            outfile=utils.uniqueFileName(outfile, type=fmt)
        else:
            print "Warning: solver: expecting file name after -out option."
            #outfile = "solution."+fmt
    if output=="":        
        basename = os.path.basename(filename)
        #basename=basename.split(".")[0]
        outfile = utils.timed_file_name("Solution-for-"+basename,fmt)
           
    #f.write("    outfile = \""+utils.uniqueFileName(outfile, type=fmt)+"\"\n")
    f.write("    outfile = \""+outfile+"\"\n")
    f.write("    of = open(outfile,\"w\")\n")
    f.write("    for i in range(len(t)):\n")
    f.write("        time = t[i]\n")
    f.write("        data = map(str,sol[i])\n")
    demarcater={"CSV":"\",\"","TSV":"\"\t\"","TXT":"\" \""}
    f.write("        data = str(time) + " + demarcater[fmt] + "+(" + demarcater[fmt]+".join(data))\n")
    f.write("        of.write(data+\"\\n\")\n")
    f.write("    of.close()\n")
       
    
    #print "argv is: ", sys.argv
     
    if ("-plot" in sys.argv) or plot==True:
        plotvars = vars
        f.write("    plotvars = "+svars+"\n")
        
        
        if ("-plot" in sys.argv):
            i = (sys.argv).index("-plot")
            if i>0:
                i += 1
                while i<len(sys.argv):
                    nextvar = sys.argv[i]
                    i += 1
                    if nextvar[0]=="-": break
                    nextvar = utils.deindex(str(nextvar))
                    if nextvar in variables:
                        plotvars.append(nextvar)
                    else:
                        print "Error: requested plot variable: "+nextvar+ " does not "\
                        + " exist in the model."
                plotvars = list (set (plotvars) ) # remove any dupes 
                f.write("    plotvars = "+str(plotvars)+"\n")
        if len(plotvars)==0:
            plotvars = variables
            f.write("    plotvars = "+str(plotvars)+"\n")  

        
        plotcols = plotcolumns
        if ("-plotcolumns" in sys.argv):
            i = (sys.argv).index("-plotcolumns") + 1
            if i < len(sys.argv):
                plotcols=int(sys.argv[i])
        plotcols = max(plotcols, 1)
        if ("-sameplot" in sys.argv) or sameplot==True:
            f.write("    cellerator.solver.samePlot(sol, t, variables, filename, plotvars)\n")
        else:  
            f.write("    cellerator.solver.plotVariables(sol, t, variables, filename,"\
            +" plotvars, columns="+str(plotcols)+")\n")
        f.write("    return\n\n")
        
    f.write("if __name__==\"__main__\":\n")
    f.write("    thesolver()\n\n")
    f.close()
    
   
        
    return os.path.abspath(pyfile)    

#*******************************************************************************  
#

    
def makesolver():
    if not ("-quiet" in sys.argv):
        d=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print "py[cellerator]: solver ("+d+")"
        
    f = generateSimulation()
    if "-quiet" in sys.argv: return f
    print "Simulation driver written to: ", f
    return f
 
#*******************************************************************************  
       
if __name__ == "__main__":
    makesolver()
    
