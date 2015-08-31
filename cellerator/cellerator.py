#!/usr/bin/env python
#
# B.E.Shapiro 11 April 2015 
# revised 30 Aug 30 for github load
#
# Convenient Cellerator Function Wrapper for iPython notebook
# usage of Python Cellerator implemenation functions
#
#
#****************************************************************************
#
#    This file is part of Cellerator
#
#    Cellerator converts reactions, expressed in a text-formatted arrow-based 
#    notation into differential equations and performs numerical simulations.
#    It is a python implementation of most the functionality of cellerator
#    Copyright (C) 2012-2015 Bruce E Shapiro.
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
import parser
import converter
import interpreter
import solver
import utils
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import random
import sbmlwrite, sbml
import time


######################################################

def getODES(model):
	"""
	getODES(file) will read a cellerator arrow text file
	run the interpreter and return the differential equations 
	as a newline delimited string.
	"""
	return interpreter.interpret(model)

######################################################

def PlotSize(width, height):
    mpl.rcParams["figure.figsize"]=width, height

######################################################

def PrintODES(model):
	""" printODES(file) will read a cellerator arrow text file
	run the interpreter to convert the list of arrows to a 
	systems of differential equations, and will print out
	the differential equations.
	
	To get the ODES as a string rather than printing the use 
	the fuction getODES(file) instead of printODES."""
	s=interpreter.interpret(model)
	print s
	
	
def PrintModel(model):
    """printModel(model) will read and print a model file
    """
    try:
	    with open(model, "r") as f:
		    text = f.readlines()
	    text = "".join(text)
	    print text
    except:
        print "Error: unable to open "+os.path.abspath(model)+"\n"\
        +"check spelling and try again"

######################################################

def Solve(inputfile, step=1, duration=100,  solverfile="",
    output="solution", mxstep=50000, scan=[], run=True, RATES={}, IC={}, 
    timer=False):
    """Solve(inputfile, step=1, duration=100, solverfile="newsolver.py",
             mxtep=50000, run=True, scan=[], RATES={}, IC={})
    
    input parameters:
    
    inputfile - pyCellerator input file containing a model
    step - (output) step size for integrator. The default integrator
		is odeint. The step size specifies the interpolated output
		values, not the actual values used by the integrator
    duration - duration of integration
    mxstep - maximum number of steps allowed by odeint. 
    solverfile - name of python program that will be automatically
		generated to solve the sytem
    run = True; if false, just generate the code, don't run simulation
    RATES={}: optional dictionary of paramaters that can be used to 
		override one or more rate values in the model
	IC={}: optional dictional of initial conditions that can be used
		to override zero or more initial conditions in the model
		
Return value: when run=True (default)
    tuple (t, variables, solutions)
	where
	t=numpy array of values at which the solution is interpolated.
	variables=list of nambed variables in the solution
	solution=numpy array of values. There is one column for each variable.
		the columns are in the same order as the variable
		names. The length of the columns is the same
		as the length of the values in t. The contents
		of each column contains the values of the corresponding
		variable at the corresponding times in t.
	when run=False, name of simulation code
    """
    
    # input: model file
    # return value: output of integrator
    
    hold=[]
    if len(scan)>0:
		if len(scan) == 4:
			scanparameter, scanstart, scanstop, scandelta = scan
			hold=[scanparameter]
			if not isinstance(scanparameter,str):
				sys.exit("Error: expecting a string for the name of the scan parameter")
		else:
			sys.exit("Error: expeciting scan=['parameter',start,stop,delta]")

    
    
    # if invoked from command line, "solver -in filename -run duration step"
    
    if timer:
        ttimer=time.time()
    (r, ic, rates, frozenvars, functions,assignments, filename)=solver.readmodel(INFILE=inputfile)
    if timer:
        tread=time.time()-ttimer
        
    #print "rates:",rates
    #print "ic:", ic
    for rate in RATES:
		if rate in rates:
			rates[rate]=RATES[rate]
		else:
			print "Warning: unexpected parameter ",rate,\
			" set in RATES not found in model was ignored."
    for variable in IC:
		if variable in ic:
			ic[variable]=IC[variable]
		else:
			print "Warning: unexpected variable ", variable,\
			" set in IC not found model was ignored."
    #print "rates:",rates
    #print "ic:", ic
		
    
    if timer:
        ttimer=time.time()
    (variables, y0, tmpdotpy) = solver.generatePythonFunction(r, rates, 
         ic, frozenvars, functions, assignments, holdrates=hold)

        
    #
    # generatePythonFunction creates the file tmp.py
    #
    f=open(tmpdotpy,"r")
    rhs=f.readlines()
    f.close()
    #
    #clean up
    try:
        os.remove(tmpdotpy)
    except:
        print "Warning: unable to remove the file "+os.path.abspath(tmpdotpy)+"\n+"\
        +"Perhaps an authorization issue? You do not need to keep this file\nand it may be deleted"
    #
    # name the solver program
    #
    pyfile = solverfile
    #print "DBG>>>solverfile:", solverfile
    #print "DBG>>> INPUT >>>", inputfile

    if pyfile=="":
        inputdrive, inputpath = os.path.splitdrive(inputfile)
        inputpath, inputfilename = os.path.split(inputpath)
        #print "DBG:drive:", inputdrive
        #print "DBG:path:", inputpath
        #print "DBG:file:", inputfilename
        pyfile=utils.timed_file_name("solver-for-"+inputfilename,"py")
        pyfile = os.path.join(inputpath, pyfile)
        #print "DBG: resolved pyfile:", pyfile
    else:
        pyfile = utils.uniqueFileName(pyfile, type="py")   
    
    File2Import=(os.path.basename(pyfile).split("."))[0]
    #
    # generate the code for the solver
    #
    f=open(pyfile,"w")
    f.write("import numpy as np\n")
    # f.write("import solver\n")
    f.write("from scipy.integrate import odeint\n")
    f.write("\n")
    #
    #  paste in the contents of the file "tmp.py" that was produced by generatePythonFunction
    #  with the right hand side of the ode defined in it ode_function_rhs
    #
    for line in rhs: f.write(line)
    #
    # write driver function to run the solver just created 
    f.write("\n\n\n")
    f.write("def thesolver():\n")
    f.write("    filename =\""+filename+"\"\n")
    svars = str(map(utils.deindex, map(str, variables)))
    f.write("    variables="+svars+"\n")
    (runtime, stepsize)=solver.getRunParameters(default_duration=duration,default_step=step)
    f.write("    runtime = "+str(runtime)+" \n")
    f.write("    stepsize = "+str(stepsize)+" \n")
    t = np.arange(0,runtime+stepsize,stepsize)
    f.write("    t = np.arange(0,runtime+stepsize,stepsize)\n")
    f.write("    y0 = "+str(y0)+"\n")   
    
    if len(scan)>0:		
        f.write("    global " + scanparameter + "\n")
        f.write("    results=[]\n")
        f.write("    "+scanparameter+"="+str(scanstart)+"\n")
        f.write("    while "+ scanparameter+" <= " + str(scanstop) + ":\n")       
        f.write("        sol = odeint(ode_function_rhs, y0, t, mxstep="+str(mxstep)+")\n") 
        f.write("        "+scanparameter+"+="+str(scandelta)+"\n")
        f.write("        res=["+ scanparameter + "] + list(sol[-1])\n")
        f.write("        results.append(res)\n")
        # f.write("        print res\n")
        f.write("    return results\n")
    else:   
		f.write("    sol = odeint(ode_function_rhs, y0, t, mxstep="+str(mxstep)+")\n")   
		f.write("    return sol\n\n")
    #
    # write main program
    #
    f.write("if __name__==\"__main__\":\n")
    f.write("    thesolver()\n\n")
    f.close()
    if timer:
        tgenerate=time.time()-ttimer
 
    #
    # stop here if all you want is the code
    #
    if not run:
        if timer:
            return (tread, tgenerate, os.path.abspath(File2Import))
        else:
            return os.path.abspath(File2Import)
    #    
    #
    # import the code just created
    #
    
    exec("import "+ File2Import)
    
    #
    # run the code just created
    #
    
    exec("temporary_solution = "+File2Import + ".thesolver()") 
    

    if len(scan)>0: 
	# returns a 2-tuple
    # names of the variables (column headers)
    # the solution - one variable per column. 
    # The first column gives the parameter scan from start to finish
    # each column j gives the value of variable j-1 at the end of the run  
	    return (variables, temporary_solution)		
    else:
	# return a tuple with a list of times (row headers)
    # names of the variables (column headers)
    # the solution - one variable per column. Each column is the time course for the corresponding variable in variables
    #  
		return (t, variables, temporary_solution)
		
###########################################################################
def SetAx(ax, scales=(), lims=(), labels=()):
    if scales!=():			 
		try:
			ax.set_xscale(scales[0])
		except:
			pass
		try:
			ax.set_yscale(scales[1])
		except:
			pass
    if lims!=():		 
		try:
			ax.set_xlim(lims[0:2])
		except:
			pass
		try:
			ax.set_ylim(lims[2:4])
		except:
			pass
    if labels!=():
		try:
			ax.set_xlabel(labels[0])
		except:
			pass
		try:
			 ax.set_ylabel(labels[1])
		except:
			pass
		try:
			ax.set_title(labels[2])
		except:
			pass

    return ax


######################################################	
def ParametricScanPlot(v,s,k,variable,marker="",color="brown", bg="antiquewhite",
	size=(), scales=(), lims=(), labels=()):
    """ParametricScanPlot(v,s,k,variable)
    input: v,s: output of Solve in scan mode
    variable: name of parameter used to do scan
    marker: optional marker to use on plot
    color: optional color to use for plot
    bg: optional background color
    
    return value: axis instance
    """
    if size!=():
		try:
			x,y=size
			PlotSize(x,y)
		except:
			pass	

    figure,ax=plt.subplots()
    try:
        i=v.index(variable)
    except:
        print "Variable ",variable,"not found in list", v
        return

    S=np.array(s)
    y=S[:,i+1]
    x=S[:,0]

    ax.plot(x,y,color=color)
    if marker != "":
        plt.scatter(x,y,color=color,marker=marker)
    ax.set_axis_bgcolor(bg)

    ax = SetAx(ax, scales=scales, lims=lims, labels=labels)

 
    return ax

######################################################
def PlotVars(t,v,s,plotvars,bg="antiquewhite", size=(), scales=(), lims=(), labels=(), ax=None):
    """PlotVars(t,v,s,variables,bg="antiquewhite")
    plots the species in the list variables on a single plot
    t,v,s are the olutput of solve
    return value: a matplotlib axis instance
    """
    if size!=():
		try:
			x,y=size
			PlotSize(x,y)
		except:
			pass	
    
    okinds=[]
    for variable in plotvars:
        try:
            i=v.index(variable)
            okinds.append(i)
        except:
            print "Error: ",variable," not found in solution. Variables are ",v
    nvars = len(okinds)       
    
    if nvars>0:
        if isinstance(ax,type(None)):
             fig,AX=plt.subplots()
        else:
            AX=ax
        for ivar in okinds:
            variable = v[ivar]
            yvals = s[:,ivar]
            AX.plot(t,yvals,label=variable)
        AX.set_axis_bgcolor(bg)
        if isinstance(ax,type(None)):
            plt.legend(loc="best")

        AX = SetAx(AX, scales=scales, lims=lims, labels=labels)
				
        return AX

######################################################
def PlotAll(t, v, s, loc="best",  bg="antiquewhite"):
    """PlotAll(t,v,s) plots the output of Solve
    t = numpy array list of times of length p
    v = list of (text) variable names of length m
    s = solution of values, p rows of length m. The jth row contains
    loc="best" optional legend location
    the p values corresponding to the jth variable in v.
    return value is an axes instance
    """
	
    nvars = len(v)
    f,ax=plt.subplots()
    for j in range(nvars):
        var=v[j]
        yvals=s[:,j]
        ax.plot(t,yvals,label=var)
    ax.legend(loc=loc)
    
    ax.set_axis_bgcolor(bg)
        
    return ax

######################################################


def PlotParametric(t,v,s,n,m, color="Red", size=(), scales=(), lims=(), labels=(),  bg="antiquewhite"):
    """PlotParametric(t,v,s,n,m,color=colorvalue)
    produces a parametric plot of two variables from a solution produced
    by Solve.
    
    Input:
    t, v, s: output from Solve
    n, m: either eithers or strings variable names. If integers, 
    they give the index into v (column of s) of the corresponding
    variable to be plotted (n for x axis, m for yaxis). If
    strings, they must be the names of the variables in v.
    color = any pyplot compatible color
    """
    if size!=():
		try:
			x,y=size
			PlotSize(x,y)
		except:
			pass	

    if isinstance(n, str):
		try:
			n=v.index(n)
		except:
			print "Requested variable "+n+" not found in ",v
			return
		xvar=s[:,n]
		xlab=v[n]
    elif isinstance(n, int):
		try:
			xvar=s[:,n]
			xlab=v[n]
		except:
			print "Index",n,"out of bounds; there are only",len(v),"variables"
			return
    else:
		print "Unable to identify variable ",n
		return
		
    if isinstance(m, str):
		try:
			m=v.index(m)
		except:
			print "Requested variable "+m+" not found in ",v
			return
		yvar=s[:,m]
		ylab=v[m]	
    elif isinstance(m,int):
		try:
			yvar=s[:,m]
			ylab=v[m]
		except:
			print "Index",m,"out of bounds; there are only",len(v),"variables"
			return
    else:
		print "Unable to identifu variable",m
		return
		
    fig,ax=plt.subplots()	
    ax.plot(xvar,yvar, color=color)
    
    thelabels=[xlab, ylab, ""]
    if len(labels)>0:
	    if labels[0]!="":  
			thelabels[0]=labels[0]
    if len(labels)>1:
		if labels[1]!="":  thelabels[1]=labels[1]
    if len(labels)>2:  thelabels[2]=labels[2]
    
    ax = SetAx(ax,scales=scales, lims=lims, labels=tuple(thelabels))
    
    ax.set_axis_bgcolor(bg)
        
    return ax

######################################################

	
def getcolornames(n):
    """getcolornames(n) returns n randomly selected color names
	from the predefined colors in matplotlib (must be n <= 140).
    """
    all = list(colors.cnames)
    random.shuffle(all)
    return all[:n]

######################################################
    
def PlotColumns(t,v,s,ncols=3, colors=[], bg="antiquewhite",
    size=(), fontsizes=()):
    """PlotColumns(t,v,s, ncols=val, colors=list, bg=color)
    plots the output of Solve as a grid of plots, with one
    variable on each plot.  
    
The size of the plots displayed can be changed by setting
	matplotlib.rcParams["figure.figsize"]=width,height
	
Input: 
    t,v,s: output of Solve
    ncols: number of columns in grid
    colors = []: a list of colors to be assigned to the curves.
    If shorter than the number of variables, random colors will be
    selected from maplotlib.colors.cnames.
    bg = desired background color for 
    size = (width, height)
    fontsizes=(xtick font size, ytick font size, title font size)

    """
    nvars=len(v)
    if size!=():
		try:
			x,y=size
			PlotSize(x,y)
		except:
			pass	
                
    if len(colors) < nvars:
        xcolors = getcolornames(nvars-len(colors))
        colors = colors+xcolors
    
    def makeplot(i, ax):
        variable = v[i]
        yvals = s[:,i]
         
        ax.plot(t,yvals, color=colors[i])
        if len(fontsizes)>0:
            ax.tick_params(axis="x",labelsize=fontsizes[0])
        if len(fontsizes)>1:
            ax.tick_params(axis="y",labelsize=fontsizes[1])
        if len(fontsizes)>2:
            ax.set_title(variable,fontsize=fontsizes[2])
        else:
            ax.set_title(variable)
        ax.set_axis_bgcolor(bg)
           
    if ncols == 1:
        fig,AX=plt.subplots(nrows=nvars,ncols=1, sharex=True)
        i=-1
        for ax in AX:   
            i=i+1
            makeplot(i, ax)
        
    elif nvars <= ncols:
        fig,AX=plt.subplots(nrows=1,ncols=nvars)
        i=-1
        for ax in AX:
            i=i+1
            makeplot(i, ax)
    else:
        nrows = nvars/ncols
        if ncols * nrows < nvars:
            nrows = nrows + 1
        fig, AX = plt.subplots(nrows=nrows, ncols=ncols)
        
        i=-1
        for row in range(nrows):
            for ax in AX[row]:
                i=i+1
                if i<nvars:
                    makeplot(i,ax)
                else:
                    ax.axis("off")
    return fig
 
 
####################################################
   
def newModel(reactions, IC="", RATES="", FUNCTIONS="", ASSIGNMENTS="", 
    output="newmodel"):
    """newModel(reactions, IC, RATES, FUNCTIONS, ASSIGNMENTS)
    creates a new model file from strings.
    """
    reactions = "$REACTIONS\n"+reactions
    ic = "$IC\n"+IC
    rates = "$RATES\n"+RATES
    funcs= "$FUNCTIONS\n"+FUNCTIONS
    assigns="$ASSIGNMENTS\n"+ASSIGNMENTS
    
    model = "\n".join([reactions, ic, assigns, funcs, rates])
    model = model.replace(";", "\n")
    model = model.replace("\n\n","\n")
    
    # clean up
    m = model.split("\n")
    newmodel=[]
    for line in m:
        line = line.strip()
        if len(line)>0:
            if line[0]=="$":
                newmodel.append(line)
            else:
                newmodel.append(" "+line)

    model = "\n".join(newmodel)
    
    outputfile = utils.uniqueFileName(output,type="model")
    f=open(outputfile,"w")
    f.write(model)
    f.close()
    
    fname=os.path.abspath(outputfile)

    return fname
####################################################
def ConvertSBML(sbmlfile, modelfile=""):
	"""Converts an SBML file to a model file.
	"""
	if modelfile=="":
		modelfile = utils.uniqueFileName("newmodel.model")
		
	f=sbml.SBMLREAD([], INFILE=sbmlfile, OUTFILE=modelfile)
	return os.path.abspath(f)
	
####################################################
def GenerateSBML(infile, output=""):
    """GenerateSBML(inputfilename, output="")
    Reads the pyCellerator model infile and generates SBML.
    
    If output is specified, the SMBL is written to the specified
    file and the absolute path name of the file is returned.
    
    if the out is not specified, the SBML is returned as a single
    string (with embedded newline characters). 
    """
    s=sbmlwrite.genmodel(infile,"")
    
    if len(output)>0:
        filename = utils.uniqueFileName(output,type="xml")
        f=open(filename,"w")
        f.write(s)
        f.close()
        p = os.path.abspath(filename)
        return(p)
####################################################
        
def PrintSBML(xml):
	"""PrintSBML(xml) prints an SBML model to the 
	standard output.
	"""
	with open(xml,"r") as f:
		sbml=f.readlines()
		print "".join(sbml)  
		    
####################################################

def ToMathematica(infile, outfile=""):
	"""ToMathematica(infile, [outfile=""])
	converts the specified model file to legacy Cellerator 
	Mathematica notebook.
	"""
	if len(outfile)==0:
		outfile="translated-model.nb"
	
	out=utils.uniqueFileName(outfile, type="nb")
	
	fout = converter.convert(infile, outfile)
    
	return fout
####################################################
                   
def PrintMathematica(nb):
	"""PrintMathematica(nb) prints a legacy cellerator
	model to the standard output.
	"""
	with open(nb,"r") as f:
		mathematica=f.readlines()
		print "".join(mathematica)
