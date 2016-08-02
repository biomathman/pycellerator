#!/usr/bin/env python
#
# Cellerator Flux Analyis Module
# B.E.Shapiro 30 Aug 2015
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
#  python flux.py
#       -in inputfilename
#       [-out outputfile]
#       [-dump]    # to dump the parser
import datetime
import os.path
import sys
from reader import readmodel
from parser import Reaction,StripComments
from utils import *
import numpy as np
from sympy import *

try:
	from pulp import *
except:
	print "Cellerator.flux unable to import pulp. Will not be able to process flux models."
	
import matplotlib.pyplot as plt
from pylab import *

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

def initialize_fluxtools(infile,  ofile="", dumpparser=False):

    DUMP = dumpparser or ("-dump" in sys.argv)
    
    if infile=="":
        if("-in" in sys.argv):
            i = (sys.argv).index("-in")+1
            if i<len(sys.argv):
                infile=sys.argv[i]
            else:
                raise SystemExit("Error: use '-in filename' to specify input file.")
        else:
            raise SystemExit("Error: use '-in filename' to specify input file.")
            
    if os.path.isfile(infile):
        datafile = infile
    else:
        raise SystemExit("Error: "+os.path.abspath(infile)+" File not found.")
        
    outfile=ofile    
    if outfile == "":
        if ("-out" in sys.argv):
            i = (sys.argv).index("-out")+1
            if i<len(sys.argv):
                outfile=sys.argv[i]
            else:
                outfile="flux-analysis.out"
        else:
            outfile="flux-analysis.out"
    outfile=os.path.abspath(uniqueFileName(outfile))
            
    lines=readmodel(INFILE=datafile)[0]
    results = invokeParser(lines, dump=DUMP)
    return (results, outfile)
 

#*******************************************************************************    
def rint(x,digs):
    if abs(round(x,digs)-int(x))<10**-digs:
        return int(x)
    return round(x,digs)
#*******************************************************************************    
def make_species_dict(db):
    # input: database of reactions
    # output: dictionary of species with an integer assigned to each species
    dict={}
    i = 0
    for r in db:
        LS = listify(r.LHS())
        RS = listify(r.RHS())
        S = list ( set(LS) | set(RS) )
        for s in S:
            if s=='Nil':
                continue
            if not(s in dict):
                dict[s]=i
                i = i+1
    return dict
#*******************************************************************************    
def species_dict_to_list(db):
    """
    input: species dictionary {species:int, species:int, ...}
    output: list [species, species, species,...]
    where the ith species has index i in the dictionary (i.e., correspoinds to 
    row i of the stoichiometry matrix
    """
    n = len(db)
    species=[0 for i in range(n)]
    for s in db:
        i = db[s]
        species[i]=s
    return species
#*******************************************************************************    

def identify_cmap(UCmap):
	maps = plt.cm.datad.keys()
	for m in maps:
		if UCmap==m.upper():
			return(m)
	return plt.cm.ocean

#*******************************************************************************    
def visualize_matrix(Nraw):
    """ add options 
	-cmap colormap
	-image filename
	-vmin value
	-vmax value
	-abs take absolute value of matrix
"""
	# revis 8/2/16
    
    #fig = plt.figure()
    #ax = fig.add_subplot(111)
    
    
    N=np.absolute(Nraw) if ("-abs" in sys.argv) else Nraw 
    big  = 1
    little = 0 if ("abs" in sys.argv) else -1
    color_map=plt.cm.ocean
    if "-cmap" in sys.argv:
        i = (sys.argv).index("-cmap")+1
        if i<len(sys.argv):
            #color_map = identify_cmap(sys.argv[i].upper())
            color_map = sys.argv[i]
    if "-vmax" in sys.argv:
		i=(sys.argv).index("-vmax")+1
		if i<len(sys.argv):
			big=float(sys.argv[i])
    if "-vmin" in sys.argv:
		i=(sys.argv).index("-vmin")+1
		if i<len(sys.argv):
			little=float(sys.argv[i]) 
				  
    plt.imshow(N,aspect='auto',interpolation='nearest', vmin=little,
		vmax=big, cmap=color_map)
    ax=plt.gca()
    ax.axes.get_yaxis().set_visible(False)
    ax.axes.get_xaxis().set_visible(False)
    plt.colorbar()
    
    fig=plt.gcf()
    fig.tight_layout()
    
    picfile="Stoichiometryy-Matrix.png"
    if "-image" in sys.argv:
		i=(sys.argv).index("-image")+1
		if i<len(sys.argv):
			picfile=sys.argv[i]
	

    fig.savefig(uniqueFileName(picfile), dpi=300)
    
    show()
    
    return
    
#*******************************************************************************    
def write_matrix(N, S, filename):
    N=N.tolist()
    Z=zip(S,N)
    Z=[[y]+map(str,row) for (y,row) in Z]
    Z=map(lambda x:','.join(x),Z)
    # fout=uniqueFileName("Stoichiometry-Matrix.csv") if filename=="" else filename
    f=open(filename,"w")
    for z in Z:
        f.write(z+"\n")
    f.close()
    return os.path.abspath(filename)

#*******************************************************************************    
def stoich_matrix(db, out=""):
    """
    input: Cellerator database of flux reactions
    output: stoichiometry matrix
    """
    i = 0
    ispecies=0
    speciesDictionary=make_species_dict(db)
    # print "speciesDictionary:", speciesDictionary
    nr = len(db)
    # print "number of reactions: ", nr
    ns = len(speciesDictionary)
    # print "number of species: ", ns    
    N = [[0 for columns in range(nr)] for rows in range(ns)]
    
    for r in db:
        LS = zip(listify(r.LHS()), \
                 map(lambda x: rint(float(x),12), list(r.LH_STOIC())))
        RS = zip(listify(r.RHS()), \
                 map(lambda x: rint(float(x),12), list(r.RH_STOIC())))
        S = list ( set([x for (x,s) in LS]) | set([x for (x,s) in RS]) )
               
        LD = list2Dict(LS)
        RD = list2Dict(RS)
  
        # print "Reaction: ", i
        # print "LS = ", LS
        # print "RS = ", RS 
        # print "S = ", S

        # get the unique index for each species and figure out its stoichiometry
        stoichDict = {}
        for s in S:
            if s == "Nil": continue
            stoich = 0
            if s in LD: stoich -= LD[s]
            if s in RD: stoich += RD[s]
            row=speciesDictionary[s]
            column = i
            N[row][i] = stoich
        i = i+1
        N = np.array(N)
    species = species_dict_to_list(speciesDictionary)
    
    filewritten=write_matrix(N, species, uniqueFileName("Stoichiometry-Matrix.csv"))
    print "Stoichiometry Matrix --> "+filewritten
    
    visualize_matrix(N)
    
    return(N, species)


#*******************************************************************************    
def GetNullspace(M):
#
#   This function is not currently in use because the sympy
#   nullspace function returns a more useful form
#
#   returns a numpy array whose rows span the nullspace of M
#
#   SVD gives M = U.S.V (normally VT, but in Python, its V)
#   where the rows of V corresponding to zero singular values span the null space of M
#   and the rows of V corresponding to the non-zero singular values span the range of M
#
#   Note that the sympy solution is much "cleaner" in that it normally consists of matrices of 1's
#   and -1's because it uses row reduction rather than SVD

    TEENY = 1e-10
    U,S,V=np.linalg.svd(M,full_matrices=True)
    # print "V=",V
    shapeOfV=V.shape
    S = S.tolist()
    nonzero = map(lambda x:abs(x)>TEENY,S)
    RangeOfM = sum(nonzero)
    # print "shape of V", shapeOfV
    # print "singular values: ", S
    # print "nonzero: ", nonzero, " Range of M = ", RangeOfM
    NS = np.array((V.tolist())[RangeOfM:])
  
    
    
    return NS
    
#*******************************************************************************    
def cons_matrix(db,out=""):
    """
    input: Cellerator database of flux-only reactions
    output: conservation matrix as np.array()
    """
    from sympy.matrices import Matrix
    N, S = stoich_matrix(db)
    # 
    # convert to a sympy matrix from a numpy array
    #
    #NT = matrices.Matrix(N.T)
    NT = Matrix(N.T)
    # determnine nullspace
    gamma = NT.nullspace()
    nullspace=[]
    for vec in gamma:
        row = []
        for v in vec:
            row.append(v)
        nullspace.append(row)
    nullspace = np.array(nullspace)    
    # print "Sympy solution:\n", nullspace
    
    # Nullspace=GetNullspace(N.T)
    # print "SVD solution:\n", Nullspace 
    
    filewritten=write_matrix(nullspace, S, uniqueFileName("Conservation-Matrix.csv"))
    print "Conservation Matrix --> "+filewritten
  
    
    return (nullspace, S)
#*******************************************************************************    
def FBA_Solve(db, out=""):
    print 'solving'
    N, S = stoich_matrix(db)   
    prob=LpProblem("Flux Balance", LpMaximize)     
    # print "N=:\n",N
    V=[]
    f=[]
    OF = 0
    for r in db:
        (lower, v, upper, fluxval, obj)= r.Rates()
        if lower=='-inf': 
            lower = None
        else:
            lower = float(lower)
        if upper=='inf': 
            upper = None
        else:
            upper = float(upper)
        variable = LpVariable(v,lowBound=lower, upBound=upper)
        V.append(variable)
        # f.append(float(obj))
        OF += float(obj)*variable

    # print V
    # objective function
    print "The objective function is ", OF
    prob += OF
    # constraints
    # -- N.v = 0
    i = 0
    for row in N:
        c = [s*rate for (s,rate) in zip(row,V)]
        constraint = sum(c)
        prob += constraint == 0, str(S[i]) # "species"+'{0:06d}'.format(i)
        # print i," ", constraint    
        i = i+1
    # write LP file
    lpfilename=uniqueFileName("fluxmodel.lp")
    prob.writeLP(lpfilename)
    print "LP file --> " + os.path.abspath(lpfilename)
    prob.solve()
    status = LpStatus[prob.status]
    logfilename = uniqueFileName("Flux-Analysis.log") if out=="" else out
    lf=open(logfilename,"w")
    lf.write(status+"\n")
    # lf.close()

    names =[]; values=[]
    for v in prob.variables():
        names.append(v.name);
        values.append(v.varValue)
        lf.write(str(v.name)+" = "+str(v.varValue)+"\n")
    lf.close()
    print "Results --> " + os.path.abspath(logfilename)
    return (names, values)
#*******************************************************************************    
def check_db(db,out=""):
    """
    returns T/F if the model in the database consists of only flux-based reactions
    """
    ok = True
    for r in db:
        if r.ArrowType() != "Flux":
            ok = False
            break
    return ok
#*******************************************************************************    

FOPTIONS = {"STOICH":stoich_matrix,
            "CONS":cons_matrix,
            "SOLVE":FBA_Solve,
            "CHECK":check_db}


def flux_tools(option):
    if option in FOPTIONS:
        data, outputfile = initialize_fluxtools("")
        ok = check_db(data)
        if not(ok):
            raise SystemExit("Error: FluxTools: Not a flux model.")
        F = FOPTIONS[option]
        result = F(data,out=outputfile)
    else:
        print "Invalid option ", option," valid options are: ",  list(FOPTIONS)
        result=[]
        
    return(result)

#*******************************************************************************    

def fat():
    if not ("-quiet" in sys.argv):
        d=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print "py[cellerator]: Flux Analysis Tool ("+d+")"
    option="STOICH"
    if "-option" in sys.argv:
        i = (sys.argv).index("-option")+1
        if i<len(sys.argv):
            option = sys.argv[i].upper()
        
    data = flux_tools(option)
    
    if "-quiet" in sys.argv: return
    print data
#*******************************************************************************    


if __name__ == "__main__":
    fat()

