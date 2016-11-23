#!/usr/bin/env python
#
#   reader.py
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

from utils import *
import sys, os.path
from sympy import *

def inputLinesToDictionary(lines):
    # input: lines of data purpotedly of the form "symbol = value"
    # output: dictionary {"symbol":value,...}
    d={}
    for line in lines:
        try:
            [var, val] = line.split("=")
        except:
            sys.exit("Error reading line '"+line+"'\nPerhaps a missing or extra '='")
        var = var.strip()
        var = deindex(var)
        val=val.rstrip("\n").strip()
        if is_number(val):
            d[var] = float(val)
        else:
            print "Error: "+line+" is not a numerical value"
            raise SystemExit("Correct Input Model File and re submit job.")
    return d
 
 #****************************************************************************
   
    
def readAssignmentFile(filename):
    if os.path.isfile(filename):
        f=open(filename)
        lines=f.readlines()
        f.close()
    else:
        raise SystemExit("Error: "+os.path.abspath(filename)+" does not exist.")
    d=inputLinesToDictionary(lines)
    return d
    
#****************************************************************************
 
def lambdify(fstring):
    f = fstring.strip()
    if "=" in fstring:
        function = f.split("=")
        if len(function)>1:
            fname = function[0]
            fdef = function[1]
            if fdef.strip()=="": 
                print "\nError: the function '"+str(fstring)+"' is missing a function "\
                    +"definition."
                fdef="0"
            if ("(" in fname) and (")" in fname):
                LP = fname.index("(")
                RP = fname.index(")")
                fn = fname[:LP]
                args=fname[LP+1:RP]
                lf = fn + " = lambda "+args+" : "+fdef
            else:
                print "\nError: the function '"+str(fstring)+"' is not formatted properly."\
                            + " Expecting f(args)=expression; no (..args..) found"
                lf =  fname + " = lambda x: x"
        else:
            print "\nError: the function '"+fstring+"' is not formatted properly."\
            + " No function definition was found. Expecting f(args)=expression"
            lf =  function + " = lambda x: x"
    else:
        print "\nError: Incorrectly formatted function '"+fstring+"'. A function definition "\
            +" should have the form f(args)=expression"
        lf = "bad_function_input = lambda x: x"

    return lf 
#**************************************************************************** 

def readmodel(INFILE=""):
    if INFILE=="":
        if ("-in" in sys.argv):
            i = (sys.argv).index("-in")+1
            if i<len(sys.argv):
                infile=sys.argv[i]
            else:
                raise SystemExit("Error: readmodel: expecting filename after -in")
        else:
            raise SystemExit("Error: readmodel: No Input File Given; use -in option.")
    else:
        infile=INFILE
        
    if not os.path.isfile(infile):
        raise SystemExit("Error: "+os.path.abspath(infile)+" File not found.")
    fullfile = os.path.abspath(infile)
    f=open(infile)
    lines = f.readlines()
    f.close()
    # convert to a single string, replace semicolons with newlines, then split back up
    singlestring = "".join(lines)       # join on newlines
    singlestring = singlestring.replace(";","\n")  # replace semicolons
    singlestring = singlestring.replace("$","\n$") # ensure $ are on first column
    singlestring = singlestring.replace("\\\n","\n") # end of line characters
    singlestring = singlestring.replace("\n\n", "\n") # ignore blank lines
    lines = singlestring.split("\n")    # split on newlines
    #print lines
    #sys.exit()
    #
    newlines=[]
    icommands=[]
    i=0
    continuation = False
    newline=""
    for line in lines:
        if len(line)<1: continue  
        previous=newline
        newline = line.strip()
        if len(newline)<1: continue
        #print newline
        if continuation:
            previous = previous[:-1]+newline
            newlines[-1]=previous
            i = i-1
        else:
            newlines.append(newline)
        if newline[0]=='$':
            icommands.append(i)
        i+=1
        continuation=(newline[-1]=="\\")
        
    # print newlines
    commands=[]
    icommands.append(len(newlines))
    newlines.append('$EOF')
    # print icommands
    for i in icommands:
        commands.append(newlines[i].upper())
    # print commands
    
    def _find_key(key):
        if key in commands:
            return(commands.index(key))
        else:
            raise SystemExit("Error: "+key+\
            " keyword not found in "+os.path.abspath(infile))
    def _find_optional_key(key):
        if key in commands:
            return(commands.index(key))
        else:
            return -1
    
    
    ireactions = _find_key('$REACTIONS')
    iic = _find_key('$IC')
    irates = _find_key('$RATES')
    ifrozen = _find_optional_key('$FROZEN')
    ifunctions = _find_optional_key('$FUNCTIONS')
    iassignments = _find_optional_key('$ASSIGNMENTS')

    # print ireactions, iic, irates
    
    reactions_start = icommands[ireactions]+1
    reactions_end = icommands[ireactions+1]
    
    ic_start = icommands[iic]+1
    ic_end = icommands[iic+1]
    
    rates_start = icommands[irates]+1
    rates_end = icommands[irates+1]
    
    if ifrozen>=0:
        frozen_start = icommands[ifrozen]+1
        frozen_end = icommands[ifrozen+1]
    if ifunctions>=0:
        functions_start = icommands[ifunctions]+1
        functions_end = icommands[ifunctions+1]
    if iassignments>=0:
        assignments_start = icommands[iassignments]+1
        assignments_end = icommands[iassignments+1]
     
    #print "reactions:", reactions_start, reactions_end
    #print "ic:", ic_start, ic_end
    #print "rates:", rates_start, rates_end
    
    reactions = newlines[reactions_start:reactions_end]
    # print "reactions:", reactions
    ic = newlines[ic_start:ic_end]
    # print "ic:",ic
    ics = inputLinesToDictionary(ic) 
    
    rates = newlines[rates_start:rates_end]
    # print "rates:", rates
    
    rates = inputLinesToDictionary(rates)

    if ifrozen>=0:
        frozenvars = newlines[frozen_start:frozen_end]
    else:
        frozenvars=[]
    # print "frozenvars = ", frozenvars
    
    # instantiate all parameter values in functions
    def ratify(expression):
        e1, e2 = expression.split(":")
        s = sympify(e2)
        for r in rates:
            s = s.subs(Symbol(r), rates[r])
        return e1+":"+str(s)
    
    if ifunctions>=0:
        functions = newlines[functions_start:functions_end]
        functions = map(lambdify, functions)
        functions = map(ratify, functions)
    else:
        functions=[]
    # print "functions=", functions
    
    if iassignments>=0:
        assignments = newlines[assignments_start:assignments_end]
    else:
        assignments = []

    # print "readmodel: rates: ", rates 
    
    return(reactions, ics, rates, frozenvars, functions, assignments, fullfile)
