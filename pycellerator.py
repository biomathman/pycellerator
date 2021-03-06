#
#****************************************************************************
#    pycellerator converts reactions, expressed in a text-formatted arrow-based 
#    notation into differential equations and performs numerical simulations.
#
#****************************************************************************
#
#    Copyright (C) 2012-15 Bruce E Shapiro
#    8/30/15 Set version number to 1.0 from alpha sequence for release
#    on github
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
import cellerator.converter
import cellerator.expander
import cellerator.interpreter
import cellerator.reader
import cellerator.solver
import cellerator.parser
import cellerator.sbml
import cellerator.flux
import sys
import subprocess

versionNumber = 1.0

#
#****************************************************************************
#
def solve():
    f = cellerator.solver.makesolver()
    if "-norun" in sys.argv: return
    subprocess.call("python "+f, shell=True)
    return
    
def version():
    print "This is py[cellerator] version "+str(versionNumber)
    return
#
#****************************************************************************
#
handlers={"EXPAND":   cellerator.expander.expander,
         "INTERPRET": cellerator.interpreter.interpreter,
         "CONVERT":   cellerator.converter.converter,
         "PARSE":     cellerator.parser.runpyx,
         "SOLVE":     solve,
         "VERSION":   version, 
         "SBML":      cellerator.sbml.runsbml,
         "FLUX":      cellerator.flux.fat
         }

#
#****************************************************************************
#
def runpycellerator(args):

    if type(args)==type(""):
        arguments = filter((lambda x:len(x)>0), args.split(" "))

        sys.argv=[]
        for arg in arguments:
            sys.argv.append(arg)
    else:
        arguments = args 
                
    n=1
    ARGS = map(lambda x:x.upper(), arguments)
    # print "ARGS = ", ARGS
    for option in handlers:
        # print " checking ",option
        if option in ARGS:
            # print "option is ", option
            f = handlers[option]
            f()
            break
        n+=1
    if n> len(handlers):
        print "Unknown or missing py[cellerator] option reqested.\n"+\
            "Valid options are:",
        for opt in handlers: print opt+",",
        
    return
#
#****************************************************************************
#
if __name__ == "__main__":

    args = sys.argv 
    runpycellerator(args)
    
