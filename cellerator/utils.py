#!/usr/bin/env python
#
# Cellerator Utilities Module
# B.E.Shapiro 1 March 2012 
# Revised 11/11/14 to allow updates in listtodict
# Varions corrections 2-8 Apr 2015
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

import os.path
import datetime


#*******************************************************************************  
def ensureList(x):
    y = x if type(x)==type([]) else [x]
    return y

#*******************************************************************************  

def uppercase(string): 
    return string.upper()

#*******************************************************************************  
#
# changes a string of the form abc(i,j,k)
# into a string of the form abc_i_j_k
#

    
def deindex(species):
    s=str(species).strip()
    if "(" in s:
        sp,index = s.split("(")
        index = (index[:-1]).replace(",","_")
        s = sp+"_"+index
    return s    

#*******************************************************************************  

def uniqueFileName(filename, type=""):
    # make sure the file name is unique (new)
    
    nameparts = filename.rsplit(".",1)
    
    # determine the base file name

    if len(nameparts)<2:
        nameparts.append(type)
    left, right = nameparts

    # generate a unique output file name
    i = 1
    fname = left+"."+right
    while os.path.isfile(fname):
        fname = left+str(i)+"."+right
        i += 1
    
    return(fname)
#*******************************************************************************  
#
# timed_file_name("filname","ftype") -->fname-yymmdd-hhmm.ftype
# then passes that through uniqueFileName
#

def timed_file_name(base,ftype):
    now=datetime.datetime.now().strftime("%y%m%d-%H%M")
    f=base+"-"+now
    f=f.replace(".","_")
    f=f.replace("-","_")
    f=f+"."+ftype
    f=f.replace("..",".")
    f=uniqueFileName(f)
    return(f)
    
#*******************************************************************************************

def atomic(x):
    return type(x) in [type("x"), type(1), type(1.0)]
    
def deatomize(x):
    if atomic(x): return([x])
    return x
#*******************************************************************************

def evenq(x):
    return ((type(x)==type(1)) and (x % 2 == 0))
    
def oddq(x):
    return not(evenq(x))
    
    

#*******************************************************************************
def bracketsmatch(s):
    i=0; l=0; r=0
    for c in s:
        if c=='[':
            i+=1
            l+=1
        elif c==']':
            i-=1
            r+=1
    return (i==0, l, r)        

#*******************************************************************************
def sumlist(z):
#   simple summing doesn't work on sympy lists
    thesum = 0
    for x in z:
        thesum = thesum + x
    return thesum
    
#*******************************************************************************

def prodlist(z):
    theprod = 1
    for x in z:
        theprod = theprod * x
    return theprod
    
#*******************************************************************************

 
def flatten(x):
    result = []
    for el in x:
        if hasattr(el, "__iter__") and not isinstance(el, basestring):
            result.extend(flatten(el))
        else:
            result.append(el)
    return result
    
#*******************************************************************************

def isanum(x):
    try:
        float(x)
        return True
    except:
        return False

def list2Dict(l):
    d={}
    for pair in l:
        (a,b) = pair
        if a in d:
            if isanum(b) and isanum(d[a]):
                d[a]+= b
            else:
                d[a]=b
        else:
            d[a]=b
    return d

#*******************************************************************************

def getlast(l,default):
    if len(l)>0:
        return l[-1]
    else:
        return default

#*******************************************************************************
        
def listify(l):
    if type(l)==type("Hello World"):
        return [l]
    else:
        return list(l)    

#*******************************************************************************

def Power(y,x):
    if y>0: return y**x
    return 0
    
#*******************************************************************************
    
def is_number(s):
    # print "DBG: is_number: ", s
    try:
        float(str(s))
        return True
    except ValueError:
        return False
#*******************************************************************************
def not_a_number(s):
    try: 
        float(s)
        return False
    except:
        return True

#*******************************************************************************

def is_one(s):
    if not is_number(s):
        return False
    return float(s)==1

#*******************************************************************************

def is_minusone(s):
    if not is_number(s):
        return False
    return float(s)==-1
    
#*******************************************************************************

def reac(lhs, rhs, k):
    # generates a simple text reaction of the form [A + B + C + ... -> P + Q + R + ..., k]
    # input: lhs = list containing the [A, B, ...]
    #        rhs = list containing the [P, Q, ...]
    #        k = rate constant
    S, P = lhs, rhs
    # print "reac:S, P, k = ", S, P, k

    if not (type(S) == type([])): S = [S]
    if not (type(P) == type([])): P = [P]
    S = flatten(S)
    P = flatten(P)
    
    S = "+".join(map(str, S))
    P = "+".join(map(str, P))
    r = "["+S+"->"+P+","+str(k)+"]"
    return r

#*******************************************************************************
    
def revreacs(lhs,rhs,k):
   [kf,kr]=k
   r1 = reac(lhs, rhs, kf)
   r2 = reac(rhs, lhs, kr)
   return [r1, r2]

    
#*******************************************************************************
    
def enzreac(S, P, E, rates):
    r = "["+str(S)+"=>"+str(P)+", mod["+str(E)+"], rates["+( ",".join(map(str,rates)))+"] ]"   
    return r    
    
