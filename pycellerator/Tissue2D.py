# alpha code not ready for use yet
# this code is not required for pycellerator -- it is developmental
# and may be changed or withdrawn. It is not used pycellerator but parts
# of it may be used in the future
#
#  Tissue2D.py is NOT USED by pycellerator
#  it is developmental code here for fposible future incorporation
#
#****************************************************************************
#    pycellerator converts reactions, expressed in a text-formatted arrow-based 
#    notation into differential equations and performs numerical simulations.
#
#****************************************************************************
#
#    Copyright (C) 2012-2015
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

from sys import argv,exit

from matplotlib.pyplot import *

#
# alpha code for multicellular models
#


####################################################################
#
#
#  
class Tissue2D(object):
    """Tissue2D T.vertices = [[x,y], [x,y], ... ]
            T.cells = [cell1, cell2, ...], each cell i is a list
                of cell vertex indices

    Initializaton: Tissue2D(vertices, cells)
                   
            Computed automatically:
            
            T.cell_edges = [cell1, cell2, ...] each cell is a list of 
                of edge indices
            T.edges = [edge1, edge2, ...] each edge is  a list of vertic
                indices
            T.edge_lengths = [length1, length2, ...]
            T.areas = [area1, area2, ...]
                
    Methods: 
            T.edge_length(i) - length of edge i
               equivalent to (T.edge_lengths)[i]
            T.cell_area(i) - area of cell i
               equivalent to (T.areas)[i]
                       
    """
    # determine the edges from the cell list
    #
    def determine_edges(self):
        nv=len(self.vertices)
        nc=len(self.cells)
        def pairedges(cell):
            second=cell[1:]
            second.append(cell[0])
            return( zip(cell, second) )           
        e=[]
        cell_edges=[]
        cell_number=-1
        for cell in self.cells:
            cell_number += 1
            cell_edge_list=[]
            edgeset=pairedges(cell)
            for edge in edgeset:
                x,y=edge
                if x>y:
                    y,x=x,y
                next_potential_edge=[x,y]
                if next_potential_edge in e:
                    k=e.index(next_potential_edge)
                else:
                    e.append(next_potential_edge)
                    k=len(e)-1
                cell_edge_list.append(k)
            cell_edges.append(cell_edge_list)               
        return(e, cell_edges)
        
    # determine length of edge i
    def edge_length(self, i):
        edgenumber=self.edges[i]
        v1,v2 = edgenumber
        v1 = self.vertices[v1]
        v2 = self.vertices[v2]
        x1,y1=v1
        x2,y2=v2
        return ((x1-x2)**2.0 + (y1-y2)**2.0)**0.5
    
    # determine lengths of all edges    
    def determine_edge_lengths(self):
        n = len(self.edges)
        lengths=[]
        for i in range(n):
            length = self.edge_length(i)
            lengths.append(length)
        return(lengths)
    
    # get cell i as list of vertex coordinates
    def list_cell_vertex_coordinates(self, i):
        cell=self.cells[i]
        vertex_list = []
        for v in cell:
            x, y = self.vertices[v]
            vertex_list.append([x,y])
        return(vertex_list)
        
    #get cell area of cell i
    # Gauss' "Shoelace Formula"
    def cell_area(self, i):
        vertex_list = self.list_cell_vertex_coordinates(i)
        first_vertex=vertex_list[0]
        vertex_list.append(first_vertex)
        area = 0
        
        x1, y1=vertex_list[0]
        for x2, y2 in vertex_list[1:]:
            area += (x1*y2 - y1*x2)  # add determinant
            x1, y1, = x2, y2
        area = abs(0.5*area)
        return area
        
    # determine areas of all cells
    def cell_areas(self):
        areas=[]
        for i in range(len(self.cells)):
            area=self.cell_area(i)
            areas.append(area)
        return(areas)
        
    # determine neighbors of a cell
    def neighbors_of(self, i):
        neighbors=[]
        me =(self.cell_edges)[i]
        ncells = len(self.cells)
        A = (self.areas)[i]
        for edge in me:
            for j in range(ncells):
                if j==i:
                    continue
                them = (self.cell_edges)[j]
                if edge in them:
                    el = (self.edge_lengths)[edge]
                    neighbors.append([j, edge, el/A])
        return neighbors
    #
    # Connection list returns a two dimensional list of integers
    # CL[i] is a list that gives the neighbors of cell i in the form
    # [[nbr1, edge1, el_1/A_i], [nbr2, edge2, el_2/A_i], [nbr3, edge3, el_3/A_i],...]
    # where each nbr is a cell number
    # and each edge is the shared edge
    # the el_i is the length of the shared edge
    # the A_i is the area of cell A
    # Diffusion is (d/dt)[X_i] = beta*(el_j/A_i)*([X_j]-[X_i])
    # where beta is wall permeability e.g. in cm/sec
    # [ x[i] <-> Nil, beta*(l_j/A_i), beta*(l_j/A_i) X[j] ]
    # 
    def connection_list(self):
        connections=[]
        ncells=len(self.cells)
        for j in range(ncells):
            nbrs=self.neighbors_of(j)
            connections.append(nbrs)
        return(connections)
        
    # show tissue with a specific color scheme
    #
    def show_cell(self, i, color):
        c=(self.cells)[i]
        xypairs=[(self.vertices)[vertex] for vertex in c]
        x,y=zip(*xypairs)
        fill(x,y, color=color)
    #
    def show(self, colors):
        nc=len(self.cells)
        rn=range(nc)
        ln=len(colors)
        figure()
        axes().set_aspect("equal")   
        for i,color in zip(rn,colors):
            self.show_cell(i,color)
        show()
         
        
    # initialize a new tissue object    
    def __init__(self, v, c):
        self.vertices=v
        self.cells=c
        
        list_of_edges, cells_as_edges=self.determine_edges()
        
        self.edges=list_of_edges
        self.cell_edges=cells_as_edges
        self.edge_lengths=self.determine_edge_lengths()
        self.areas=self.cell_areas()
        
##########################################################

def readGeometryFile(gfile):
    f=open(gfile,"r")
    geometry=f.readlines()
    f.close()
    dimensions = geometry[0]
    n_vertices =int((geometry[1].strip().split(","))[0])
    n_faces =int((geometry[2].strip().split(","))[0])
    geometry=geometry[3:]
    vertices=geometry[:n_vertices]
    v=[]
    for xypair in vertices:
        v.append(map(float, xypair.strip().split(",")))
    cells=[]
    faces=geometry[n_vertices:]
    
    if "STARTSAT1" in argv:
        for face in faces:
             vertex_numbers = map(int,face.strip().split(","))
             cells.append([p-1 for p in vertex_numbers])
       
    else:   
        for face in faces:
            cells.append(map(int,face.strip().split(",")))
    
    T=Tissue2D(v, cells)
    
    return T

#####################################################################

def RGBColor(value, scalemin, scalemax):
    if scalemin > scalemax:
        scalemin, scalemax = scalemax, scalemin
    width = scalemax - scalemin
    color = (value - scalemin)/width
    color = max(0, color)
    color = min(1, color)
    R=1; G=0; B=0; 
    if color<.33333:
        R=color/.33333
    elif color <.66666:
        DR=color-.33333
        R=1-3*DR
        G=DR/.33333
    else:
        DG=color-.66666
        R=0.0
        G=1-3*DG
        B=DG/.333
    return (R, G, B)
 
#####################################################################
   
def RGBInterpolate(value, scalemin, scalemax, RGBMin, RGBMax):
    if scalemin > scalemax:
        scalemin, scalemax = scalemax, scalemin
    width = scalemax - scalemin
    color = (value - scalemin)/width
    color = max(0, color)
    color = min(1, color)
    (R1,G1,B1)=RGBMin
    (R2,G2,B2)=RGBMax
    R=R1+color*(R2-R1)
    G=G1+color*(G2-G1)
    B=B1+color*(B2-B1)
    return (R,G,B)
    
       
def display_geometry(T):
    vertices=T.vertices
    cells=T.cells


    if "-RAINBOW" in argv:
        def center(c):
            vc = [vertices[j] for j in c]
            x,y=zip(*vc)
            n = len(x)*1.0
            return( sum(x)/n, sum(y)/n )
        xvals,yvals = zip(*vertices)
        xmin = min(xvals)
        xmax = max(xvals)
        xs = [center(c)[0] for c in cells]
        colors = [RGBColor(x, xmin, xmax) for x in xs]
        T.show(colors)
        return()
            
    nv = len(vertices)
    nc = len(cells)
    
    fill_polys="-FILL" in argv
 
    
    def drawcell(c):
        xcell=[]; ycell=[];
        for vertex in c:
            x,y=vertices[vertex]
            xcell.append(x)
            ycell.append(y)
        if fill_polys:
            fill(xcell, ycell)
        else:
            plot(xcell, ycell, "b")
 
    axes().set_aspect("equal")   
    for cell in cells:
        drawcell(cell)
    show()
    
##################################################################  

#
#
# -GEOMETRY filename [-DISPLAY [-NOFILL]]
#    
if __name__ == "__main__":
    print "*** DEBUGGING MODE -- CALL FROM MAIN PROGRAM ***"
    args=argv
    upper_args = [s.upper() for s in args]
    if "-TISSUE" in upper_args:
        i = upper_args.index("-TISSUE")
        if i+1<len(args):
            gfile = args[i+1]
            T = readGeometryFile(gfile)

            if "-DISPLAY" in upper_args:
                display_geometry(T)
            if "-LISTEDGES" in upper_args:
                k=0
                for x in T.edges:
                    print k, x; k+=1
            if "-LISTCELLEDGES" in upper_args:
                k=0
                for x in T.cell_edges:
                    print k, x; k+=1
            if "-LISTEDGELENGTHS" in upper_args:
                k=0
                for x in T.edge_lengths: 
                    print k, x; k+=1
            if "-LISTAREAS" in upper_args: 
                k=0; AT=0
                for x in T.areas:
                    print k, x
                    AT += x; k+= 1
                print "total area = ", AT
            if "-CONNECTIONS" in upper_args:
                k=0; 
                CL = T.connection_list()
                for x in CL:
                   print x
                    
        else:
            exit("-GEOMETRY file_name expected")
    else:
        exit("Huh?")
        
            

