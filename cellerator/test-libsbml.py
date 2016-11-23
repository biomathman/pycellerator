
from libsbml import *
reader = SBMLReader()
SBML = reader.readSBML("../../BMDB/BM0001.xml")
nerrors = SBML.getNumErrors()
if nerrors > 0:
    errorLog = SBML.getErrorLog()
    print errorLog
    raise SystemExit("Bad SBML File")

level = SBML.getLevel()
version = SBML.getVersion()

print "Level ",level, " Version ", version

model = SBML.getModel()
print model

LOC = model.getListOfCompartments()
n = len(LOC)
print n, " compartments" 
for i in range (n):
    ctype = "?"
    
    c = model.getCompartment(i)
    ctype = c.getCompartmentType()
    print c.getId()," ", c.getName()," ", c.getSize(), " ", c.getVolume(),\
          " ", c.getUnits(),\
          c.getSpatialDimensions(), " '", ctype,"'"
