import os.path
import os
import shutil


def make_install():
     
    ipathbase="Install-pycellerator-v-1.1"
    i=0
    ipath=ipathbase
    while os.path.exists(ipath):
        i+=1
        ipath=ipathbase+"-"+str(i)
     
    os.makedirs(ipath)
    print "folder",os.path.abspath(ipath),"created"
    f="how-to-install-pycellerator.txt"
    src=os.path.abspath(f)
    tgt=os.path.join(ipath,f)
    with open(tgt,"w") as F:
        F.write("Installation Instructions\n")
        F.write("-------------------------\n")
        F.write("Copy the entire folder `pycellerator' to your hard drive\n")
    #shutil.copy(src,tgt)    
    
    # create destination folders
    
    cpath=os.path.join(ipath,"pycellerator", "cellerator")
    mpath=os.path.join(ipath,"pycellerator","models")
    top=os.path.join(ipath,"pycellerator")
    os.makedirs(cpath)
    os.makedirs(mpath)  
    
    print "folder",os.path.abspath(cpath),"created"
    print "folder",os.path.abspath(mpath),"created"
    
    # copy documentation
   
    doc = os.path.join("doc","latex", "pycellerator.pdf")
    src = os.path.abspath(doc)
    tgt = os.path.join(ipath,"pycellerator","pycellerator.pdf")
    shutil.copy(src,tgt)
    print src,"--->",tgt
    
    # COPY models
   
    modelsourcepath=os.path.abspath("models")
    modelfiles = os.listdir(modelsourcepath)
    for modelfilename in modelfiles:
        if modelfilename[-6:]==".model":
            src=os.path.join(modelsourcepath, modelfilename)
            tgt=os.path.join(os.path.abspath(mpath),modelfilename)
            shutil.copy(src,tgt)
            print src,"--->",tgt
    
    # copy python code
    
    pypath=os.path.abspath("cellerator")
    pyfiles = os.listdir(pypath)
    for pyfile in pyfiles:
        if pyfile in ["Tissue2D.py","cellzillafy.py"]:
            continue # ignore these files
        if pyfile[-3:]==".py":
            src=os.path.join(pypath, pyfile)
            tgt=os.path.join(os.path.abspath(cpath),pyfile)
            shutil.copy(src,tgt)
            print src,"--->",tgt
    # copy top level files
            
    topfiles=["readme.txt","pycellerator.py","LICENSE","Gold1.model","demo.ipynb"]
    for f in topfiles:
        src=os.path.abspath(f)
        tgt=os.path.join(top,f)
        if os.path.exists(src):
            shutil.copy(src,tgt)
            print src,"--->",tgt
        else:
            print src,"not found"
  
            


if __name__ == "__main__":

    make_install()
