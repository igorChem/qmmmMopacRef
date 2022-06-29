import sys                                                                                
sys.path.append("/home/igorchem/VisMol/easyhybrid/pDynamoMethods") 
import pymp
#-----------------------------------------------------------------------
import os, glob
from commonFunctions import *
from CoreInterface import *
import SimulationsPreset

from pBabel                    import *                                     
from pCore                     import *                                     
from pMolecule                 import *                              
from pMolecule.MMModel         import *
from pMolecule.NBModel         import *                                     
from pMolecule.QCModel         import *
from pScientific               import *                                
from pSimulation               import *

#put the path of the amber tools bin
_path   = "/home/igorchem/programs/amber20/bin"
Reduce  = os.path.join(_path,"reduce")
tleap   = os.path.join(_path,"tleap" )
parmchk = os.path.join(_path,"parmchk2" )
antech  = os.path.join(_path,"antechamber")

#=======================================================================
def TleapPars(pdb):
	'''
	Generate force field parameters using AMBER tleap by passing the pdb file path
	The pdb cannot have residues different from the ones previuos known  by the force field
	Thus, this function works only for the protein without residue. 	
	'''
	#Creating tleap input to save ligand library
	
	pdb2 = pdb[:-4] + "_wh.pdb"	
	
	os.system("/home/igorchem/programs/amber20/bin/reduce -Trim {}".format(pdb) +" > "+pdb2 )
	tleap_in =  "source oldff/leaprc.ff99SB \n"	
	tleap_in += "source leaprc.water.tip3p \n"
	tleap_in += "protein = loadPdb {} \n".format(pdb2)
	tleap_in += "savePdb protein {}\n".format(pdb2)
	tleap_in += "saveamberparm protein " +pdb[:-4]+".top "+ pdb[:-4] +".crd\n"
	tleap_in += "quit"

	tleap_file = open('tleap_in','w')
	tleap_file.write(tleap_in)
	tleap_file.close()
	
	os.system("/home/igorchem/programs/amber20/bin/tleap -f tleap_in")
	
#-----------------------------------------------------------------------
def ParametrizeLig(pdb_file,lig_name):
	'''
	Fast parametrization of organic molecules using the antechamber.
	To do this, we need the pdb with only the ligand molecule and the three-letter residue code.
	'''
	
	command = antech +" -i "+pdb_file+" -fi pdb -o "+lig_name+".mol2 -fo mol2 -c bcc -nc 0 -m 1"
	os.system(command)
	command = parmchk+" -i "+lig_name+".mol2 -f mol2 -o " +lig_name+".frcmod"
	
	input()
	tleap_in =  "source leaprc.gaff2 \n"	
	tleap_in += "lig = loadmol2 {}\n".format(mol_file)
	tleap_in += "savePdb lig {}.pdb\n".format(lig_name)
	tleap_in += "saveoff lig {}.lib.".format(lig_name)
	tleap_in += "quit"

	tleap_file = open('tleap_in','w')
	tleap_file.write(tleap_in)
	tleap_file.close()
	
	os.system( tleap + " -f tleap_in")
	
#-----------------------------------------------------------------------
def TleapPars_Lig(pdb,lig):
	'''
	Generate force field parameters
	'''
	#Creating tleap input to save ligand library
	
	pdb2 = pdb[:-4] + "_wh.pdb"		
	
	os.system("/home/igorchem/programs/amber20/bin/reduce -Trim {}".format(pdb) +" > "+pdb2 )
	tleap_in =  "source oldff/leaprc.ff99SB \n"	
	tleap_in =  "source leaprc.gaff2 \n"
	tleap_in += "source leaprc.water.tip3p \n"
	tleap_in += "protein = loadPdb {} \n".format(pdb2)
	tleap_in += "savePdb protein {}\n".format(pdb2)
	tleap_in += "saveamberparm protein " +pdb[:-4]+".top "+ pdb[:-4] +".crd\n"
	tleap_in += "quit"

	tleap_file = open('tleap_in','w')
	tleap_file.write(tleap_in)
	tleap_file.close()
	
	os.system("/home/igorchem/programs/amber20/bin/tleap -f tleap_in")
#-----------------------------------------------------------------------
def load_system(top,crd):
	'''
	'''
	system = ImportSystem(top)
	system.coordinates3 = ImportCoordinates3(crd, log=None)
	nbmodel = NBModelCutOff.WithDefaults( )
	system.DefineNBModel(nbmodel)
	system.Energy()    
	Pickle(top[:-4]+".pkl",system)	

#-----------------------------------------------------------------------
def mopac_input(pkl,center,size):
	'''
	'''
	proj=SimulationProject( "mopac_refine") 	
	proj.LoadSystemFromSavedProject(pkl)
	methods = ["pm6","pm7","am1","rm1","pm3"]
	
	#setting reaction coordinates for ploting labels
	
	parameters = { "xnbins":1			               ,
				   "ynbins":0			               ,
				   "change_qc_region":True             ,
				   "radius":size                       ,
				   "center":center                     , 
				   "folder":pkl[:-4]+"_results"        ,
				   "charge":0		                   ,
				   "multiplicity":1 	               ,
				   "methods_lists":methods             ,					   
				   "NmaxThreads":1                     ,
				   #"mopac_keywords":["grad qmmm"]     ,
				   "simulation_type":"Energy_Refinement",				   
				   "Software":"mopac"	}

	proj.RunSimulation(parameters)
#=======================================================================
if __name__ == "__main__":
	if sys.argv[1] == "-gp":
		print("generating force field parameters with amber force field")			
		TleapPars(sys.argv[2])
	elif sys.argv[1] == "-lig":
		print("generating force field parameters with amber force field")			
		ParametrizeLig(sys.argv[2],sys.argv[3])
	elif sys.argv[1] == "-ls":
		print("Testing the parameter loading")
		load_system(sys.argv[2],sys.argv[3])
	elif sys.argv[1] == "-mop":
		print("Testing Mopac Refinement")
		center = [ float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]) ] 
		mopac_input(sys.argv[2],center,float(sys.argv[6]))
	
