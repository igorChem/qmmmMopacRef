import sys                                                                                
sys.path.append("/home/igorchem/pDynamo3_scripts") 
import pymp
#-----------------------------------------------------------------------
import os, glob, math
from SimulationProject import *
from Simulation import * 

from pBabel                    import *                                     
from pCore                     import *                                     
from pMolecule                 import *                              
from pMolecule.MMModel         import *
from pMolecule.NBModel         import *                                     
from pMolecule.QCModel         import *
from pScientific               import *                                
from pSimulation               import *

#put the path of the amber tools bin
_path   = "/home/igorchem/programs/amber_20/bin"
Reduce  = os.path.join(_path,"reduce")
tleap   = os.path.join(_path,"tleap" )
parmchk = os.path.join(_path,"parmchk2" )
antech  = os.path.join(_path,"antechamber")


#=======================================================================
def SS_bonds(pdb):
	'''
	'''
	text_lines = ""
	new_pdb    = "" 
	pdb2       = pdb[:-4] + "_cyx.pdb"
	ss_bonds   = False
	#-------------------------------------------------------------------
	pdb_file = open(pdb,"r")
	lines = pdb_file.readlines()
	cyx_residues = []
	res_number   = []
	j_coord_list = []
	# ------------------------------------------------------------------
	for line in lines:
		if line.__contains__('SG  CYS'):
			ss_bonds = True
			resnumber = int(line[22:30])
			x = float(line[31:38])
			y = float(line[39:46])
			z = float(line[47:54])
			cyx_residues.append([resnumber,x,y,z])
			res_number.append(resnumber)
		else:
			new_pdb += line
	#-------------------------------------------------------------------
	for i in cyx_residues:
		for j in cyx_residues:
			# ---
			if i[0] not in j_coord_list:
				dst = math.sqrt((i[1]-j[1])**2 + (i[2]-j[2])**2 + (i[3]-j[3])**2)
			else:
				continue
			# ---
			if dst > 0 and dst < 2.5 and i[0]:
				j_coord_list.append(j[0])
				text_lines+= "bond protein.{}.SG protein.{}.SG\n".format(i[0],j[0])
	#-------------------------------------------------------------------
	return(text_lines,ss_bonds,res_number)

#=======================================================================
def Prepare_PDB(pdb):
	'''
	'''
	pdb2 = pdb[:-4] + "_wh.pdb" # pdb without hydrogens
	pdb3 = pdb[:-4] + "_am.pdb" # pdb treated by amber 
	pdb4 = pdb[:-4] + "_cyx.pdb" # pdb com cyx residues
	os.system("/home/igorchem/programs/amber_20/bin/reduce -Trim {}".format(pdb) +" > "+pdb2 )
	
	tleap_in =  "source leaprc.protein.ff14SB \n"
	tleap_in += "source leaprc.water.tip3p \n"
	tleap_in += "protein = loadPdb {} \n".format(pdb2)
	tleap_in += "savePdb protein {} \n".format(pdb3)
	tleap_in += "quit"
	
	tleap_file = open('tleap_in','w')
	tleap_file.write(tleap_in)
	tleap_file.close()
	
	os.system("/home/igorchem/programs/amber_20/bin/tleap -f tleap_in")
	
	txt_cyx = SS_bonds(pdb3)
	
	if txt_cyx[1]:
		cyx_residues = txt_cyx[2]
		new_pdb      = ""
		pdb_file     = open(pdb3,'r')
		for line in pdb_file:
			if line.__contains__("CYS"):
				if int(line[22:30]) in cyx_residues:
					nline = line.replace("CYS","CYX")
					new_pdb+=nline					
				else:
					new_pdb+=line
			else:
				new_pdb+=line
		new_pdb_file = open(pdb4,"w")
		new_pdb_file.write(new_pdb)
		new_pdb_file.close()
		return(txt_cyx[0])
	else:
		print(txt_cyx[1])
		
#=======================================================================
def TleapPars(pdb):
	'''
	Generate force field parameters using AMBER tleap by passing the pdb file path
	The pdb cannot have residues different from the ones previuos known  by the force field
	Thus, this function works only for the protein without residue. 	
	'''
	#Creating tleap input to save ligand library
	
		
	pdb2 = pdb[:-4] + "_cyx.pdb"
	pdb3 = pdb[:-4] + "_wh2.pdb"
	pdb4 = pdb[:-4] + "_tleap.pdb"
	pdb5 = pdb[:-4] + "_am.pdb"
	info = Prepare_PDB(pdb)

	tleap_in = ""
	tleap_in += "source leaprc.protein.ff14SB \n"
	tleap_in += "source leaprc.water.tip3p \n"
	if info:
		os.system("/home/igorchem/programs/amber_20/bin/reduce -Trim {}".format(pdb2) +" > "+pdb3 )
	
		tleap_in += "protein = loadPdb {} \n".format(pdb3)
		tleap_in += info
		tleap_in += "savePdb protein {}\n".format(pdb4)
		tleap_in += "saveamberparm protein " +pdb[:-4]+".top "+ pdb[:-4] +".crd\n"
		tleap_in += "quit"
	else:		
		tleap_in += "protein = loadPdb {} \n".format(pdb5)
		tleap_in += "savePdb protein {}\n".format(pdb4)
		tleap_in += "saveamberparm protein " +pdb[:-4]+".top "+ pdb[:-4] +".crd\n"
		tleap_in += "quit"
		

	tleap_file = open('tleap_in','w')
	tleap_file.write(tleap_in)
	tleap_file.close()
	
	os.system("/home/igorchem/programs/amber_20/bin/tleap -f tleap_in")
	
#-----------------------------------------------------------------------
def ParametrizeLig(pdb_file,lig_name):
	'''
	Fast parametrization of organic molecules using the antechamber.
	To do this, we need the pdb with only the ligand molecule and the three-letter residue code.
	'''
	
	command = antech +" -i "+pdb_file+" -fi pdb -o "+lig_name+".mol2 -fo mol2 -c bcc -nc 0 -m 1"
	os.system(command)
	command = parmchk+" -i "+lig_name+".mol2 -f mol2 -o " +lig_name+".frcmod"
	
	
	tleap_in =  "source leaprc.gaff2 \n"	
	tleap_in += "lig = loadmol2 {}.mol2\n".format(lig_name)
	tleap_in += "savePdb lig {}.pdb\n".format(lig_name)
	tleap_in += "saveoff lig {}.lib\n".format(lig_name)
	tleap_in += "quit"

	tleap_file = open('tleap_in','w')
	tleap_file.write(tleap_in)
	tleap_file.close()
	
	os.system( tleap + " -f tleap_in")
	
#-----------------------------------------------------------------------
def TleapPars_Lig(pdb,lig=None):
	'''
	Generate force field parameters
	'''
	#Creating tleap input to save ligand library
	
		
	tleap_in =  "source leaprc.ff19SB \n"
	tleap_in =  "source leaprc.gaff2 \n"
	tleap_in += "source leaprc.water.tip3p \n"
	tleap_in += "{} = loadmol2 {}.lib\n".format(lig,lig)
	tleap_in += "complex = loadPdb {} \n".format(pdb)
	
	tleap_in += "saveamberparm complex " +pdb[:-4]+".top "+ pdb[:-4] +".crd\n"
	tleap_in += "quit"

	tleap_file = open('tleap_in','w')
	tleap_file.write(tleap_in)
	tleap_file.close()
	
	os.system("/home/igorchem/programs/amber_20/bin/tleap -f tleap_in")
	
#-----------------------------------------------------------------------
def mopac_input(top,crd,center,size):
	'''
	'''
	proj = SimulationProject.From_Force_Field(top,crd)
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

#-----------------------------------------------------------------------
def Treat_PDB_bind(path):
	'''
	Make all process of the functions above to test the methods for PDB_Bind data set
	'''
	path2    = path + "/*/*_protein.pdb"
	pdb_list = glob.glob(path2)
	
	error_list = [] 
	for pdb_ in pdb_list:
		try:
			top_ = pdb_[:-4] + ".top"
			crd_ = pdb_[:-4] + ".crd"
			#pkl_ = pdb_[:-4] + ".pkl"
			TleapPars(pdb_)
			load_system(top_,crd_)
			#mopac_input(pkl_)
		except:
			error_list.append(pdb_)
			pass

#-----------------------------------------------------------------------
def Treat_PDB_bind_ligand(path):
	'''
	Make all process of the functions above to test the methods for PDB_Bind data set
	with ligand
	'''
	path2    = path + "/*/*_protein.pdb"
	pdb_list = glob.glob(path2)
	lig_list = glob.glob(path3)
	
	error_list = [] 
	for pdb_ in pdb_list:
		try:
			top_ = pdb_[:-4] + ".top"
			crd_ = pdb_[:-4] + ".crd"
			TleapPars_Lig(pdb_,lig_)
			mopac_input()
		except:
			error_list.append(pdb_)
			pass


#=======================================================================
if __name__ == "__main__":
	if 	 sys.argv[1] == "-pPDB":
		print("Preparing PDB files")
		Prepare_PDB(sys.argv[2])
	elif sys.argv[1] == "-gp":
		print("Generating force field parameters with amber force field")
		TleapPars(sys.argv[2])
	elif sys.argv[1] == "-lig":
		print("Generating force field parameters with amber force field")
		ParametrizeLig(sys.argv[2],sys.argv[3])
	elif sys.argv[1] == "-ls":
		print("Testing the parameter loading")
		load_system(sys.argv[2],sys.argv[3])
	elif sys.argv[1] == "-mop":
		print("Testing Mopac Refinement")
		center = [ float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]) ] 
		mopac_input(sys.argv[2],center,float(sys.argv[6]))
	elif sys.argv[1] == "-PDBbind":
		print("Testing functions on PDB Binding!")
		Treat_PDB_bind(sys.argv[2])
	
