#!/usr/bin/env python
# coding: utf-8

from rdkit import Chem
from rdkit.Chem import AllChem
import os,sys,string,argparse
from optparse import OptionParser

parser = argparse.ArgumentParser(description="create docking configure file for glide grid.\n")
parser.add_argument('ref_mol_file',metavar='<reference SDF file>',help="Mol2/SDF format reference molecule to define the binding site center.")
parser.add_argument('receptor',metavar='<receptor MAE file>',help="MAE format receptor file.")

args = parser.parse_args()
ref_mol_file = args.ref_mol_file
receptor = args.receptor

if not os.path.exists(ref_mol_file):
   #message = "Sorry, cannot find the "%s" file."
   print("Sorry, cannot find the %s file" % ref_mol_file)
   sys.exit()

dir,ifile = os.path.split(ref_mol_file)
format = ifile.split('.')[-1]
if format == 'mol2' :
   mol = Chem.rdmolfiles.MolFromMol2File(ref_mol_file,sanitize=True,removeHs=False)
elif format == 'sdf' :
   suppl = Chem.SDMolSupplier(ref_mol_file, removeHs=False)
   mol = suppl[0]

n = mol.GetNumAtoms()
x = []
y = []
z = []

for i in range(0,n):
   pos = mol.GetConformer().GetAtomPosition(i)
   x.append(round(pos.x,2))
   y.append(round(pos.y,2))
   z.append(round(pos.z,2))
center_x = (min(x) + max(x))/2
center_y = (min(y) + max(y))/2
center_z = (min(z) + max(z))/2
size_x = 30
size_y = 30
size_z =30
print('USECOMPMAE YES')
print('INNERBOX 10, 10, 10')
print('ACTXRANGE 30.000000')
print('ACTYRANGE 30.000000')
print('ACTZRANGE 30.000000')
print('GRID_CENTER %s,%s,%s' %(center_x,center_y,center_z))
print('OUTERBOX 30.000000, 30.000000, 30.000000')
print('GRIDFILE grid.zip')
print('RECEP_FILE',receptor)
