#!/usr/bin/env python
# coding: utf-8

import os,sys,string,argparse
from optparse import OptionParser

parser = argparse.ArgumentParser(description="create docking configure file.\n")
parser.add_argument('database',metavar='<database SDF file>',help="SDF format database file.")
parser.add_argument('grid',metavar='<docking grid file>',help="GLIDE docking grid file.")
parser.add_argument('method',metavar='<Docking Method>',help="Docking method: confgen, rigid, mininplace, inplace.")
parser.add_argument('precision',metavar='<Docking Precision>',help="Docking mode: SP, Normal, Accurate, HTVS and XP.")
parser.add_argument('npose',metavar='<Number of poses>',help="Maximum number of poses to write per ligand.")

args = parser.parse_args()
dbase = args.database
grid = args.grid
method = args.method
mode = args.precision
n_pose = args.npose

if not os.path.exists(dbase):
   #message = "Sorry, cannot find the "%s" file."
   print("Sorry, cannot find the %s file" % dbase)
   sys.exit()

print('WRITEREPT YES')
print('USECOMPMAE YES')
print('POSTDOCK_NPOSE 10')
print('POSES_PER_LIG %s' % n_pose)
print('MAXREF 800')
print('RINGCONFCUT 2.500000')
print('GRIDFILE %s' % grid)
print('LIGANDFILE %s' % dbase)
print('DOCKING_METHOD %s' %method)
print('PRECISION %s' %mode)
print('POSE_OUTTYPE ligandlib_sd')
