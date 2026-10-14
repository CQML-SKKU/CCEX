#!/usr/bin/env python
import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../tools/SpinTensorReader"))
from spHamilReader.vasp.extract_spinHamils_vasp import extSpinHamilTensors

if len(sys.argv) != 3:
	print("[1] OUTCARfilename")
	print("[2] savefilename")
	sys.exit("Check above")

foutcar = sys.argv[1]
fsavefile = sys.argv[2]

Dtensordict = {}

print(f"OUTCAR :  {foutcar}")
dat = extSpinHamilTensors(foutcar)
Dtensordict["st"] = dat.Ddata


f = open(fsavefile,'w')

for k,v in Dtensordict.items():
	f.write(f"---------------------------------\n")
	f.write(f"D tensor (MHz) = \n")
	f.write(f"   {v['ZFS'][0,0]:>10.3f}    {v['ZFS'][0,1]:>10.3f}    {v['ZFS'][0,2]:>10.3f}\n")
	f.write(f"   {v['ZFS'][1,0]:>10.3f}    {v['ZFS'][1,1]:>10.3f}    {v['ZFS'][1,2]:>10.3f}\n")
	f.write(f"   {v['ZFS'][2,0]:>10.3f}    {v['ZFS'][2,1]:>10.3f}    {v['ZFS'][2,2]:>10.3f}\n")
	f.write(f"\n")
	f.write(f"D tensor (diag) = \n")
	f.write(f"   {v['D_diag'][0]:>10.3f}    {v['D_diag'][1]:>10.3f}    {v['D_diag'][2]:>10.3f}\n")
	f.write(f"\n")
	f.write(f"Eigenvectors = \n")
	f.write(f"   {v['D_eigvec'][0][0]:>10.3f}    {v['D_eigvec'][0][1]:>10.3f}    {v['D_eigvec'][0][2]:>10.3f}\n")
	f.write(f"   {v['D_eigvec'][1][0]:>10.3f}    {v['D_eigvec'][1][1]:>10.3f}    {v['D_eigvec'][1][2]:>10.3f}\n")
	f.write(f"   {v['D_eigvec'][2][0]:>10.3f}    {v['D_eigvec'][2][1]:>10.3f}    {v['D_eigvec'][2][2]:>10.3f}\n")
	f.write(f"\n")
	f.write(f"ZFS parameters = \n")
	f.write(f"   D = {v['D']:>10.3f}\n")
	f.write(f"   E = {v['E']:>10.3f}\n")
	f.write(f"\n")
f.write(f"---------------------------------\n")
f.close()
