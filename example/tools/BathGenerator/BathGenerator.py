#!/usr/bin/env python
import argparse
import os
import sys
import time
import psutil

import re
import copy
import math
import numpy as np

import numpy.linalg as lin
import random as r
from collections import Counter
from itertools import product

from mpi4py import MPI

###################
#####Functions#####
# 1 : Make the pure configure
# 2 : Make the configure using defect cell
# 3 : Make the lattice deformated configure
# 4 : Convert the Configure to POSCAR file (to check)
###################

##########################
#0. Set the options
##########################


####################################################
####################################################

##########################
#1. Read the argument
##########################
def parse():
    #Create the parser
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                    description="\
Make the Configure(Bath) using the POSCAR file\n \
            ([unit] size : angstrom, prob : %)\n\n \
 Example of simple usage : \n\
  1) Make the pure configure \n\
    BathGenerator.py -o {output-file} -p {prob file} -i {pure POSCAR} --size {size} --Defect {defect} \n\
\n\
  2) Make the configure with defect cell \n\
    BathGenerator.py -o {output-file} -p {prob file} -i {pure POSCAR} -I {defect POSCAR} --size {size} --Defect {defect} \n\
\n\
  3) Make the confgure with lattice deformation \n\
    BathGenerator.py -o {output-file} -p {prob file} -i {pure POSCAR} -I {defect POSCAR} -D {deformated POSCARs} --axis {deformated axiss}  --size {size} --Defect {defect} \n\
\n\
  4) Convert the Configure to POSCAR file (only for check) \n\
    BathGenerator.py -o {output-POSCAR} -i {input Configure} --convert\n\n ")
    
    #outputfile name (txt)
    parser.add_argument('-o','--output', dest='output',type=str,required=True, help='Output-file name of configure (Bath)') #result name

    #probability file name (txt)
    parser.add_argument('-p','--prob',dest='prob',type=str, help='Probability file name') # probfile : [name of iso, probability]

    #inputfile name (POSCAR or Configure)
    parser.add_argument('-i','--input',dest='input',type=str,required=True, help='Input-file of pure POSCAR or Configure (POSCAR or Configure)') 
    
    #defect poscar name (POSCAR)
    parser.add_argument('-I','--Insert',dest='Insert',type=str, help='Input-file of relaxed POSCAR file with defect (relaxed POSCAR)') 

    #lattice deformated poscar name (POSCAR)
    parser.add_argument('-D','--Deformation',dest='Deformation',nargs='*',type=str, help='Input-file of relaxed POSCAR file with defect (relaxed POSCAR)') 
    parser.add_argument('--axis',dest='Deformation_axis',nargs='*',default=None,choices=['x','y','z','xy','xz','yz'], type=str, help='Axis of lattice deformation [{x,y,z} : plane-shape, {xy,xz,yz} : line-shape] (default : None)') #default : z
    #parser.add_argument('--Deformation-axis',dest='Deformation_axis',default="z",choices=['x','y','z','xy','xz','yz'], type=str, help='Axis of lattice deformation [{x,y,z} : plane-shape, {xy,xz,yz} : line-shape] (default : Z-axis line)') #default : z

    #extra options    
    parser.add_argument('--size', dest='size', type=int, help='Number of making configure (default : None --> only 1)') #unit : angstrom
    parser.add_argument('--Defect',dest='Defect',type=str, help='Name of defect position')  #defect atom name
    parser.add_argument('--Defect_center',dest='Defect_center',type=float, default=[0.5,0.5,0.5],nargs=3, help='pivot point of fractional defect position')  #defect atom name
    

    parser.add_argument('--range', dest='range', default=None, type=int, help='Number of making configure (default : None --> only 1)') #default : None
    parser.add_argument('--opt2D', dest='opt2D',action='store_true', help='Flag of making the only 2d-bath (default : False)') #default : False
    parser.add_argument('--test', dest='test', action='store_true', help='Flag of testing the configure using poscar (default : False)') #default : False
    parser.add_argument('--shape', dest='shape', action='store_false', help='Flag of creating the Bath shape with cubic shape (default : True)') #default : True
    parser.add_argument('--convert', dest='convert', action='store_true', help='Flag of converting the Configure to POSCAR for test (default : False)') #default : False
    parser.add_argument('--verbosity', dest='verbosity', action='store_true', help='Flag of printing detail options (default : False)') #default : False
    parser.add_argument('--minbond', dest='minbond', action='store_true', help='Flag of printing minimum bonding of configure (default : False)') #default : False
    
    #analysis the parser
    args = parser.parse_args()
    return args


####################################################
####################################################

##########################
#2.1. Read the inputfile (POSCAR)
##########################
class AtomData:
    def __init__(self):
        self.inputfile=None

        self.LatMult=1          #multiplicity of lattice parameter
        self.Cell_parameter=[]  #poscar cell parameter
        self.Atom=[]        #atom specises
        self.Atom_num=[]    #number of atoms
    
        self.Atom_data={}   #dict of atom data ; {atom name : atom number} 
        self.Atomic_info=[] #total atomic data ; {atom name, atomic position}
        
        self.Defect_info=[] #defect atomic data ; {defect name, atomic position}

    #multitple unit_cell_parameter
    def Multi_cell(self):
        for i in range(len(self.Cell_parameter)):
            for j in range(3):
                self.Cell_parameter[i][j]=self.Cell_parameter[i][j]*self.LatMult

def ReadPOSCAR(inputfile,opt):
    #check the exiting the inputfile
    if not (os.path.isfile(inputfile)) :
        print("\t---------------------------------------------------------------------------------")
        print("\tThere is no file!! : "+str(inputfile))
        print("\tCheck the Option : ",opt)
        print("\t---------------------------------------------------------------------------------\n")
        sys.exit(1)

    #create the atom data from poscar
    POSCARdata=AtomData()
    POSCARdata.inputfile=inputfile

    #read the poscar file
    with open(inputfile, 'r') as f:

        CheckAtomic=False
        temp_atomic=[]
        for i, line in enumerate(f):
            inputdata=line.strip().split()
            if i == 1: #multiple parameter
                POSCARdata.LatMult=list(map(float,inputdata))[0]
                continue
            elif i >= 2  and i <= 4: #cell parameter
                POSCARdata.Cell_parameter.append(list(map(float,inputdata)))
                continue
            elif i == 5: #atom specise
                POSCARdata.Atom.append(inputdata)
                continue
            elif i == 6: #num of atoms
                POSCARdata.Atom_num.append(list(map(int,inputdata)))
                continue

            if 'Direct' in line or 'direct' in line :
                CheckAtomic=True
                continue
            if (CheckAtomic): #atomic position (fracation corrodination)
                temp_atomic.append(list(map(float,inputdata)))
                
    f.close()

    #multiply the lattice multiple parameter
    if POSCARdata.LatMult != 1:
        POSCARdata.Multi_cell()
    
    #form 2D array to 1D array
    POSCARdata.Atom=sum(POSCARdata.Atom,[])
    POSCARdata.Atom_num=sum(POSCARdata.Atom_num,[])

    #check the errors
    if len(POSCARdata.Atom) != len(POSCARdata.Atom_num):
        print("\t---------------------------------------------------------------------------------")
        print("\tAtom :",POSCARdata.Atom)
        print("\tAtom num :",POSCARdata.Atom_num)
        print("\tPlz, check the atom species and number!!")
        print("\t---------------------------------------------------------------------------------\n")
        sys.exit(1)
    elif len(POSCARdata.Cell_parameter) != 3:
        print("\t---------------------------------------------------------------------------------")
        print("\tCell_parameter :",POSCARdata.Cell_parameter)
        print("\tPlz, check the cell parameter!!")
        print("\t---------------------------------------------------------------------------------\n")
        sys.exit(1)
    elif len(temp_atomic) != sum(POSCARdata.Atom_num):
        print("\t---------------------------------------------------------------------------------")
        print("\tLength of Atomic position :",len(temp_atomic))
        print("\tTotal Atom num :",sum(POSCARdata.Atom_num))
        print("\tPlz, check the cell informatioin number!!")
        print("\t---------------------------------------------------------------------------------\n")
        sys.exit(1)

    #create the atomic positions
    line=0
    for i , count in enumerate(POSCARdata.Atom_num):
        Atom=POSCARdata.Atom[i]
        for j in range(count):
            temp_list=sum([[Atom],temp_atomic[line]],[])
            POSCARdata.Atomic_info.append(temp_list)
            line+=1

    #for j, Atom in enumerate(POSCARdata.Atom):
    #    for i in range(0, POSCARdata.Atom_num[j]):
    #        temp_list=sum([[Atom], temp_atomic[i]],[])
    #        POSCARdata.Atomic_info.append(temp_list)
    #        print(temp_list)
    
    #making dict of atom data
    for i,name in enumerate(POSCARdata.Atom):
        POSCARdata.Atom_data[name]=POSCARdata.Atom_num[i]

    return POSCARdata


##########################
#2.2. Read the probability file
##########################
class ProbData:
    def __init__(self):
        self.inputfile=None

        self.prob_name=[]   #Atom name
        self.prob_num=[]    #Nuclear spin of atom
        self.prob_spin=[]   #prob of nuclear spin (unit : %)

        self.prob={}        #dict of prob {atom name : atom prob}

        self.iso_name = {}  #dict of atom and iso  {atom name : iso name} (__Nope__ --> fail)
        self.iso_prob = {}  #dict of atom and prob {atom anem : prob of iso}

def ReadProb(inputfile, opt):
    #check the exiting the inputfile
    if not (os.path.isfile(inputfile)) :
        print("\t---------------------------------------------------------------------------------")
        print("\tThere is no file!! : "+str(inputfile))
        print("\tCheck the Option : ",opt)
        print("\t---------------------------------------------------------------------------------\n")
        sys.exit(1)

    #create the prob class
    Prob=ProbData()
    Prob.inputfile=inputfile

    #read the prob file
    with open(inputfile,'r') as f:
        for m, line in enumerate(f):
            data=line.split()
            if m % 3 == 0: #Atom name
                Prob.prob_name.append(data)
            elif m % 3 == 1: #Nuclear spin of atom
                Prob.prob_spin.append(data)
            elif m % 3 == 2: #prob of nuclear spin (unit : %)
                Prob.prob_num.append(list(map(float,data)))
        if(len(Prob.prob_num)+len(Prob.prob_spin)+len(Prob.prob_name))%3 != 0:
            print("\tPlz, check the probability fiel and number error !!\n")
    f.close()
   
    #form 2D array to 1D array
    Prob.prob_name = sum(Prob.prob_name,[])
    Prob.prob_spin = sum(Prob.prob_spin,[]) #doesn't use it
    Prob.prob_num  = sum(Prob.prob_num,[])

    #making dict of probability
    for i,name in enumerate(Prob.prob_name):
        Prob.prob[name]=Prob.prob_num[i]
    

    #making dict of prob and iso
    atom_name = [re.sub(r'\d', '', prob_name) for prob_name in Prob.prob_name] #remove the number in iso
    atom_name = sorted(list(set(atom_name))) #remove the overlap atom
    #print(atom_name)
    for temp_atom in atom_name:  #make empty list
        Prob.iso_name[temp_atom] = []
        Prob.iso_prob[temp_atom] = []

    
    #add the iso and prob
    for i, iso in enumerate(Prob.prob_name):
        temp_atom = re.sub(r'\d', '', iso)
        Prob.iso_name[temp_atom].append(iso) # add iso name
        Prob.iso_prob[temp_atom].append(Prob.prob_num[i]) # add iso prob (%)
    
    #add the fail to dict
    for temp_atom in atom_name:
        Sum_prob = sum(Prob.iso_prob[temp_atom]) 
        if (Sum_prob < 100): #add fail
            Prob.iso_name[temp_atom].append("__Nope__") # add fail name
            Prob.iso_prob[temp_atom].append(100-Sum_prob) # add fail prob
        elif Sum_prob > 100 : #err, the prob is over 100 %!!
            print("\tError, The probability of isotope is over 100%!!")
            print(Prob.iso_name)
            print(Prob.iso_prob)    
            sys.exit(1)

        Prob.iso_prob[temp_atom] = list(np.array(Prob.iso_prob[temp_atom])/100) #sum goes to 1

    return Prob

def CheckProb(ListPOSCAR, probData, args, myrank):
    AtomList=[];     Check=[];

    #Get all atom into list
    for POSCARdata in ListPOSCAR:
        for Atom in POSCARdata.Atom_data.keys():
            AtomList.append(Atom)

    #Check the Atoms with POSCAR
    AtomList=list(set(AtomList))
    for Iso in probData.prob_name:
        for Atom in AtomList:
            if IsoCheck(Iso,Atom):
                Check.append(Atom)

    #Check the Atoms with Defect
    Check.append(args.Defect)
    Check=list(set(Check))

    #Make the relative complement
    Err_Atom=list(set(AtomList)-set(Check))

    if len(Err_Atom) != 0 and myrank==0:
        print("\t---------------------------------------------------------------------------------")
        print("\t>> Error in Atom Probability!! Check the probability file!!")
        print("\t   Some Atoms does not have probability data!!",Err_Atom)
        print("\t---------------------------------------------------------------------------------\n")
        sys.exit(1)

    #print the prob information
    if myrank == 0:
        print("\t< Probability information > -------------------------------------")
        print("\t (Iso)\t (Prob)\t (Spin)\t")
        for i in range(len(probData.prob_name)):
            print("\t {}\t {}\t {}".format(probData.prob_name[i],probData.prob_num[i],probData.prob_spin[i]))
        print(flush=True)

def IsoCheck(iso,Atom):
    #Example :
    #########################
    #   iso = '10Ba'        #
    #   Atom = 'Ba'         #
    #   --> return True     #
    #########################
    #   iso = '10Ba'        #
    #   Atom = 'BA'         #
    #   --> return False    #
    #########################

    p=re.compile(Atom)
    if p.search(iso) != None:
        return True
    else:
        return False
    

##########################
#2.3. Check the lattice informatoin
##########################
def CheckLattice(ListPOSCAR, Deformation_axis,opt2D, myrank):

    Error=False

    #No consider the just one POSCAR data
    if len(ListPOSCAR) == 1:
        pass

    #Only consider the pure and defect POSCAR
    elif len(ListPOSCAR) == 2:
        #check the lattice parameter
        if (ListPOSCAR[0].Cell_parameter != ListPOSCAR[1].Cell_parameter):
            Error=True #2D or bulk material

        if (Error) and myrank==0:
            print("\t---------------------------------------------------------------------------------")
            print("\t>> Error in Cell Parameter, Need to Check it!!\n")
            print("\t   pure POSCAR Cell parameter   : ",ListPOSCAR[0].inputfile)
            print("\t   ",ListPOSCAR[0].Cell_parameter)
            print()
            print("\t   defect POSCAR Cell parameter : ",ListPOSCAR[1].inputfile)
            print("\t   ",ListPOSCAR[1].Cell_parameter)
            print("\t---------------------------------------------------------------------------------\n")
            sys.exit(1)

    #Consider the deformated POSCARs
    else:
        for num, axis in  enumerate(Deformation_axis):
            Axis=[0,1,2]
            if 'z' in axis: #c-axis
                del Axis[2]
            if 'y' in axis: #b-axis
                del Axis[1]
            if 'x' in axis: #a-axis
                del Axis[0]

            #Check the lattice parameter mismatch
            for axis_num in [0,1,2]:
                #pure == deform (except axis)
                if axis_num in Axis: #except Deformation_axis
                    if (ListPOSCAR[0].Cell_parameter[axis_num] != ListPOSCAR[num+2].Cell_parameter[axis_num]):
                        Error=True

                #defect == deform (only axis)
                else: #Deformation_axis
                    if (ListPOSCAR[1].Cell_parameter[axis_num] != ListPOSCAR[num+2].Cell_parameter[axis_num]):
                        Error=True

            if (Error) and myrank==0:
                print("\t---------------------------------------------------------------------------------")
                print("\t>> Error in Cell Parameter, Need to Check it!!\n")
                print("\t   pure POSCAR Cell parameter   : ",ListPOSCAR[0].inputfile)
                print("\t   ",ListPOSCAR[0].Cell_parameter)
                print()
                print("\t   defect POSCAR Cell parameter : ",ListPOSCAR[1].inputfile)
                print("\t   ",ListPOSCAR[1].Cell_parameter)
                print()
                print("\t   deformated POSCAR Cell parameter : ",ListPOSCAR[num+2].inputfile)
                print("\t   ",ListPOSCAR[num+2].Cell_parameter)
                print()
                print("\t   deformated Axis parameter : ",axis)
                print("\t---------------------------------------------------------------------------------\n")
                sys.exit(1)

####################################################
####################################################

##########################
#3.0. prepare the functions
##########################

#3.0.1. distance between two positions
############################
def dist(pos1,pos2):    #array type or list type
    tmp = (np.array(pos2) - np.array(pos1))
    return np.sqrt((tmp[0])**2 + (tmp[1])**2 + (tmp[2])**2)

def Length(listA):
    tmp=0
    for i in listA:
        tmp+=i**2
    return tmp**0.5
    

def AngsPosition(latticeParameter, fractional):
    Lattice=np.array(latticeParameter).T
    Fract=np.array(fractional)

    #matrixl multiplcity
    RealPosi=Lattice@Fract     

    return RealPosi


#3.0.2.judge the boundary of Bath shape
############################
def BDCut(Realpoint,vertex):
    #print("Relapoint :",Realpoint)
    #print("vertex :",vertex)

    #Check the real point is in the BD using the vector between 8-points and defect

    #judge the realpoint is in the bath cell shape
    for i, point in enumerate(vertex):
        tmp=(np.array(Realpoint) - np.array(point))
        if i == 0:   #0. [0,0,0]
            if (tmp[0]>=0 and tmp[1]>=0 and tmp[2]>=0) == False:
                return False
        elif i == 1: #1. [x,0,0]
            if (tmp[0]<=0 and tmp[1]>=0 and tmp[2]>=0) == False:
                return False
        elif i == 2: #2. [0,y,0]
            if (tmp[0]>=0 and tmp[1]<=0 and tmp[2]>=0) == False:
                return False
        elif i == 3: #3. [x,y,0]
            if (tmp[0]<=0 and tmp[1]<=0 and tmp[2]>=0) == False:
                return False
        elif i == 4: #4. [0,0,z]
            if (tmp[0]>=0 and tmp[1]>=0 and tmp[2]<=0) == False:
                return False
        elif i == 5: #5. [x,0,z]
            if (tmp[0]<=0 and tmp[1]>=0 and tmp[2]<=0) == False:
                return False
        elif i == 6: #6. [0,y,z]
            if (tmp[0]>=0 and tmp[1]<=0 and tmp[2]<=0) == False:
                return False
        elif i == 7: #7. [x,y,z]
            if (tmp[0]<=0 and tmp[1]<=0 and tmp[2]<=0) == False:
                return False
    #all cases pass (== include the real point)
    return True 

def DefineCellShape(TargetPOSCAR, args):
    ####################################
    #Shape is same to POSCAR Cell Shape
    #POSCARdata.Cell_parameter
    
    #[a1, a2, a3] --> [a1, a2, a3] * sizeX/(a1**2 + a2**2 + a3**2)**0.5
    #[b1, b2, b3] --> [b1, b2, b3] * sizeY/(b1**2 + b2**2 + c3**2)**0.5
    #[c1, c2, c3] --> [c1, c2, c3] * sizeZ/(c1**2 + c2**2 + b3**2)**0.5
    ####################################

    #Boundary Cell Shape
    if args.shape: #cubic Bath shape
        if args.opt2D == False:
            CellShape=[[args.size,0,0],[0,args.size,0],[0,0,args.size]]
        else:
            RealZ=np.array(TargetPOSCAR.Cell_parameter).T@np.array([0,0,1])
            CellShape=[[args.size,0,0],[0,args.size,0],[0,0,RealZ[2]]]

    else : #defect POSCAR Bath shape
        CellShape=[]
        for i, line in enumerate(TargetPOSCAR.Cell_parameter):
            if args.opt2D == False:
                CellShape.append(line*args.size/Length(line))
            else:
                CellShape.append(line)

    return CellShape

def Findvertex3D(CellShape):
    #CellShape --> create the 8-point vertex of 3D Bath Cell
    #
    #there are samples of vertex
    #0. [0, 0, 0]
    #1. [x, 0, 0]
    #2. [0, y, 0]
    #3. [x, y, 0]
    #4. [0, 0, z]
    #5. [x, 0, z]
    #6. [0, y, z]
    #7. [x, y, z]
   
    vertex=[]
    for z in  [0,1]:
        for y in  [0,1]:
            for x in [0,1]:
                vertex.append(np.array(CellShape).T@np.array([x,y,z]))
    return vertex


def PointCut(Realpoint, minCut, maxCut):
    if Realpoint <= maxCut and Realpoint >= minCut:
        return True
    else :
        return False


##########################
#3.1. Set a Defect position
##########################
def FindDefect(ListPOSCAR, args, myrank):
    #Choose the POSCAR to get the defect
    if len(ListPOSCAR) == 1 : #pure POSCAR
        #TargetPOSCAR=copy.deepcopy(ListPOSCAR[0]) 
        TargetPOSCAR=(ListPOSCAR[0])
    else : #defect POSCAR
        #TargetPOSCAR=copy.deepcopy(ListPOSCAR[1])
        TargetPOSCAR=(ListPOSCAR[1])
    
    #check the defect is in the POSCAR
    if args.Defect not in TargetPOSCAR.Atom_data.keys() and myrank ==0:
        print("\t---------------------------------------------------------------------------------")
        print("\t>> Error the defect (%s) is not in POSCAR(%s)!!"%(args.Defect,TargetPOSCAR.inputfile))
        print("\t   Need to Check the Defect option (--Defect %s)"%(args.Defect))
        print("\t---------------------------------------------------------------------------------\n")
        sys.exit(1)

    #check the defect center have range 
    for i in args.Defect_center:
        if ((i<0) or (i>1)) and myrank==0:
            print("\t---------------------------------------------------------------------------------")
            print("\t>> Error the Defect center should be in (0 ~ 1)!! Check the --Defect_center !!")
            print("\t   Now Defect center :",args.Defect_center)
            print("\t---------------------------------------------------------------------------------\n")
            sys.exit(1)

    
    #Choose the Defect close to center (0.5,0.5,0.5)
    Defect=None
    distance=np.inf #tmp value
    center=AngsPosition(TargetPOSCAR.Cell_parameter,args.Defect_center)
    for i, atom in enumerate(TargetPOSCAR.Atomic_info):
        if args.Defect == atom[0]:
            posi=AngsPosition(TargetPOSCAR.Cell_parameter, [atom[1],atom[2],atom[3]])
            tmpDist=dist(posi, center)
            if distance > tmpDist:
                Defect=[atom[0],atom[1],atom[2],atom[3]]
                distance=copy.deepcopy(tmpDist)

    if myrank ==0:
        print("\t< Defect information > ------------------------------------------")
        print("\t Defect Atom Name        :",Defect[0])
        print("\t Dist from Defect_center :",distance," [Angstrom] (from ",args.Defect_center,")")
        print("\t Defect fractional data  :",Defect[1], Defect[2], Defect[3], flush=True)

    #Insert the Defect information into ListPOSCAR
    TargetPOSCAR.Defect_info=Defect 
 
##########################
#3.2. Set a pivot list for configure
##########################
class PivotData:
    def __init__(self):
        self.vertex=None        #total Bath vertex (real coordi)
        self.CellShape=None     #total Bath cell shape (unit : angstrom)
        self.Center=None        #total Bath center (real coordi)

        self.CellIndex=[]       #possible indexs to generate atoms

        self.DefectPivot=None   #vertex of defect poscar  (real coordi)
        self.Defect_info=None   #defect information (fractional coordi)
    

def SetPivotList(ListPOSCAR, args, myrank):
    #just find the minimum and maximum number index of lattice

    #Choose the POSCAR at the center
    if len(ListPOSCAR) == 1 :   #pure POSCAR
        TargetPOSCAR=copy.deepcopy(ListPOSCAR[0])
    else :                      #defect POSCAR
        TargetPOSCAR=copy.deepcopy(ListPOSCAR[1])

    #Create the pivot data class
    Pivot=PivotData()
    Pivot.Defect_info=TargetPOSCAR.Defect_info

    #Define the total Bath cell shape
    Pivot.CellShape = DefineCellShape(TargetPOSCAR,args)
    #Find the vertex of total Bath cell shape
    Pivot.vertex = Findvertex3D(Pivot.CellShape)

    #Find the defect's 8-vertex in real-position
    Pivot.Center=np.array(Pivot.CellShape).T@np.array([0.5,0.5,0.5]) #total Bath's center
    Pivot.DefectPivot=FindDefectPivot(TargetPOSCAR,Pivot.Center,args)

    #Find the minimum cell parameter for using the lattice BDcut ruler
    min_length=np.inf
    Ref_Cell_parameter=[]
    for i in range(3):
        ref_cell=None
        for POSCARdata in ListPOSCAR:
            temp_length=Length(POSCARdata.Cell_parameter[i])
            if Length(POSCARdata.Cell_parameter[i]) <= temp_length:
                min_length=temp_length
                ref_cell=POSCARdata.Cell_parameter[i]
        Ref_Cell_parameter.append(ref_cell)

    #print(Pivot.__dict__)

    #set the center index into CellIndex
    XIndex=[0,0];    YIndex=[0,0];    ZIndex=[0,0];
    Pivot.CellIndex=[XIndex,YIndex,ZIndex]

    #Find the possible lattice indexes
    for i, pivot in enumerate(Pivot.DefectPivot):
        add_x=0; add_y=0; add_z=0;

        Vector=np.array(pivot)-np.array(Pivot.Center) #unit : angstrom

        signlist=[1,1,1]
        if Vector[0] < 0:  signlist[0]=-1 #sign a-vector
        if Vector[1] < 0:  signlist[1]=-1 #sign b-vector
        if Vector[2] < 0:  signlist[2]=-1 #sign c-vector

        FindBDPivot(pivot, Ref_Cell_parameter, signlist, Pivot,args)
        
    #Increase the Cell index for extra
    for i, index in enumerate(Pivot.CellIndex):
        if i==2 and args.opt2D == True:
            Pivot.CellIndex[i][1] += 1
            continue
        else:
            Pivot.CellIndex[i][0] -= 3
            Pivot.CellIndex[i][1] += 3

    if myrank == 0:
        print()
        print("\t< Bath Range candidate information > ----------------------------")
        print("\t Bath Cell Shape          : ",Pivot.CellShape)
        print("\t Range of the cell indexs : ",Pivot.CellIndex)

    return Pivot

def FindBDPivot(pivot, Ref_Cell_parameter, signlist, Pivot,args):
    
    #Make the 7-direction from defect pivot
    # Example of BD pivot
    # 1.  (signlist[0],0,0) 
    # 2.  (0,signlist[1],0)  
    # 3.  (signlist[0],signlist[1],0)
    # 4.  (0,0,signlist[2])
    # 5.  (signlist[0],0,signlist[2]) 
    # 6.  (0,signlist[1],signlist[2])
    # 7.  (signlist[0],signlist[1],signlist[2])
    tmpDirect=[]
    for zdir in [0,signlist[2]]:
        if args.opt2D == True and zdir != 0:
            continue
        for ydir in [0,signlist[1]]:
            for xdir in [0,signlist[0]]:
                tmpDirect.append([xdir,ydir,zdir])
    del tmpDirect[0]


    #Check the BD pivot using 7-direction
    for d, Direct in enumerate(tmpDirect):
        #set the temporary Ruler vector
        tmpRuler=np.array(Ref_Cell_parameter).T@np.array(Direct)

        #set the multiplicity for using Ruleer
        multi=np.inf
        IndexSign=[]
        for i, cellindex in enumerate(Direct):
            tmpmul=np.inf
            if cellindex >0:  
                IndexSign.append(1) #positive direct
                tmpmul=copy.deepcopy(np.abs(Pivot.CellIndex[i][1]))
            elif cellindex < 0:
                IndexSign.append(0) #negative direct
                tmpmul=copy.deepcopy(np.abs(Pivot.CellIndex[i][0]))
            else:
                continue
            if multi > tmpmul:
                multi=tmpmul

        if multi == 0:  multi=1

        # Find the maximum multiplicity 
        while True:
            test=(BDCut(np.array(Pivot.DefectPivot[d])+multi*tmpRuler, Pivot.vertex))
            if (test == False):
                # Renewal the Maximum or Minimum Cell Index for generating bath
                for i, cellindex in enumerate(Direct):
                    #get the reference multiplicity
                    if cellindex >0:
                        Ref_multi=copy.deepcopy(Pivot.CellIndex[i][1])
                    elif cellindex <0 :
                        Ref_multi=copy.deepcopy(Pivot.CellIndex[i][0])
                    else : #newal ==0
                        #skip the this time
                        continue 

                    #renewal the multiplcity into Pivot class
                    if (np.abs(multi) > np.abs(Ref_multi)):
                        if cellindex >0:
                            Pivot.CellIndex[i][1] = copy.deepcopy(multi)
                        elif cellindex <0:
                            Pivot.CellIndex[i][0] = -copy.deepcopy(multi)
                            
                break
            multi+=1


def FindDefectPivot(defectPOSCAR,Center, args):
    Pivot=[]
    #create the real position pivot (8-point)
    for z in [-args.Defect_center[2], 1-args.Defect_center[2]]:
        for y in [-args.Defect_center[1], 1-args.Defect_center[1]]:
            for x in [-args.Defect_center[0], 1-args.Defect_center[0]]:
                posi=AngsPosition(defectPOSCAR.Cell_parameter, [x,y,z])
                Pivot.append([Center[0]+posi[0],Center[1]+posi[1],Center[2]+posi[2]]) 

    #Pivot have 8 pivot-point; (unit : angstrom)
    #0. [center-x, center-y, center-z]
    #1. [center+x, center-y, center-z]
    #2. [center-x, center+y, center-z]
    #3. [center+x, center+y, center-z]
    #4. [center-x, center-y, center+z]
    #5. [center+x, center-y, center+z]
    #6. [center-x, center+y, center+z]
    #7. [center+x, center+y, center+z]

    return Pivot

    
##########################
#3.3. Make a configure atoms
##########################
class ConfigData:
    def __init__(self):
        self.outputfile=None    #args.output
        self.defectfile=None    #args.output + '_defect'
        self.poscarfile=None    #'POSCAR_' + args.output 

        self.bath=[]            #bath info    : {iso, x_real, y_real, z_real}
        self.defect=[]          #defect info  : {iso, x_real, y_real, z_real}

        self.tot_iso_Num={}     #tot isotopic number info (created iso and In BD)
        self.Atom_data={}       #dict of atom data ; {atom name : atom number} 
        self.TableAtom2Iso={}   #dict key : Atom & value : Iso
        self.TableIso2Atom={}   #dict key : Atom & value : Iso
        
        self.limit={}           #limit r_bath : {bottom, up, left, right, back, front)
                                #using bath vertex


def MakeConfigure(ListPOSCAR,probData,Pivot,args, myrank, world_size):
    
    #Create the Configure data class
    Config=ConfigData()
    Config.outputfile=args.output
    Config.defectfile=args.output+'_defect'
    Config.poscarfile=Add_string(args.output,'POSCAR_','/')
    
    #define pure poscar and defect poscar 
    if len(ListPOSCAR) == 1: #pure poscar
        purePOSCAR=copy.deepcopy(ListPOSCAR[0])    #pure POSCAR
        defectPOSCAR=copy.deepcopy(ListPOSCAR[0])  #defect POSCAR
        #Cell_parameter = np.array(ListPOSCAR[0].Cell_parameter)

    elif len(ListPOSCAR) == 2: #insert defect poscar
        purePOSCAR=copy.deepcopy(ListPOSCAR[0])    #pure POSCAR
        defectPOSCAR=copy.deepcopy(ListPOSCAR[1])  #defect POSCAR
        #Cell_parameter = np.array(ListPOSCAR[0].Cell_parameter)

    else : #deformated POSCAR
        purePOSCAR=copy.deepcopy(ListPOSCAR[0])    #pure POSCAR
        defectPOSCAR=copy.deepcopy(ListPOSCAR[1])  #defect POSCAR
        #Cell_parameter = np.array(ListPOSCAR[0].Cell_parameter)
    
    #set the zero point in real position 
    #                  and cell parameter
    Zeropoint=np.array(Pivot.DefectPivot[0])

    #print the atom probability using verbosity
    if myrank == 0:
        print()
        print("\t< Now Generate random isotope > ---------------------------------", flush=True)
        if args.verbosity == False:
            print("\t .....")
            print("\t  Skip printing of process of isotope detail ")
            print("\t  if you want to see detail, turn on the option : use '--verbosity'")
            print("\t .....", flush=True)

    #Generage the xyz index list, previous
    Xrange=range(Pivot.CellIndex[0][0],Pivot.CellIndex[0][1])
    Yrange=range(Pivot.CellIndex[1][0],Pivot.CellIndex[1][1])
    Zrange=range(Pivot.CellIndex[2][0],Pivot.CellIndex[2][1])
    
    xyzInd_list = list(product(Xrange, Yrange, Zrange))

    #for mpi process, distribute the index
    index_MPI = FindMPIindex(myrank, len(xyzInd_list), world_size)

    mpi_xyzInd_list = []
    for i_mpi in index_MPI:
        mpi_xyzInd_list.append(xyzInd_list[i_mpi])

    #to reduce meory
    del index_MPI
    del xyzInd_list 
    del Xrange, Yrange, Zrange
    
    #############################################################
    # Main process to generate bath
    #############################################################
    #Generage the nuclear spin along CellIndex
    #for index_eachmpi in index_MPI:
    while len(mpi_xyzInd_list) != 0:
        cell = mpi_xyzInd_list.pop(0) #get and remove

        xInd = cell[0]
        yInd = cell[1]
        zInd = cell[2]

        ####################################
        #Check and Set the defect atomic position information
        ####################################
        
        #Related defect and pure poscar
        if xInd == 0 and yInd == 0 and zInd ==0:
            Atomic_info = defectPOSCAR.Atomic_info
            Cell_parameter = np.array(defectPOSCAR.Cell_parameter)
        else :
            Atomic_info = purePOSCAR.Atomic_info
            Cell_parameter = np.array(purePOSCAR.Cell_parameter)

        #Related deformated poscar
        if len(ListPOSCAR) > 2:
            if (args.opt2D == False):
                #example of axis : ['x','y','z','xy','xz','yz']
                #about axis 'x'
                if xInd == 0 and yInd != 0 and zInd != 0 and 'x' in args.Deformation_axis:
                    index = 2 + args.Deformation_axis.index('x')
                    Atomic_info = ListPOSCAR[index].Atomic_info
                    Cell_parameter = np.array(ListPOSCAR[index].Cell_parameter)

                #about axis 'y'
                elif xInd != 0 and yInd == 0 and zInd != 0 and 'y' in args.Deformation_axis:
                    index = 2 + args.Deformation_axis.index('y')
                    Atomic_info = ListPOSCAR[index].Atomic_info
                    Cell_parameter = np.array(ListPOSCAR[index].Cell_parameter)
                    
                #about axis 'z'
                elif xInd != 0 and yInd != 0 and zInd == 0 and 'z' in args.Deformation_axis:
                    index = 2 + args.Deformation_axis.index('z')
                    Atomic_info = ListPOSCAR[index].Atomic_info
                    Cell_parameter = np.array(ListPOSCAR[index].Cell_parameter)

                #about axis 'xy'
                elif xInd == 0 and yInd == 0 and zInd != 0 and 'xy' in args.Deformation_axis:
                    index = 2 + args.Deformation_axis.index('xy')
                    Atomic_info = ListPOSCAR[index].Atomic_info
                    Cell_parameter = np.array(dwiListPOSCAR[index].Cell_parameter)

                #about axis 'yz'
                elif xInd != 0 and yInd == 0 and zInd == 0 and 'yz' in args.Deformation_axis:
                    index = 2 + args.Deformation_axis.index('yz')
                    Atomic_info = ListPOSCAR[index].Atomic_info
                    Cell_parameter = np.array(ListPOSCAR[index].Cell_parameter)
                    
                #about axis 'xz'
                elif xInd == 0 and yInd != 0 and zInd == 0 and 'xz' in args.Deformation_axis:
                    index = 2 + args.Deformation_axis.index('xz')
                    Atomic_info = ListPOSCAR[index].Atomic_info
                    Cell_parameter = np.array(ListPOSCAR[index].Cell_parameter)
            else:
                #example of axis : ['x','y','xy']
                #about axis 'x'
                if xInd == 0 and yInd != 0 and 'x' in args.Deformation_axis:
                    index = 2 + args.Deformation_axis.index('x')
                    Atomic_info = ListPOSCAR[index].Atomic_info
                    Cell_parameter = np.array(ListPOSCAR[index].Cell_parameter)

                #about axis 'y'
                elif xInd != 0 and yInd == 0 and 'y' in args.Deformation_axis:
                    index = 2 + args.Deformation_axis.index('y')
                    Atomic_info = ListPOSCAR[index].Atomic_info
                    Cell_parameter = np.array(ListPOSCAR[index].Cell_parameter)

                #about axis 'xy'
                elif xInd == 0 and yInd == 0 and 'xy' in args.Deformation_axis:
                    index = 2 + args.Deformation_axis.index('xy')
                    Atomic_info = ListPOSCAR[index].Atomic_info
                    Cell_parameter = np.array(dwiListPOSCAR[index].Cell_parameter)
        ##########################################################################

        #check the defect atom 
        Defect_index = None
        if xInd == 0 and yInd == 0 and zInd == 0 :
            for i, atomic in enumerate(Atomic_info):
                if atomic == defectPOSCAR.Defect_info :
                    Realpoint = Zeropoint + Cell_parameter.T@np.array([xInd,yInd,zInd]) \
                                + Cell_parameter.T @ np.array([atomic[1],atomic[2],atomic[3]])
                    defect_info=Realpoint.tolist()
                    defect_info.insert(0,defectPOSCAR.Defect_info[0]) #{atom, x_real, y_real, z_real}
                    Config.defect.append(defect_info)
                    Config.defect=sum(Config.defect,[])

                    #we don't need to generate the spin at defect position
                    Atomic_info.pop(i)
                    continue
        ##########################################################################

        #extract the only atom name
        Atom_name_data = list(list(zip(*Atomic_info))[0])

        #0. judge the spin is generated and isotope name
        Atom_spin_iso_info = CheckGeneratedSpin(Atom_name_data, probData) #generate isotope name
        Checking_spin_board = np.where(Atom_spin_iso_info != "__Nope__")[0] #only get the index of generated spin

        #1. check the spin is generated
        if len(Checking_spin_board) == 0: #there is no spin
            continue #skip
        #print("Checking_spin_board :", Checking_spin_board)
        
        #2. consider the spin which is generated (not '__Nope__')
        atomic = np.array(Atomic_info)
        atomic = atomic[Checking_spin_board,1:].astype(float) #only get the position information if it has spin

        #3. compute the real position of generated the spin
        Realpoints = Zeropoint + Cell_parameter.T @ np.array([xInd, yInd, zInd])
        Realpoints = np.expand_dims(Realpoints, axis=0) #expaned dimension
        Realpoints = np.tile(Realpoints, len(atomic)).reshape(len(atomic),3) #fill empty array
        
        tmp_x = np.expand_dims(atomic[:,0], axis=1)
        tmp_y = np.expand_dims(atomic[:,1], axis=1)
        tmp_z = np.expand_dims(atomic[:,2], axis=1)

        Realpoints += (Cell_parameter[0] * tmp_x + Cell_parameter[1] * tmp_y + Cell_parameter[2] * tmp_z)
        #print(Realpoints)
        

        #4. check the boundary condition
        Checking_boundary_in = np.array(CheckBoundaryCut(Realpoints, Pivot.vertex)) #if spin is boundary-in, index is called

        if len(Checking_boundary_in) == 0: #all spin is boundary out
            continue
        else: #will remove BD out index
            Checking_spin_board = Checking_spin_board[Checking_boundary_in] #spin index in BD
            Realpoints = Realpoints[Checking_boundary_in] #real position of spin in BD
            Atom_spin_iso_info = Atom_spin_iso_info[Checking_spin_board] #iso name of spin in BD
        
        #5. add the spin to list and total spin number info
        for i, tmp_iso in enumerate(Atom_spin_iso_info):
            Config.bath.append([tmp_iso, Realpoints[i][0], Realpoints[i][1], Realpoints[i][2]]) 
            GatherAtomInfo(tmp_iso, Config.tot_iso_Num)
        ##########################################################################
        
        ################################################################################
        # --> too slow
        ##make the real positioni using index and atomic information
        #for atomic in Atomic_info:
        #    Realpoint = Zeropoint + Cell_parameter.T@np.array([xInd,yInd,zInd]) \
        #                + Cell_parameter.T @ np.array([atomic[1],atomic[2],atomic[3]])
        #    
        #    #check the defect atom 
        #    if xInd == 0 and yInd == 0 and zInd == 0 and atomic == defectPOSCAR.Defect_info :
        #        defect_info=Realpoint.tolist()
        #        defect_info.insert(0,defectPOSCAR.Defect_info[0]) #{atom, x_real, y_real, z_real}
        #        Config.defect.append(defect_info)
        #        Config.defect=sum(Config.defect,[])
        #        continue
    
        #    if args.verbosity:
        #        print("\t Indexs         :",xInd,yInd,zInd)
        #        print("\t Atom Name      :",atomic[0])
        #        print("\t Atomic Posi    :",atomic[1],atomic[2],atomic[3])
        #        print("\t Real Posi      :",Realpoint[0],Realpoint[1],Realpoint[2])
        #   
        #    #atom is in-side BD
        #    if BDCut(Realpoint,Pivot.vertex) == True:
        #        GatherAtomInfo(atomic[0], Config.tot_iso_Num)

        #        #Find the nuclear spin using prob file
        #        tmp_iso = CheckSpin(atomic, probData, args)
        #        
        #        if ((tmp_iso) != str(False)): #yes, generate spin
        #            #{iso, x_real, y_real, z_real}
        #            Config.bath.append([tmp_iso, Realpoint[0], Realpoint[1], Realpoint[2]]) 
        #            GatherAtomInfo(tmp_iso, Config.tot_iso_Num)

        #            if args.verbosity == True: 
        #                print("\t O (rand, {})\n".format(tmp_iso))

        #        else: #no, no spin

        #            if args.verbosity == True:
        #                print("\t X (rand)\n")

        #    #atom is out of BD
        #    else:
        #        if args.verbosity == True:
        #            print("\t X (out of Boundary!)\n")
    #############################################################

    #Gather the config data into myrank==0
    #   after this line, we only use the myrank == 0
    #   Config.defect, Config.bath, Config.tot_iso_Num
    Config.defect = comm.gather(Config.defect, root=0) #it's dim is 2D
    Config.bath = comm.gather(Config.bath, root=0) #it's dim is 3D
    tmp_tot_iso_Num = comm.gather(Config.tot_iso_Num, root=0) #list format
    ##########################################################################


    if myrank == 0:
        Config.defect = sum(Config.defect, []) #it should be 1D
        Config.bath = sum(Config.bath, [])  #it should be 2D
        Config.bath = list(filter(None, Config.bath)) #remove the empty list in 2D list

        Config.tot_iso_Num = {} #initialize the dict
        for tmp in tmp_tot_iso_Num:
            keys = list(tmp.keys())
            for k in keys:
                GatherAtomInfo(k, Config.tot_iso_Num, tmp[k]) #re calculate the dicts
                        
        #check the defect information                 
        if len(Config.defect) == 0:
            print("\t---------------------------------------------------------------------------------")
            print("\t>> Error in defect information in Configure class")
            print("\t   Need to check the defect of configure")
            print("\t---------------------------------------------------------------------------------\n")
            sys.exit(1)

        #find the radius limitation to cal CCE
        FindLimit(Config,Pivot)

        #Get all atom into list
        atomlist=[]
        for POSCARdata in ListPOSCAR:
            for atom in POSCARdata.Atom_data.keys():
                atomlist.append(atom)
        atomlist=list(set(atomlist))

        #print the created atom probability
        PrintIsotopInfo(Config,probData,atomlist)
    ##########################################################################

    return Config
##########################################################################

def CheckGeneratedSpin(atom_name_data, probData):
    #
    # atom_name_ata : list
    # probData : class
    #

    #check the atom_name_data is in prob file
    for tmp_atom in list(set(atom_name_data)):
        if tmp_atom not in list(probData.iso_name.keys()):
            print("\tThe Atom in bath configure is not in probfile!!")
            print("\tTarget atom in Configure : %s"%(tmp_atom))
            print("\tAtom name in probfile :",probData.iso_name.keys())
            sys.exit(1)
    
    iso_list = []
    for tmp_atom in atom_name_data:
        iso_list.append(np.random.choice(probData.iso_name[tmp_atom], 1, p=probData.iso_prob[tmp_atom])[0])
    
    return np.array(iso_list)
##########################################################################

def CheckBoundaryCut(Realpoints, vertex):
    #Check the real point is in the BD using the vector between 8-points 
    # if spin is inside-in BD, the index is returned
    # return : list
    #
    # Realpoints    : np.array (# of atom spin, 3)
    # vertex        : list --> np.array
    #                 there are samples of vertex
    #                 0. [0, 0, 0]
    #                 1. [x, 0, 0]
    #                 2. [0, y, 0]
    #                 3. [x, y, 0]
    #                 4. [0, 0, z]
    #                 5. [x, 0, z]
    #                 6. [0, y, z]
    #                 7. [x, y, z]
    #

    #0. [0, 0, 0]
    condi_0 = np.where( (Realpoints[:,0] >= vertex[0][0]) \
                      & (Realpoints[:,1] >= vertex[0][1]) \
                      & (Realpoints[:,2] >= vertex[0][2]) )[0]
    if len(condi_0) == 0: return []
    tmp_index = set(condi_0)
    

    #1. [x, 0, 0]
    condi_1 = np.where( (Realpoints[:,0] <= vertex[1][0]) \
                      & (Realpoints[:,1] >= vertex[1][1]) \
                      & (Realpoints[:,2] >= vertex[1][2]) )[0]
    if len(condi_1) == 0: return condi_1
    tmp_index = tmp_index & set(condi_1)
    if len(tmp_index) == 0: return []

    #2. [0, y, 0]
    condi_2 = np.where( (Realpoints[:,0] >= vertex[2][0]) \
                      & (Realpoints[:,1] <= vertex[2][1]) \
                      & (Realpoints[:,2] >= vertex[2][2]) )[0]
    if len(condi_2) == 0: return condi_2
    tmp_index = tmp_index & set(condi_2)
    if len(tmp_index) == 0: return []

    #3. [x, y, 0]
    condi_3 = np.where( (Realpoints[:,0] <= vertex[3][0]) \
                      & (Realpoints[:,1] <= vertex[3][1]) \
                      & (Realpoints[:,2] >= vertex[3][2]) )[0]
    if len(condi_3) == 0: return condi_3
    tmp_index = tmp_index & set(condi_3)
    if len(tmp_index) == 0: return []

    #4. [0, 0, z]
    condi_4 = np.where( (Realpoints[:,0] >= vertex[4][0]) \
                      & (Realpoints[:,1] >= vertex[4][1]) \
                      & (Realpoints[:,2] <= vertex[4][2]) )[0]
    if len(condi_4) == 0: return condi_4
    tmp_index = tmp_index & set(condi_4)
    if len(tmp_index) == 0: return []

    #5. [x, 0, z]
    condi_5 = np.where( (Realpoints[:,0] <= vertex[5][0]) \
                      & (Realpoints[:,1] >= vertex[5][1]) \
                      & (Realpoints[:,2] <= vertex[5][2]) )[0]
    if len(condi_5) == 0: return condi_5
    tmp_index = tmp_index & set(condi_5)
    if len(tmp_index) == 0: return []

    #6. [0, y, z]
    condi_6 = np.where( (Realpoints[:,0] >= vertex[6][0]) \
                      & (Realpoints[:,1] <= vertex[6][1]) \
                      & (Realpoints[:,2] <= vertex[6][2]) )[0]
    if len(condi_6) == 0: return []
    tmp_index = tmp_index & set(condi_6)
    if len(tmp_index) == 0: return []

    #7. [x, y, z]
    condi_7 = np.where( (Realpoints[:,0] <= vertex[7][0]) \
                      & (Realpoints[:,1] <= vertex[7][1]) \
                      & (Realpoints[:,2] <= vertex[7][2]) )[0]
    if len(condi_7) == 0: return []
    tmp_index = tmp_index & set(condi_7)
    if len(tmp_index) == 0: return []


    #return overlap index which is in BD
    return list(tmp_index)
    

def FindMPIindex(myrank, length_data, world_size):
    #Distribute the index for myrank
    index_list =[]
    for i in range(length_data):
        j = i*world_size + myrank
        if j >= length_data:
            break
        index_list.append(j)
    return index_list

def FindLimit(Config,Pivot):

    #The vertex sample :
    #0. [0, 0, 0]    #1. [x, 0, 0]
    #2. [0, y, 0]    #3. [x, y, 0]
    #4. [0, 0, z]    #5. [x, 0, z]
    #6. [0, y, z]    #7. [x, y, z]
    
    #The plain point sample
    #[0,2,4,6] (left), [1,3,5,7] (right),
    #[2,3,6,7] (back),  [0,1,4,5] (front),
    #[0,1,2,3] (bottom), [4,5,6,7] (up)

    #limit r_bath : {bottom, up, left, right, back, front)

    #make the plain point
    plain_list={'Left   (-x)':[0,2,4,6],'Right  (+x)':[1,3,5,7],
                'Back   (-y)':[2,3,6,7],'Front  (+y)':[0,1,4,5],
                'Bottom (-z)':[0,1,2,3],'Up     (+z)':[4,5,6,7]}

    point=[Config.defect[1],Config.defect[2],Config.defect[3]]

    #add the defect point and Cell plain
    for direct, plain in plain_list.items():
        plain_points=[]
        for pp in plain:
            plain_points.append(Pivot.vertex[pp])

        #calculate the length between point and plain
        tmp=dist_p2p(point, plain_points)

        #add the length between point and plain
        Config.limit[direct] = tmp


def Add_string(output,string, splitter):
    tmp=output.split(splitter)
    tmp[-1] = string + tmp[-1]

    newString=splitter.join(tmp)

    return newString

def dist_p2p(point,plain_points):
    #array type or list type
    #point = (x1, y1, z1)

    #Make the plain eqs using plain_points
    #in this case, we just use the three point for making plain eqs

    #Ax + By+ Cz + D =0
    P1=plain_points[0]
    P2=plain_points[1]
    P3=plain_points[2]
    
    A=P1[1]*(P2[2]-P3[2]) + P2[1]*(P3[2]-P1[2]) + P3[1]*(P1[2]-P2[2])
    B=P1[2]*(P2[0]-P3[0]) + P2[2]*(P3[0]-P1[0]) + P3[2]*(P1[0]-P2[0])    
    C=P1[0]*(P2[1]-P3[1]) + P2[0]*(P3[1]-P1[1]) + P3[0]*(P1[1]-P2[1])    
    D = P1[0]*(P2[1]*P3[2]-P3[1]*P2[2]) \
        + P2[0]*(P3[1]*P1[2]-P1[1]*P3[2]) \
        + P3[0]*(P1[1]*P2[2]-P2[1]*P1[2])
    D=-D

    #Calculate the distance using formular
    l = np.abs(A*point[0] + B*point[1] + C*point[2] + D )/math.sqrt(A**2 + B**2 + C**2)

    return l



def PrintIsotopInfo(Config,probData,atomlist):

    #Check the Atoms with POSCAR
    for Iso in probData.prob_name:
        for Atom in atomlist:
            if IsoCheck(Iso,Atom):
                count=Config.Atom_data.get(Atom, 0)
                try:
                    add = Config.tot_iso_Num[Iso]
                except:
                    add = Config.tot_iso_Num[Iso]=0    

                Config.Atom_data[Atom] = count + add
                
                Config.TableIso2Atom[Iso]=Atom

                conutlist = Config.TableAtom2Iso.get(Atom,[])
                if conutlist == []:
                    Config.TableAtom2Iso[Atom]=[]
                Config.TableAtom2Iso[Atom].append(Iso)

    #print the created atom probability
    totAtom=0;    totIso=0;
    for iso, num in sorted(Config.tot_iso_Num.items()):
        #check the isotope in bath (generated)
        if iso in list(probData.prob.keys()):
            totIso+=num
        #check the atom in cell (in BD)
        if iso in atomlist:
            totAtom+=num
        

    #print the isotope information
    print()
    print("\t< Created Isotope information > ---------------------------------")
    print("\t Total Atom # in Boundary ({}) (except for Defect)".format(totAtom))
    #example    'B : 10'
    for iso, num in sorted(Config.tot_iso_Num.items()):
        if iso in atomlist:
            print("\t  {}\t: {}".format(iso, num)) 
    print()
    
    print("\t Spin percentage[%](#) in configure (except for Defect)")
    #example    '10B : 8'
    #           '11B : 2'
    for iso, num in sorted(Config.tot_iso_Num.items()):
        prob_of_iso=0
        if iso in list(probData.prob.keys()):
            for Atom, tot_num in sorted(Config.Atom_data.items()):
                print(f"iso : {iso}, Atom : {Atom}, num : {num} ,tot_num : {tot_num}")
                if IsoCheck(iso,Atom):
                    if tot_num != 0:
                        prob_of_iso = num / tot_num * 100
                    elif tot_num == 0 and num ==0:
                        prob_of_iso = 0
                    else:
                        sys.exit("Error!! tot_num is zero")
            print("\t  {}  : {:>8.4f}% ({})".format(iso, prob_of_iso, Config.tot_iso_Num[iso]))
    print("\t --> All created Spins (#) : {}".format(totIso))
    print()

def GatherAtomInfo(atom, tot_iso_Num, add_value=None):
    count=tot_iso_Num.get(atom, 0)
    if add_value == None: add_value = 1
    tot_iso_Num[atom] = count + add_value

def CheckSpin(atomic, probData, args):

    #gather all probability
    problist=[]; probnamelist=[];

    #check the prob
    for iso in list(probData.prob.keys()):
        if (IsoCheck(iso,atomic[0]) == True):
            problist.append(probData.prob[iso])
            probnamelist.append(iso)
    
    #error of no match the iso and atom
    if len(problist)==0:
        print("\t---------------------------------------------------------------------------------")
        print("\t>> Error in matching the prob of iso and atom in making bath!")
        print("\t   Need to check the Isotopic data and Atom ({} & {})".format(list(probData.prob.keys()),atomic[0]))
        print("\t---------------------------------------------------------------------------------\n")
        sys.exit(1)
    
    #error in prob number (unit : %)
    if sum(problist) > 100:
        print("\t---------------------------------------------------------------------------------")
        print("\t>> Error in probabiltiy number in making bath!")
        print("\t   Summation of probability ({}) is over 100%!".format(sum(problist)))
        print("\t---------------------------------------------------------------------------------\n")
        sys.exit(1)
    
    #making the nuclear spin using ramdom number
    problist.append(100-sum(problist))
    probnamelist.append(False)
    
    if args.verbosity:
        print("\t Iso Name       :",probnamelist)
        print("\t Prob List      :",problist)
    
    problist=np.array(problist)/100
    
    return str(np.random.choice(probnamelist ,1,p=problist)[0])


####################################################
####################################################

##########################
#4. Write the configure
##########################
def WriteConfigure(Config, args):
    print("\t< Now Save the Bath information > -------------------------------")
    print("\t Save the Bath info into output file   :",Config.outputfile)

    with open(Config.outputfile,'w') as f:
        f.write("{:>10.5f}\t{:>10.5f}\t{:>10.5f}\t{:>10.5f}\n".format(len(Config.bath)+1, 0,0,0))
        for bath in Config.bath:
            f.write("{:>10.5f}\t{:>10.5f}\t{:>10.5f}\t{:>10}\n".format(bath[1],bath[2],bath[3],bath[0]))
    print("\t Total Bath info : {}".format(len(Config.bath)))
    print()

def WriteDefect(Config, args):
    print("\t< Now Save the Defect information > -----------------------------")
    print("\t Save the Defect info into defect file     :",Config.defectfile)

    with open(Config.defectfile,'w') as f:
        #write the Defect info
        f.write("{:>10.5f}\t{:>10.5f}\t{:>10.5f}\n".format(Config.defect[1],Config.defect[2],Config.defect[3]))
        f.write("\t(Defect atom name : {})\n".format(Config.defect[0]))
        print("\n\t Defect atom name : {}".format(Config.defect[0]))
        print("\t {:>10.5f}\t{:>10.5f}\t{:>10.5f}".format(Config.defect[1],Config.defect[2],Config.defect[3]))
        
        #write the rBath limit
        rBathCut=np.inf
        f.write("\n\tLength between Bath cell and defect point [Angstrom]\n")
        print("\n\t Length between Bath cell and defect point [Angstrom]")
        for direct, leng_limt in Config.limit.items():
            f.write("\t {:11}-plain \t : {:>10.5f}\n".format(direct, leng_limt))
            print("\t  {:11}-plain \t : {:>10.5f}".format(direct, leng_limt))
            if rBathCut > leng_limt:
                if (args.opt2D) and (direct == 'Bottom (-z)' or direct == 'Up     (+z)') :
                        continue
                rBathCut = leng_limt

        f.write("\n")
        print()
        if (args.opt2D):
            f.write("\tCaution of using Bath size in 2D!\n")
            print("\t Caution of using Bath size in 2D!")
        f.write(">> Do not jump the rBath limitation!! (limit :{:>6.5})\n".format(rBathCut))
        print("\t Do not jump the rBath limitation!! (limit :{:>6.5})".format(rBathCut))

        f.write("\n\tMade time : "+time.strftime('%c',time.localtime(time.time())))
    print()


def WritePOSCAR(Config,Pivot, args):

    #prepare for lattice parameter
    lattice_para=Pivot.CellShape
    rev_lattice=np.linalg.inv(np.array(lattice_para))

    #prepare for poscar atomic position
    Atom_list=[];    Bath_pos=[];
    for bath in Config.bath:
        #bath[0] #iso
        Atom=Config.TableIso2Atom[bath[0]]
        Atom_list.append(Atom)

        #bath[1], bath[2], bath[3] #real position
        tmp=np.array(rev_lattice)@np.array([bath[1],bath[2],bath[3]])
        tmp=tmp.tolist()
        #Bath_pos.append(sum([[Atom],[tmp[0],tmp[1],tmp[2]],[]))
        Bath_pos.append([Atom,tmp[0],tmp[1],tmp[2]])
            

    #sort the Bath along atom name and atom name
    Bath_pos.sort(key=lambda x:x[0])
    Atom_list=list(set(Atom_list))
    Atom_list.sort()


    #write the poscar file using bath
    print("\t< Now Save the POSCAR with Bath information > -------------------")
    print("\t Save the  info into POSCAR     :",Config.poscarfile)

    with open(Config.poscarfile, 'w') as f:
        f.write("POSCAR_of_{}\n".format(args.output))
        f.write("1.0\n")

        #lattice parameter
        f.write(" {:>10.5f}\t {:>10.5f}\t {:>10.5f}\n".format(lattice_para[0][0],lattice_para[0][1],lattice_para[0][2]))
        f.write(" {:>10.5f}\t {:>10.5f}\t {:>10.5f}\n".format(lattice_para[1][0],lattice_para[1][1],lattice_para[1][2]))
        f.write(" {:>10.5f}\t {:>10.5f}\t {:>10.5f}\n".format(lattice_para[2][0],lattice_para[2][1],lattice_para[2][2]))
        
        #For atom name and number
        tmp_atom=[];        tmp_num=[];
        for Atom in Atom_list:
            tmp_atom.append(Atom)
            tmp_num.append(Config.Atom_data[Atom])

        #Atom name
        for atom in tmp_atom:
            f.write(" {}  ".format(atom))
        f.write("\n")
        #Atom number
        for num in tmp_num:
            f.write(" {}  ".format(num))
        f.write("\n")
        
        f.write("Direct\n")

        ##atomic position
        for atomic in Bath_pos:
            f.write("{:>12.9f}\t{:>12.9f}\t{:>12.9f}\n".format(atomic[1],atomic[2],atomic[3]))


####################################################
####################################################

##########################
#5. Convert the configure to poscar
##########################
def ReadConfig(inputConfig,CheckMinBond, myrank):
    #check the inputfile 
    if not (os.path.isfile(inputConfig)) :
        print("\t---------------------------------------------------------------------------------")
        print("\tThere is no file!! : "+str(inputConfig))
        print("\tCheck the Option : --convert")
        print("\t---------------------------------------------------------------------------------\n")
        sys.exit(1)
    
    #create the configure class
    Configdata=ConfigData()

    #read the bath info
    with open(inputConfig, 'r') as f:
        for i, line in enumerate(f):
            inputdata=line.strip().split()
            if i == 0: #line
                temp_line=list(map(float,inputdata))[0]
                temp_line=int(temp_line)-1
                continue
            else: #bath info
                temp_bath=list(map(float,inputdata[0:3]))
                temp_bath.insert(0,inputdata[3])
                Configdata.bath.append(temp_bath)
    f.close()

    #check the error in number of bath info
    if temp_line != len(Configdata.bath) and myrank==0:
        print("\t---------------------------------------------------------------------------------")
        print("\t>> Error of number of bath info while converting configure to poscar!")
        print("\t   Length of Configure : {}, total line : {}".format(len(Configdata.bat),temp_line))
        print("\t   Check the Configure file ({})".format(inputConfig))
        print("\t---------------------------------------------------------------------------------\n")
        sys.exit(1)

    #check the minimum length from configure
    if CheckMinBond == True:
        MinBonding=np.inf
        for i in range(len(Configdata.bath)):
            for j in range(i+1,len(Configdata.bath)):
                tmp = dist(np.array(Configdata.bath[i][1:]), np.array(Configdata.bath[j][1:]))

                if tmp < MinBonding:
                    MinBonding = tmp
                    index=[i,j]

                    if MinBonding < 1 and myrank ==0:
                        print("\t>> Becareful! Some bonding is lower than 1A!")
                        print("\t   Minimum length of Bonding :",MinBonding)
                        print("\t   Atom positions related with minimum are,")
                        print("\t   Atom [{}] :\t".format(index[0]), Configdata.bath[index[0]])
                        print("\t   Atom [{}] :\t".format(index[1]), Configdata.bath[index[1]])
                        print()
                    
        if myrank == 0:
            print("\t>> Minimum length of Bonding :",MinBonding)
            print("\t   Atom positions related with minimum are,")
            print("\t   Atom [{}] :\t".format(index[0]), Configdata.bath[index[0]])
            print("\t   Atom [{}] :\t".format(index[1]), Configdata.bath[index[1]])
            print()
    
    
    #Gather isotope information
    for bath_info in Configdata.bath:
        GatherAtomInfo(bath_info[0],Configdata.tot_iso_Num)

    #convert the isotope info to Atom info
    for iso, tot_num in Configdata.tot_iso_Num.items():
        atom="".join(re.findall("[a-zA-Z]+",iso))
        count=Configdata.Atom_data.get(atom, 0)
        Configdata.Atom_data[atom] = count + tot_num

        Configdata.TableIso2Atom[iso]=atom

    return Configdata

def ConvertConfig(outputPoscar,Configdata):
    #convert the info to POSCAR info
    limX=[0,0];    limY=[0,0];    limZ=[0,0];
    for bath_info in Configdata.bath:
        #lattice parameter (cubic shape)
        #X-axis
        if bath_info[1] > limX[1]:      limX[1] = bath_info[1]
        elif bath_info[1] < limX[0]:    limX[0] = bath_info[1]
        #Y-axis
        if bath_info[2] > limY[1]:      limY[1] = bath_info[2]
        elif bath_info[2] < limY[0]:    limY[0] = bath_info[2]
        #Z-axis
        if bath_info[3] > limZ[1]:      limZ[1] = bath_info[3]
        elif bath_info[3] < limZ[0]:    limZ[0] = bath_info[3]

    #prepare for lattice parameter
    lattice_para=[[limX[1]-limX[0],0,0],[0,limY[1]-limY[0],0],[0,0,limZ[1]-limZ[0]]]
    lattice_para=np.array(lattice_para)
    rev_lattice=np.linalg.inv(np.array(lattice_para))

    #prepare for poscar atomic position
    Atom_list=[];    Bath_pos=[];
    for bath in Config.bath:
        #bath[0] #iso
        Atom=Config.TableIso2Atom[bath[0]]
        Atom_list.append(Atom)

        #bath[1], bath[2], bath[3] #real position
        tmp=np.array(rev_lattice)@np.array([bath[1],bath[2],bath[3]])
        tmp=tmp.tolist()
        #Bath_pos.append(sum([[Atom],[tmp[0],tmp[1],tmp[2]],[]))
        Bath_pos.append([Atom,tmp[0],tmp[1],tmp[2]])
            

    #sort the Bath along atom name and atom name
    Bath_pos.sort(key=lambda x:x[0])
    Atom_list=list(set(Atom_list))
    Atom_list.sort()


    #write the poscar file using bath
    print("\t< Now Save the POSCAR with Bath information > -------------------")
    print("\t Save the  info into POSCAR     :",outputPoscar)

    with open(outputPoscar, 'w') as f:
        f.write("POSCAR_of_{}\n".format(args.input))
        f.write("1.0\n")

        #lattice parameter
        f.write(" {:>10.5f}\t {:>10.5f}\t {:>10.5f}\n".format(lattice_para[0][0],lattice_para[0][1],lattice_para[0][2]))
        f.write(" {:>10.5f}\t {:>10.5f}\t {:>10.5f}\n".format(lattice_para[1][0],lattice_para[1][1],lattice_para[1][2]))
        f.write(" {:>10.5f}\t {:>10.5f}\t {:>10.5f}\n".format(lattice_para[2][0],lattice_para[2][1],lattice_para[2][2]))
        
        #For atom name and number
        tmp_atom=[];        tmp_num=[];
        for Atom in Atom_list:
            tmp_atom.append(Atom)
            tmp_num.append(Config.Atom_data[Atom])

        #Atom name
        for atom in tmp_atom:
            f.write(" {}  ".format(atom))
        f.write("\n")
        #Atom number
        for num in tmp_num:
            f.write(" {}  ".format(num))
        f.write("\n")
        
        f.write("Direct\n")

        ##atomic position
        for atomic in Bath_pos:
            f.write("{:>12.9f}\t{:>12.9f}\t{:>12.9f}\n".format(atomic[1],atomic[2],atomic[3]))

    ###############################


####################################################
####################################################

##########################
#6. Print the End-Phase
##########################
def PrintFinal(args):
    print()
    print("\t       !! Success creating the configure !!")
    print("\t    Now Making time : "+time.strftime('%c',time.localtime(time.time())))
    print()
    print("\t___________________________________________________\n")


####################################################
####################################################


####################################################
####################################################
#main process
####################################################
####################################################
if __name__ == "__main__":
    
    #0. Prepare the MPI option (a caller running this file with runpy may pass its own communicator as bathgen_comm,
    #   e.g. MPI.COMM_SELF to make one bath per process; example/tutorial/make_natab_baths.py does)
    comm = globals().get("bathgen_comm", MPI.COMM_WORLD)
    myrank = comm.Get_rank()
    world_size = comm.Get_size()
    
    #1. read arg
    args = parse()
    if myrank == 0:
        print()
        print("\t___________________________________________________\n")
        print("\t   Now create the Bath configure using POSCAR!!\n")
        print("\t< Your input options >  -----------------------------------------")
        print("\t",args.__dict__)
        print(flush=True)
        print("\tYou will use the MPI process (cpu = %d)\n"%(world_size))
        print(flush=True)

    #####################################
    #2.0. convert the configure to poscar
    if args.convert:
        Config=ReadConfig(args.input, args.minbond, myrank) #read the configure file
        ConvertConfig(args.output,Config) #convert to poscar
        if myrank == 0: PrintFinal(args)
        sys.exit() #done
    
    #list of poscar data
    ListPOSCAR=[]
    #2.1. read the pure POSCAR (pureData)
    ListPOSCAR.append(ReadPOSCAR(args.input,"-i {pure POSCAR}"))

    #2.2. read the defect POSCAR (defectData)
    if (args.Insert != None):
        ListPOSCAR.append(ReadPOSCAR(args.Insert,"-I {defect POSCAR}"))

    #2.3. read the deformated POSCAR (deformData)
    if (args.Deformation != None):
        if (args.Deformation_axis==None and myrank == 0):
            print("\t---------------------------------------------------------------------------------")
            print("\t>> Error of Deformation Axis information!!")
            print("\t   Check the Deformation Axis!!")
            print("\t---------------------------------------------------------------------------------\n")
            sys.exit(1) 
        if (len(args.Deformation_axis) != len(list(set(args.Deformation_axis))) and myrank == 0):
            print("\t---------------------------------------------------------------------------------")
            print("\t>> Error in Deformation axis!! It is Overlapped!!")
            print("\t   Now Deformatio axis : ", args.Deformation_axis)
            print("\t---------------------------------------------------------------------------------\n")
            sys.exit(1)

        if (args.Insert != None):
            #check the defomated axis
            if len(args.Deformation) != len(args.Deformation_axis) and myrank ==0:
                print("\t---------------------------------------------------------------------------------")
                print("\t>> Error of length of deformation POSCAR information!!")
                print("\t   Check the Deformation-POSCAR file and Deformation axis!!")
                print("\t   len(Deformation POSCAR) : ",len(args.Deformation))
                print("\t   len(Deformation Axis)   : ",len(args.Deformation_axis))
                print("\t---------------------------------------------------------------------------------\n")
                sys.exit(1)

            #read the deformated POSCAR
            for deformPOSCAR in args.Deformation:
                ListPOSCAR.append(ReadPOSCAR(deformPOSCAR,"--Deformation {deformated POSCAR}"))

        else:
            if myrank == 0:
                #print err (need defect cell)
                print("\t  To use the deformated POSCAR, we need the defect POSCAR!!")
                print("\t  Use -I {defect POSCAR}")
                print()
                sys.exit(1) 

    #2.4. check the lattice parameter for making configure
    CheckLattice(ListPOSCAR, args.Deformation_axis, args.opt2D, myrank)

    #2.5. read the probability
    probData = ReadProb(args.prob,"-p {prob file}")
    CheckProb(ListPOSCAR,probData, args, myrank)

    #####################################
    #3.1. Set a Defect position
    FindDefect(ListPOSCAR, args, myrank)
 
    #3.2. Set a pivot list for configure
    Pivot=SetPivotList(ListPOSCAR, args, myrank)
    
    #3.3. Make a configure atoms
    comm.Barrier()
    Configure = MakeConfigure(ListPOSCAR, probData, Pivot, args, myrank, world_size)
    comm.Barrier()

    #####################################
    #4. Write the configure
    if myrank == 0:
        WriteConfigure(Configure, args)  #make the configure
        WriteDefect(Configure, args) #make the defect

        if (args.test):
            WritePOSCAR(Configure,Pivot, args) #make the poscarfile

    #####################################
    #Final print code
    if myrank == 0: PrintFinal(args) #print the End-Phase

####################################################
####################################################
####################################################
####################################################
