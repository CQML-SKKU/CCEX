import numpy as np
import sys
import os
from spHamilReader.spin_database.constants import *
from spHamilReader.base.cnvt_coord import frac2cart
from spHamilReader.base.Calgeometry import distance
from spHamilReader.vasp.extract_spinHamils_vasp import printform, extSpinHamilTensors
from spHamilReader.base.ioPoscar import read_POSCAR
from spHamilReader.base.dictproc import print_dict
from spHamilReader.base.rwfiles import provenanceline, check_spin
from spHamilReader.base import pyform

###############################################################
#
#        Find vertex
#
###############################################################
def findvertex3D(CellShape,defect):
    Def = np.array(CellShape).T@np.array(defect)
    
    vertex=[]
    for z in [0,1]:
        for y in [0,1]:
            for x in [0,1]:
                temp = np.array(CellShape).T@np.array([x,y,z])
                temp -= Def
                vertex.append(temp)

    return vertex

###############################################################
#
#        FILE FORMAT
#
###############################################################

##########################
# comm
##########################

def divline():
    return  "{:->170}\n".format("")

##########################
# A tensor file
##########################

def Afilehead():
    return   "{:>10}\t{:>10}\t{:>10}\t{:>10}\t\
{:>10}\t{:>10}\t{:>10}\t{:>10}\t\
{:>10}\t{:>10}\t{:>10}\t\
{:>10}\t{:>10}\t{:>10}\n"\
.format("AtomName","r_1i_X","r_1i_Y","r_1i_Z",\
        "iso[MHz]", "Axx", "Axy", "Axz",\
                    "Ayx", "Ayy", "Ayz",\
                      "Azx", "Azy", "Azz")                 

def Afileinfo():
    atm = "{:5}"
    gfactorline = "g-factor___: {:>10f} [No-unit]\n"
    etcline     = "etc____:\t{:10}\t{:10}\t{:10}\t {:10.5f}\n"\
                              .format("---","---","---",0.00)
    minline = "MinDif[A]:\t {:10.5f}\t {:10.5f}\t {:10.5f}\n"
    maxline = "MaxDif[A]:\t {:10.5f}\t {:10.5f}\t {:10.5f}\n"
    vertex= "\
v1       :\t {:10.5f}\t {:10.5f}\t {:10.5f}\n\
v2       :\t {:10.5f}\t {:10.5f}\t {:10.5f}\n\
v3       :\t {:10.5f}\t {:10.5f}\t {:10.5f}\n\
v4       :\t {:10.5f}\t {:10.5f}\t {:10.5f}\n\
v5       :\t {:10.5f}\t {:10.5f}\t {:10.5f}\n\
v6       :\t {:10.5f}\t {:10.5f}\t {:10.5f}\n\
v7       :\t {:10.5f}\t {:10.5f}\t {:10.5f}\n\
v8       :\t {:10.5f}\t {:10.5f}\t {:10.5f}\n"


    Dict = {'gline'    : atm+gfactorline,
            'etcline'  : atm+etcline,
            'minline'  : minline,
            'maxline'  : maxline,
            'vertex'   : vertex}

    return Dict

def Atensorform():
    return "{:10}\t {:10.5f}\t {:10.5f}\t {:10.5f}\t\
 {:10.5f}\t {:10.5f}\t {:10.5f}\t {:10.5f}\t\
 {:10.5f}\t {:10.5f}\t {:10.5f}\t\
 {:10.5f}\t {:10.5f}\t {:10.5f}\n"


##########################
# Q tensor file
##########################

def Qfilehead():
    return   "{:>10}\t{:>10}\t{:>10}\t{:>10}\t\
{:>10}\t{:>10}\t{:>10}\t\
{:>10}\t{:>10}\t{:>10}\t\
{:>10}\t{:>10}\t{:>10}\t\
[Hatree/Bohr_radius^2]\n"\
.format("AtomName","r_1i_X","r_1i_Y","r_1i_Z"\
                      ,"EFGxx","EFGxy","EFGxz"\
                      ,"EFGyx","EFGyy","EFGyz"\
                      ,"EFGzx","EFGzy","EFGzz")    

def Qfileinfo():
    atm = "{:5}"
    eQline         = "eQ___:\t {:>10f}\t[Q/milibarn]*(10^-1)\n"
    etcline     = "etc____:\t{:10}\t{:10}\t{:10}\t"\
                              .format("---","---","---")
    etclinetail= " {:10.5f}\t {:10.5f}\t {:10.5f}\t\
 {:10.5f}\t {:10.5f}\t {:10.5f}\t\
 {:10.5f}\t {:10.5f}\t {:10.5f}\n"
    minline = "MinDif[A]:\t {:10.5f}\t {:10.5f}\t {:10.5f}\n"
    maxline = "MaxDif[A]:\t {:10.5f}\t {:10.5f}\t {:10.5f}\n"
    vertex= "\
v1       :\t {:10.5f}\t {:10.5f}\t {:10.5f}\n\
v2       :\t {:10.5f}\t {:10.5f}\t {:10.5f}\n\
v3       :\t {:10.5f}\t {:10.5f}\t {:10.5f}\n\
v4       :\t {:10.5f}\t {:10.5f}\t {:10.5f}\n\
v5       :\t {:10.5f}\t {:10.5f}\t {:10.5f}\n\
v6       :\t {:10.5f}\t {:10.5f}\t {:10.5f}\n\
v7       :\t {:10.5f}\t {:10.5f}\t {:10.5f}\n\
v8       :\t {:10.5f}\t {:10.5f}\t {:10.5f}\n"


    Dict = {'eQline'    : atm+eQline,
            'etcline'   : atm+etcline+etclinetail,
            'minline'   : minline,
            'maxline'   : maxline,
            'vertex'   : vertex}

    return Dict

def Qtensorform():
    return "{:10}\t {:10.5f}\t {:10.5f}\t {:10.5f}\t\
 {:10.5f}\t {:10.5f}\t {:10.5f}\t\
 {:10.5f}\t {:10.5f}\t {:10.5f}\t\
 {:10.5f}\t {:10.5f}\t {:10.5f}\n"

def addelement(_list,elem):
    if elem not in _list:
        _list.append(elem)
    return _list

def addlist(_list,sublist):
    if sublist not in _list:
        _list.append(sublist)
    return _list

###############################################################
#
#        Check poscar if the file include central defect
#
###############################################################
def doesItHaveDefect(dposcar,dfctname):
    if dfctname in dposcar["atom_data"]:
        return  dposcar["atom_data"].index(dfctname)
    else:
        return None

def mkDefectKeyInPOSCAR(dposcar,dfctname,defectposition_frac=None):

    """ To simplify, make a key for only central defect in poscar data"""

    idx = doesItHaveDefect(dposcar,dfctname)

    #print('idx: ',idx)
    #print('dfctname :',dfctname)

    if idx != None: 
        # remove defect atomic position in unit cell information
        defectindex = sum(dposcar['atom_data_num'][:idx])
        defectposition = dposcar['unit_cell_information'][defectindex]
        del dposcar['unit_cell_information'][defectindex]
        del dposcar['atom_data_num'][idx]
        del dposcar['atom_data'][idx]

        # add defect atomic position as a new key
        dposcar['central_defect'] = defectposition    
    else:
        if defectposition_frac==None:
            sys.exit("\t"+"="*100 +"\n\n\tError!!! \
you should add defect position in argument of code \
or add defect position in file"+"\n\n\t"+"="*100)
        else:
            dposcar['central_defect'] = [dfctname] + defectposition_frac

    return dposcar    

###############################################################
#
#         Relative Distance 
#
###############################################################

def relativeDistance(latt_constant,unit_cell_parameter,frac_unitCellInfo,frac_centralDefect,_rangeformat='vertex'):
    
    relativepos = { 'r_1i': [],
                    'MinDif[A]': [], 
                    'MaxDif[A]': [],
                    'farthest_atoms_idx' : [],
                    'farthest_atoms_info' : [],
                    'vertex' : []} 

    cellparam = np.asarray(unit_cell_parameter)*latt_constant
    cart_centralDefect = frac2cart(cellparam,\
                            np.array(frac_centralDefect[1:4]))

    r_1i_dist = []

    atm_name=frac_unitCellInfo[0][0]

    for i,frac_atm in enumerate(frac_unitCellInfo):

        cart_atm = frac2cart(cellparam,\
                        np.array(frac_atm[1:4]))

        # Find farthest atom to etc line ##########################
        r_1i_dist.append([i,distance(cart_atm,cart_centralDefect)])
        if atm_name != frac_atm[0] or i == len(frac_unitCellInfo)-1:
            
            if i != len(frac_centralDefect)-1:
                del r_1i_dist[-1]

            zip_r_1i_dist = list(zip(*r_1i_dist))
            idx = zip_r_1i_dist[1].index(max(zip_r_1i_dist[1]))
            farthest_atom_idx  = r_1i_dist[idx][0]
            farthest_atom_info = relativepos['r_1i'][farthest_atom_idx]
            relativepos['farthest_atoms_idx'].append(farthest_atom_idx)
            relativepos['farthest_atoms_info'].append(farthest_atom_info)
            r_1i_dist = []
            r_1i_dist.append([i,distance(cart_atm,cart_centralDefect)])

        atm_name= frac_atm[0]
        ###########################################################

        r_1i = cart_atm - cart_centralDefect

        relativepos['r_1i'].append( [frac_atm[0]] + list(r_1i) )

        if frac_atm[0] == frac_centralDefect[0]:# atomic name..
            sys.exit("Error!!! \
                        There is central defect in unit cell information")

    ri_zip = list(zip(*relativepos['r_1i']))

    relativepos['MinDif[A]'] = [min(ri_zip[1]), min(ri_zip[2]), min(ri_zip[3])]
    relativepos['MaxDif[A]'] = [max(ri_zip[1]), max(ri_zip[2]), max(ri_zip[3])]

    if _rangeformat == 'vertex': 
        relativepos['vertex'] = findvertex3D(cellparam,frac_centralDefect[1:])

    return relativepos 

###############################################################
#
#        MAKE A TENSOR FILE
#
###############################################################

def Afile(relativepos,doutcar_Adata,doutcar_isodata,fout,_spin,_rangeformat='vertex'):

    """
        Args:

            relativepos(dict)
             - Keys 
                 'r_1i'(2d list)     : relative atomic position 
                                     between atom and central defect 
                                     e.g. [Xi,Yi,Zi] - [X_df,Y_df,Z_df]
                'MinDif[A]'(list)    : minimum r_1i
                'MaxDif[A]'(list)    : maximum r_1i
                'vertex' (list)      : vertex

            doutcar_Adata(dict) 
             - Keys 
                'Aiso'            :[]        # ion# x Fermi contact(MHz) 
                'Adip'            :[]        # ion# x [Axx, Axy, Axz, Ayx .., Azz] (MHz)
                'Adiag'           :[]        # ion# x [Axx(D), Ayy(D), Azz(D)] (MHz) for total A
                'givenr'          :{}        # ion# x given gamma (MHz/T) in incar 
                'cAiso'            :[]        # ion# x corrected Fermi contact(MHz) 
                'cAdip'            :[]        # ion# x corrected Adip(9) (MHz)
                'cAdiag'           :[]        # ion# x corrected Adiag(9) (MHz)
                                              --> cAdiag : diagonalize cAiso + cAdip

            doutcar_isodata(dict)
             - Keys
                key(str) : atom name (ex. B)
                value(myspindict) : isotope information for each atom (ex. Isodata['B'] =  {'10B': (.. , .., ) })
                ==> ex. : Isodata['B'] =  {'10B': (3.0, 2.8747, 0.08450, 0.1990) }
                isodata['B']['10B'].s          : spin         (no unit)
                isodata['B']['11B'].gyro      : gyro..     (rad/ms/G)
                isodata['N']['14N'].q          : eQ         (barn = 10^-28 m^2)
                isodata['N']['15N'].conc      : concent.. (%/100)

            fout(str)
    """

    dgfactor=[]
    dmindiff=[]
    dmaxdiff=[]
    dreldist=[]
    dAtensor=[]
    
    # g factor 
    # NOTE     : nuclear magnetic moment(\mu) = g * nuclear magneton(\mu_N) * I / h_bar
    #        : g = ( h_bar * \mu ) / ( I * \mu_N )
    #        : gyro magnetic ratio(\gamma_n) = g * nuclear magneton(\mu_N) / h_bar
    #        : g = ( h_bar * \gamma_nn ) / \mu_N

    dgfactor = []
    for atm, isotmp in doutcar_isodata.items():
        for iso, dat in isotmp.items(): 
            g = dat.gyro * 1e7 * HBAR_SI / NUCLEAR_MAGNETON
            dgfactor = addlist(dgfactor,[iso,g])

    #print(dgfactor)
    Dict = Afileinfo()
    # Make A tensor file
    with open(fout, 'w') as f:
        # HEAD: the tag and the provenance line CCEX reads (format v3, 2026-10-13).
        # VASP divides the hyperfine by <S_z> of the total spin already, so the values are those of the spin-S defect.
        f.write("VASP\n")
        f.write(provenanceline("VASP OUTCAR",_spin))
        f.write(divline())
        f.write(Afilehead())
        # ISOTOPE DATA
        f.write(divline())
        for iso_g in dgfactor:
            f.write(Dict['gline'].format(*iso_g))
        for iso_g in dgfactor:
            f.write(Dict['etcline'].format(iso_g[0]))
        # CELL RANGE 
        if _rangeformat=='coord':
            f.write(Dict['minline'].format(*relativepos['MinDif[A]']))
            f.write(Dict['maxline'].format(*relativepos['MaxDif[A]']))
        elif _rangeformat=='vertex':
            f.write(Dict['vertex'].format(*relativepos['vertex'][0],\
                                          *relativepos['vertex'][1],\
                                          *relativepos['vertex'][2],\
                                          *relativepos['vertex'][3],\
                                          *relativepos['vertex'][4],\
                                          *relativepos['vertex'][5],\
                                          *relativepos['vertex'][6],\
                                          *relativepos['vertex'][7]))
        # MAIN DATA
        doutcar_Adata['cAiso']=[]
        doutcar_Adata['cAdip']=[]
        doutcar_Adata['cAdiag']=[]
        f.write(divline())
        for i,atm in enumerate(relativepos['r_1i']):
            # NOTE : We extract the A tensor into MHz unit.
            #         We should multipy them by γ_N/g_N (If inputted γ_N=1.0)
            #         in which the value would be "μ_N/h" = "7.62259 MHz/T"
            #         In Summary:
            #            1) γ_N   = 1.0 case is considred in this code
            #            input γ_N  =1.0 , then New_Aiso = Aiso * 7.62259
            #            2) γ_N  != 1.0 case is NOT considred in this code
            #            input γ_N !=1.0 , then New_Aiso = Aiso * 7.62259 / γ_N 
            corrected_Aiso = np.asarray(doutcar_Adata['Aiso'][i])/doutcar_Adata['givenr'][atm[0]] * NUCLEAR_MAGNETON_DIV_PLANKCONST
            corrected_Adip = np.asarray(doutcar_Adata['Adip'][i])/doutcar_Adata['givenr'][atm[0]] * NUCLEAR_MAGNETON_DIV_PLANKCONST
            f.write(Atensorform().format(*atm,corrected_Aiso,*corrected_Adip))

            doutcar_Adata['cAiso'].append(corrected_Aiso)
            doutcar_Adata['cAdip'].append(corrected_Adip)
            w,v = np.linalg.eig(corrected_Aiso * np.eye(3) + corrected_Adip.reshape(3,3))
            doutcar_Adata['cAdiag'].append(w)

    return doutcar_Adata 

def _Adata(relativepos,doutcar_Adata,doutcar_isodata):

    Adata = { 'relatpos' : [],
              'cAiso'    : [],
              'cAdip'    : [],
              'cAdiag'   : [],
              'dgfactor' : {} }

    dgfactor = {}
    for atm, isotmp in doutcar_isodata.items():
        for iso, dat in isotmp.items(): 
            g = dat.gyro * 1e7 * HBAR_SI / NUCLEAR_MAGNETON
            dgfactor[iso] = g

    Adata['relatpos'] = relativepos['r_1i']
    Adata['dgfactor'] = dgfactor
    Adata['cAiso']    = doutcar_Adata['cAiso']
    Adata['cAdip']    = doutcar_Adata['cAdip']
    Adata['cAdiag']   = doutcar_Adata['cAdiag']

    return Adata

###############################################################
#
#        MAKE Q TENSOR FILE
#
###############################################################

def Qfile(relativepos,doutcar_Qdata,doutcar_isodata,fout,ignoreatm=None,_rangeformat='vertex'):

    # V/A^2 -> Hatree/(Bohr_radius^2)
    # Converting factor : 
    """
        here, the V means eV

        EFG(V/A^2) = EFG * eV_TO_HARTREE / (ANGSTROM_TO_BOHR^2)
                   = EFG'(Hartree/Bohr_radius^2)

    """
    VpAsq_TO_HARTREEpBOHRsq = eV_TO_HARTREE / pow(ANGSTROM_TO_BOHR,2)

    # eQ [unit : barn = 10^-28 m^2 ]

    deQ = []
    for atm, isotmp in doutcar_isodata.items():
        for iso, dat in isotmp.items(): 
            # barn to "milibarn(10^-31 m^2) * 10^1" (10^-30 m^2)
            eQ = dat.q * pow(10,2) # * 10^-2 (10^-28 m^2) = (10^-30 m^2)
            deQ = addlist(deQ, [iso,eQ])

    # Farthest atoms EFG tensor
    detcEFG = []
    i=0
    for atm, isotmp in doutcar_isodata.items():
        if atm != ignoreatm :
            idx = relativepos['farthest_atoms_idx'][i]
            EFG = np.array(doutcar_Qdata['EFG'][idx]) * VpAsq_TO_HARTREEpBOHRsq
            i += 1
            for iso, dat in isotmp.items():
                detcEFG = addlist(detcEFG,[iso, EFG])

    Dict = Qfileinfo() 
    # Make Q tensor file
    with open(fout, 'w') as f:
        # HEAD
        f.write(divline())
        f.write(Qfilehead())
        # ISOTOPE DATA
        f.write(divline())
        for iso_eQ in deQ:
            f.write(Dict['eQline'].format(*iso_eQ))
        for iso_etcEFC in detcEFG:
            f.write(Dict['etcline'].format(iso_etcEFC[0],*iso_etcEFC[1]))

        # CELL RANGE 
        if _rangeformat=='coord':
            f.write(Dict['minline'].format(*relativepos['MinDif[A]']))
            f.write(Dict['maxline'].format(*relativepos['MaxDif[A]']))
        elif _rangeformat=='vertex':
            f.write(Dict['vertex'].format(*relativepos['vertex'][0],\
                                          *relativepos['vertex'][1],\
                                          *relativepos['vertex'][2],\
                                          *relativepos['vertex'][3],\
                                          *relativepos['vertex'][4],\
                                          *relativepos['vertex'][5],\
                                          *relativepos['vertex'][6],\
                                          *relativepos['vertex'][7]))

        # MAIN DATA
        doutcar_Qdata['cEFG']=[]
        f.write(divline())
        for i,atm in enumerate(relativepos['r_1i']):
            corrected_EFG = np.array(doutcar_Qdata['EFG'][i])*VpAsq_TO_HARTREEpBOHRsq
            f.write(Qtensorform().format(*atm,*corrected_EFG))

            doutcar_Qdata['cEFG'].append(corrected_EFG)

    return doutcar_Qdata

def _Qdata(relativepos,doutcar_Qdata,doutcar_isodata):

    # V/A^2 -> Hatree/(Bohr_radius^2)
    # Converting factor : 
    """
        here, the V means eV

        cEFG(V/A^2) = EFG * eV_TO_HARTREE / (ANGSTROM_TO_BOHR^2)
                    = EFG'(Hartree/Bohr_radius^2)

    """
    
    Qdata = { 'relatpos': [],
              'cEFG' : [],
              'ceQ'  : {} }

    # eQ [unit : barn = 10^-28 m^2 ]
    deQ = {}
    for atm, isotmp in doutcar_isodata.items():
        for iso, dat in isotmp.items(): 
            # barn to "milibarn(10^-31 m^2) * 10^1" (10^-30 m^2)

            eQ = dat.q * pow(10,2) # * 10^-2 (10^-28 m^2) = (10^-30 m^2)
            deQ[iso] = eQ

    Qdata['ceQ'] = deQ
    Qdata['relatpos'] = relativepos['r_1i']
    Qdata['cEFG'] = doutcar_Qdata['cEFG'] 

    return Qdata


def spin_from_outcar(foutcar):
    """S = |total magnetization| / 2 from the last 'number of electron ... magnetization' line of the OUTCAR."""
    mag = None
    with open(foutcar) as f:
        for line in f:
            if line.startswith(" number of electron") and "magnetization" in line:
                mag = float(line.split("magnetization")[1].split()[0])
    if mag is None:
        pyform._error(f"{foutcar}: no 'magnetization' line, the defect spin cannot be read; is it a spin-polarized run?")
    S = round(abs(mag)) / 2
    if S == 0 or abs(abs(mag) - 2*S) > 0.05:
        pyform._error(f"{foutcar}: total magnetization {mag} is not 2S for a half-integer S > 0")
    return S

class do_mkAQDdata:

    def __init__(self,foutcar,fposcar,dfname,dfposition_frac=None,_ignoreatm=None,_addhfcore=True,_spin=None,_rangeformat='vertex'):
        self.foutcar = foutcar
        self.fposcar = fposcar
        self.dfname = dfname
        self._ignoreatm = _ignoreatm
        self.dfposition_frac = dfposition_frac
        # The defect spin S is read for the A file alone (mkAdata): a Q file does not depend on it, and an EFG run need
        # not be spin polarized (review 20261013-180049 T1).
        self._spin_given = _spin
        self._spin = None
        self._addhfcore = _addhfcore
        self._rangeformat = _rangeformat
        print("\n\n\t=========================================\n")
        print(f"\t\t\t{type(self).__name__} : \n")
        print("\t=========================================\n\n")

        #print("\t",self.foutcar)
        #print("\t",self.fposcar)
        #print("\t",self.dfname)
        #print("\t",self.dfposition_frac)
        #print("\n")
        # df position : fractional ..
    
        ###########################################################
        
        # Read outcar informations (A,Q,D tensor,,, isotopes... etc)
    
        self.outcar_spHam = extSpinHamilTensors(foutcar,print_io=False,addhfcore=_addhfcore)
        self.readA = (self.outcar_spHam).readA
        self.readQ = (self.outcar_spHam).readQ
        self.readD = (self.outcar_spHam).readD

        self.Qdata = (self.outcar_spHam).Qdata
        self.Adata = (self.outcar_spHam).Adata
        self.Ddata = (self.outcar_spHam).Ddata
        self.isodata = (self.outcar_spHam).isodata

        ###########################################################
    
        # Calculate relative position 
    
        self.poscar = read_POSCAR(fposcar)
        self.poscarmod = mkDefectKeyInPOSCAR(self.poscar,dfname,defectposition_frac=dfposition_frac)
        #if you want to make central defect
        #    or the defect information is not in poscar data
        #     , you should add defect position 
        #To obtain relative position, defect position is needed
        self.relativepos  = relativeDistance(self.poscarmod["latt_constant"],\
                                             self.poscarmod["unit_cell_parameter"],\
                                             self.poscarmod["unit_cell_information"],\
                                             self.poscarmod['central_defect'],_rangeformat=self._rangeformat)

        ###########################################################
        #    
        #   # Make A tensor file
        #
        #   if 'T' in self.readA :
        #       Afile(self.relativepos, self.outcar_spHam.Adata, self.outcar_spHam.isodata,self.fout_Afile)
        #   if 'T' in self.readQ : 
        #       Qfile(self.relativepos, self.outcar_spHam.Qdata, self.outcar_spHam.isodata,self.fout_Qfile,ignoreatm=_ignoreatm)
        #
        ############################################################

    def mkAdata(self,fout_Afile):
        self.fout_Afile = fout_Afile
        if 'T' not in self.readA :
            sys.exit("\t Error, outcar doesn't include A tensor information")
        # The defect spin S, written in the provenance line of the A file: from the OUTCAR's total magnetization;
        # a --spin given on the command line must agree with it (the values themselves do not depend on it).
        self._spin = spin_from_outcar(self.foutcar)
        if self._spin_given is not None and abs(check_spin(self._spin_given) - self._spin) > 1e-9:
            pyform._error(f"--spin {self._spin_given:g} contradicts the OUTCAR {self.foutcar}: its total magnetization "
                          f"gives S = {self._spin:g}")
        self.outcar_spHam.Adata = Afile(self.relativepos, self.outcar_spHam.Adata, self.outcar_spHam.isodata,fout_Afile,self._spin,_rangeformat=self._rangeformat)
    def mkQdata(self,fout_Qfile):
        self.fout_Qfile = fout_Qfile
        if 'T' not in self.readQ : 
            sys.exit("\t Error, outcar doesn't include Q tensor information")
        self.outcar_spHam.Qdata = Qfile(self.relativepos, self.outcar_spHam.Qdata, self.outcar_spHam.isodata,fout_Qfile,ignoreatm=self._ignoreatm,_rangeformat=self._rangeformat)

    def getAdata(self):
        return _Adata(self.relativepos,self.outcar_spHam.Adata,self.outcar_spHam.isodata)
    def getQdata(self):
        return _Qdata(self.relativepos,self.outcar_spHam.Qdata,self.outcar_spHam.isodata)
    def getDdata(self):
        return self.outcar_spHam.Ddata

    def __repr__(self,read):
        message = ""
        message += "\n"

        message+="\n\n\t=========================================\n"
        message+="\t READ OUTCAR (KEY : outcar_spHam) \n"
        message+="\t=========================================\n\n"
        message+=self.outcar_spHam.__repr__()
        message+="\n\n\t=========================================\n"
        message+="\t READ POSCAR (KEY : poscar or poscarmod) \n"
        message+="\t=========================================\n\n"
        message+=self.prt_poscar()
        message+="\n\n\t=========================================\n"
        message+="\t          MAKE A,Q FILE                     \n"
        message+="\t=========================================\n\n"
        message+=self.prt_relativepos()
        message+="\n\t--- Save files --- \n"
        if 'T' in (self.outcar_spHam).readA and read[0] :
            message+=f"\tAfile : {self.fout_Afile}\n"
        else:
            message+=f"\tAfile : No read (see incar)\n"
            
        if 'T' in (self.outcar_spHam).readQ and read[1] : 
            message+=f"\tQfile : {self.fout_Qfile}\n"
        else:
            message+=f"\tfile : No read (see incar)\n"
        message+="\n\n\t=========================================\n"
        message+="\t                 DONE                     \n"
        message+="\t========================================="

        return message

    def prt_poscar(self):
        message = "" 
        message += f"\n\t POSCAR from {self.fposcar} (KEY : poscar)\n"
        message += "\t" + "-"*30 + "\n"
        message += f"\t *atom_data : {self.poscar['atom_data']}\n" 
        message += f"\t *atom_data_num : {self.poscar['atom_data_num']}\n"
        message += f"\t *unit_cell_information :  \n\n"
        message += printform(3,["frac_x","frac_y","frac_z"],\
                                [v[1:] for v in self.poscar["unit_cell_information"]])
        return message

    def prt_relativepos(self):
        message = "" 
        message += f"\n\t Relative positions from central defect (KEY : relativepos)\n"
        message += "\t" + "-"*30 + "\n\n"
        message += f"\t *MinDif[A] : {np.round(np.asarray(self.relativepos['MinDif[A]']),3)}\n"
        message += f"\t *MaxDif[A] : {np.round(np.asarray(self.relativepos['MaxDif[A]']),3)}\n"
        message += f"\t *farthest_atoms_idx : {self.relativepos['farthest_atoms_idx']}\n"
        for i in range(len(self.relativepos['farthest_atoms_info'])):
            if i==0:
                l=self.relativepos['farthest_atoms_info'][i]
                message += f"\t *farthest_atoms_info : {l[0]} , {np.round(np.asarray(l[1:]),3)}\n"
            else:
                message += f"\t                      : {l[0]} , {np.round(np.asarray(l[1:]),3)}\n"

        message += f"\t *r_1i (cart) :  \n\n"
        message += printform(4,["atom","r_1ix","r_1iy","r_1iz"],\
                                [v[:] for v in self.relativepos['r_1i']])
        return message

if __name__ =='__main__':

#    ###########################################################
#    
#    # Read outcar informations (A,Q,D tensor,,, isotopes... etc)
#
#    from extract_spinHamils import extSpinHamilTensors
#    outcar_spHam = extSpinHamilTensors("./src/OUTCAR_flat",print_io=False)
#    #print(outcar_spHam)
#    ###########################################################
#
#    # Calculate relative position 
#
#    from ioPoscar import read_POSCAR
#    from dictproc import print_dict
#    poscar = read_POSCAR("./src/POSCAR_flat")
#    #if you want to make central defect
#    #    or the defect information is not in poscar data
#    #     , you should add defect position 
#    #To obtain relative position, defect position is needed
#    poscarmod = mkDefectKeyInPOSCAR(poscar,'VB',defectposition_frac=None)
#    
#    #relative positions are obtained as cartesitan coord.
#    relativepos  = relativeDistance(poscarmod["latt_constant"],\
#                                    poscarmod["unit_cell_parameter"],\
#                                    poscarmod["unit_cell_information"],\
#                                    poscarmod['central_defect'])
#    #print(relativepos)
#    #print_dict(relativepos)
#
#    ###########################################################
#
#    # Make A tensor file
#    Afile(relativepos, outcar_spHam.Adata, outcar_spHam.isodata,'./AtensorTest')
#    Qfile(relativepos, outcar_spHam.Qdata, outcar_spHam.isodata,'./QtensorTest')
#
    ###########################################################


#    ###########################################################
#    # Example : 
#    ###########################################################
#    foutcar = "./src/OUTCAR_g90_VB_S1"
#    fposcar = "./src/POSCAR_g90_VB_S1_rlx_mod_vb_mvatoms"
#    dfname  = 'VB'
#    fout_Afile = "./AtensorTest"
#    fout_Qfile = "./QtensorTest"
#    a = do_mkAQDdata(foutcar,fposcar,dfname,fout_Afile,fout_Qfile,dfposition_frac=None)
#    print(a.__repr__())
#    
#
    ###########################################################
    # Main : 
    ###########################################################
    import sys
    if len(sys.argv)!=7:
        print("[1] foutcar")
        print("[2] fposcar")
        print("[3] dfname")
        print("[4] fout_Afile")
        print("[5] fout_Qfile")
        print("[6] fout_Qfile_ignoreatm")
        sys.exit('Check the arguments')

    foutcar = sys.argv[1]
    fposcar = sys.argv[2]
    dfname  = sys.argv[3]
    fout_Afile = sys.argv[4]
    fout_Qfile = sys.argv[5]
    Qfile_ignoreatm = sys.argv[6]

    a = do_mkAQDdata(foutcar,fposcar,dfname,dfposition_frac=None,_ignoreatm=Qfile_ignoreatm)
    print(a)

    
