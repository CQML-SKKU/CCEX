# -*- coding: utf-8 -*-
"""
Created on Tue Jul 18 00:41:14 2023

@author: skku
"""

import sys
import numpy as np
import copy
from scipy.optimize import curve_fit, OptimizeWarning
from CoherencePlot_v2 import EasyPlot, EasyPlotCoherence, EasyPlotCoherence_Power, EasyPlotCoherence_T2
from CoherencePlot_v2 import IntensityPlot
from scipy.stats import norm
from coherence import writefile, read_coherence, error_correct, ensemble_average

def function_analysis_coherence(L):
    return -np.log(np.abs(L)) #-ln(L)

class DataSpace:
    
    """
        ---------------------------------------------------------------
        
        filenames_keys : "VAR1 | VAR2 | ..."
        filenames_dict : {"VAR1 | VAR2 | VAR3 .." : fileloacation }
    
        ---------------------------------------------------------------
        
        variables_keys : [ VAR1, VAR2, .. VARN ]
    
        variables_list : [ [ "VAR1" , [0,1,2,3,4] ],
                           [ "VAR2" , [0,1,2,3,4] ],
                               :
                           [ "VARN" , [0,1,2,3,4] ] ]
            
        variables_dict : { "VAR1" : [0,1,2,3,4],
                           "VAR2" : [0,1,2,3,4],
                             :
                           "VARN" : [0,1,2,3,4] }
            
        ---------------------------------------------------------------
        
        data_keys : "VAR1 | VAR2 | ..."
        data_dict : {"VAR1 | VAR2 | VAR3 .." : data }
        
        ---------------------------------------------------------------
            
    """
    
    def __init__(self):
        self.filenames_dict = {}
        self.filenames_keys = []
        
        self.variables_list = []
        self.variables_keys = []
        self.variables_dict = {}

        self.data_dict = {}        
        self.data_keys = self.filenames_keys
        
        self.key_string_tmp = ""
        
        self.fig_dict = {}
        self.line_dict = {}
        self.fig_t2 = {}
        self.fig_power = {}
        self.line_t2 = {}
        self.line_power = {}

        self.verbosity=False
    
    def get_filenames_format(self,directory,file):
        self.format_directory = directory
        self.format_file = file
        
    def get_variable(self,variable):
        variable_key = copy.copy(variable[0])
        variable_values = copy.copy(variable[1])
        self.variables_list.append(variable)
        self.variables_keys.append(variable_key)
        self.variables_dict.update({variable_key : variable_values})
        
        if self.filenames_keys == []:
            self.filenames_keys = variable_values  
        else:
            filenames_keys_copy = copy.deepcopy(self.filenames_keys)
            self.filenames_keys = [] 
            for variable_value_old in filenames_keys_copy:
                for variable_value_new in variable_values:
                    variable_value_update = variable_value_old + self.spliter() + variable_value_new    
                    self.filenames_keys.append(variable_value_update)
        self.data_keys = self.filenames_keys
        
    def get_variable_all(self,variables):
        for i,v in enumerate(variables):
            self.get_variable(v)
                        
    def update_filenames(self):
        
        for filenames_key in self.filenames_keys:            
            for i, variable_key in enumerate(self.variables_keys):                
                variable_value = self._get_variables_from_key(filenames_key)[i]
                format_file = None
                
                if i==0:
                    format_file = copy.deepcopy(self.format_file)
                    format_directory = copy.deepcopy(self.format_directory)
                else:
                    format_file = copy.deepcopy(self.filenames_dict[filenames_key])
                    format_directory = ""
                
                self.filenames_dict[filenames_key] = self._replace_filename_with_variable\
                                                (
                                                   format_directory, \
                                                   format_file, \
                                                   variable_key,\
                                                   variable_value\
                                                )
                                                    
    """______________________________________________________________________
    
        Read t2 file
    ______________________________________________________________________ """
 
    def readtry_t2file(self,timeunit_converter=1):

        counts={}
        for key, filename in self.filenames_dict.items():
            
            try:
                self.data_dict[key] = np.array(self._read_file(filename)).astype(float)
                self.data_dict[key][:,2] = self.data_dict[key][:,2] * timeunit_converter
                
                rows_to_delete = []
                count=0
                for i, vlist in enumerate(self.data_dict[key]):
                    c = vlist[0]
                    s = vlist[1]
                    t2 = vlist[2]
                    n = vlist[3]
                    if t2==0.0 and n==0.0:
                        rows_to_delete.append(i)
                    else:
                        count+=1
                
                self.data_dict[key] = np.delete(self.data_dict[key],rows_to_delete,axis=0)
                counts[key] = count
            except:
                pass;
        return counts
                    


    """______________________________________________________________________
    
        Read coherence result files (real value)
    ______________________________________________________________________ """
 
    def read_file_coherence(self,timeunit_converter=1.0):
        
        for key, filename in self.filenames_dict.items():
            self.data_dict[key] = np.array(self._read_file(filename)).astype(complex)
            self.data_dict[key][:,0] = self.data_dict[key][:,0] * timeunit_converter
    
    def do_fitting_coherence(self,space=1,init_fitPower=None):
        
        self.t2_dict = {}
        self.power_dict = {}
        
        for key, data in self.data_dict.items():
            # print(key,data)
            ht,hL,fL,T2,n = self._do_fitting_coherence(data[:,0], data[:,1],space=space,init_fitPower=init_fitPower)          
            self.t2_dict[key] = T2
            self.power_dict[key] = n
    
    """______________________________________________________________________
    
        Read coherence result files (Complex value)
    ______________________________________________________________________ """
    
    def read_file_coherence_complex(self,timeunit_converter=1.0, interval_converter=1):
        
        for key, filename in self.filenames_dict.items():
            data = np.array(self._read_file(filename)).astype(complex)
            t = data[:,0]
            coherence_real = data[:,1]
            coherence_imag = data[:,2]
            coherence_complex = coherence_real + 1j * coherence_imag
            self.data_dict[key] = np.empty((len(t),2),dtype=complex)
            
            # convert time unit
            self.data_dict[key][:,0] = t.astype(complex) * timeunit_converter
            self.data_dict[key][:,1] = coherence_complex
            
            # convert interval
            self.data_dict[key][:,0] = self.data_dict[key][:,0][::interval_converter]
            self.data_dict[key][:,1] = self.data_dict[key][:,1][::interval_converter]
            
    def readtry_file_coherence_complex(self,timeunit_converter=1.0, interval_converter=1):
        for key, filename in self.filenames_dict.items():
            try:
                data = np.array(self._read_file(filename)).astype(complex)
                t = data[:,0]
                coherence_real = data[:,1]
                coherence_imag = data[:,2]
                coherence_complex = coherence_real + 1j * coherence_imag
                self.data_dict[key] = np.empty((len(t),2),dtype=complex)
                
                # convert time unit
                self.data_dict[key][:,0] = t.astype(complex) * timeunit_converter
                self.data_dict[key][:,1] = coherence_complex
                
                # convert interval
                self.data_dict[key][:,0] = self.data_dict[key][:,0][::interval_converter]
                self.data_dict[key][:,1] = self.data_dict[key][:,1][::interval_converter]
                
            except: pass;

    def do_fitting_coherence_complex(self,ista=0,iend=None,space=1,init_fitPower=None):
        
        self.t2_dict = {}
        self.power_dict = {}
        
        for key, data in self.data_dict.items():
            if iend ==None:
                ht,hL,fL,T2,n = self._do_fitting_coherence(data[:,0][ista:].real, data[:,1][ista:].real,space=space,init_fitPower=init_fitPower)
            else:
                ht,hL,fL,T2,n = self._do_fitting_coherence(data[:,0][ista:iend].real, data[:,1][ista:iend].real,space=space,init_fitPower=init_fitPower)
            self.t2_dict[key] = T2
            self.power_dict[key] = n
    
    def do_fitting_ramsey_complex(self,space=1,yconverter=1.0,init_fitPower=None):
        
        self.t2_dict = {}
        self.power_dict = {}
        
        for key, data in self.data_dict.items():
            data[:,1] = data[:,1] * yconverter
            ht,hL,fL,T2,n = self._do_fitting_ramsey(data[:,0].real, data[:,1].real,space=space,init_fitPower=init_fitPower)
            # res = self._do_fitting_ramsey(data[:,0].real, data[:,1].real,space=space,init_fitPower=init_fitPower)
            # print(res)
            self.t2_dict[key] = T2
            self.power_dict[key] = n    
    
    """______________________________________________________________________
    
        Read coherence result files (Complex value when they are : "a+bi" or "a b")
    ______________________________________________________________________ """
    
    def read_file_coherence_complex_v2(self,timeunit_converter=1.0, interval_converter=1):

        count = 0
        for key, filename in self.filenames_dict.items():
            try:
                t, L = read_coherence(filename, timeunit_converter, interval_converter)
            except Exception:
                print(f"\t\t ! Fail to read : {filename}")
                continue
            self.data_dict[key] = np.empty((len(t),2),dtype=complex)
            self.data_dict[key][:,0] = t
            self.data_dict[key][:,1] = L
            count += 1
            if self.verbosity:
                print(f"\t\tRead done : {filename}")

        return count

    """______________________________________________________________________
    
        Read density matrix result files (Complex value when they are : "a+bi")
    ______________________________________________________________________ """
    
    # format : t rho_00 rho_01 rho_02   4x4
    #            rho_03 rho_10 rho_11
    #            rho_12 rho_13 rho_20
    #            rho_21 rho_22 rho_23
    #            rho_30 rho_31 rho_32
    #            rho_33
    def read_file_dm_complex_v2(self,timeunit_converter=1.0, interval_converter=1):
        
        for key, filename in self.filenames_dict.items():
            
            data = None
            # Check if it is possible to read
            try:
                data = self._read_file(filename)
                
                if len(data[0]) == 4: # then density matrix form
                    pass;
                else:
                    sys.exit("Error, data[0] length (=%d) is not 4 \n ",len(data[0]))
                    
            except:
                pass;
            
            if data != None:
                # read
                dim = 0
                dm_ti = []
                dm_ts = []
                ts = []  
                for i, d in enumerate(data):
                    
                    
                    # Update dm_ts
                    if ((i!=0) and (len(d)==4)):
                        dim = np.sqrt(len(dm_ti)).astype(int) # dm dimension
                        if dim*dim != len(dm_ti):
                            sys.exit("Error, The dimension is different.")
                        dm_ts.append(np.asarray(dm_ti).reshape(dim,dim).astype(complex))
                        dm_ti = []
                    
                    # Get data
                    if len(d) == 4: # t dm_00 dm_01 dm_02
                        ts.append(d[0])
                        dm_ti.append(d[1])
                        dm_ti.append(d[2])
                        dm_ti.append(d[3])
                    else:
                        for j in range(len(d)):
                            dm_ti.append(d[j])
                    
                    # Get last data
                    if i==(len(data)-1):
                        dim = np.sqrt(len(dm_ti)).astype(int) # dm dimension
                        if dim*dim != len(dm_ti):
                            sys.exit("Error, The dimension is different.")
                        dm_ts.append(np.asarray(dm_ti).reshape(dim,dim).astype(complex))
                        dm_ti = []
                
                            
                self.data_dict[key] = (np.asarray(ts).astype(float),np.asarray(dm_ts))
                self._mask_dm(key)
                
    
    def _mask_dm(self,key):
        ts,dm_ts = self.data_dict[key]
        for i,dm_ti in enumerate(dm_ts):
            if i==0:
                idx_mask = np.where((np.nan_to_num(dm_ts[0]) == dm_ts[0])==False)
            self.data_dict[key][1][i][idx_mask] = 0.0+0.0j

    """______________________________________________________________________
    
        Refine result
    ______________________________________________________________________ """
    
    def refine_nan2num_coherence_all(self):
        
        for key, data in self.data_dict.items():
            self.data_dict[key][:,1] = np.nan_to_num(data[:,1])
            
    def refine_onlyimag_coherence_all(self):
         
        for key, data in self.data_dict.items():
            self.data_dict[key][:,1] = np.imag(data[:,1])
            
    def refine_onlyreal_coherence_all(self):
             
        for key, data in self.data_dict.items():
            self.data_dict[key][:,1] = np.real(data[:,1])

    def refine_error_coherence(self,key, key_ref, tol=1.0):
        self.data_dict[key][:,1] = self._error_correction_coherence(self.data_dict[key][:,1],self.data_dict[key_ref][:,1],tol)
    
    def refine_error_coherence_all(self,key_ref, tol=1.0):
        for key, data in self.data_dict.items():
            self.refine_error_coherence(key, key_ref, tol)
    
    def _error_correction_coherence(self,L_wD,L_nD,tolerance=2.0): #L_nD = L_ref
        return error_correct(L_wD, L_nD, tolerance)

    def refine_result_with_function(self,func):
        for key, data in self.data_dict.items():
            self.data_dict[key][:,1] = func(data[:,1])

    """______________________________________________________________________
    
        Read bath file
    ______________________________________________________________________ """

    def read_file_bath(self):
        
        for key, filename in self.filenames_dict.items():
            
            data = self._read_file(filename)
            data_new = []
            for i,line in enumerate(data):
                try:
                    line_new = [float(line[0]),float(line[1]),float(line[2]),line[3]]
                except:
                    # Defect file case
                    line_new = [float(line[0]),float(line[1]),float(line[2])]
                data_new.append(line_new)                    
                            
            self.data_dict[key] = data_new

    """______________________________________________________________________
     
         Read file general
    ______________________________________________________________________ """
       
    def read_file(self):
        
        for key, filename in self.filenames_dict.items():
            self.data_dict[key] = self._read_file(filename)
             
    #########################################################################
    
    # Get functions : if you give variables n#, then result would be returned
    
    def get_filename(self,*args):
        
        key_string = args[0]
        
        for arg in list(args)[1:]:
            key_string+= self.spliter() + str(arg)
        
        return self.filenames_dict[key_string]
    
    def get_data(self,*args):
        
        key_string = args[0]
        
        for arg in list(args)[1:]:
            key_string+= self.spliter() + str(arg)
        
        return self.data_dict[key_string]
    
    def get_key(self,*args):
        
        key_string = args[0]
        
        for arg in list(args)[1:]:
            key_string+= self.spliter() + str(arg)
        
        return key_string
    
    def get_key_in_order(self,variable_name,variable_value):
        
        
        if self.key_string_tmp == "":
            self._initialize_key_string_tmp()

        i_variable = list(self.variables_dict.keys()).index(variable_name)
        self.key_string_tmp = self.key_string_tmp.replace(f"VAR{i_variable}",variable_value)
                
        return self.key_string_tmp
    
    def get_keyformat(self):
        
        KEYFORMAT = ""
        nkey = len(self.variables_keys)
        for i,k in enumerate(self.variables_keys):
            KEYFORMAT += k
            if i!=nkey-1:
                KEYFORMAT += self.spliter()
                
        return KEYFORMAT
                        
    
    def get_key_from_keyformat(self,variable_name,variable_value,base=None):
        
        if not variable_name in self.variables_keys:
            sys.exit(f"Error! input key : {variable_name} , possible key : {self.variable_keys}")
        
        if not str(variable_value) in self.variables_dict[variable_name]:
            sys.exit(f"Error! input key, value :  {variable_name}:{variable_value}, possible value : {self.variables_dict[variable_name]}")
        
        keyformat=""
        
        if base==None or base=="":
            keyformat = self.get_keyformat()
        else:
            keyformat = copy.deepcopy(base)
        
        return keyformat.replace(variable_name,variable_value)
    
    
    """______________________________________________________________________
     
         get data or its T2,n for one variable
    ______________________________________________________________________ """
    
    # args : 
    #    x_variable_name : the variable name that you want to get
    #    other_variable_info : the names(keys) and values of variable_dict except x_variable_name
    #                  : (dict) other_variable_info = {key : value , key: value ...}
    #    
    
    def get_data_for_one_variable(self,x_variable_name,other_variable_info, passio=False):
      
        x_variable_values = self._get_variable_values(x_variable_name)
        keybase = self._get_keybase(other_variable_info)
        
        x_variable_values_list = []
        datalist = [] 
        
        for i,v in enumerate(x_variable_values):
            if not passio:
                key = self.get_key_from_keyformat(x_variable_name, v, base=keybase)
                datalist.append(self.data_dict[key])
                x_variable_values_list.append(v)
            else:
                try:
                    key = self.get_key_from_keyformat(x_variable_name, v, base=keybase)
                    datalist.append(self.data_dict[key])
                    x_variable_values_list.append(v)
                except KeyError:
                    pass;
            
        return x_variable_values_list, datalist
    
    def get_params_for_one_variable(self,x_variable_name,other_variable_info):
        
        x_variable_values = self._get_variable_values(x_variable_name)
        keybase = self._get_keybase(other_variable_info)
            
        t2list = [] 
        nlist = []
      
        for i,v in enumerate(x_variable_values):
            key = self.get_key_from_keyformat(x_variable_name, v, base=keybase)
            T2 = self.t2_dict[key]
            n = self.power_dict[key]
            t2list.append(T2)
            nlist.append(n)
            
        return x_variable_values,t2list,nlist
    
    """______________________________________________________________________
     
         2d Plot
     ______________________________________________________________________ """
    
    # args : 
    #    x_variable_name : the variable name that you want to get
    #    other_variable_info : the names(keys) and values of variable_dict except x_variable_name
    #                  : (dict) other_variable_info = {key : value , key: value ...}
    #    
    
    def plot_t2(self,x_variable_name,other_variable_info,xlabel=None,ylabel=r"$T_2$",xunit="",yunit=r"$\mu s$",plotlabel="",xlistrep=[]):

        #get values       
        xlist, t2list,nlist = self.get_params_for_one_variable(x_variable_name, other_variable_info)
        
        if xlabel==None:
            xlabel = x_variable_name

        try:
            xlist = xlist.astype(float)
        except:
            pass
        
        if xlistrep != []:
            xlist = copy.copy(xlistrep)
        
        self.fig_t2 = EasyPlotCoherence_T2(xlabel,xunit=xunit,yunit=yunit,ylabel=ylabel)
        self.line_t2 = self.fig_t2.do(xlist,t2list,variable=plotlabel)
        
        return (np.array(xlist),t2list)
            
    def plot_power(self,x_variable_name,other_variable_info,xlabel=None,xunit="",yunit=r"$\mu s$",plotlabel="",xlistrep=[]):

        #get values       
        xlist, t2list,nlist = self.get_params_for_one_variable(x_variable_name, other_variable_info)
        
        if xlabel==None:
            xlabel = x_variable_name

        try:
            xlist = xlist.astype(float)
        except:
            pass
        
        if xlistrep != []:
            xlist = copy.copy(xlistrep)
        
        self.fig_power = EasyPlotCoherence_Power(xlabel,xunit=xunit,yunit=yunit)
        self.line_power = self.fig_power.do(xlist,nlist,variable=plotlabel)    
        
        return (xlist,nlist)
    
    def addplot_t2(self,x_variable_name,other_variable_info,plotlabel=""):
        
        #get values       
        xlist, t2list,nlist = self.get_params_for_one_variable(x_variable_name, other_variable_info)
        try:
            xlist = xlist.astype(float)
        except:
            pass
        
        self.line_t2 = self.fig_t2.do(xlist,t2list,variable=plotlabel)
        
    def addplot_power(self,x_variable_name,other_variable_info,plotlabel=""):
        
        #get values       
        xlist, t2list,nlist = self.get_params_for_one_variable(x_variable_name, other_variable_info)
        try:
            xlist = xlist.astype(float)
        except:
            pass
        
        self.line_power = self.fig_power.do(xlist,nlist,variable=plotlabel)

        
    def plot_coherence_all(self,x_variable_name,other_variable_info, passio= False, variable_name=None,xunit=r"$\mu$s", colors=[], alpha=None, markers=None, markersizes=None, linestyles=None, figuresize=(7,5)):
    
        #get values       
        xlist, datalist = self.get_data_for_one_variable(x_variable_name, other_variable_info, passio = passio)
        
        if variable_name==None:
            #variable_name=xlist
            pass;
        
        self.fig_coherence = EasyPlotCoherence(unit=xunit,figsize=figuresize)
        
        for i,data in enumerate(datalist):
            if len(colors) != 0:
                self.fig_coherence.Fig.color= colors[i]
            if markers != None:
                self.fig_coherence.Fig.marker = markers[i]
                if markersizes != None:
                    self.fig_coherence.Fig.markersize = markersizes[i]
            if linestyles != None:
                self.fig_coherence.Fig.ls = linestyles[i]
                
            if alpha !=None:
                self.fig_coherence.Fig.alpha = alpha

            if variable_name!=None:
                self.fig_coherence.do(data[:,0],data[:,1],variable=variable_name[i])
            else:
                self.fig_coherence.do(data[:,0],data[:,1])
                
            
    def plot_coherence_all_imag(self,x_variable_name,other_variable_info, passio= False, variable_name=None,xunit=r"$\mu$s", colors=[], alpha=None, markers=None, markersizes=None, linestyles=None, figuresize=(7,5)):
    
        #get values       
        xlist, datalist = self.get_data_for_one_variable(x_variable_name, other_variable_info, passio = passio)
        
        if variable_name==None:
            #variable_name=xlist
            pass;
        
        self.fig_coherence = EasyPlotCoherence(unit=xunit,figsize=figuresize)
        
        for i,data in enumerate(datalist):
            if len(colors) != 0:
                self.fig_coherence.Fig.color= colors[i]
            if markers != None:
                self.fig_coherence.Fig.marker = markers[i]
                if markersizes != None:
                    self.fig_coherence.Fig.markersize = markersizes[i]
            if linestyles != None:
                self.fig_coherence.Fig.ls = linestyles[i]
                
            if alpha !=None:
                self.fig_coherence.Fig.alpha = alpha

            if variable_name!=None:
                self.fig_coherence.do(np.real(data[:,0]),np.imag(data[:,1]),variable=variable_name[i])
            else:
                self.fig_coherence.do(np.real(data[:,0]),np.imag(data[:,1]))
            
    def addplot_coherence_all(self,figbase,x_variable_name,other_variable_info,variable_name=None,xunit=r"$\mu$s", colors=[], alpha=None, markers=None, markersizes=None, linestyles=None, figuresize=(7,5)):

        #get values       
        xlist, datalist = self.get_data_for_one_variable(x_variable_name, other_variable_info)
        
        if variable_name==None:
            #variable_name=xlist
            pass;
        
        for i,data in enumerate(datalist):
            if len(colors) != 0:
                figbase.Fig.color= colors[i]
            if markers != None:
                figbase.Fig.marker = markers[i]
                if markersizes != None:
                    figbase.Fig.markersize = markersizes[i]
            if linestyles != None:
                figbase.Fig.ls = linestyles[i]
                
            if alpha !=None:
                figbase.Fig.alpha = alpha

            if variable_name!=None:
                figbase.do(data[:,0],data[:,1],variable=variable_name[i])
            else:
                figbase.do(data[:,0],data[:,1])

        return figbase
         
    def plot_coherence_one(self,x_variable_name,x_variable_value,other_variable_info,variable_name=None,xunit=r"$\mu$s",showWfitting=False,iend=None,space=1,init_fitPower=None, adjust_x=0.2,adjust_y=0.7,adjust_space=0.0, colors=[], alpha=None, markers=None, markersizes=None, linestyles=None, figuresize=(7,5)):

        key = None

        if (x_variable_name != "None"):
            #get key
            keybase = self._get_keybase(other_variable_info)
            key = self.get_key_from_keyformat(x_variable_name, x_variable_value,base=keybase)

        else:
            if x_variable_value == "avg":
                key = x_variable_value
            else:
                sys.exit("Error! key is not \"avg\"")

        data = self.data_dict[key]
        
        self.fig_dict[key] = EasyPlotCoherence(unit=xunit,figsize=figuresize)
        
        if len(colors) != 0:
            self.fig_dict[key].Fig.color= colors[0]
        if markers != None:
            self.fig_dict[key].Fig.marker = markers[0]
            if markersizes != None:
                self.fig_dict[key].Fig.markersize = markersizes[0]
        if linestyles != None:
            self.fig_dict[key].Fig.ls = linestyles[0]
        if alpha !=None:
            self.fig_dict[key].Fig.alpha = alpha
            
        self.line_dict[key] = self.fig_dict[key].do(data[:,0].real,data[:,1].real,variable=variable_name)
        
        if showWfitting:
            ht,hL,fL,T2,n = (None, None, None, None, None) 
            if iend==None:
                ht,hL,fL,T2,n = self._do_fitting_coherence(data[:,0].real, data[:,1].real,space=space,init_fitPower=init_fitPower)
            else:
                ht,hL,fL,T2,n = self._do_fitting_coherence(data[:,0][:iend].real, data[:,1][:iend].real,space=space,init_fitPower=init_fitPower)
           
            self.fig_dict[key].Fig.ls='--'
            self.fig_dict[key].Fig.color='r'

            if T2!=None:

                if iend==None:
                    self.fig_dict[key].do(data[:,0].real,fL,variable=f"{variable_name} , fitting")
                else:
                    fL = self._fitting_function(data[:,0].real, T2, n)
                    self.fig_dict[key].do(data[:,0].real,fL,variable=f"{variable_name} , fitting")
                xmax = data[:,0][-1].real
                ymax = max(data[:,1].real)
                self.fig_dict[key].show_params(T2,n,showT2=True,showPower=True,adjust_x=xmax*adjust_x,adjust_y=ymax*adjust_y,adjust_space=adjust_space)
                
                self.fig_dict[key].Fig.ls=''
                self.fig_dict[key].Fig.marker='o'
                self.fig_dict[key].Fig.markersize=10
                self.fig_dict[key].Fig.color='k'
                
                self.fig_dict[key].scatter(ht,hL,variable=f"{variable_name} , highest points")
                
                self.fig_dict[key].Fig.marker=""
                self.fig_dict[key].Fig.markersize=2
                self.fig_dict[key].Fig.color=None
            else:
                print(f"\n\t\tFitting error {self.filenames_dict[key]}")

        return self.fig_dict[key]
            
    def addplot_coherence_one(self,figbase,x_variable_name,x_variable_value,other_variable_info,variable_name=None,xunit=r"$\mu$s",showWfitting=False,space=1,init_fitPower=None):
        #get key
        keybase = self._get_keybase(other_variable_info)
        key = self.get_key_from_keyformat(x_variable_name, x_variable_value,base=keybase)
        data = self.data_dict[key]
        
        self.line_dict[key] = figbase.do(data[:,0].real,data[:,1].real,variable=variable_name)
        
        if showWfitting:
            ht,hL,fL,T2,n = self._do_fitting_coherence(data[:,0], data[:,1],space=space,init_fitPower=init_fitPower)
            figbase.do(data[:,0],fL,variable=f"{variable_name} , fitting")
            figbase.show_params(T2,n,showT2=True,showPower=True,adjust_x=0.0,adjust_y=0.0,adjust_space=0.0)
            figbase.Fig.marker='o'
            figbase.Fig.markersize=10
            figbase.Fig.color='k'
            figbase.scatter(ht,hL,variable=f"{variable_name} , highest points")
            figbase.Fig.marker=""
            figbase.Fig.markersize=2
            figbase.Fig.color=None
        return figbase

    def plot_coherence_one_certainvalues(self,x_variable_name,x_variable_values,other_variable_info,variable_names=None,xunit=r"$\mu$s",showWfitting=False,space=1,init_fitPower=None):
        for i,v in enumerate(x_variable_values):

            if variable_names==None:
                self.plot_coherence_one(x_variable_name, v, other_variable_info,variable_name=variable_names,xunit=xunit,showWfitting=showWfitting,space=space,init_fitPower=init_fitPower)
            else:
                self.plot_coherence_one(x_variable_name, v, other_variable_info,variable_name=variable_names[i],xunit=xunit,showWfitting=showWfitting,space=space,init_fitPower=init_fitPower)
                
    def plot_coherence_all_certainvalues(self,x_variable_name,x_variable_values,other_variable_info,variable_names=None,xunit=r"$\mu$s",showWfitting=False,space=1,init_fitPower=None):
        figbase=None
        for i,v in enumerate(x_variable_values):
            
            if i==0:
                if variable_names==None:
                    figbase = self.plot_coherence_one(x_variable_name, v, other_variable_info,variable_name=variable_names,xunit=xunit,showWfitting=showWfitting,space=space,init_fitPower=init_fitPower)
                    
                else:
                    figbase = self.plot_coherence_one(x_variable_name, v, other_variable_info,variable_name=variable_names[i],xunit=xunit,showWfitting=showWfitting,space=space,init_fitPower=init_fitPower)
            else:
                if variable_names==None:
                    self.addplot_coherence_one(figbase,x_variable_name, v, other_variable_info,variable_name=variable_names,xunit=xunit,showWfitting=showWfitting,space=space,init_fitPower=init_fitPower)
                    
                else:
                    self.addplot_coherence_one(figbase,x_variable_name, v, other_variable_info,variable_name=variable_names[i],xunit=xunit,showWfitting=showWfitting,space=space,init_fitPower=init_fitPower)
        return figbase
                
      
    """______________________________________________________________________
     
         T2 histogram Plot
     ______________________________________________________________________ """    
     
                   

    """______________________________________________________________________
     
         Intensity Plot
     ______________________________________________________________________ """    
     
    def plot_coherence_intensity(self,x_variable_name,x_variable_values,other_variable_info,interpolate=(),xunit=r"$\mu$s",ylabel="B",yunit="G",cmap="RdYlBu_r",clim=[0,1]):
        
        variable_array = np.asarray(x_variable_values).astype(float)
        coherence_array = None
        time_array = None
        
        initial_key = None
        for i,x_variable_value in enumerate(x_variable_values):
            
            keybase = self._get_keybase(other_variable_info)
            key = self.get_key_from_keyformat(x_variable_name, x_variable_value,base=keybase)
            data = self.data_dict[key]
            
            t = data[:,0]
            L = data[:,1]
            
            if i==0:
                coherence_array = copy.deepcopy(L)
                time_array = copy.deepcopy(t)
                initial_key = copy.deepcopy(key)
            else:
                if (t != time_array).all():
                    sys.exit(f"{key} has different t value from {initial_key}" )
                coherence_array = np.vstack((coherence_array, L))
                
        print("coherence array shape : ",coherence_array.shape)
        print("variable array shape : ",variable_array.shape)
        print("time array shape : ",time_array.shape)
        
        Iplot = IntensityPlot()
        Iplot.clim=clim
        Iplot.cmap=cmap
        
        if interpolate == ():
            Iplot.intensityplot(time_array.real,variable_array,coherence_array.real,xunit=xunit,ylabel=ylabel,yunit=yunit)
        else:
            coherence_array = Iplot.interpolate_array(coherence_array, interpolate)
            variable_array = Iplot.interpolate_1d_array(variable_array, interpolate[0])
            time_array = Iplot.interpolate_1d_array(time_array, interpolate[1])
            Iplot.intensityplot(time_array.real,variable_array,coherence_array.real,xunit=xunit,ylabel=ylabel,yunit=yunit)
        
        return Iplot
    
    """______________________________________________________________________
    
         Coherence average
    ______________________________________________________________________ """  

    def average_coherence_all(self):
        for key in self.data_dict:
            if self.verbosity:
                print(f"\t\tAverage of variable \"{key}\"")
        try:
            return ensemble_average((d[:,0], d[:,1]) for d in self.data_dict.values())
        except ValueError as e:
            sys.exit(f"Error. average_coherence_all : {e}")

    #########################################################################
    def _replace_filename_with_variable(self,directory,file,key,value):
        string = directory + file
        return string.replace(key,value)
    
    def _get_variables_from_key(self,key):
        return key.split(self.spliter())
        
    def spliter(self):
        return " | "
    
    def _read_file(self,filename):
        
        data = []
        try:
            with open(filename,'r') as f:
                for i,line in enumerate(f):
                    data.append(line.strip().split())   
        except FileNotFoundError:
            print(f"Error ! {filename} file not found !")
            
        return data
    
    def _initialize_key_string_tmp(self):
        
        self.key_string_tmp = ""
        n_variable = len(list(self.variables_dict.keys()))    
        for i in range(n_variable):
            if i!=n_variable-1:
                self.key_string_tmp += f"VAR{i}" + self.spliter()
            else:
                self.key_string_tmp += f"VAR{i}"
    
    def _do_fitting_coherence(self,t,L,space=1,init_fitPower=None,showT2=True,showPower=True):
        
        if init_fitPower!=None:
            ht,hL,fL,T2,n =  self._fitting(t,L,space,init_fitPower=init_fitPower,function=self._fitting_function)
        else:
            ht,hL,fL,T2,n =  self._fitting(t,L,space)
        
        return ht,hL,fL,T2,n
    
    def _do_fitting_ramsey(self,t,L,space=1,init_fitPower=None,showT2=True,showPower=True):
        
        if init_fitPower!=None:
            ht,hL,fL,T2,n =  self._fitting(t,L,space,init_fitPower=init_fitPower,function=self._fitting_function_ramsey)
        else:
            ht,hL,fL,T2,n =  self._fitting(t,L,space)
            
        #print(T2,n)
        return ht,hL,fL,T2,n
    
    def _fitting_function(self,t,T2,n):
        return np.exp(-(t/T2)**n)

    def _fitting_function_ramsey(self,t,T2sta,n,w):
        return np.exp(-(t/T2sta)**n)*np.sin(w*t)

    def _get_keybase(self,other_variable_info):
        
        keybase=""
        for k,v in other_variable_info.items():
            variable_name = k
            variable_value = v
            keybase = self.get_key_from_keyformat(variable_name, variable_value, base=keybase)

        return keybase

    def _get_variable_values(self,x_variable_name):
        
        x_variable_values = self.variables_dict[x_variable_name]
        
        xlist = None
        try:
            xlist = np.asarray(x_variable_values).astype(float)
        except:
            xlist = np.asarray(x_variable_values)
            
        return np.asarray(x_variable_values)
    
    """______________________________________________________________________
     
         Fitting
     ______________________________________________________________________ """    
     
    def _fitting(self,t,coherence,space,function=None,init_fitPower=None):
        
        if (function == None):
            function = self._fitting_function
        
        high_t, high_coherence=self._pick_high_data(t,coherence,space)
        popt, pcov = self._find_best_initial_condition(high_t,high_coherence,init_fitPower=init_fitPower)
        if np.all(popt == None) and np.all(pcov == None):
            return (None,None,None,None,None)

        fit_func = np.abs(function(t,*popt))

        return (high_t,high_coherence,fit_func,popt[0],popt[1])
    
    def _find_best_initial_condition(self,high_t,high_coherence,function=None,init_fitPower=None):
        
        if (function == None):
            function = self._fitting_function
            
        arr = np.zeros((len(high_t),3))
        for i,y in enumerate(high_coherence):
            arr[i][0] = high_t[i]
            arr[i][1] = y 
            arr[i][2] = abs(1/np.e - y)

        sortedarr = np.asarray(sorted(arr, key=lambda x: x[2]))
        Notdone=True
        idx=0
        while(Notdone):

            init_fitT = sortedarr[idx][0] #t_free
            if init_fitPower == None:
                init_fitPower = 2    
            else:pass;

            try:
                popt,pcov = curve_fit(function,high_t,high_coherence,p0=[init_fitT,init_fitPower])
                T2=popt[0]
                n=popt[1]
                if n <= 0.1 or n > 4.0: pass;
                elif T2 > high_t[-1] or T2 <= 0: pass;
                else:
                    Notdone=False
                    break;
            except OptimizeWarning:pass;
            except RuntimeWarning:pass;
            except RuntimeError:pass;

            if idx == (len(sortedarr)-1):
                #sys.exit("Error! Couldn't get the T2!!")
                return (None, None)

            idx = idx + 1

        return (popt,pcov)
    
    def _BD_range(self,inital, max_range, spaceBD):
        if (int(inital)+int(spaceBD)) <= max_range:
            return range(inital, inital+int(spaceBD))
        else:
            return range(inital, max_range)
        
    def _pick_high_data(self,t,cce,space):
        high_point=[]
        high_point+=[0]
        index_num=0
        
        while index_num < len(t)-1:
            #among(index_num, index_num+1)
            min_diff=(cce[index_num+1]-cce[index_num])/(t[index_num+1]-t[index_num])
            diff_num=index_num+1

            #among(index_num, i)
            for i in self._BD_range(index_num+1, len(t), space):
                diff=(cce[i]-cce[index_num])/(t[i]-t[index_num])
                if diff > min_diff:
                    min_diff=diff
                    diff_num=i
            if diff_num == index_num+1:
                index_num = index_num+1
            else:
                index_num = diff_num
            high_point+=[index_num]

        sort_high=[]
        for v in high_point:
            if v not in sort_high:
                sort_high.append(v)
        sort_high.sort()
        high_point=sort_high
        
        #select high point
        high_time=[]
        high_cce=[]
        for i in high_point:
            high_time+=[t[i]]
            high_cce+=[cce[i]]

        return high_time, high_cce

    def _get_pdf_T2(self,keys=None):
        if keys==None:
            keys=self.data_dict.keys()

        if self.t2_dict == {}: 
            sys.exit("t2_dict is None")

        t2s = np.empty(0,dtype=float)

        for i, k in enumerate(keys):
            t2s = np.append(t2s,self.t2_dict[k])

        mu, std = np.mean(t2s), np.std(t2s)
        x = np.linspace(mu - 3*std, mu + 3*std, 100)
        y = norm.pdf(x, mu, std)
        return x, y

    def _get_key_frequentT2_from_pdf(self,probability,pdfx,pdfy,keys):
        idx = np.where(pdfy>probability)
        xmin = np.min(pdfx[idx])
        xmax = np.max(pdfx[idx])

        t2frequent_key = []
        for i, k in enumerate(keys):
            value = self.t2_dict[k]
            if value >= xmin and value <= xmax:
                t2frequent_key.append(k)

        return t2frequent_key

    
    #%%        
if __name__ == "__main__":
    
    timescale={'ns':1e-9, 'us':1e-6, 'ms':1e-3,'s':1}
    ####################################################################################################
    import argparse
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                    description="\
    ./DataSpace.py  -d  <path>/cal/pycce/result/STRAIN/BFIELDG/rawdata/ \
                    -f  gCCE2_hBNVB_HE_BFIELD_Conv10 \
                    -vn STRAIN BFIELD  \
                    -v  g95 g96 g97 flat t110  // Variable values for STRAIN \
                    -v  0 2 4 6 8  // Variable values for BFIELD")

    parser.add_argument('-d','--input_dir',dest='inputfile_dir',type=str,help="The directory of inputfile")
    parser.add_argument('-fi','--input',dest='inputfile',type=str,help="The inputfile name")
    parser.add_argument('-vn','--varnames',dest='variable_names',type=str,nargs='*',help="Variable names")
    parser.add_argument('-v','--vars',dest='variable_values',action='append',type=str,nargs='*',help="Variable values for each names")
    parser.add_argument('-ea','--ens_avg',dest='ensemble_average',action='store_true' ,help="Do ensemble average for all variables")
    parser.add_argument('-uc','--unit_conv',dest='unit_conversion',default="ms", choices=['ms','us','ns','s'] ,help="Time unit conversion")
    parser.add_argument('-pl','--plot',dest='do_plot',action='store_true', help="Plot for all data")
    parser.add_argument('-fit','--fitting',dest='do_fitting',action='store_true', help="Show fitting on plot")
    parser.add_argument('-fo','--fout',dest='outfile',nargs='*',type=str, help="output file")
    
    args = parser.parse_args()
    print("\t",args,end='\n\n')
    ####################################################################################################

    def is_2d(lst):
        if not lst:
            return False
        return isinstance(lst[0], list)

    def group_variables(args):
        groups = []
        current_group = []
        
        for arg in args:
            if arg in ['-v', '--var']:
                if current_group:
                    groups.append(current_group)
                    current_group = []
                else:
                    current_group.append(arg)
        if current_group:
            groups.append(current_group)

        return groups

    #if args.variable_values:
    #    grouped_variable_values = group_variables(sys.argv[1:])
    #    print("\tvariable_values : ",group_variables,end='\n\n')

    try:
        if (len(args.variable_names) == len(args.variable_values)):
            pass;
        else:
            print("variable_names : ",args.variable_names,end='\n\n') 
            print("variable_values : ",args.variable_values,end='\n\n') 
            sys.exit("variable_values length is different from variable_names")

    except TypeError:
        print("variable_names : ",args.variable_names,end='\n\n') 
        print("variable_values : ",args.variable_values,end='\n\n') 
        sys.exit("TypeError")
    ####################################################################################################

    variables = []
    for i,varname in enumerate(args.variable_names):
        variables.append([args.variable_names[i],args.variable_values[i]])


    print("\t directory : ", args.inputfile_dir,end='\n')
    print("\t file : ", args.inputfile,end='\n')
    print("\t variable list : ",variables,end='\n\n')

    Dsp = DataSpace()
    Dsp.get_filenames_format(args.inputfile_dir,args.inputfile)
    Dsp.get_variable_all(variables)

    # Change the variable_name in the file/directory into the actual value
    Dsp.update_filenames()

    # Read all files
    Dsp.read_file_coherence_complex_v2(timeunit_converter=timescale['ms']/timescale[args.unit_conversion])

    ####################################################################################################

    # 
    if (args.ensemble_average):

        # Write
        ts, L_ens = Dsp.average_coherence_all()
        writefile(f"{args.outfile[0]}",['f','c'],[ts.real,L_ens])


    if (args.do_plot):

        if len(args.variable_names) == 1:
            for i,varname in enumerate(args.variable_names):
                for j,varval in enumerate(args.variable_values[i]):
                    Fig = Dsp.plot_coherence_one(varname,varval,{} ,variable_name=varval,xunit=args.unit_conversion,showWfitting=args.do_fitting)
                    Fig.Fig.save(f"{args.outfile[0]}")

    #Dsp.do_fitting_coherence_complex()
    #data = DATASPACE.data_dict
    #t2_dict = DATASPACE.t2_dict


