# coding: utf-8

# """
#    
#   NAME:
#     camb4py
# 
#   PURPOSE:
#     Wrapper for fortran90 CAMB, which you should have installed
#     somewhere (the env variable CAMB_DIR should point there). Designed
#     particularly to handle multiple transfer function/power spectrum redshifts
#     and handle the parameter file accordingly.
# 
#   USE:
#     result=camb4py(paramfile='input_params.ini',
#                     outparamfile='output_params.ini',
#                     runcamb=True,[camb_keys...])
# 
#   INPUT:
#     
#   OPTIONAL INPUT:
#     paramfile - an input CAMB parameter file.  Defaults to
#                 'default_params.ini'
#     outparamfile - the paramfile that will be written out and passed
#                    to camb.  Identical to input if no keywords set.
#     output_root, etc. - inputs to CAMB.  See paramfile.
# 
#   KEYWORDS:
#     runcamb - Will call CAMB and run on output param file.  Must have
#               environment variable CAMB_DIR defined to location of CAMB executable.
# 
#   OUTPUT:
#     IF runcamb is set, outputs all of the CAMB files as set in param
#     file. Returns a dict with all of this information, and more.
# 
#   NOTES: 
#     Param file should have spaces before/after '=' signs (this is
#     required for accurate parsing)
# 
#     Does not currently support multiple scalar amplitudes, spectral
#     indices, etc.
# 
#     DOES support an array of redshifts for the transfer
#     function/matter power spectrum, and the code will handle naming
#     them and writing them to files (in fact this functionality is
#     largely why I wrote it!).  
# 
#     You shouldn't put a list of redshifts/files in the default param
#     file by hand (even if just two) - let the code do it for you.
# 
#     Note that if you set
#     transfer_redshift to an array with n_elements > 1, the code will automatically
#     set transfer_num_redshifts.  If not supplied with an equivalent set of file
#     names for the transfer function/matter power spectrum, will default to nameing them
#     "transfer_" and "matterpower_" with the z value appended.
#     
#     CAMB will only do up to 150 redshifts at a time.
# 
#   HISTORY:
#     6-2-16 - Written - MAD (Dartmouth)
# 
# """

import subprocess
import os
import numpy as np
def camb4py(paramfile='default_params.ini', 
            outparamfile='my_params.ini',runcamb=False,**kwargs):
    #Read in default parms
    defaults = read_defaults(paramfile)
    #Read in list of keywords
    params = [line.rstrip('\n') for 
             line in open('camb_keywords.txt')]
    
    #If given an array of redshifts, but transfer_redshift_num not
    #set or list of out files not given, fix that
    if 'transfer_redshift' in kwargs:
        if type(kwargs['transfer_redshift']) is list:
            if len(kwargs['transfer_redshift']) > 1:
                kwargs['transfer_num_redshifts'] = len(kwargs['transfer_redshift'])
                if 'transfer_filename' not in kwargs:
                    transnames=[]
                    for z in kwargs['transfer_redshift']:
                        transnames.append('transfer_' + 
                                          str(z) + '.dat')
                    kwargs['transfer_filename'] = transnames
                if 'transfer_matterpower' not in kwargs:
                    mpnames=[]
                    for z in kwargs['transfer_redshift']:
                        mpnames.append('matterpower_' + 
                                      str(z) + '.dat')
                    kwargs['transfer_matterpower'] = mpnames
                    
    
    #Initialize string describing what camb did
    done=[]
    
    #Setup dict of output files
    outfiles={}
    
    #Open output paramfile for writing
    f = open(outparamfile,'w')
    
    #Loop over default paramfile values
    for default in defaults:
        #Check if line has keyword - if not print it as is,
        #or set the kewyord and print
        result = get_keyword(default,params)
        if result == 'notkey':
            f.write(default+'\n')
        else:
            if result == 'output_root':
                if 'output_root' in kwargs:
                    outroot=kwargs['output_root']
                else:
                    outroot=get_value(default)
            if result == 'transfer_num_redshifts':
                if 'transfer_num_redshifts' in kwargs:
                    nz=kwargs['transfer_num_redshifts']
                else:
                    nz=float(get_value(default))
        
            #Some keywords require special treatment - check it
            special=check_special(result)
        
            #Variable not set in call, use default
            if result not in kwargs:
                f.write(default+'\n')
                
                if special == 3 or special == 5:
                    filestruc = get_outfile(result,get_value(default),outroot)
                    outfiles.update(filestruc)
                if special == 4:
                    did = what_did_camb_do(result,get_value(default))
                    done.append(did)
            #Variable set in call, use it
            if result in kwargs:
                tmp = kwargs[result]
                if special == 3 or special == 5:
                    filestruc = get_outfile(result,kwargs[result],outroot)
                    outfiles.update(filestruc)
                if special == 4:
                    did=what_did_camb_do(result,kwargs[result])
                    done.append(did)
                
                if special != 1 and special != 2 and special !=5:
                    f.write(result + ' = ' + str(kwargs[result]) + '\n')
                if special == 1:
                    f.write(result + '(1) = ' + str(kwargs[result]) + '\n')
                if special == 2 or special == 5:
                    for i in range(len(kwargs[result])):
                        f.write(result + '(' + str(i+1) + ') = ' + 
                                str(kwargs[result][i]) + '\n')
    f.close()
    
    #If nothing set to be done but runcamb=True, you must have
    #set it by mistake
    flag = 0
    if 'S' not in done and 'V' not in done and 'T' not in done and runcamb:
        print 'There''s nothing to do'
        flag = 1
        return -1
    
    if runcamb and not flag:
        highlfile = os.path.join(os.path.abspath(os.environ['CAMB_DIR']),
                                 'HighLExtrapTemplate_lenspotentialCls.dat')
        runcp = subprocess.call(['cp',highlfile,'.'])
        
        cmd = [os.path.join(os.path.abspath(os.environ['CAMB_DIR']),'camb'),outparamfile]
        print 'Calling CAMB with ' + cmd[0] +  ' ' + cmd[1]
        camboutput = subprocess.check_output(cmd)
        cambres = parse_cambout(camboutput,done,nz)

        cambout = read_camb_output(done,outfiles)
        cambres.update(cambout)
    
    return cambres        


# """  
#   NAME:
#     check_special
# 
#   PURPOSE:
#     Identifies keywords that require special action in the main code.    
# 
#   USE:
#     result=check_special(keyword)
# 
#   INPUT:
#     keyword - the (string) name of the keyword to check against the list
#     
#   OUTPUT:
#     Returns 0 (not special), 1 (needs additional string), 2 (related to redshift
#     list), 3 (contains name of output file), 4 (contains info on
#     what calculations CAMB did), or 5 (both redshift and output file)
# 
#   HISTORY:
#     5-25-16 - Written - MAD (Dartmouth)
# 
# """

def check_special(keyword):

    mod_list=['scalar_amp',
             'scalar_spectral_index',
             'scalar_nrun'
             'scalar_nrunrun',
             'tensor_spectral_index',
             'tensor_nrun',
             'initial_ratio',
             'tensor_amp']
    
    redshift_list=['transfer_redshift']

    output_list=['scalar_output_file', 
                   'vector_output_file', 
                   'tensor_output_file', 
                   'total_output_file', 
                   'lensed_output_file', 
                   'lensed_total_output_file', 
                   'lens_potential_output_file']

    do_list=['get_scalar_cls', 
             'get_vector_cls', 
             'get_tensor_cls', 
             'get_transfer', 
             'do_lensing']

    z_and_out_list=['transfer_filename', 
                    'transfer_matterpower']

    out=0

    if keyword in mod_list:
        out=1
    if keyword in redshift_list:
        out=2
    if keyword in output_list:
        out=3
    if keyword in do_list:
        out=4
    if keyword in z_and_out_list:
        out=5
        
    return out


# """"
# 
# NAME:
#     get_keyword
# 
#   PURPOSE:
#     Splits up a string from the parameter file to pull out a
#     potential keyword.  Identifies as a keyword or not.
# 
#   USE:
#     result=get_keyword(string,list)
# 
#   INPUT:
#     string - the string to parse
#     list - the list of possible keywords to check against
#     
#   OUTPUT:
#     Returns string - 'notkey' or the name of the keyword.
# 
#   HISTORY:
#     5-25-16 - Written - MAD (Dartmouth)
#     
# """

def get_keyword(line,list):
    if len(line) > 0:
        split1=line.split()
        split2=split1[0].split('(')
    
        if split2[0] in list:
            out=split2[0]
        else:
            out='notkey'
    else:
        out='notkey'
            
    return out


# """
# 
# NAME:
#     get_outfile
# 
#   PURPOSE:
#     builds a structure of output file names for use in read in parser
# 
#   USE:
#     result=get_outfile(key,val,root)
# 
#   INPUT:
#     key - the keyword of the output file
#     val - the name of the output file
#     root - the root of the output files
#     
#   OUTPUT:
#     Returns structure with all the possible out file names (not all
#     necesarily used, depends on your CAMB settings)
# 
#   HISTORY:
#     5-6-16 - Written - MAD (Dartmouth)
#     
# """

def get_outfile(key,val,root):
    output_list=['transfer_filename', 
                 'transfer_matterpower', 
                 'scalar_output_file', 
                 'vector_output_file', 
                 'tensor_output_file', 
                 'total_output_file', 
                 'lensed_output_file', 
                 'lensed_total_output_file', 
                 'lens_potential_output_file']

    x = [s for s in output_list if key in s]
#    t = Table([[root+'_'+val]],names=(x[0],),dtype=('str',))

    if type(val) is str:
        t = {x[0]:root+'_'+val}
    else:
        out=[]
        for v in val:
            out.append(root+'_'+v)
        t = {x[0]:out}
    return t


# """
# 
#   NAME:
#     get_value
# 
#   PURPOSE:
#     Some cases need to save the value of a parameter, even if not set
#     in call (i.e. need the default).  This gets it.
# 
#   USE:
#     result=get_value(line)
# 
#   INPUT:
#     line - the parameter file line to read
#     
#   OUTPUT:
#     Returns value of parameter as a string
# 
#   HISTORY:
#     5-26-16 - Written - MAD (Dartmouth)
#     
# """

def get_value(line):
    split1=line.split()
    return split1[-1]


def parse_cambout(inp,done,numz):
    t = {'raw_output':inp}
    return t


# """
# 
#   NAME:
#     read_camb_output
# 
#   PURPOSE:
#     Read in files that CAMB writes out
# 
#   USE:
#     result=read_camb_output(used,file_struct)
# 
#   INPUT:
#     used - string from what_did_camb_do listing what was calculated
#     file_struct - structure containing possible file names
#     
#   OUTPUT:
#     Returns a structure with all the CAMB data
# 
#   HISTORY:
#     5-26-16 - Written - MAD (Dartmouth)
# 
# """

def read_camb_output(used,files):
    out = dict()
    if 'S' in list(used):
        out['scalar_cls'] = read_matrix(files['scalar_output_file'])
    if 'V' in list(used):
        out['vector_cls'] = read_matrix(files['vector_output_file'])
    if 'T' in list(used):
        out['tensor_cls'] = read_matrix(files['tensor_output_file'])
    if 'S' in list(used) and 'T' in list(used):
        out['total_cls'] = read_matrix(files['total_output_file'])
    if 'l' in list(used):
        out['lensed_cls'] = read_matrix(files['lensed_output_file'])
    if 'S' in list(used) and 'T' in list(used) and 'l' in list(used):
        out['lensed_total_cls'] = read_matrix(files['lensed_total_output_file'])
    if 'tr' in list(used):
        if type(files['transfer_filename']) is list:
            i=1
            for fname in files['transfer_filename']:
                tagname = 'transfer_func_'+str(i)
                out[tagname] = read_matrix(fname)
                i=i+1
        else:
            out['transfer_func'] = read_matrix(files['transfer_filename'])
        
        if (type(files['transfer_matterpower'])) is list:
            i=1
            for fname in files['transfer_matterpower']:
                tagname = 'matter_power_'+str(i)
                out[tagname] = read_matrix(fname)
                i=i+1
        else:
            out['transfer_matterpower'] = read_matrix(files['transfer_matterpower'])
    return out


# """
# 
#   NAME:
#     read_matrix
#   
#   PURPOSE:
#     read in a matrix from a text file into a numpy array
# 
#   USE:
#     matrix = read_matrix('filename.txt')           
# 
#   INPUT:
#     file - string name of file to read
# 
#   RETURNS:
#     matrix - the array read from the file
# 
#   HISTORY:
#     6-1-16 - Written - MAD (Dartmouth)
#  
#  """

def read_matrix(filename):
    strarr=[line.rstrip('\n') for line in open(filename)]
    width=len(strarr[0].split())
    length=len(strarr)
    outar=np.zeros((length,width))
    for i in range(length):
        outar[i]=np.fromstring(strarr[i], 
                               dtype=(float), sep=' ')
        
    return outar


# """
# 
#   NAME:
#     read_defaults
# 
#   PURPOSE:
#     Read in initial parameter file, split into string array
# 
#   USE:
#     result=read_defaults(paramfile)
# 
#   INPUT:
#     file - string name of parameter file
#     
#   OUTPUT:
#     Returns string array, each element a line of the parameter file
# 
#   HISTORY:
#     5-26-16 - Written - MAD (Darmtouth)
#     
#  """

def read_defaults(file):
    lines=[line.rstrip('\n') for 
           line in open(file)]
    return lines


# """
# 
#   NAME:
#     what_did_camb_do
# 
#   PURPOSE:
#     Determines what your CAMB run did, based on 'T' and 'F'
#     statements for things like get_scalar_cls, get_tensor_cls,
#     do_lensing, etc.
# 
#   USE:
#     result=what_did_camb_do(key,val)
# 
#   INPUT:
#     key - the keyword you're checking
#     val - the value of the keyword (T or F)
#     
#   OUTPUT:
#     Returns 5 element string (scalar, vector, tensor, transfer,
#     lensing), set to N if not done, and S/V/T/tr/l if done.
# 
#   HISTORY:
#     5-26-16 - Written - MAD (Dartmouth)
#     
#  """

def what_did_camb_do(key,val):
    if key == 'get_scalar_cls':
        out='S'
    if key == 'get_vector_cls':
        out='V'
    if key == 'get_tensor_cls':
        out='T'
    if key == 'get_transfer':
        out='tr'
    if key == 'do_lensing':
        out='l'
    if val == 'F':
        out='N'
    
    return out

