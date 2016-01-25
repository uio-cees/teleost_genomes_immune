"""

Classes and functions for manipulation of files, including program execution.

"""

import os, sys, types, subprocess
from seqtools import *

def getlines(filename):
  """Open and read lines from a file.
  
  Parameters:
    filename -- name of the file from which lines will be read.
    
  Return:
    List containing lines from the file, trailing whitespace removed.
    
  """
  
  file = open(filename)
  lines = [line.rstrip() for line in file]
  
  file.close()

  return lines
  
def getpath(executable):
  """Find the path for an executable (must be in shell's PATH env. variable).

  If the executable is not found, raise RunError.
  
  Parameters:
    executable -- name of the executable to look for.
    
  Return:
    Full path of the executable (as a string).

  """

  path = ''
  args = ('which', executable)
  process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
  
  out, err = process.communicate()
  if err:
   sys.stderr.write('Should not see any stderr output:\n')
   sys.stderr.write(err)
  
  if out:
    path = out.strip()

  else:
    raise RunError('No path to %s found!' % executable)

  
  if path == '' or path[0] != '/':
    raise RunError('Executable %s not found in PATH.' % executable)
  
  return path


def removefiles(files):
  """Delete all files in a list or tuple.
  
  Parameters:
    files -- a list or tuple containing names of files to be deleted.

  """
  
  for f in files:
    os.remove(f)
    
    
def renamefiles(filedict, reverse=False):
  """Rename files listed in filesdict.

  The keys in filesdict are first sorted, then files are renamed in order from
  "smallest" to "largest" (alphabetical order), unless reverse is set to True.

  Parameters:
    filedict -- a dictionary with old filenames as keys, new names as values.
    reverse -- boolean that determines order in which files should be renamed.

"""
  keys = filedict.keys()
  keys.sort()
  
  if reverse:
    keys.reverse()
  
  for k in keys:
    try:
      os.rename(k, filedict[k])
    except OSError:
      sys.stderr.write('Failed to rename file %s. Files to rename:\n' % k)
      for k in filedict:
        sys.stderr.write('%s: %s\n' % (k, filedict[k]))
      raise
    
def runbin(binary, cmds, nice = 0, debug = False):
  """Run an executable with a list of input commands.  

  Designed for executables with necessary user interaction, such as proml, but 
  will work with other programs as well if an empty tuple (or list) is supplied 
  for the cmds parameter.

  Paramters:
    binary -- a string containing the name of the executable to run
    cmds -- a list (or tuple) of commands to be passed to the executable
    nice -- nice level to use with the binary (must be an int) [default = 0]
    debug -- if True, output and error streams from the executed program
             will be written to standard error

  """

  args = ('nice', '-n', '%d' % nice)
  args += tuple(binary.split(' '))  
  cmdstr = '\n'.join(cmds)

  if debug:
    sys.stderr.write('Will run command: %s\n' %  ' '.join(args))

  process = subprocess.Popen(args, stderr=subprocess.PIPE, stdout=subprocess.PIPE, stdin=subprocess.PIPE)  
  
  out, err = process.communicate(cmdstr)

  if err:
    if debug:
      sys.stderr.write(err)
  if out:
    if debug:
      sys.stderr.write(out)
  
def getraxmlversion(raxmlpath):
  args = (raxmlpath, '-v')
  version = -1
  versionre = re.compile(r'This is RAxML version (\d+\.\d+)')
  process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
  
  out, err = process.communicate()
  
  if err:
    sys.stderr.write('Strange, nothing should be sent to stderr:\n')
    sys.stderr.write(err)
    
  if out:
    versionstr = out.strip()
    result = versionre.search(versionstr)
    if result:
      try:
        version = float(result.group(1))
      except IndexError:
        version = -1
        
  return version
  
def runraxml(raxmlpath, alignname, suffix, model, version=7.0, debug=False):
  """Run RAxML with the specified parameters.
  
  Parameters:
    raxmlpath -- the path to the RAxML executable
    alignname -- The file for RAxML to build a tree with
    model -- the model to use, alid models are 'DAYHOFF', 'DCMUT', 'JTT',
            'MTREV', 'WAG', 'RTREV', 'CPREV', 'VT', 'BLOSUM62', 'MTMAM', 'GTR')
  """
  
  if version == 7.0:
    modelclass = 'MIX'
  elif version == 7.2:
    modelclass = 'CAT'
  else:
    raise RunError('RAxML version %s is not supported (yet).' % str(version))  
  if model == 'GTR':
    binstr = '%s -n %s -m GTR%s -s %s' % (raxmlpath, suffix, modelclass, alignname)
  else:  
    binstr = '%s -n %s -m PROT%s%s -s %s' % (raxmlpath, suffix, modelclass, model, alignname)
  runbin(binstr, (), debug=debug)

def runphyml(phymlpath, aligns = None, alignname = 'infile', datatype = 'aa', \
             format = 's', nsets = 1, model = 'JTT', kappa = 'e', invar = 0.0, \
             ncats = 4, alpha = 'e', trees = None, treefile = None, \
             opttop = True, optbl = True, debug=False):
  """Run PHYML with the specified parameters.
  
  Parameters:
    phymlpath -- full path of the PHYML binary to be used.
    aligns -- a list of Alignment objects with which PHYML will be run.
              If None, nsets must be specified, and alignname is assumed
              to already have the sets in it.  As well, format must be 
              specified. [default = None]
    alignname -- name to give to the sequence file that PHYML will analyse.
                 [default = 'infile']
    datatype -- indicates amino acid ('aa') or nucleotide ('nt'). 
                [default = 'aa']
    format -- indicates whether input alignment file is interleaved ('i') or 
              sequential ('s') [default = 's']
    nsets -- number of alignments in input file (alignname).  If aligns is not
             None, this value is ignored, and the number of Alignment objects
             in aligns is used instead. [default = 1]
    model -- indicates which model to use (JC69|K2P|F81|HKY|F84|TN93|GTR (DNA)
             JTT|MtREV|Dayhoff|WAG (Amino Acids) ). [default = 'JTT']
    kappa -- transition/transversion ratio, ignored for amino acid.  For ML 
             estimate, use 'e'; otherwise, should be a float. [default = 'e']
    invar -- fraction of invariable sites.  For ML estimate, use 'e'; otherwise,
             should be a float.  [default = 0.0]
    ncats -- number of rate categories to use with gamma distribution.  For 
             uniform distribution of rates across sites, use 1. [default = 4]
    alpha -- shape parameter for the gamma distribution model of variation of
             rates across sites. For ML estimate, use 'e'; otherwise, should be
             a float. [default = 'e']
    trees -- a list (or tuple) containing the tree strings to be analysed.  If
             trees are specified, but no treefile is given, the trees will be
             written to the file 'intree'. [default = None]
    treefile -- name of the file to which input trees will be written.  If
                no trees are specified, this file is assumed to already exist
                and contain input treefiles.  If None, a BIONJ tree will
                be used as a starting tree.  [default = None]
    opttop -- if True, topology will be optimised. [default = True]
    optbl -- if True, branch lengths will be optimised. [default = True]
    
  """
  
  if format == 's':
    interleaved = False
  elif format == 'i':
    interleaved = True
  else:
    raise RunError('Invalid format parameter!')
  
  if aligns:
    writephy(aligns[0], alignname, interleaved=interleaved)
    for align in aligns[1:]:
      writephy(align, alignname, append=True, interleaved=interleaved)

    nsets = len(aligns)
    
  ntmodels = ('JC69', 'K2P', 'F81', 'HKY', 'F84', 'TN93', 'GTR')
  aamodels = ('JTT', 'MtRev', 'Dayhoff', 'WAG')
  
  for param in (kappa, invar, alpha):
    if not (param == 'e' or isinstance(param, types.FloatType)):
      raise RunError('Unsupported value for parameter: %s. ' % str(param)+ \
                     'Must be a float or \'e\'')
    
  if datatype == 'aa':
    if model not in aamodels:
      raise RunError('Model %s not supported for datatype aa!' % model)
    datatype = 1 # acceptable format for PHYML
    kappa = ''
    
  elif datatype == 'nt':
    if model not in ntmodels:    
      raise RunError('Model %s not supported for datatype nt!' % model)
    datatype = 0 # acceptable format for PHYML
    
  else:
    raise RunError('Unsupported datatype: %s' % str(datatype))
  
  if ncats < 1:
    raise RunError('Must have at least 1 rate category!')
  
  if trees:
    if not treefile:
      treefile = 'intree'
    tf = open(treefile, 'w')
    for tree in trees:
      tf.write('%s\n' % tree)
    tf.close()
    
  if not treefile:
    treefile = 'BIONJ'
    
  if opttop:
    opttop = 'y'
  else:
    opttop = 'n'

  if optbl:
    optbl = 'y'
  else:
    optbl = 'n'
    
  binstr = '%s %s %d %s %d 0 %s %s %f %d %s %s %s %s' \
         % (phymlpath, alignname, datatype, format, nsets, model, str(kappa), \
            invar, ncats, str(alpha), treefile, opttop, optbl) 
  
  runbin(binstr, (), debug = debug)
  
  
class RunError(Exception):
  """
  Exception designed for errors related to running executables.
  
  """

  def __init__(self, value = ''):
    """Constructor for RunError exceptions."""
    self.value = value
      
  def __str__(self):
    """Returns a string representation of the exception."""
    return repr(self.value)
    
  
