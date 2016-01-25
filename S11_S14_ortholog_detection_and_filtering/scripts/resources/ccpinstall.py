#!/usr/bin/python
"""
ccpinstall.py

Modifies the paths to RAxML in the file concaterpillar.py

Usage: python ccpinstall.py [options]

  Options:
    -h,                          print this message
    -r,                          path to raxml
    
  The program will look for RAxML using the which command.
  If paths are provided as arguments, they will overwrite the paths
  found using which
    
"""

import filetools
import fileinput
import sys
import re
import os



def getoptdict(arglist):
  """Gets arguments from the list arglist
  
  Return: 
    a dictionary mapping flags to values
  
  """
  singleopts = ('-h', )
  multiopts = ('-r',)
  x = 0
  optdict = {}
  while x < len(arglist):
    arg = arglist[x]
    if arg in multiopts:
      optdict[arg] = arglist[x+1]
      x += 1
    elif arg in singleopts:
      optdict[arg] = None
    x += 1
  return optdict
  


def selectpath(path, program):
  """Prompt the user to ask where the executable program is located . 
  
  Parameters:
    path - the path the os for program, or none if none was found
    program - the name of the executable we want to find
  
  Return: 
    string
  
  """
  
  while True:
    if path:
      print "%s was found at %s , is this the correct path? y/n" % (program, path)
      answer = raw_input().strip().lower()
      while not (answer == 'n' or answer == 'y'):
        print "please answer yes or no, y/n"
        answer = raw_input().strip().lower()
      if answer == 'y':
        return path
      else:
        path = None
        
    else:
      print "No path for %s is specified, please enter a path, q to quit" % program
      path = raw_input().strip()
      if path.lower() == 'q':
        print "cancelling installation..."
        sys.exit()
      else:
        verified = os.access(path, os.X_OK)
        if verified:
          print "path verified..."
          return path
        else:
          print "could not find %s at %s, please specify a different path" % (program, path)
          path = None
      


if __name__ == '__main__':
  
  optdict = getoptdict(sys.argv[1:])
  if '-h' in optdict:
    sys.stderr.write(__doc__)
    sys.exit(0)
    
  
  print 'Welcome to the Concaterpillar installation script!'
  print 'This script will help you find where RAxML is located.'
  print '\n\
                     o  o\
\n __  __  __  __  __  _\\/\
\n/  \\/  \\/  \\/  \\/  \\/ .\\\
\n\\__/\\__/\\__/\\__/\\__/\\_~/\
\n     ""  ""  ""  ""\n' 
  
  file = 'concaterpillar.py'
  
  # look for paths
  try:
    raxmlpath = filetools.getpath('raxml')
  except filetools.RunError:
    try:
      raxmlpath = filetools.getpath('raxmlHPC')
    except filetools.RunError:
      raxmlpath = None
  
  # get paths from arguments
  if '-r' in optdict:
    raxmlpath = optdict['-r']
    verified = os.access(raxmlpath, os.X_OK)
    if not verified:
      print "the given path %s was not correct" % raxmlpath
      raxmlpath = None
      
  # ask the user if the paths are correct
  raxmlpath = selectpath(raxmlpath, 'raxml')
  
  raxmlRE = re.compile(r"#?[ \t]*raxmlpath[ \t]*=")
  
  # modify concaterpillar.py
  for line in fileinput.input(file, inplace=1):
    raxmlmatch = raxmlRE.match(line)
    if raxmlmatch:
      sys.stdout.write("raxmlpath = '%s'\n" % raxmlpath)
    else:
      sys.stdout.write(line)
  fileinput.close()
  print "installation complete."
  
