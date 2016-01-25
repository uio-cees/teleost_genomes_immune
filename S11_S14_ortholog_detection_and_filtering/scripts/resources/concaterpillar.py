#!/usr/bin/python
"""
concaterpillar.py

Assessment of congruence for protein alignment concatenation.

Usage: python concaterpillar.py [options]

  Options:
    -h, --help                print this message
    -t, --topology            begin topological congruence test
    -i, --interrupt           continue interrupted topological congruence test
    -b, --blength             begin branch length compatibility assessment
    -c ncpu, --cpus=ncpu      use specified number of processors (ncpu must be 
                              integer greater than 0, default = 1)
    -p cutoff, --pval=cutoff  use specified p-value as a cutoff for tests 
                              (must be a number between 0 and 1, default = 0.05)
    -m model, --model=model   The substitution model to be used by RAxML
                              (Valid models are: 'DAYHOFF', 'DCMUT', 'JTT', 'MTREV', 
                              'WAG', 'RTREV', 'CPREV', 'VT', 'BLOSUM62', 'MTMAM', 
                              'GTR', default = WAG)
    -s nrounds --save=nrounds save a backup once it finishes the first run of phyml
                              and puzzle, and then every nrounds round of toptest.
                              the files will be saved in the file savepoint.tbz,
                              which will be overwritten each time.
                              (must be integer greater than 0)
    -d, --debug               print all output from calls to other programs to stderr
    -n, --noask               cleanup leftover RAxML_* and boot* files from previous
                              without asking
    
  Notes:
    At least one option of -t, -i, and -b must be specified.  More
    than one may be used.  Note also that if you run with one of these switches
    at a time, -a must be run before -r, and -r before -t.  -i should only be
    used if trees have already been build with PHYML (-a) and evaluated with
    TREE-PUZZLE (-r), and the topology test (-t) has begun, but has been
    interrupted prior to completion for some reason.  In theory, -b can be used
    independently of all other options, in which case branch length 
    compatibility of all alignments will be evaluated.  Otherwise, if congruence
    test has been completed, the branch length compatibility test will be
    performed only on the alignments in the largest cluster.

"""


import cPickle, getopt, glob, shutil, sys
from ccptools import *

raxmlpath = './raxmlHPC'
# if you want to time the run
timing = True
# if you want to limit the memory usage
memorylimit = False

memorysize = 1500000000 # size in bytes program will take at most

version = '1.7.2'

if __name__ == '__main__':
  
  # see if we want to limit resource use
  if memorylimit:
    import resource
    resource.setrlimit(resource.RLIMIT_AS, (memorysize, memorysize))
  
  # see if mpi is used
  mpipresent = False
  try:
    import mpi
    if mpi.size > 1:
      mpipresent = True
      if mpi.rank > 0:
        processdispatcher()
  except ImportError:
    pass

  if not mpipresent or mpi.rank == 0:
    if os.access('_mpiusagemap_', os.R_OK):
      print "file _mpiusagemap_ detected, will be used during the process"  
    
    if timing:
      timing = [time.time()]
    
    # get and process command line options
    optdict = getoptdict(sys.argv[1:])

    if '-h' in optdict:
      usage(0)
      
    coreopts = ('-t', '-i', '-b')#('-a', '-r', '-t', '-i', '-b')
    noopt = True
    
    for opt in coreopts:
      if opt in optdict:
        noopt = False
        break
    
    if noopt:
      usage(1)
    
    # default substituion model
    model = 'WAG'
    
    try:
      if '-c' in optdict:
        nprocs = int(optdict['-c'])
      else:
        nprocs = 1
        
      if nprocs < 1:
        raise ValueError
        
      if '-p' in optdict:
        cutoff = float(optdict['-p'])
      else:
        cutoff = 0.05
      
      if cutoff > 1 or cutoff < 0:
        raise ValueError
        
      if '-s' in optdict:
        saverounds = int(optdict['-s'])
      else:
        saverounds = None
      
      if '-d' in optdict:
        debug = True
      else:
        debug = False
        
      if '-n' in optdict:
        ask = False
      else:
        ask = True
      
      if saverounds:
        if saverounds < 1:
          raise ValueError
      
      goodmodels = ('DAYHOFF', 'DCMUT', 'JTT', 'MTREV', 'WAG', 'RTREV',
                    'CPREV', 'VT', 'BLOSUM62', 'MTMAM', 'GTR')
      if '-m' in optdict:
        model = optdict['-m'].upper()
      if not model in goodmodels:
        raise ValueError
    
    except ValueError:
      usage(1)
            
    if '-i' in optdict and '-t' in optdict:
        usage(1)
    
    cleanup(ask=ask)
    seqlist = getseqlist()
    if raxmlversion != 7.0 and raxmlversion != 7.2 and raxmlversion != 7.3:
      sys.stderr.write('RAxML version %f is not supported (yet).\n\n' % raxmlversion)
      sys.stderr.write('If you are using the latest version of Concaterpillar, please email\n')
      sys.stderr.write('jessica.w.leigh@gmail.com to request support of this version.\n')
      sys.exit(1)

    # resume interrupted congruence test  
    if '-i' in optdict:  
      try:
        catsets = recover(model, nprocs, rebuild=True, debug=debug) 
      except IOError:
        sys.stderr.write('Could not recover concatenated sets.  Sorry.\n')
        shutdown(0)
      
      likes = getalllikes(len(catsets))
      
      toptest(catsets, likes, cutoff, model, nprocs, saverounds, debug=debug)
    
    # start a new congruence test
    # build trees with RAxML
    if '-t' in optdict:
      try:
        sets = preparesets(seqlist)
      except ConcaterpillarError, ce:
        sys.stderr.write('%s\n' % ce.value)
        sys.stderr.write('An exception was received when parsing alignments.\n')
        sys.stderr.write('Normally, this means that a pair of alignments has fewer than 4 taxa in\n')
        sys.stderr.write('common. Keep in mind that taxa must have the same names in all alignments\n')
        sys.stderr.write('or they will be parsed as different taxa!\n')
        
        sys.exit(1)
        
      buildmultitrees(sets, model, nprocs, debug=debug)
      
      if timing:
        timing.append(time.time() - timing[0])
        print "Time to finish building trees:", timing[-1]

      # evaluate likelihoods with RAxML
      naligns = len(seqlist)
      runmultipuzzle(model, naligns, nprocs = nprocs, debug=debug)
      
      if timing:
        timing.append(time.time() - timing[0])
        print "Time to finish reanalysis:", timing[-1]
    
      # start new congruence test
      catsets = [[set] for set in seqlist]
      print '\nStarting compatibility analysis.'
      resfile = open('results.ccp', 'w')
      resfile.write('Concaterpillar results:\n\n')
      resfile.close()
      
      # create savepoint if nessecary
      if saverounds:
        createsavepoint()
      
      likes = getalllikes(len(seqlist))
      try:
        toptest(catsets, likes, cutoff, model, nprocs, saverounds, debug=debug)
      except IOError: 
        sys.stderr.write('Invalid starting conditions for congruence test.\n')
        raise
      
      if timing:
        timing.append(time.time() - timing[0])
        print "Time to finish congruence test:", timing[-1]
      
    # test branch length compatibility
    if '-b' in optdict:
      allcatsetslist = brlensetup(model, debug=debug)
      if len(allcatsetslist) == 1 and 'BLtest-set000' not in os.listdir('.'):
        bltestall = True
      else:
        bltestall = False
      
      
      setcount = 0
      
      for catsets in allcatsetslist:
        catsets = [[set] for set in catsets]
        if bltestall:
          setdir = ''
        else:
          setdir = 'BLtest-set%03d/' % setcount
          os.chdir(setdir)
          print 'Starting branch-length congruence test on set%03d.' % setcount
        
        for i in range(len(catsets)):
          writephy(getseqs(catsets[i][0]), 'set%03d.seq' % i)
          
        fixtreefile('global.tre', roottree = True, clean = True)        
        brlentest(catsets, setdir, cutoff, model, nprocs, debug=True)
    
        if not bltestall:
          os.chdir('..')
          setcount += 1
        
      # if MPI is running, tell processes to terminate
    if timing:  
      timing.append(time.time() - timing[0])
      print "total time taken:", timing[-1]
      
    shutdown(0)
  
