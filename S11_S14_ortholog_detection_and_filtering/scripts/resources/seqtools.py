"""
A reusable module of sequence manipulation tools.  Originally included in
Concaterpillar.

"""

import os, random, re, string, sys, time

version = '1.0b'

def bootstrap(align, nsites = -1):
  """Return a bootstrap alignment containing the specified number of sites.
  
  Sites are resampled with replacement.  Normally, the same number of sites
  sampled should be the same as in the original alignment, but if nsites is 
  specified, the number of sites can be controlled.
  
  Parameters:
    align -- an Alignment object from which sites will be resampled
    nsites -- number of sites to sample from seqs (default = len(seqs[0][1],
              the number of sites in the seqs alignment)
  
  Return:
    The bootstrap alignment, a new Alignment object.
  """
  if nsites == -1:
    nsites = align.length
    
  if nsites <= 0:
    raise SeqError, 'Invalid number of sites to bootstrap: %d' % nsites
  
  while True:
    emptyseqs = False
    bootseqs = []
    siteids = []
  
    for i in range(nsites):
      siteids.append(random.randint(0, align.length - 1))
    
    siteids.sort()
    
    for seq in align:
      bootseqs.append(seq.resample(siteids))
      bootseqs[-1].setseqtype()
      if bootseqs[-1].type == 'empty':
        emptyseqs = True
        break
      
    if emptyseqs:
      continue
    
    return Alignment(bootseqs)

def repartition(align, partlengths, random=True):
  """(Randomly) partition sites in alignments of given length.
  
  The lengths should sum to the length of the alignment, but no error will be
  raised if they sum to less than the alignment length.  SeqError will be raised
  if the sum exceeds alignment length.

  Parameters:
    align -- Alignment object to repartition
    partlengths -- list or tuple containing lengths of partition.
    random -- determines whether partitioning should be done randomly, or
              whether alignment should just be split (default = random)
  
  Return:
    List of (randomly) partitioned Alignment objects.
  
  """
  
  totallength = sum(partlengths)
  if totallength > align.length:
    raise SeqError, 'Sum of partition lengths exceeds alignment length: %d' \
          % totallength
  
  # produce a list of alignment site indices in random order
  if random:
    siteids = random.sample(xrange(align.length), align.length)
  else:
    siteids = range(align.length)
  
  partitions = []
  start = 0
  for length in partlengths:
    end = start + length
    partseqs = []
    
    for seq in align:
      partseqs.append(seq.resample(siteids[start:end]))
    
    partitions.append(Alignment(partseqs))
    start = end
    
  return partitions

def jackknife(align, nsites = -1):
  """Produce a resampled (without replacement) Alignment object.
  
  Parameters:
    align -- Alignment object from which sites are sampled.
    nsites -- number of sites to sample from align (default = 50%)
    
  Return:
    Alignment object with resampled sites.

  """
  
  if nsites == -1:
    nsites = align.length / 2
    
  if nsites <= 0:
    raise SeqError, 'Invalid number of sites to jackknife: %d' % nsites
  
  siteids = random.sample(xrange(align.length), nsites)
  jackseqs = []
  
  for seq in align:
    jackseqs.append(seq.resample(siteids))
    
  return Alignment(jackseqs)

def jackstrap(align, nsites):
  
  if nsites <= 0 or nsites > align.length:
    raise SeqError, 'Invalid number of sites to jackstrap: %d' % nsites
  
  jacksiteids = random.sample(xrange(align.length), nsites)
  siteids = []
  
  for i in range(nsites):
    siteids.append(jacksiteids[random.randint(0, nsites - 1)])
    
  jstrapseqs = []
  
  for seq in align:
    jstrapseqs.append(seq.resample(siteids))
  
  return Alignment(jstrapseqs)
  
def catseq(alignlist, logfilename=None):
  """Concatenate a list of alignments.

  Sequences need not be in the same order in input alignments, but must have
  EXACTLY the same name in all alignments.
  
  Parameters:
    alignlist -- a list of Alignment objects
    logfilename -- filename to which various statistics about the alignment are
                   written (number of missing genes/sequence, etc.)
                   
  Return:
    Alignment containing concatenated alignments.
    
  """
  charcount = 0
  seqdict = {}
  seqcounter = {}
  charsets = []
  typelist = []
  
  for alignment in alignlist:
  
    charsets.append((alignment.name, alignment.seqlist[0].type, charcount + 1, charcount + len(alignment)))
    charcount += len(alignment)
    
    for seq in alignment:
      newseq = Sequence(seq.name, seq.sequence, seq.code, seq.chksum, seq.type)
      if len(seq.sequence) < alignment.length:
        newseq.pad(alignment.length)
      try:
        seqdict[seq.name] += newseq
        seqcounter[seq.name] += 1

      except KeyError:
        newseq.fpad(charcount)
        seqdict[seq.name] = newseq
        seqcounter[seq.name] = 1
    

    for name in seqdict:
      if len(seqdict[name].sequence) < charcount:
        seqdict[name].pad(charcount)
    
  alignment = Alignment(seqdict.values(), charsets=charsets)

  # log output if asked to do so
  if logfilename:
    logfile = open(logfilename, 'w')
    totalgenes = len(alignlist)
    totalchars = len(alignment)
    logfile.write('SUMMARY:\n')
    logfile.write('Total genes: %d\n' % totalgenes)
    logfile.write('%-20s%10s%10s%10s%10s\n' % \
                  ('Taxon', '# Genes', '# Missing', '% Missing', '% Gaps'))
      
    seqnames = seqcounter.keys()
    seqnames = [(name.lower(), name) for name in seqnames]
    seqnames.sort()
    seqnames = [name[1] for name in seqnames]
    for seqname in seqnames:
      ngenes = seqcounter[seqname]
      nchar = seqdict[seqname].countnongapchars()
      gapchars = (100.0 * (totalchars - nchar)) / totalchars
      logfile.write('%-20s%10d%10d%10.2f%10.2f\n' % \
                    (seqname[:20], ngenes, totalgenes - ngenes, \
                     (100.0 * (totalgenes - ngenes)) / totalgenes, gapchars))
      
    logfile.write('\nDETAILS:\n')
    logfile.write('%-20s ' % 'Taxon')
    for a in alignlist:
      logfile.write('%8s ' % a.name.rsplit('.', 1)[0][:8])
    logfile.write('\n')
    
    for name in seqnames:
      logfile.write('%-20s ' % name[:20])
      for a in alignlist:
        if name in a.taxlist:
          logfile.write('%8d ' % 1)
        else:
          logfile.write('%8d ' % 0)
      logfile.write('\n')
    
    logfile.close()
      
  return alignment


def fixnames(align):
  """Removes genbank-type codes from the front of names of Sequence objects.

  For example, ZP_0299987Geobacter_metallired will become Geobacter_metallired.
  This function will only remove codes that contain 1 to 3 letters, an optional
  underscore, then at least 4 digits.  If the name of the sequence actually
  starts with a number, then this will be removed.  

  The code is removed from sequence names, but added to the Sequence object's
  code field.  
  """
  
  regex = re.compile(r'(^[a-zA-Z]{1,3}_?\d{4,})') # re to match Genbank code
  goodseqlist = []

  for seq in align:
    # because re must be at the beginning, there won't be more than 3 arguments
    # however, there will be an empty first argument because the re splits
    # nothing from the rest of the string at the beginning
    spl = regex.split(seq.name)
    if len(spl) < 3:
      sys.stdout.write('Alignment %s: no change made to name %s: '\
                       % (align.name, seq.name))
      sys.stdout.write('please change this name manually.\n')

    else:
      empty, code, name = spl
      seq.code = code
      seq.name = name
      
    goodseqlist.append(seq)

  align.seqlist = goodseqlist
 
def getaligns(seqfiles, seqtype = 'auto', verbose = False):
  """Read an alignment object from each of the given alignment files.
  
  Parameters:
    seqfiles -- list of alignment filenames
    seqtype -- format of alignment files (default = automatic)
    verbose -- print more garbage to the screen if True
    
  Return:
    List of Alignment objects read from the files in seqfiles.
    
  """
  
  aligns = []
  for sf in seqfiles:
    aligns.append(getseqs(sf))
    aligns[-1].name = sf
    
  return aligns
    
def getseqs(sfname, seqtype = 'auto', verbose = False):
  """Read sequences from the given alignment file.
  
  Parameters:
    sfname -- name of the alignment file.
    seqtype -- format of the alignment file (default = automatic)
    verbose -- print more garbage to the screen if True
  
  Return:
    Alignment object with sequences from input file.
    
  """
  seqs = []  
  seqfile = open(sfname, 'U')
  
  if seqtype == 'auto':
    line = seqfile.readline()  
    while line != '':
      line = line.strip()
      if line == '':
        line = seqfile.readline()
      else:
        break

    line = line.lower()    
    
    if line == '#nexus':
      seqtype = 'nexus'
    elif line == '#mega':
      seqtype = 'mega'
    elif '#=au' in line.split():
      # old HMMER format
      seqtype = 'selex' 
    # comments in fasta-related files
    elif line[0] == '#' or line[0] == ';': 
      while line[0] == '#' or line[0] == ';' or len(line.strip()) == 0:
        line = seqfile.readline()
      if line[0] == '>':
        seqtype = 'fasta'
    elif line == '{':
      seqtype = 'gde'
    elif line == 'pileup' or line == '!!aa_multiple_alignment' or line == '!!na_multiple_alignment':
      seqtype = 'msf'
    elif line[0] == '>':
      seqtype = 'fasta'
    else:
      wordlist = line.split()
      if 'msf:' in wordlist:
        seqtype = 'msf'
      elif 'clustal' in wordlist:
        seqtype = 'clustal'
      else:
        try:
          ntaxa, nchar = int(wordlist[0]), int(wordlist[1])
        except (ValueError, IndexError):
          raise SeqError('This file is of unknown format!')
        seqtype = 'phylip'
        
    seqfile.seek(0, 0)
      
  if verbose:
    print 'File %s is in %s format.' % (sfname, seqtype)
  
  if seqtype == 'nexus':
    for seq in readnex(seqfile):
      # as long as seq isn't None, save it
      if seq: 
        seqs.append(seq)
  elif seqtype == 'mega':
    for seq in readmega(seqfile):
      if seq:
        seqs.append(seq)
  elif seqtype == 'selex':
    for seq in readselex(seqfile):
      if seq:
        seqs.append(seq)
  elif seqtype == 'gde':
    for seq in readgde(seqfile):
      if seq:
        seqs.append(seq)
  # note that this reads pir, must, or fasta files
  elif seqtype == 'fasta': 
    for seq in readfasta(seqfile): 
      if seq:
        seqs.append(seq)
  elif seqtype == 'phylip':
    for seq in readphy(seqfile):
      if seq:
        seqs.append(seq)
  elif seqtype == 'msf':
    for seq in readmsf(seqfile):
      if seq:
        seqs.append(seq)
  elif seqtype == 'clustal':
    for seq in readclustal(seqfile):
      if seq:
        seqs.append(seq)
  else:
    raise SeqError('This file is of unknown format!')
      
  seqfile.close()
  
  return Alignment(seqs, sfname)

def getsequnion(align1, align2):
  """Produce Alignments containing only sequences of taxa in both alignments.

  Parameters:
    align1, align2 -- Alignment objects

  Return:
    a tuple containing copies of the input Alignments, minus any Sequence
    objects with names not in the taxlists of both Alignments
  """

  seqlist1 = []
  seqlist2 = []
  for seq in align1.seqlist:
    if seq.name in align2.taxlist:
      seqlist1.append(seq)
  for seq in align2.seqlist:
    if seq.name in align1.taxlist:
      seqlist2.append(seq)

  return (Alignment(seqlist1), Alignment(seqlist2))

def readclustal(seqfile):
  """Iterates over sequences in an open clustal alignment.

  Parameters:
    seqfile -- an open clustal format alignment
    
  Return:
    tuple containing sequence as (seqname, seq)
    
  """
  
  def isqualline(line):
    """Returns True iff string appears to be a clustal quality line."""
    qualchars = '.:*'

    for c in line:
      if c not in qualchars:
        return False
    
    return True      
  
  seqdict = {}
  
  for line in seqfile:
    wordlist = line.split()
    if len(wordlist) != 2 or isqualline(wordlist[1]):
      continue
    try:
      seqdict[wordlist[0]].extend(wordlist[1])
    except KeyError:
      seqdict[wordlist[0]] = Sequence(wordlist[0], wordlist[1])
  
  for name in seqdict:
    if seqdict[name].type == 'empty':
      yield None
    else:
      yield seqdict[name]
  
   

def readfasta(seqfile):
  """Iterate over sequences in an open fasta alignment.
  
  Parameters:
  seqfile -- an open fasta format alignment
  
  Yield:
    Sequence objects.
    
  """
  ismust = ispir = False
  seq = name = ''
  pirtypes = ('p1', 'f1', 'dl', 'dc', 'rl', 'rc', 'xx')
  missingchar = ''
  skipline = False
  code = None
  mask = ''
  goodmask = False
  
  seqs = []
  
      
  for line in seqfile:
    if skipline:
      skipline = False
      continue
    line = line[:-1]
    if len(line) == 0:
      continue
    if line[0] == '#':
      ismust = True
      continue
    if line[0] == ';':
      continue
    if line[0] == '>':
      if '@' in line:
        ismust = True
      if seq != '':
        if name != 'mask' and seq[-1] == '*':
          seq = seq[:-1]
        if name == 'mask':
          mask = seq
        else:
          seqs.append(Sequence(name, seq, code))
        seq = ''
        code = None
      if ispir:
        wordlist = line[1:].split(';')
        name = wordlist[1].strip()
        skipline = True
      elif ismust:
        if '@' in line:
          wordlist = line[1:].split('@')
          name = wordlist[0].strip()
          code = wordlist[1].strip()
        else:
          name = line[1:].strip()
      else:
        wordlist = line[1:].split(';')
        if wordlist[0].lower() in pirtypes:
          ispir = True
          name = wordlist[1].strip()
          skipline = True
        else:
          name = wordlist[0].strip()
    else:
      if ismust:
        seq += line.replace(' ', '-').replace('*', '-')
      else:
        seq += line.replace(' ', '')

  if name != mask and seq[-1] == '*':
    seq = seq[:-1]
  if name == 'mask':
    mask = seq
  else:
    seqs.append(Sequence(name, seq, code))

  # remove masked positions
  if mask:
    goodsites = []
    goodmask = True
    for i in xrange(len(mask)):
      if mask[i] == '*':
        goodsites.append(i)
      elif mask[i] != '-':
        print 'Invalid mask in file %s. Ignoring.' % seqfile.name
        goodmask = False
        break

  for seq in seqs:
    if seq.type == 'empty':
      yield None
    else:
      if seq.type == 'dna':
        missingchar = 'N'
      else:
        missingchar = 'X'

      seq.sequence = seq.sequence.replace('?', missingchar)
      if goodmask:
        seq = seq.resample(goodsites)
      yield seq
  
  
def readgde(seqfile):
  """Iterate over sequences in an open gde alignment.

  Parameters:
    seqfile -- an open gde format alignment
  
  Yield:
    Sequence objects.
    
  """

  inseq = False
  seqs = []
  mask = ''
  seqtype = 'dna'
  offset = 0

  for line in seqfile:
    line = string.strip(line)
    if len(line) == 0:
      continue
    if inseq:
      if '"' in line:
        inseq = False
        line = line[:-1]
      if seqtype == 'mask':
        mask += line
      elif seqtype == 'dna':
        seqs[-1].extend(line.replace('?', 'N'))
      elif seqtype == 'prot':
        seqs[-1].extend(line.replace('?', 'X'))
      continue
      
    wordlist = string.split(line)
    if wordlist[0] == 'name':
      seqname = wordlist[1][1:-1]
      offset = 0 # reset offset
    elif wordlist[0] == 'sequence-ID':
      try:
        seqcode = wordlist[1][1:-1]
      except IndexError:
        seqcode = 'NO_CODE'
    elif wordlist[0] == 'sequence':
      inseq = True
      if seqtype == 'mask':
        mask = '0' * offset # if more than 1 mask specified, use last one
      elif seqtype != 'text': # neither mask nor text
        seqs.append(Sequence(seqname, sequence = ('-' * offset), type = seqtype, code = seqcode, 
                             chksum = -1)) # temporary chksum
      wordlist = line.split('"') # check to see if sequence data on this line
      if len(wordlist) > 1 and wordlist[1] != '':
        if line[-1] == '"': # if the whole sequence is on 1 line
          inseq = False
        if seqtype == 'mask':
          mask += wordlist[1]
        elif seqtype == 'dna':
          seqs[-1].extend(wordlist[1].replace('?', 'N'))
        elif seqtype == 'prot':
          seqs[-1].extend(wordlist[1].replace('?', 'X'))
      
    elif wordlist[0] == 'type':
      seqtype = wordlist[1][1:-1].lower()
      if seqtype == 'protein':
        seqtype = 'prot'
      elif seqtype == 'mask' or seqtype == 'dna' or seqtype == 'text':
        pass
      else:
        raise SeqError('GDE sequence of unrecognised type: %s' % seqtype.upper())

    elif wordlist[0] == 'offset':
      offset = int(wordlist[-1])
      
  if mask: # if mask exists
    sitelist = []
    for i in range(len(mask)):
      if mask[i] == '1':
        sitelist.append(i)
   
    for seq in seqs:
      seq.sequence = seq.resample(sitelist).sequence
        
  for seq in seqs:
    if seq.type == 'empty':
      yield None
    else:
      seq.setchecksum()
    yield seq
    

def readmega(seqfile):
  """Iterate over sequences in an open mega alignment.

  Parameters:
    seqfile -- an open mega format alignment
    
  Yield:
    Sequence objects.

  """
  seqdict = {}
  
  line = seqfile.readline().strip().lower()
  if line != '#mega':
    raise SeqError('Not a mega-format file!')
  
  ntaxa = nchar = 0
  missing = unknown = '-'
  identical = '.'
  statetype = ''
  
  firstname = ''
  
  incomment = instatement = False
  
  for line in seqfile:
    line = line.strip()
    if len(line) == 0:
      continue
    if line[0] == '[':
      incomment = True
    elif line[0] == '!':
      instatement = True
      statetype = line[1:].split()[0].lower()
    if incomment:
      if line[-1] == ']':
        incomment = False
      continue
    elif instatement:
      if line[-1] == ';':
        instatement = False
      if statetype == 'format':
        wordlist = line.lower().split()
        for word in wordlist:
          if word[:8] == 'datatype':
            seqtype = word[9:]
            if seqtype == 'nucleotide' or seqtype == 'rna' or seqtype == 'dna':
              seqtype = 'dna'
            elif seqtype == 'protein':
              seqtype = 'prot'
            else:
              raise SeqError('Sequence is neither nucleic acid nor protein!')
          if word[:5] == 'nseqs' or word[:5] == 'ntaxa':
            ntaxa = int(word[6:])
          elif word[:6] == 'nsites':
            nchar = int(word[7:])
          elif word[:9] == 'identical' or word[:9] == 'matchchar':
            identical = word[10:]
          elif word[:7] == 'missing':
            missing = word[8:]
          elif word[:5] == 'indel':
            indel = word[6:]
    else:
      if line[0] == '#': # this is a sequence
        wordlist = line[1:].split()
         # if there is a taxon group appended to the name, get rid of it
        groupidx = wordlist[0].find('{')
        
        if groupidx != -1:
          if wordlist[0][groupidx - 1] == '_':
            wordlist[0] = wordlist[0][:groupidx - 1]
          else:
            wordlist[0] = wordlist[0][:groupidx]
        if firstname == '':
          firstname = wordlist[0]
        try:
          seqdict[wordlist[0]].extend(''.join(wordlist[1:]))
        except KeyError:
          seqdict[wordlist[0]] = Sequence(wordlist[0], string.join(wordlist[1:], ''), 
                                          type = seqtype)
  
  for seq in seqdict.values():
    if seq.type == 'empty':
      yield None
    else:
      if seq.type == 'dna':
        missingchar = 'N'
      else: # if seq.type is prot
        missingchar = 'X'
      seq.sequence = seq.sequence.replace(missing, missingchar)
      if identical in seq.sequence:
        seq.homogenize(seqdict[firstname], identical)
      yield seq

def readmsf(seqfile):
  """Iterate over sequences in an open msf alignment.
  
  Parameters:
    seqfile -- an open msf format alignment
    
  Yield:
    Sequence object.
    
  """
  seqdict = {}
  missingchar = ''
  gapchar = '-'
  
  # read header
  for line in seqfile:
    line = line.strip()
    if len(line) == 0:
      continue
    if line == '//':
      break
    wordlist = line.split()
    if wordlist[0].lower() == 'name:':
      seqdict[wordlist[1]] = Sequence(wordlist[1], type = 'none')
      
  for line in seqfile:
    line = line.strip()
    if len(line) == 0:
      continue
    wordlist = string.split(line)
    try:
      seqdict[wordlist[0]].extend(string.join(wordlist[1:], ''))
    except KeyError:
      sys.stderr.write('Sequence %s has no corresponding header entry, and will be ignored!\n'
                       % wordlist[0])

  for seq in seqdict.values():
    if seq.type == 'empty':
      yield None
    else:
      seq.setchecksum()
      seq.setseqtype()
      if seq.type == 'dna':
        missingchar = 'N'
      else:# seq.type is 'prot':
        missingchar = 'X'
      
      yield seq
      
    
def readnex(seqfile):
  """Iterate over sequences in an open NEXUS alignment.
  
  Parameters:
    seqfile -- an open nexus format alignment
    
  Yield:
    Sequence object.
    
  """

  seqdict = {}
  
  line = seqfile.readline().strip().lower()
  if line != '#nexus':
    raise SeqError('Not a nexus-format file!')
      
  isnexus = False
  ntaxa = nchar = 0
  missing = '?'
  gap = '*'
  gapchar = '-'
  incomment = False
       
  # read header part of file
  for line in seqfile:
    line = line.strip()
    
    # loop to remove any comments
    while line:
      if incomment:
        commentend = line.find(']')
        if commentend == -1:
          line = ''
          break
        line = line[commentend + 1:].strip()
        incomment = False
      else:
        commentstart = line.find('[')
        
        # no more comments on this line
        if commentstart == -1:
          break
        commentend = line.find(']')
        if commentend == -1:
          line = line[:commentstart].strip()
          incomment = True
          break
        line = line[:commentstart] + line[commentend + 1:]
        line = line.strip()
          
    wordlist = line.lower().split()

    # ignored blank lines
    if len(wordlist) == 0:
      continue
    
    if wordlist[0] == 'format':
      for word in wordlist:
        if word[-1] == ';':
          word = word[:-1]
        if word[:8] == 'datatype':
          seqtype = 'dna'#word[9:]
        elif word[:7] == 'missing':
          missing = word[8:]
        elif word[:3] == 'gap':
          gap = word[4:]
    if wordlist[0] == 'dimensions':
      for word in wordlist:
        if word[-1] == ';':
          word = word[:-1]
        if word[:4] == 'ntax':
          ntaxa = int(word[5:])
        elif word[:5] == 'nchar':
          nchar = int(word[6:])
    elif wordlist[0] == 'matrix':
      isnexus = True
      break
  
  if seqtype == 'dna' or seqtype == 'nucleotide':
    missingchar = 'N'
  elif seqtype == 'protein':
    seqtype = 'prot'
    missingchar = 'X'
  else:
    raise SeqError('Sequence is neither nucleic acid nor protein!')
      
  if not (isnexus and ntaxa and nchar):
    raise SeqError('Not a nexus-format file!')
  
  # read matrix
  for line in seqfile:
    line = line.strip()

    if not line:
      continue
    
    while line:
      if incomment:
        commentend = line.find(']')
        if commentend == -1:
          line = ''
          break
        line = line[commentend + 1:].strip()
        incomment = False
        
      else:
        if line[0] == '[':
          commentend = line.find(']')
          # whole line is an unfinished comment
          if commentend == -1:
            line = ''
            incomment = True
            break
          line = line[:commentstart] + line[commentend + 1:]
          line = line.strip()
          
        # no more comments at the beginning of this line
        else:
          break
        if commentstart == -1:
          break
        

    # after comments removed, if there's nothing left of this line, skip it     
    if not line:
      continue

    # end of block
    if ';' in line:
      break

    namefound = False
    if line[0] == '\'':
      nameend = line.find('\'', 1)
      if nameend != -1:
        namefound = True
    
    while line:
      if namefound:
        commentstart = line.find('[', nameend + 1)
      else:
        commentstart = line.find('[')
        
      # no more commments on this line
      if commentstart == -1:
        break

      commentend = line.find(']', commentstart)
      if commentend == -1:
        incomment = True
        line = line[:commentstart].strip()
        break
      
      line = line[:commentstart] + line[commentend + 1:]
      line = line.strip()
    
    # after that loop, if there's nothing left of line, skip it
    if not line:
      continue
              
    if namefound:
      name = line[1: nameend]
      seq = line[nameend + 1:].strip()
    else:
      wordlist = line.split()
      if len(wordlist) == 0:
        continue
      name = wordlist[0]
      seq = ''.join(wordlist[1:])
    
    try:
      seqdict[name].extend(seq)
    except KeyError:
      seqdict[name] = Sequence(name, seq, type = seqtype)

    if len(seqdict[name]) > nchar:
      errmsg = 'Warning! The number of characters in sequence '
      errmsg += '%s exceeds nchar (%d). ' % (name, nchar)
      errmsg += 'Normally this is because two sequences in your NEXUS '
      errmsg += 'alignment have the same name. If this is the case, and your '
      errmsg += 'alignment is interleaved. there\'s going to be a big problem! '
      errmsg += 'In any case, look into it.'
      
      errmsglst = errmsg.split()
      
      errmsg = ''
      linestart = 0

      for word in errmsglst:
        linelength = len(errmsg) - linestart
        if linelength + len(word) >= 80:
          linestart = len(errmsg) + 1
          errmsg = '%s\n%s' % (errmsg, word)
        else:
          errmsg = '%s %s' % (errmsg, word)
      sys.stderr.write('%s\n' % errmsg)
      seqdict[name].truncate(nchar)

  inexset = False
  rangere = re.compile('\d+-\d+')

  exsets = []
  for line in seqfile:
    line = line.lower()
    wordlist = line.split()
    if wordlist and wordlist[0] == 'exset':
      inexset = True
    if inexset:
      ranges = [r.split('-') for r in rangere.findall(line)]
      for r in ranges:
        exsets.append((int(r[0]) - 1, int(r[1])))
      if ';' in line:
        break
      
  sitelist = []
  sitecounter = 0
  rangeptr = 0

  #print 'exsets', exsets,
  
  nexsites = 0
  for set in exsets:
    nexsites += (set[1] - set[0])
    
  #print '(%d sites excluded)' % nexsites
  
  if exsets:
    while sitecounter < nchar:
      #print 'rangeptr:', rangeptr, 'exsets[rangeptr]:', exsets[rangeptr]
      if sitecounter < exsets[rangeptr][0]:
        sitelist.extend(range(sitecounter, exsets[rangeptr][0]))
      
      sitecounter = exsets[rangeptr][1]
      rangeptr += 1
      
      if rangeptr >= len(exsets):
        sitelist.extend(range(sitecounter, nchar))
        sitecounter = nchar
      
    #print sitelist
  
  # Need to check for empty sequences!!!!!!!!!!
  for seq in seqdict.values():
    seq.sequence = seq.sequence.replace(missing, missingchar).replace(gap, gapchar)
    
    if exsets:
      seq.sequence = seq.resample(sitelist).sequence
    seqtype = seq.type
    seq.setseqtype()
    if seq.type == 'empty':
      yield None
    
    else:
      if seq.type != seqtype:
        sys.stderr.write('WARNING: identified type (%s) different from declared type (%s)\n' % (seq.type, seqtype))
        sys.stderr.write('for sequence %s. Keeping declared type.\n' % seq.name)
        seq.type = seqtype
      yield seq

def readphy(seqfile):
  """Iterate over sequences in an open PHYLIP alignment.
  
  Both interleaved and sequential-format alignments are read by this function.

  Parameters:
    seqfile -- an open phylip format alignment
    
  Yield:
    Sequence object
    
  """
  
  ntaxa, nchar = seqfile.readline().split()[:2]
  
  ntaxa = int(ntaxa)
  nchar = int(nchar)
  
  seqs = []
  issequential = False
  isrelaxed = False
  linecount = 0
  charsperline = 0
  
  while True:
    for line in seqfile:
      if isrelaxed:
        wordlist = line.split()
        if len(wordlist) != 2 or len(wordlist[1]) != nchar:
          sys.stderr.write('Error! Detected relaxed PHYLIP format in file %s, but seq length is not right!\n' % seqfile.name)
          sys.exit(1)
        seqs.append(Sequence(wordlist[0], wordlist[1], type = 'none', chksum = -1))
      else:
        wordlist = [line[:10].strip(), line[10:].strip().replace(' ', '')] 
        if charsperline == 0:
          charsperline = len(wordlist[1])
        if len(wordlist[1]) < charsperline:
          wordlist[1] = ''.join([wordlist[0], wordlist[1]])
          wordlist[0] = ''
        if len(wordlist[0]) == 0:
          if len(wordlist[1]) == 0: 
            continue
          if linecount < ntaxa:
            issequential = True
          # seq too long: indicates relaxed PHYLIP format
          if len(wordlist[1]) > nchar:
            #isrelaxed = True
            #seqfile.seek(0)
            #seqfile.readline()
            #seqs = []
            #linecount = 0
            break
          if issequential:
              seqs[-1].extend(wordlist[1])
          else:
            seqs[linecount % ntaxa].extend(wordlist[1]) # if interleaved, extend appropriate seq
        else:
          seqs.append(Sequence(wordlist[0].replace(' ', '_'), wordlist[1], type = 'none', chksum = -1))
      
      linecount += 1
    if isrelaxed:
      break
    if linecount < ntaxa:
      isrelaxed = True
      seqs = []
      linecount = 0
      seqfile.seek(0)
      seqfile.readline()
    else:
      for seq in seqs:
        # if any seqs are too short, probably relaxed PHYLIP format
        if len(seq.sequence) < nchar:
          isrelaxed = True
          seqs = []
          linecount = 0
          seqfile.seek(0)
          seqfile.readline()
          break
      if not isrelaxed:
        break
 
  
  for seq in seqs:
    seq.setseqtype()
    if seq.type == 'empty':
      yield None
    else:
      seq.setchecksum()
      if seq.type == 'dna':
        missing = 'N'
      else: # if seq.type is prot
        missing = 'X'
      seq.sequence = seq.sequence.replace('?', missing)
      yield seq
      
def readselex(seqfile):
  """Iterate over sequences in an open SELEX (old HMMER) format alignment file.
  
  Parameters:
    seqfile -- an open phylip format alignment
    
  Yield:
    Sequence object
    
  """
  
  seqdict = {}
  
  for line in seqfile:
    line = line.strip() # I'm ignoring the option to use spaces as gaps    
    
    if line == '': # skip blank lines
      continue 
    
    if line[0] == '%' or line[0] == '#': # skip comments
      continue
    
    name, seq = line.split()
    seq = seq.replace('.', '-').replace('_', '-').upper() # ignore "space" gaps
    
    try:
      seqdict[name].extend(seq)
    except KeyError:
      seqdict[name] = Sequence(name=name, sequence=seq)

  for seq in seqdict.values():
    seq.setseqtype()
    if seq.type == 'empty':
      yield None
    else:
      seq.setchecksum()
      yield seq
      
def writeclustal(align, outfile, verbose=False):
  """Write sequences to outfile in CLUSTAL format.
  
  Parameters:
    align -- Alignment object to be written
    outfile -- name of the file to which align will be written
    verbose -- write more garbage to the screen
  
  """

  outfile = open(outfile, 'w')
  
  outfile.write('CLUSTAL multiple sequence alignment\n\n')
  charsperline = 60
  
  seqlist = []
  
  for seq in align:
    seqlen = len(seq.sequence)
    if seqlen < align.length:
      seqlist.append(Sequence(seq.name, seq.sequence, seq.code, seq.chksum, seq.type))
      seqlist[-1].pad(align.length)
    else:
      seqlist.append(seq)
    
  for i in range(0, align.length, charsperline):
    for seq in seqlist:
      outfile.write('%-*s %s\n' 
                    % (align.namelen, seq.name, seq.sequence[i:i + charsperline]))
    outfile.write('\n')
    
  outfile.close()    
  

def writefasta(align, outfile, variation=None, verbose=False):
  """Write sequences to outfile in FASTA format.

  The variation option would be set to either 'must', 'div',  or 'pir'.   
  Anything else will raise an exception.
  
  Parameters:
    align -- Alignment object to be written
    outfile -- name of the file to which align will be written
    variation -- string indicating subtype of FASTA format should be written
    verbose -- write more garbage to the screen
  
  """
  
  if variation != None and variation != 'must' and variation != 'pir' and \
     variation != 'div':
    raise SeqError('Unsupported sequence type: FASTA derivative %s.' \
                   % variation)

  outfile = open(outfile, 'w')
       
  if variation == 'must':
    outfile.write('#%s\n#Produced by readseq.py, version %s\n' 
                  % (outfile.name, globals()['version']))
  else:
    tag = ''
    
  if variation == 'div':
    outfile.write('%d %d\n' % (len(align.taxlist), align.length))
    seqtype = align.seqlist[0].type
    if seqtype == 'dna':
      missing = 'N'
    else: # seqtype is protein
      missing = 'X'
     
   
  for seq in align:
    if variation == 'pir':
      if seq.type == 'dna':
        outfile.write('>DL;%s\n%s %s\n' % (seq.name, seq.name, seq.code))
      else:#if seq.type is prot
        outfile.write('>P1;%s\n%s %s\n' % (seq.name, seq.name, seq.code))
      for i in range(0, len(seq.sequence), 60):
        outfile.write('%s\n' % seq.sequence[i:i + 60])
      outfile.write('*\n\n')
    elif variation == 'div':
      outfile.write('>\n%s\n%s*\n' % (seq.name, seq.sequence.replace('-', '?').\
                                      replace(missing, '?')))

    else:
      outfile.write('>%s' % seq.name)
      if variation == 'must':
        outfile.write('@%s\n' % seq.code)
        outfile.write('%s\n' % seq.sequence.replace('-', '*'))
      #elif variation == 'div':
      #  outfile.write('\n%s*\n' % seq.sequence.replace('-', '?'))
      else:
        outfile.write('\n')
        for i in range(0, len(seq.sequence), 60):
          outfile.write('%s\n' % seq.sequence[i:i + 60])
        outfile.write('\n')
      
  outfile.close()
  
def writegde(align, outfile, verbose=False):
  """Write sequences to outfile in GDE format.
  
  Parameters:
    align -- Alignment object to be written
    outfile -- name of the file to which align will be written
    verbose -- write more garbage to the screen
  
  """
  
  outfile = open(outfile, 'w')
  
  for seq in align:
    outfile.write('{\n%-16s"%s"\n' % ('name', seq.name))
    seqlen = len(seq.sequence)

    if seq.type == 'dna':
      outfile.write('%-16s"%s"\n' % ('type', 'DNA'))
    else: #if seq.type == 'prot':
      outfile.write('%-16s"%s"\n' % ('type', 'PROTEIN'))

    outfile.write('%-16s"%s"\n' % ('sequence-ID', seq.code))
    outfile.write('%-16s"%s , %d bases, %X checksum."\n' 
                  % ('descrip', seq.name, seqlen, seq.chksum))
    outfile.write('%-16s%s\n' % ('creation-date', time.strftime('%x %X')))
                  
    outfile.write('%-16s"\n' % ('sequence'))
    for i in range(0, seqlen, 60):
      outfile.write('%s\n' % (seq.sequence[i:i + 60]))
    outfile.write('"\n}\n')
  
  outfile.close()
  
def writemega(align, outfile, verbose=False):
  """Write sequences to outfile in MEGA format.  
  
  Parameters:
    align -- Alignment object to be written
    outfile -- name of the file to which align will be written
    verbose -- write more garbage to the screen
  
    """
  
  seqtype = align.seqlist[0].type
  
  if seqtype == 'dna':
    seqtype = 'Nucleotide'
    missing = 'N'
  else:#if seqtype == 'prot':
    seqtype = 'Protein'
    missing = 'X'
  
  outfile = open(outfile, 'w')
  
  outfile.write('#MEGA\n!Title %s, written using readseq.py, version %s;\n' \
                % (outfile.name, globals()['version']))
  outfile.write('Format\n   Datatype=%s\n   NSeqs=%d NSites=%d\n' \
                % (seqtype, len(align.seqlist), align.length))
  outfile.write('   Identical=. Missing=%s Indel=-;\n\n' % missing)
  
  charsperline = 60
  
  for i in range(0, align.length, charsperline):
    for seq in align:
      outfile.write('#%-*s %s\n' \
                    % (align.namelen, seq.name, seq.sequence[i:i + charsperline]))
    outfile.write('\n')
    
  outfile.close()    

def writemsf(align, outfile, verbose=False):
  """Writes sequences to outfile in GCG/MSF format.  
  
  Parameters:
    align -- Alignment object to be written
    outfile -- name of the file to which align will be written
    verbose -- write more garbage to the screen
  
  """
    
  outfile = open(outfile, 'w')
  
  seqtype = align.seqlist[0].type # assume all seqs have same type
  if seqtype == 'dna':
    seqtype = 'N'
    outfile.write('!!NA_MULTIPLE_ALIGNMENT\n\n')
  else:# seqtype == 'prot':
    seqtype = 'P'
    outfile.write('!!AA_MULTIPLE_ALIGNMENT\n\n')
  outfile.write(' %s  MSF: %d  Type: %s  Check: %d ..\n\n' 
                % (outfile.name, align.length, seqtype, align.chksum))

  seqs = []
  
  for seq in align:
    seqlen = len(seq.sequence)
    if seqlen < align.length:
      if verbose:
        print 'Sequence %s is too short, padding to fit!' % seq.name
      seqs.append(Sequence(seq.name, seq.sequence, seq.code, seq.chksum, seq.type))
      seqs[-1].pad(align.length, '.')
    else:
      seqs.append(seq)
  
    outfile.write(' Name: %-16s Len: %5d  Check: %5d  Weight:  1.00\n' \
                  % (seq.name[:16], align.length, seq.chksum))
    
  outfile.write('\n//\n\n')

  charsperline = 60
  
  for i in range(0, align.length, charsperline):
    for seq in seqs:
      outfile.write('%16s %s\n' \
                    % (seq.name[:16], seq.sequence[i:i + charsperline].replace('-', '.')))
    outfile.write('\n')
  
  outfile.close()
  
def writenex(align, outfile, interleaved=False, verbose=False):
  """Writes sequences to outfile in nexus format (optionally interleaved).  
  
  Parameters:
    align -- Alignment object to be written
    outfile -- name of the file to which align will be written
    interleaved -- if True, sequences will be interleaved in output
    verbose -- write more garbage to the screen
  
  """

  ntaxa = len(align.seqlist)

  seqtype = align.seqlist[0].type
  if seqtype == 'dna':
    missingchar = 'N'
    seqtype = 'DNA'
  else: # seqtype is'prot':
    seqtype = 'Protein'
    missingchar = 'X'

  outfile = open(outfile, 'w') 
  outfile.write('#NEXUS\nBegin data;\n')
  outfile.write('    Dimensions ntax=%d nchar=%d;\n' % (ntaxa, align.length))

  mixedtype = False

  if align.charsets:
    type = None
    for set in align.charsets:
      if type and set[1] != type:
        mixedtype = True
        break
      else:
        type = set[1]
    
  if mixedtype:
    formatstr = '    Format datatype=mixed('
    for set in align.charsets:
      formatstr += '%s:%d-%d,' % (set[1], set[2], set[3])
      
    formatstr = formatstr[:-1] # cut off trailing ','
    formatstr += ') '
      
    outfile.write(formatstr)
    
  else:
    outfile.write('    Format datatype=%s ' % seqtype)

  if interleaved:
    outfile.write('interleave missing=%s gap=-;\n' % missingchar)
    outfile.write('    Matrix\n')
    charsperline = 60
  else:
    outfile.write('missing=%s gap=-;\n' % missingchar)
    outfile.write('    Matrix\n')
    charsperline = align.length

  seqs = []
  for seq in align:
    seqlen = len(seq.sequence)
    if seqlen < align.length:
      if verbose:
        print 'Sequence %s is too short, padding to fit!' % seq.name
      seqs.append(Sequence(seq.name, seq.sequence, seq.code, seq.chksum, seq.type))
      seqs[-1].pad(align.length)
    else:
      seqs.append(seq)
      

  for i in range(0, align.length, charsperline):
    for seq in seqs:
      outfile.write('%*s %s\n' % (align.namelen, seq.name, seq.sequence[i:i + charsperline]))
    outfile.write('\n')

  outfile.write(';\nEnd;\n\n')
    
  if align.charsets:
    outfile.write('Begin assumptions;\n\n')
    
    for set in align.charsets:
      outfile.write('  charset %s=%d-%d;\n' % (set[0], set[2], set[3]))
    
    outfile.write('\nEnd;\n\n')

  outfile.close()
    

def writephy(align, outfile, interleaved=False, relaxed=True, verbose=False, append=False):
  """Writes sequences to outfile in phylip format (optionally interleaved).  
  
  Parameters:
    align -- Alignment object to be written
    outfile -- name of the file to which align will be written
    interleaved -- if True, sequences will be interleaved in output
    verbose -- write more garbage to the screen
  
"""
  ofname = outfile
  ntaxa = len(align.seqlist)
  namelist = []
  
  if append:
    outfile = open(ofname, 'a')
  else:
    outfile = open(ofname, 'w')
  outfile.write(' %d %d\n' % (ntaxa, align.length))
    
  if interleaved:
    charsperline = 60
  else:
    charsperline = align.length

  seqs = []
  for seq in align:
    if relaxed:
      seqname = seq.name
      namecount = 2
      while seqname in namelist:
        seqname = seq.name + str(namecount)
        namecount += 1
    else:
      seqname = seq.name[:10]
      namecount = 2
      while seqname in namelist:
        countstr = str(namecount)
        seqname = seqname[:10 - len(countstr)] + countstr
        namecount += 1
    namelist.append(seqname)
    
    seqs.append(Sequence(seqname, seq.sequence, seq.code, seq.chksum, seq.type))

    seqlen = len(seq.sequence)
    if seqlen < align.length:
      if verbose:
        print 'Sequence %s is too short, padding to fit!' % seq.name
      seqs[-1].pad(align.length)
    
  for i in range(0, align.length, charsperline):
    for seq in seqs:
      if i == 0:
        if relaxed:
          outfile.write('%s %s\n' % (seq.name.replace(' ', '_'), seq.sequence[i:i + charsperline]))
        else:
          outfile.write('%-10s %s\n' % (seq.name, seq.sequence[i:i + charsperline]))
      else:
        outfile.write('%s %s\n' % (' ' * 10, seq.sequence[i:i + charsperline]))
    #outfile.write('\n')
  
  outfile.close()
  if align.charsets:
    modelfile = open('%s.model' % ofname, 'w')
    
    for set in align.charsets:
      if set[1] == 'dna':
        modelfile.write('DNA, ')
      else:
        modelfile.write('WAG, ')
      genename = '.'.join(set[0].split('.')[:-1])
      if not genename:
        genename = set[0]
      modelfile.write('%s=%d-%d\n' % (genename, set[2], set[3]))
    
    modelfile.close()

class Alignment:
  """Defines a collection of sequences."""
  
  def __init__(self, seqlist, name = None, charsets = None):
    self.name = name
    self.seqlist = seqlist
    self.taxlist = []
    self.length = 0
    self.chksum = 0
    self.namelen = 0
    self.charsets = charsets
    
    namecounters = {}
    
    for seq in seqlist:    
      if seq.name in namecounters:
        namecounters[seq.name] += 1
        seq.setname('%s%d' % (seq.name, namecounters[seq.name]))
      else:
        namecounters[seq.name] = 1
      self.taxlist.append(seq.name)
      seqlen = len(seq.sequence)
      if seqlen > self.length:
        self.length = seqlen
      self.chksum += seq.chksum
      namelen = len(seq.name)
      if namelen > self.namelen:
        self.namelen = namelen
        
    self.taxlist.sort()
        
  def __iter__(self):
    """Iterate over sequences in the Alignment.
    
    Yield:
      Sequence object.
      
    """
    
    for seq in self.seqlist:
      yield seq
      
  def __len__(self):
    """Determine the length of the Alignment.

    Return:
      Integer length of the alignement.
      
    """
    
    return self.length
  
  def __getitem__(self, key):
    """Retrieve the 'key'th column of the alignment.
    
    Raise a SeqError if the key is not an integer.
    
    Parameters:
      key -- an integer, representing the column of the alignment to retrieve.
      
    Return:
      A list containing the characters found in column 'key', in the order of
      self.taxlist.
    
    """
    
    try:
      column = []
      for s in self:
        column.append(s[key])
      return column
        
    except TypeError:
      raise SeqError('Invalid key: %s.  Integer required.' % str(key))

class Sequence:
  
  """
  Defines a sequence object.
  
  """
  
  def __init__(self, name='', sequence='', code=None, chksum=0, type=None):
    """Sequence contructor.
    
    Parameters:
      name -- string, name to associate with the Sequence (default: empty)
      sequence -- string, the sequence itself (default: empty)
      code -- usually a string, probably a database code (e.g. GI number; 
              default: None)
      chksum -- required by certain alignment formats (default: calculate)
      type -- type of sequence data.  Normally 'dna', 'prot', 'rna' (default:
              automatic)
    
    """
    
    self.name = name
    self.sequence = sequence
    self.code = code
    if chksum == 0:
      self.setchecksum()
    else:
      self.chksum = chksum
    if type == None:
      self.setseqtype()
    else:
      self.type = type
      
  def __str__(self):
    """'To-string' method.
    
    Return:
      Sequence as 'name: x positions', where x is the lenght of the sequence
      
    """
    
    return '%s: %d positions' % (self.name, len(self.sequence))
  
  def __add__(self, other):
    """Overload addition operator.
    
    The second sequence is concatenated to the first.  This method modifies 
    the Sequence object in-place, but also returns the original Sequence.
    
    Parameters:
      other -- another Sequence object
      
    Return:
      The original Sequence object, with the other sequence appended.
    
    """

    try:
      self.sequence += other.sequence
    except AttributeError:
      raise SeqError('This is not a valid Sequence object!') 
    return self
    
  def __len__(self):
    """Return the length of the sequence."""
    return len(self.sequence)
  
  def __getitem__(self, key):
    """Retrieve the character in the Sequence with index 'key'. 
    
    Raise TypeError if key is not an integer.
    
    Parameters:
      key -- index of the character to be retrieved, should be an integer.
      
    Return:
      The character with index 'key'.
    
    """
    return self.sequence[key]
  
  def homogenize(self, other, identical):
    """Replace all 'identical' chars in self with matching char in other.
    
    Parameters:
      other -- another Sequence
      identical -- character that indicates that this position of Sequence
                   object is identical to the character at that index in 'other'
    
    """
    
    if len(other.sequence) != len(self.sequence):
      raise SeqError('Cannot homogenize sequences of different length!')
    
    for i in range(len(self.sequence)):
      if self.sequence[i] == identical:
        self.sequence = self.sequence[:i] + other.sequence[i] + self.sequence[i + 1:] 

  def setchecksum(self):
    """Set self.chksum to the ascii-value checksum for the sequence."""
    
    chksum = 0
    for c in self.sequence:
      chksum += ord(c)
      
    self.chksum = chksum
  
  def countnongapchars(self):
    """Count the number of characters in Sequence that are not 'missing data'.
    
    Return:
      The number of non-gap characters in the Sequence.
    
    """
    if self.type == 'dna':
      missing = ('N', '-')
  
    else:
      missing = ('X', '-')
    
    nongapcount = 0
    for char in self.sequence:
      if char not in missing:
        nongapcount += 1
        
    return nongapcount
     
  def extend(self, seq):
    """Append a string to the Sequence.
    
    Parameters:
      seq -- string to be appended to the Sequence
      
    """
    
    self.sequence += seq
    
  def setname(self, newname):
    """Change the name of the Sequence.
    
    Parameters:
      newname -- string to which the Sequence's name will be changed
    
    """
    self.name = newname

  def truncate(self, length):
    """Remove all but the first 'length' characters from the Sequence.
    
    If 'length' exceeds the length of the Sequence, it will remain unchanged.
    
    Parameters:
      length -- length to which the Sequence will be truncated.
    
    """
    
    self.sequence = self.sequence[:length]
    
  def pad(self, length, padchar = '-'):
    """Add to the length of the sequence by padding the end.
    
    Raise SeqError if the length of the Sequence exceeds length,
    
    Parameters:
      length -- length to which the Sequence will be increased
      padchar -- character used to pad the sequence.  This argument must be a
                 string.  If it is not a single character, error messages will
                 be printed, but no Exception will be raised.
    
    """
    
    if len(padchar) > 1:
      sys.stderr.write('%s is too long for a padding character.  Using first character only.\n' 
                       % padchar)
      padchar = padchar[0]
    elif len(padchar) < 1:
      sys.stderr.write('Padding character must not be an empty string! Setting to \'-\'.\n')
      padchar = '-'
    
    seqlen = len(self.sequence)
    if seqlen > length:
      raise SeqError('Length of sequence \'%s\' already exceeds %d!' 
                     % (self.name, length))
    self.sequence += (length - len(self.sequence)) * padchar
  
  def fpad(self, length, padchar = '-'):
    """Add to the length of a Sequence by padding the beginning.
    
    Raise SeqError if the length of the Sequence exceeds length,
    
    Parameters:
      length -- length to which the Sequence will be increased
      padchar -- character used to pad the sequence.  This argument must be a
                 string.  If it is not a single character, error messages will
                 be printed, but no Exception will be raised.
                 
    """
    
    if len(padchar) > 1:
      sys.stderr.write('%s is too long for a padding character.  Using first character only.\n' 
                       % padchar)
      padchar = padchar[0]
    elif len(padchar) < 1:
      sys.stderr.write('Padding character must not be an empty string! Setting to \'-\'.\n')
      padchar = '-'
    
    seqlen = len(self.sequence)
    if seqlen > length:
      raise SeqError('Length of sequence \'%s\' already exceeds %d!'
                     % (self.name, length))
    self.sequence = (length - len(self.sequence)) * '-' + self.sequence
    
  def cleangaps(self):
    """Remove all gap characters from the Sequence."""
    
    
    if self.type == 'prot':
      gapchars = 'x-?'
    elif self.type == 'dna':
      gapchars = 'n-?' 
    else:
      raise SeqError('Can only remove gaps from DNA and protein sequences!')
    
    seqlist = []
    for char in self.sequence:
      if char.lower() not in gapchars:
        seqlist.append(char)
    
    self.sequence = ''.join(seqlist)
  
  def setseqtype(self):
    """Set Sequence's type field to guessed type.
    
    This function is not foolproof!  Its decision is based on the frequency that 
    letters other than a, c, g, t, and u are observed in the sequence.  The type
    will be set to 'dna', 'prot', or 'empty'.

    """
    
    alphacount = dnacount = numcount = missingcount = 0
    chardict = {}
    dnachars = ('a', 'c', 'g', 't', 'u')
    missingchars = ('-', '?', '*', 'x', 'n')

    seq = self.sequence.lower()
    for char in seq:
      if char in missingchars:
        missingcount += 1
      else:
        if char in dnachars:
          dnacount += 1
        if char in string.letters:
          alphacount += 1
        elif char in string.digits:
          numcount += 1
      try:
        chardict[char] += 1
      except KeyError:
        chardict[char] = 1

    if dnacount == 0 and alphacount == 0:
      if missingcount == 0:
        raise SeqError('Sequence is neither nucleic acid nor protein!')
      else:
        self.type = 'empty'
        return

    # fewer than 90% of characters are nucleic acid
    if (float(dnacount) / alphacount) < 0.9:
      self.type = 'prot'
    else:
      self.type = 'dna'

  def resample(self, sites):
    """Resample specified sites from Sequence.
    
    Raise TypeError if sites can't be iterated over, or doesn't contain
    integers.

    Parameters:
      sites -- a list of indices of sites to sample (should be integers).
    
    Return:
      A new Sequence object containing only the sites specified.
      
    """
    
    bootseq = []
    for i in sites:
      try:
        bootseq.append(self.sequence[i])
      except IndexError:
        break
    
    return Sequence(self.name, ''.join(bootseq), self.code, type = self.type)
  
class SeqError(Exception):
  """Exception class designed for sequence manipulation errors."""
  def __init__(self, value=''):
    self.value = value
  def __str__(self):
    return repr(self.value)
  
