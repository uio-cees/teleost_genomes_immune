"""
A reusable module with of tree manipulation functions. Originally included in 
Concaterpillar.

"""

import os, re, string

def aresametree(treeA, treeB):
  """Determine whether two tree strings are topologically equivalent.
  
  Parameters:
    treeA, treeB -- tree strings in Newick format, may be rooted or unrooted.
    
  Return:
    True iff treeA and treeB are equivalent, False otherwise
    
  """

  partsA = getpartitions(treeA)
  
  return havesamepartitions(partsA, treeB)

def buildalltrees(taxa, verbose=False):
  """Iterate over all possible rooted trees for a given set of taxa.
  
  This function is fairly simple, and not terribly useful.  The script 
  scramble.py can be used instead to make all possible (unrooted) resolutions
  of an unresolved tree.  In other words, it produces an exhaustive list of
  unrooted topologies with a given list of constraints.
  
  Parameters:
    taxa -- a list of labels for the terminal nodes (vertices with 1 edge)
    verbose -- if True, a "status" periodically indicating how many trees have
               been built will be printed to the terminal.
           
  Yield:
    Rooted tree strings in Newick format, without the final semi-colon.
    
  """
  
  ntrees = [1, 1] 
  
  for i in xrange(3, len(taxa) + 1):
    ntrees.append(ntrees[-1] * ((i * 2) - 3))

  # ntrees now contains the number of rooted trees possible for (index + 1) taxa
  
  totaltrees = ntrees[-1]
  
  root = TreeNode(name = 'root')
    
  for i in xrange(totaltrees):
    
    # arbitrarily add 2 taxa to the tree
    root.addchild(TreeNode(taxa[0]))
    root.addchild(TreeNode(taxa[1]))
    
    # add each taxon as a sister to each other taxon already in the tree
    for j in xrange(2, len(taxa)):
      allnodes = root.getnodes()
      sisterindex = (i / (totaltrees / ntrees[j])) % ((j * 2) - 1)
      sister = allnodes[sisterindex]
      
      if sister.isroot():
        root = TreeNode('root')
        sister.name = ''
        root.addchild(sister)
        root.addchild(TreeNode(taxa[j]))
      else:
        sister.addsister(TreeNode(taxa[j]))
        
    yield str(root)[:-1] # yield the tree string w/o the semi-colon
    root.clear() # remove all children, start again
    
    if verbose:
      if i % (totaltrees / 10) == 0:
        print 'Built %d trees.' % (i + 1)
  
  

def buildstartree(taxa):
  """Produce a star tree string containing the names specified in taxa.
  
  Parameters:
    taxa -- a list of labels for terminal nodes in the start tree.

  Return:
    A star tree string in Newick format.
    
  """
  
  tree = '(' + taxa [0]

  for i in taxa[1:]:
    tree = '%s,%s' % (tree, i)
  tree = tree + ');'

  return tree

def cleantree(tree):
  """Remove branch lengths and pos int bootstrap values from a tree string.
  
  Bootstrap values are often represented as internal node labels, usually as
  integers,  If they appear as floating-point numbers (e.g. probabilities), this
  function will not work.
  
  Parameters:
    tree -- a tree string in Newick format.
  
  Return:
    Newick-format string without branch lenghts/bootstrap values.
  
  """
  tree = string.strip(tree)
  regexp = re.compile(r':\d+.\d+')
  tree = regexp.sub('', tree)
  
  regexp = re.compile(r'\)\d+') # internal node support values
  tree = regexp.sub(')', tree)
  
  return tree
        
def completeparts(partitions, taxlst):
  """Produce a sorted list of bipartitions from a list of half-bipartitions.
  
  This function should normally not be called directly.  If you wish to get a
  list of bipartitions in a tree, it's much easier to call "getpartitions".
  
  Parameters:
    partitions -- a list containing "half-bipartitions".  That is, only the taxa
                  on one side of the partition will be represented (such a list
                  is produced by the partition function).
    taxlst -- a list containing the names of all taxa in the tree to which the
              list of partitions belongs.
              
  Return:
    A list of bipartitions (pairs of lists of sequence names that appear on
    opposite sides of a bipartition).
  
  """

  completeparts = []
  taxlst.sort()

  for i in partitions:
    i.sort()
    outgroup = []
    for j in taxlst:
      if j not in i:
        outgroup.append(j)

    # keep the list sorted... comparison between first elements in lists
    if i < outgroup:
      completeparts.append([i, outgroup])
    else:
      completeparts.append([outgroup, i])

  completeparts.sort()
  return completeparts

def fixtreefile(filename, roottree = None, clean = False):
  """Remove linebreaks from within trees, given the name of a treefile.

  If the optional root argument is given, and is True, trees will be rooted.
  If it is False, all trees will be unrooted.  If it is NOT given, trees
  will be left as is.  In any case, this function works "in place": original
  tree file is modified.
  
  Parameters:
    filename -- name of a (closed) tree file to be modified.
    roottree -- optional argument that determines whether trees in the modified 
                tree file should be rooted, unrooted, or left as is (see above).
    clean -- optional argument that determines whether output trees should be
             "clean" (passed through cleantree function).
             
  """
  
  treefile = open(filename, 'U')
  newtf = open('tempfile', 'w')
  tree = ''

  for line in treefile:
    tree += line.strip()
    if tree[-1] == ';':
      if roottree == None:
        pass
      elif roottree == False:
        tree = unroot(tree)
      elif roottree == True:
        tree = root(tree)
        
      if clean == True:
        tree = cleantree(tree)
        
      newtf.write('%s\n' % tree)
      tree = ''
    
  newtf.close()
  treefile.close()
  
  os.rename('tempfile', filename)
    
def getpartitions(tree):
  """Return a list of bipartitions, given a Newick tree string.

  Parameters:
    tree -- a tree string, rooted or unrooted.  Must be in Newick format.
  
  Return:
        A list of bipartitions (pairs of lists of sequence names that appear on
    opposite sides of a bipartition).
    
  """
  
  tree = unroot(tree)
  partitions = partition(tree)

  return completeparts(partitions, gettaxa(tree))
  
def getsubtrees(tree):
  """Identify left and right subtrees from a Newick tree string.

  Parameters:
    tree -- a tree string, semicolon optional.  May be rooted, unrooted, or 
            multifurcating.
            
  Return:
    A list containing (left, right) subtrees (as Newick tree strings)
  
  """

  # remove trailing ';' and superfluous parentheses
  if tree[-1] == ';':
    tree = tree[1:-2]
  else:
    tree = tree[1:-1]

  # if "tree" has only 1 internal node, return a list of taxa
  if '(' not in tree: 
    return string.split(tree, ',')

  parentheses = 0
  curlies = 0
  startidx = 0
  subtrees = []

  for i in xrange(len(tree)):
    if tree[i] == ',':
      if parentheses == 0 and curlies == 0:
        subtrees.append(tree[startidx:i])
        startidx = i + 1
    if tree[i] == '(':
      parentheses += 1
    elif tree[i] == ')':
      parentheses -= 1
    elif tree[i] == '{':
      curlies += 1
    elif tree[i] == '}':
      curlies -= 1

  subtrees.append(tree[startidx:])

  return subtrees

def gettaxa(tree):
  """Extract taxon names from a Newick-format tree string.

  Parameters:
    tree -- a Newick-format tree string.
    
  Return:
    A list of names at leaves.
    
  """
  tree = cleantree(tree)
  tree = string.replace(tree, '(', ' ')
  tree = string.replace(tree, ',', ' ')
  tree = string.replace(tree, ')', ' ')
  tree = string.replace(tree, ';', ' ')

  return string.split(tree)

def havesamepartitions(partsA, treeB):
  """Determine whether treeB is equivalent to the tree with bipartitions partsA.

  Note that two trees are topologically equivalent iff they share the same 
  bipartition list.  This function is designed for comparison of one tree (A) 
  against a whole batch of others.  Note that if all trees are to be compared 
  multiple times each, this is not the most efficient comparison function.  
  Instead, a partition list should be constructed for each, and these should be 
  compared directly (using the comparison operator ('=='))

  Parameters:
    partsA -- a list of bipartitions for some tree (produced with getpartitions)
    treeB -- a Newick-format tree string (rooted or unrooted).

  Return True iff treeB can be described by bipartitions partsA, False otherwise
    
  """
  partsB = getpartitions(treeB)

  if partsB == partsA:
    return True

  return False

def partition(tree):
  """Produce a list of "half-bipartitions" of a tree.

  Recursive function that produces a list of bipartitions of a Newick format 
  tree.  The list contains only half of each of the bipartitions that should be 
  completed using the completeparts function.
  
  Normally neither of these functions should be called directly.  See the
  getpartitions function instead.
  
  Parameters:
    tree -- a Newick-format tree string (must be unrooted).
    
  Return:
    List of bipartitions, represented as lists, but containing only the taxa on
    one side of the bipartition.
  
  """  
  subtrees = getsubtrees(cleantree(tree))
  
  partitions = []
  for i in subtrees:
    taxa = gettaxa(i)
    if len(taxa) != 1:
      partitions.append(taxa)
      parts = (partition(i))
      if len(parts) != 0:
        partitions.extend(parts)
        
  return partitions

def prunetree(tree, taxlst):
  """Produce a Newick tree string with taxa not appering in taxlst removed.

  Any taxa appearing in taxlst, but not in tree, will be ignored.
  
  Parameters:
    tree -- Newick-format tree string
    taxlist -- list of taxa whose names should be retained
    
  Return:
    Newick tree string.
    
  """
  
  treetaxa = gettaxa(tree)
  t = string2tree(tree)
  
  for taxon in treetaxa:
    if taxon not in taxlst:
      t.removenode(taxon)
      
  return str(t)

def readtrees(treefile):
  """Iterate over over trees in a Newick tree file.
  
  If puzzle likelihoods are present, they will be removed.
  
  Parameters:
    treefile -- name of a Newick-format tree file
    
  Yield:
    Tree string read from the file.
  
  """
  
  tf = open(treefile)

  tree = ''
  for line in tf:
    line = line.strip()
    tree += line
    if tree[-1] == ';':
      if tree[0] == '[': # get rid of puzzle likelihood
        tree = tree.split(']')[1]
      yield tree
      
      tree = ''
    
  tf.close()
  
def root(tree):
  """Produce a rooted copy of a tree string.
  
  Parameters:
    tree -- a Newick-format tree string, either rooted or unrooted.
  
  Return:
    Rooted tree string.
    
  """
  
  subtrees = getsubtrees(tree)
  if len(subtrees) < 3:
    return tree
  
  if len(subtrees) > 3:
    raise TreeError('Unable to root trees with more than 3 subtrees.')
  
  if ':' not in tree:
    return '((%s,%s),%s);' % (subtrees[0], subtrees[1], subtrees[2])
  
  brstart = subtrees[2].rfind(':')
  brlen = float(subtrees[2][brstart + 1:])/2
  
  return '((%s,%s):%.4f,%s:%.4f);' % (subtrees[0], subtrees[1], brlen, \
                                     subtrees[2][:brstart], brlen)
  

def string2tree(treestr):
  """Build a tree data structure from a Newick string representation.
  
  Raise TreeError if tree string seems not to be a valid tree.
  
  Parameters:
    treestr -- Newick-format tree string, must be rooted.
    
  Return:
    TreeNode pointer to the root of the tree.
  
  """
  if treestr[-1] != ';':
    raise TreeError('This string is not a valid tree!')

  root = TreeNode('root')
  nodeptr = root

  i = 0
  
  treechars = (':', ';', ',', '(', ')')
  
  while i < len(treestr):
    if treestr[i] == ';':
      break
    if treestr[i] =='(' or treestr[i] == ',':
      if treestr[i] == ',':
        nodeptr = nodeptr.parent
      if treestr[i + 1] == '(':
        nodeptr.addchild(TreeNode())
        i += 1
      else:
        j = i + 1
        while treestr[j] not in treechars:
          j += 1
        if treestr[j] == ':':
          k = j + 1
          while treestr[k] in string.digits or treestr[k] == '.':
            k += 1
          nodeptr.addchild(TreeNode(treestr[i + 1:j], float(treestr[j + 1:k])))
          i = k
        else:
          nodeptr.addchild(TreeNode(treestr[i + 1:j]))
          i = j
      if nodeptr.right: # if there's a right child to move to, we've just added it
        nodeptr = nodeptr.right
        
      else: # if not, we need to move left
        nodeptr = nodeptr.left
   
    elif treestr[i] == ')':
      nodeptr = nodeptr.parent
      i += 1
      if treestr[i] == ':':
        j = i + 1
        while treestr[j] == '.' or treestr[j] in string.digits:
          j += 1
        
        nodeptr.setbrlen(float(treestr[i + 1:j]))
        i = j
        
    else:
      raise TreeError('This string is not a valid tree!')
    
  return root

def string2unresolvedtree(treestr):
  """Build a tree data structure from a Newick string representation.
  
  Raise TreeError if tree string seems not to be a valid tree (multifurcations
  are allowed, though).

  Parameters:
    treestr -- Newick-format tree string

  Return:
    UnresolvedNode pointer to the root.
    
  """
  
  intcounter = 0
  if treestr[-1] != ';':
    raise TreeError('This string is not a valid tree!')

  root = UnresolvedNode('root')
  nodeptr = root

  i = 0
  
  treechars = (':', ';', ',', '(', ')', '{', '}')
  
  while i < len(treestr):
    if treestr[i] == ';':
      break
    if treestr[i] =='(' or treestr[i] == '{' or treestr[i] == ',':
      if treestr[i] == ',':
        nodeptr = nodeptr.parent
      if treestr[i + 1] == '(' or treestr[i + 1] == '{':
        nodeptr.addchild(UnresolvedNode('int%d' % intcounter))
        intcounter += 1
        i += 1
      else:
        j = i + 1
        while treestr[j] not in treechars:
          j += 1
        if treestr[j] == ':':
          k = j + 1
          while treestr[k] in string.digits or treestr[k] == '.':
            k += 1
          nodeptr.addchild(UnresolvedNode(treestr[i + 1:j], float(treestr[j + 1:k])))
          i = k
        else:          
          nodeptr.addchild(UnresolvedNode(treestr[i + 1:j]))          
          i = j
      nodeptr = nodeptr.lastchild()
      
    elif treestr[i] == ')' or treestr[i] == '}':
      nodeptr = nodeptr.parent
      j = i + 1
      if treestr[j] not in treechars:
        while treestr[j] not in treechars:
          j += 1
        try:
          nodeptr.bootstrap = float(treestr[i + 2: j - 1])
        except:
          nodeptr.bootstrap = float(treestr[i+1:j])
      i = j
      if treestr[j] == ':':
        j += 1
        while treestr[j] == '.' or treestr[j] in string.digits:
          j += 1
        nodeptr.setbrlen(float(treestr[i + 1:j]))
        i = j
        
    else:
      print 'treestr[i]:', treestr[i]
      raise TreeError('This string is not a valid tree!')
    
  return root

def unroot(tree):
  """Returns an unrooted copy of a tree.

  Tree passed may be rooted or unrooted, must be in Newick format.
  
  Parameters:
    tree -- Newick string representation of a tree.
    
  Return:
    Unrooted Newick tree string.
  
  """
    
  subtrees = getsubtrees(tree)

  if len(subtrees) == 3:
    return tree
  
  if ':' not in tree:
    if '(' in subtrees[0]:
      return '(%s,%s);' % (subtrees[0][1:-1], subtrees[1])
    if '(' in subtrees[1]:
      return '(%s,%s);' % (subtrees[0], subtrees[1][1:-1])
    else:
      return tree  
  
  br0start = subtrees[0].rfind(':')
  brlen0 = float(subtrees[0][br0start + 1:])
  
  br1start = subtrees[1].rfind(':')
  brlen1 = float(subtrees[1][br1start + 1:])

  if '(' in subtrees[0]:
    return '%s,%s:%.4f);' % (subtrees[0][:br0start - 1], subtrees[1][:br1start], \
                             brlen0 + brlen1)                               
  if '(' in subtrees[1]:
    return '(%s:%.4f,%s;' % (subtrees[0][:br0start], brlen0 + brlen1, \
                             subtrees[1][1:br1start])
  return tree # if there are only 2 species

class TreeNode:
  """Specification for nodes in a bifurcating tree data structure."""
  
  def __init__(self, name = '', brlen = 0.0):
    """TreeNode constructor."""
    self.name = name
    self.brlen = brlen
    self.parent = self.left = self.right = None
    
  def __str__(self):
    """Return a Newick format string representation of a tree."""
    if self == None:
      return
    
    if not self.isinternal():
      if self.brlen:
        return '%s:%.4f' % (self.name, self.brlen)
      else:
        return self.name
    
    if self.isroot():
      return '(%s,%s);' % (str(self.left), str(self.right))
    
    if self.brlen:
      return '(%s,%s):%.4f' % (str(self.left), str(self.right), self.brlen)
    else:
      return '(%s,%s)' % (str(self.left), str(self.right))
    
  def findnode(self, name):
    """Find the node in the tree with the given name field.
    
    Parameters:
      name -- string containing the name to be found in the tree.
    
    Return: 
      TreeNode object (None, if name is not found).
    
    """
    
    if self.name == name:
      return self

    found = None
    
    if self.left != None:
      found = self.left.findnode(name)
        
    if found == None and self.right != None:
      found = self.right.findnode(name)
      
    return found
  
  def isroot(self):
    """Return True iff TreeNode has no parent."""
    
    if self.parent == None:
      return True
    
    return False
  
  def isinternal(self):
    """Return True iff TreeNode has at least one child."""
    
    if self.left or self.right: # if either child is not None
      return True
    
    return False
  
  def isleft(self):
    """Return True iff node is the left child of its parent."""
    
    try:
      return self.parent.left is self
    except AttributeError:
      return False
    
  def setbrlen(self, brlen):
    """Change the node's branch length.
    
    Parameters:
      brlen -- new branch length (should be a float, but this is not tested).
    
    """
    
    self.brlen = brlen
  
  def removenode(self, name):
    """Remove the node with the given name field from the tree below self.

    Should be called on a pointer to the root of a tree, which is searched for a
    node with the given name field.  If the node is not found, no action is 
    taken.  If the node is found, it is removed, as is its parent node.  If the 
    tree has branch lengths, the node's "sibling's" branch length is increased 
    by the length of the branch leading to the parent that was removed.  Slight 
    modifications to this procedure are made in the case that the node to be 
    removed is a direct child of the root node.
    
    If the name is not found in the tree, no node is removed [This may be
    changed at some point, an exception should probably be raised].
    
    Parameters:
      name -- the name (string) of the node to be removed
      
    """
    
    node = self.findnode(name)
    
    if node == None:
      return
    
    # special cases: node is root, node's parent is root
    if node.isroot():
      raise TreeError('Can\'t remove root node!')
    
    parent = node.parent
    
    if parent.isroot():
      if node.isleft():
        if parent.right.isinternal():
          # remove the right child as well, moving its pointers to root
          parent.left = parent.right.left
          parent.left.parent = parent
          parent.right = parent.right.right
          parent.right.parent = parent
          
        else: # if parent's right is external
          parent.left = None
          
      else: # node is right child
        if parent.left.isinternal():
          parent.right = parent.left.right
          parent.right.parent = parent
          parent.left = parent.left.left
          parent.left.parent = parent
        else:
          parent.right = None
          
      return
    
    # end of special cases
  
    if node.isleft():
      parent.right.parent = parent.parent
      parent.right.brlen += parent.brlen
      if parent.isleft():
        parent.parent.left = parent.right
      else:
        parent.parent.right = parent.right
        
    else:
      parent.left.parent = parent.parent
      parent.left.brlen += parent.brlen
      if parent.isleft():
        parent.parent.left = parent.left
      else:
        parent.parent.right = parent.left
      
  def addchild(self, child):
    """Add a child node to a node.

    If the node has no children, the child is added as left.  If the node has a 
    left child already, the child is added as right.  If the node already has 2
    children, an exception is raised.
    
    Parameters:
      child -- a TreeNode object to be added.
    
    """
    
    if self.left == None:
      self.left = child
    elif self.right == None:
      self.right = child
    else:
      raise TreeError('Attempted to add child to parent with two children.')
    
    child.parent = self
  
  def addsister(self, sister):
    """Add a node to a tree as a sister to a node.
    
    Raise TreeError if attempting to add a sister to the root.  A new parent
    node is created.
    
    Parameters:
      sister -- a TreeNode object
    
    """
    
    if self.isroot():
      raise TreeError('Root node can\'t have a sister!')
    
    if self.isleft():
      self.parent.left = None
    else:
      self.parent.right = None
       
    nodeptr = TreeNode()
    self.parent.addchild(nodeptr)
    nodeptr.addchild(self)
    nodeptr.addchild(sister)
    
  def clear(self):
    """Set right and left subtrees to None."""
    self.left = self.right = None
    
  def getnodes(self):
    """Return a list containing pointers to nodes in self's subtree."""
   
    nodelist = [self]
    if self.left: # has a left subtree
      nodelist.extend(self.left.getnodes())
     
    if self.right: # has a right subtree
      nodelist.extend(self.right.getnodes())
     
    return nodelist
  
class UnresolvedNode:
  """
  Node structure for an optionally multifurcating tree.
  
  This structure does NOT use Felsenstein's cyclic representation for internal
  nodes, but a list representation of children.

  """
  
  def __init__(self, name, bootstrap = 0, brlen = None):
    """Constructor for UnresolvedNode.

    Parameters:
      name -- string to be assigned as a label.
      bootstrap -- support value for the node, only relevant for internal nodes
      brlen -- length of the branch leading to the node (should be a float)
    
    """
    
    self.parent = None
    self.children = []
    self.name = name
    self.size = None
    self.bootstrap = bootstrap
    self.brlen = brlen
    
  def __len__(self):
    """Set the number of branches leading down from the node.
    
    This value corresponds to the number of children.
    
    Return:
      Number of children
      
    """
    
    if self.size == None:
      self.size = 1
      for child in self.children:
        self.size += len(child)
        
    return self.size
  
  def __iter__(self):
    """Iterate over the nodes in the tree starting with this node.
    
    Iteration is by preorder traversal of the tree.
    
    Yield:
      UnresolvedNode
    
    """

    yield self
    
    for child in self.children:
      for grandchild in child:
        yield grandchild
        
  def __str__(self):
    """Produce a Newick format string representation of a tree.
    
    Return:
      String representation.
    
    """

    if self == None:
      return
    
    if not self.isinternal():
      return self.name
    
    # This code was to indicate unresolved nodes with "{}"... it was irritating.
    #if len(self.children) > 2:
    #  obrace = '{'
    #  cbrace = '}'
    #else:
    obrace = '('
    cbrace = ')'
      
    treestr = obrace
    for child in self.children:
      treestr += str(child) + ','
    
    treestr = treestr[:-1] + cbrace
    

    if self.isroot():
      return treestr + ';'
      
    return treestr
        
  def addchild(self, child):
    """Add a child node to a node.

    Parameters:
      child -- a TreeNode object to be added.
    
    """

    self.children.append(child)
    child.parent = self
  
  def findnode(self, name):
    """Returns a TreeNode object with the given name field."""
    if self.name == name:
      return self

    found = None
    
    for node in self.children:
      found = node.findnode(name)
      if found:
        break
      
    return found

  def rank(self):
    """Find the name of the descendent leaf that comes first alphabetically.
    
    This method is used to sort the tree so that leaves are in alphabetical
    order (or as close to it as the topology allows) in the Newick string.
    
    Return:
      String
    
    """
    
    if self.children:
      return min([child.rank() for child in self.children])
   
    return self.name

  def sort(self):
    """Sort descendents such that names of leaves are in alphabetical order.
    
    Order will be as close to alphabetical as the topology allows.
    
    """
    
    if self.children:
      self.children = [(child.rank(), child) for child in self.children]
      self.children.sort()
      self.children = [child for (rank, child) in self.children]

      for child in self.children:
        child.sort()

  def removechild(self, child):
    """Remove the child node from the list of children.

    Raise TreeError if child is not among this node's children.
    
    Parameters:
      child -- node to remove
    
    """
    
    try:
      self.children.remove(child)
      child.parent = None
    except ValueError:
      raise TreeError('%s is not a child of %s!' % (child.name, self.name))
        
  def removenode(self, name):
    """Remove the node with the given name field from the tree below self.

    Should be called on a pointer to the root of a tree, which is searched for a
    node with the given name field.  If the node is not found, no action is
    taken.  If the node is found, it is removed, as is its parent node, if only 
    one sibling node remains.  This function does NOT deal appropriately with 
    branch lengths.
    
    Parameters:
      name -- name of the node to be removed (a string).
    
    """
    
    node = self.findnode(name)
    
    if node == None:
      return
    
    # special cases: node is root, only one child remains
    if node.isroot():
      raise TreeError('Can\'t remove root node!')
    
    parent = node.parent
    
    for child in node.children:
      parent.addchild(child)
      parent.children.append(child)
      child.parent = parent

    parent.removechild(node)
    
    # Only one child left? Don't need parent anymore (unless parent is root)
    if len(parent.children) == 1:
      if not parent.isroot():
        parent.parent.addchild(parent.children[0])
        parent.parent.removechild(parent)
    
  def collapse(self, bootlimit):
    """Collapse any internal nodes below self with low bootstrap values.
    
    Parameters:
      bootlimit -- bootstrap value below which a node is collapsed.
    
    """
    
    if not self.isinternal():
      return
    
    for child in self.children[:]:
      child.collapse(bootlimit)
      
    if not self.isroot() and self.bootstrap < bootlimit:
      for child in self.children:
        self.parent.addchild(child)
        
      self.parent.removechild(self)

  def setbrlen(self, brlen):
    """Change the node's branch length.
    
    Parameters:
      brlen -- new branch length (should be a float, but this is not tested).
    
    """
    
    self.brlen = brlen
        
  def davinprint(self, level = 0):
    """Print an indented list of nodes to indicate the structure of the tree.
    
    Parameters:
      level -- used in recursive calls to this function to indicated level of
               indentation (call with the default value!)
    
    """
    
    print ' ' * level * 4 + self.name
    
    for child in self.children:
      child.davinprint(level + 1)
      
  def isinternal(self):
    """Return True iff this node is not a leaf."""
    
    return len(self.children) > 0
  
  def isresolved(self):
    """Return True if this node has exactly 2 children."""
    return len(self.children) == 2
  
  def isroot(self):
    """Return True iff this node has no parent."""
    
    if self.parent:
      return False
    
    return True
  
  def reroot(self):
    """Reroot the tree such that this node is a child of the root."""
    
    if self.isroot():
      return
    
    self.parent.rotate(self)
    
  def resolutions(self):
    """Calculate number of bifurcating trees congruent with subtree below self.
    
    Return:
      int (number of trees)
    
    """
    
    tops = 1
    for i in xrange(2, len(self.children) + 1):
      tops *= (2 * i) - 3

    for child in self.children:
      tops *= child.resolutions()
            
    return tops
  
  def rotate(self, child):
    """Make self a child of the given child.
    
    Raise TreeError if child is not among this node's children.
    
    This function is used for re-rooting the tree.
    
    Parameters:
      child -- the node to be made parent of self (must be a child of self
               currently!)
    
    """
    
    if child not in self.children:
      raise TreeError('Can\'t rotate with a node that isn\'t a child!')

    if self.isroot():
      if len(self.children) > 2:
        child.addchild(self)
        self.removechild(child)
      else:
        for rootchild in self.children:
          if rootchild != child:
            child.addchild(rootchild)
        child.addchild(self)
        self.removechild(child)
          
    else:
      self.parent.rotate(self)
      child.addchild(self)
      self.removechild(child)
      
  def lastchild(self):
    """Return the child of self that appears last in list of children."""
    
    if len(self.children) == 0:
      raise TreeError('Node has no children!')
    
    return self.children[-1]
    
    
class TreeError(Exception):
  """Exception designed for tree-related errors."""
  
  def __init__(self, value = ''):
    """Constructor for TreeError exceptions."""
    self.value = value
      
  def __str__(self):
    """Returns a string representation of the exception."""
    return repr(self.value)
    
if __name__ == '__main__':
  
  # Some random test code follows. Ignore. Or add more.
  
  t = TreeNode()
  print t.isleft()
  t = '(A,(B,C));'
  print unroot(t)
  
  t = '((A:1.0,D:1.0):2.0,(B:1.0,C:1.0):1.0);'
  print t
  t2 = unroot(t)
  print t2
  print 'rooted:', root(t2)
  #t = '((A,B),(C,(D,E)),(F,G))'
  print getpartitions(t)
  
  
