#!/lusr/bin/python
'''
Created on Jun 3, 2013

@author: smirarab
'''
import dendropy
import sys
import os
import copy

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
    
if __name__ == '__main__':

    if (len(sys.argv) < 2):
        print "USAGE: %s tree_file [threshold - default 75] [outfile name; - uses default] [-strip-internal|-strip-bl|strip-both|-nostrip; default: nostrip]" %sys.argv[0]

    treeName = sys.argv[1]            
    t = 75 if len (sys.argv) < 3 else float(sys.argv[2])
    resultsFile="%s.%d" % (treeName,t) if len (sys.argv) < 4 or sys.argv[3]=="-" else sys.argv[3]
    #print "outputting to", resultsFile    
    strip_internal=True if len (sys.argv) > 4 and ( sys.argv[4]=="-strip-internal" or sys.argv[4]=="-strip-both" ) else False 
    strip_bl=True if len (sys.argv) > 4 and ( sys.argv[4]=="-strip-bl" or sys.argv[4]=="-strip-both" ) else False
    
    trees = dendropy.TreeList.get_from_path(treeName, 'newick')
    filt = lambda edge: False if (edge.label is None or (is_number(edge.label) and float(edge.label) >= t)) else True
    for tree in trees:
        for n in tree.internal_nodes():
            if n.label is not None:
                n.label = float (n.label)
                n.edge.label = n.label
        edges = tree.get_edge_set(filt)
        for e in edges:
            e.collapse()
        if strip_internal:
            for n in tree.internal_nodes():
                n.label = None
        if strip_bl:
            for e in tree.get_edge_set():
                e.length = None
        
    trees.write(open(resultsFile,'w'),'newick',write_rooting=False)
