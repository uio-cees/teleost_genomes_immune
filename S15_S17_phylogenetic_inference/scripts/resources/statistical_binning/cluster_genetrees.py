#!/usr/bin/env python

import sys
import os

global rep

GRAPHCOLORINGJAVA="../../../../scripts/resources/statistical_binning/graphColoring/code"

rep=""

def locate(pattern, root=os.curdir):
    '''Locate all files matching supplied filename pattern in and below
    supplied root directory.'''
    if not os.path.exists(root):
       raise RuntimeError ("path not found: %s" % root)    
    for path, dirs, files in os.walk(os.path.abspath(root)):
       for filename in fnmatch.filter(files, pattern):
           yield os.path.join(path, filename)

def gene_name(filename):
    return filename.split(".")[0]

def stat_filename(g,ths):
    return "%s%s.%d" %(rep,g,ths)

def subsets(genes,ths):
    genetoi = {}
    itogene = {}
    grows=[]
    for i, gn in enumerate(genes):
        genetoi[gn] = i
        itogene[i] = gn
    edges = 0
    for i in range(0,len(genes)):
       gn = itogene[i]
       row = ["0"]*len(genes) 
       #row = [] 
       f = stat_filename(gn,ths)
       nb = []
       for line in open(f,'r'):
           #print line
           r = line.split()
           if genetoi.has_key(r[1]):
               herid = genetoi[r[1]]
               incompatible = grows[herid][i] == "1" if herid < len(grows) else False
               if len(r) < 3:
                  print r
               incompatible = incompatible or int(r[2]) != 0 or int(r[3]) != 0
               row[herid] = "1" if incompatible else "0"
               edges += 1 if incompatible else 0
               if herid < len(grows):
                   grows[herid][i] = row[herid]
       grows.append(row)
    print "Total number of edges: ",edges
    if edges == 0:
        raise Exception("All genes are fully compatible. No binning required. Just run concatenation.")
    gr="p edge %d %d\n%s" %(len(genes),edges,'\n'.join(filter(lambda y: y!="",("\n".join(("e %d %d"%(j+1,i+1) for i,x in enumerate(r) if x!="0")) for j,r in enumerate(grows)))))
    #print gr
    from subprocess import Popen, PIPE, STDOUT
    p = Popen(['java','-Xmx2000M','-classpath',GRAPHCOLORINGJAVA,'GraphColoring'], stdout=PIPE, stdin=PIPE, stderr=None)
    stdout = p.communicate(input=gr)[0]
    print stdout
    res = [[itogene[int(y)] for y in x.split()] for x in stdout[1:-1].replace("||","|").split("|") if x.strip() !=""]
    return res

if __name__ == '__main__':
    thresholds = [int(sys.argv[2])]
    allgenenames=sys.argv[1]
    rep=sys.argv[3] if len(sys.argv) > 3 else ""
    maxsize = 1
    genesets = [filter(lambda x: x is not None and x.strip() != '', open(allgenenames).read().split('\n'))]
    for t in thresholds:
        newgenesets = []
        for gs in genesets:
            if len(gs) > maxsize:
                s = subsets(gs,t)
                print "Number of subsets:", len(s)
                newgenesets.extend(s)
            else:
                newgenesets.append(gs)
        genesets = newgenesets
        print "%d:\n%s" %(t,",".join((str(gs) for gs in genesets)))
    for i,gs in enumerate(genesets):
        fo = open(os.path.join(os.path.dirname(allgenenames),'bin.%d.txt' %i),'w')
        fo.write('\n'.join(gs))
        fo.write('\n')
        fo.close()
