# -*- coding: utf-8 -*-
"""
Created on Mon Jan 26 21:45:04 2015
@author: bfetler

Synteny Blocks in Genetics:

Genomic sequencing and mapping have enabled comparison of the
general structures of genomes of many different species. The 
general finding is that organisms of relatively recent divergence 
show similar blocks of genes in the same relative positions in
the genome. This situation is called synteny, translated roughly 
as possessing common chromosome sequences. For example, many of 
the genes of humans are syntenic with those of other mammals — 
not only apes but also cows, mice, and so on. Study of synteny 
can show how the genome is cut and pasted in the course of 
evolution.
                --from https://en.wikipedia.org/wiki/Synteny

The basic question is how to map synteny blocks between
two species.  By connecting the 5' and 3' ends of the gene, 
we create a circular gene representation for each species.  
By comparing the regions that swapped during mutations, we 
identify large blocks of swapped regions.  We can create a 
Graph between two species, identifying connected components, 
and therefore possible swaps that ocurred during evolution.  

The exercise here is to create an algorithm to identify
the gene graph and connected components.

"""

import re    # regex
import datetime
import sys

class GeneGraph:
    '''Gene sequence graph, from DiGraph, using gene synteny blocks'''
#    vn = 0  # number of vertices
#    en = 0  # number of edges
#    adj = []  # adjacency list of lists
    
    def __init__(self, v, nam=''):
        if v < 0:
            raise TypeError, 'cannot be less than zero'
        self.vn = v
        self.en = 0
        self.name = nam
        self.adj = []
#        self.adj = [[]] * v
        for i in range(v):
            self.adj.append([])

    def validateVertex(self, v):
        if (v < 0 or v >= self.vn):
            raise TypeError, 'invalid vertex number'
    
    def addEdge(self, v, w):
#        self.validateVertex(self, v)
        self.validateVertex(v)
        self.validateVertex(w)
        self.en += 1
        self.adj[v].append(w)
        self.adj[w].append(v)  # digraph or not
        
    def degree(self, v):
        self.validateVertex(v)
        return len(self.adj[v])
        
    def printConnections(self):
        print self.name, 'connections:',
        for i, a in enumerate(self.adj):
            for b in a:
                print i, '->', b, ',',
        print ''
        
class ConnComp:    # connected components
    '''connected components of a gene graph'''
    def __init__(self, gGraph):
        self.count = 0
        self.name = gGraph.name
        sz = gGraph.vn
        self.marked = [False] * sz
        self.id = [-1] * sz
        self.size = [0] * sz
#        self.marked = []
#        self.id = []
#        self.size = []
#        for v in range(sz):
#            self.marked.append(False)
#            self.id.append(-1)
#            self.size.append(0)
        for v in range(sz):
            if (not self.marked[v]):
                self.dfs(gGraph, v)
                self.count += 1
    
    def dfs(self, gGraph, v):
        '''depth first search'''
#        print 'called dfs'
        self.marked[v] = True
        self.id[v] = self.count
        self.size[self.count] += 1
        for w in gGraph.adj[v]:
            if (not self.marked[w]):
                self.dfs(gGraph, w)
    
    def printConn(self):
        print self.name, 'connComp count', self.count,' marked'
        for i, a in enumerate(self.marked):
            print i, ':', a, ',',
        print ''
        print 'connComp id',
        for i, a in enumerate(self.id):
            print i, ':', a, ',',
        print ''

class ConnComp2:    # connected components
    '''connected components of a gene graph'''
    def __init__(self, gGraph):
        self.count = 0
        self.name = gGraph.name
        sz = gGraph.vn
        self.marked = [False] * sz
        self.id = [-1] * sz
        self.size = [0] * sz
        for v in range(sz):
            if (not self.marked[v]):
                vlist = [v]
                while (vlist != []):
#                    print 'init processing', vlist, 'count', self.count
                    vlist = self.dfs2(gGraph, vlist)
                self.count += 1

# how to do dfs() without recursion?  while loop, for loop?
# dfs() returns node or nodes?
    
    def dfs2(self, gGraph, vlist):
        '''depth first search'''
#        print 'called dfs'
        wlist = []
        for v in vlist:
            self.marked[v] = True
            self.id[v] = self.count
            self.size[self.count] += 1
            for w in gGraph.adj[v]:
                if (not self.marked[w]):
                    wlist.append(w)
        return wlist
    
    def printConn(self):
        print self.name, 'connComp count', self.count,' marked'
        for i, a in enumerate(self.marked):
            print i, ':', a, ',',
        print ''
        print 'connComp id',
        for i, a in enumerate(self.id):
            print i, ':', a, ',',
        print ''
        print self.name, 'connComp count', self.count
    
    def printSummary(self):
        print self.name, 'connComp count', self.count

def readfile(filename):
    '''read dna file, return string'''
    f = open(filename,'r')
    s1 = f.readline()
    s2 = f.readline()
    return s1, s2
    
def rmnull(x):
    return x <> ''

def getsz(r1):
    sz = 0
    for r in r1:
        for s in r.split():
            sz += 1
    return sz

def addEdges(gene, ra):
    for r in ra:
        rs = r.split()
        q1 = []
        q2 = []
        for s in rs:
            ii = int(s)
            if ii > 0:
                d1 = ii * 2 - 2
                d2 = ii * 2 - 1
            else:
                d1 = -ii * 2 - 1
                d2 = -ii * 2 - 2
            q1.append(d1)
            q2.append(d2)
        q0 = q1.pop(0)  # rotate q1[0] from start to end
        q1.append(q0)
#        q0 = q2.pop(0)  # rotate q2[0] from start to end
#        q2.append(q0)
#        for qq in zip(q1, q2):
#            print qq[0], qq[1], ';',
#        print ''
        
        for qq in zip(q1, q2):
            gene.addEdge(qq[1], qq[0])


dtime = datetime.datetime.now()
print 'start time', dtime

geneA = GeneGraph(12, 'geneA')
# add connections
geneA.addEdge(1, 2)  # index-1
geneA.addEdge(3, 4)
geneA.addEdge(5, 6)
geneA.addEdge(7, 8)
geneA.addEdge(9, 10)
geneA.addEdge(11, 0)
# standard blocks
geneA.addEdge(0, 1)
geneA.addEdge(2, 3)
geneA.addEdge(4, 5)
geneA.addEdge(6, 7)
geneA.addEdge(8, 9)
geneA.addEdge(10, 11)
geneA.printConnections()
compA = ConnComp(geneA)
compA.printConn()

geneB = GeneGraph(12, 'geneB')
# add connections
geneB.addEdge(1, 5)
geneB.addEdge(4, 11)
geneB.addEdge(10, 9)
geneB.addEdge(8, 0)
geneB.addEdge(3, 7)
geneB.addEdge(6, 2)
# standard blocks
geneB.addEdge(0, 1)
geneB.addEdge(2, 3)
geneB.addEdge(5, 4)
geneB.addEdge(7, 6)
geneB.addEdge(9, 8)
geneB.addEdge(11, 10)
geneB.printConnections()
compB = ConnComp(geneB)
compB.printConn()

# breakpoint graph of geneA, geneB
geneAB = GeneGraph(12, 'geneAB')
# add geneA connections
geneAB.addEdge(1, 2)
geneAB.addEdge(3, 4)
geneAB.addEdge(5, 6)
geneAB.addEdge(7, 8)
geneAB.addEdge(9, 10)
geneAB.addEdge(11, 0)
# add geneB connections
geneAB.addEdge(1, 5)
geneAB.addEdge(4, 11)
geneAB.addEdge(10, 9)
geneAB.addEdge(8, 0)
geneAB.addEdge(3, 7)
geneAB.addEdge(6, 2)
geneAB.printConnections()
compAB = ConnComp(geneAB)
compAB.printConn()

# how to do this from file read?
# e.g. file 1st line   (+1 -3 -6 -5)(+2 -4)
# compare with 2nd line (+1 +2 +3 +4 +5 +6)

rx = re.compile('[()\n]')

print 'starting gene read'
# s1, s2 = readfile('/Users/bfetler/Downloads/data_gene1.txt')
s1, s2 = readfile('/Users/bfetler/Downloads/dataset_288_4.txt')
r1 = filter(rmnull, rx.split(s1))
r2 = filter(rmnull, rx.split(s2))
sz1=getsz(r1)
sz2=getsz(r2)

# print 's1', s1, r1, 'sz', sz1
# print 's2', s2, r2, 'sz', sz2
print 'sizes:', sz1, sz2

geneC = GeneGraph(2 * sz2, 'geneC')
addEdges(geneC, r1)
addEdges(geneC, r2)
# geneC.printConnections()
compC = ConnComp2(geneC)  # max recursion depth exceeded
compC.printSummary()

etime = datetime.datetime.now()
print 'end time', etime
print 'time delta', etime-dtime
print 'long loop'
cc = 0
for c in [0,1000]:
    s1, s2 = readfile('/Users/bfetler/Downloads/dataset_288_4.txt')
    cc += 1
ltime = datetime.datetime.now()
print 'loop time', ltime-etime
print 'python version', sys.version
