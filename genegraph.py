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

Ref: Bioinformatics Algorithms Volume 1, P. Compeau and
P. Pevzner, Ch. 6 "Are There Fragile Regions in the
Human Genome?", Active Learning Publishers, 2015.  

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
        '''create graph w/ number of vertices v'''
        if v < 0:
            raise TypeError, 'cannot be less than zero'
        self.vn = v
        self.en = 0
        self.name = nam
        self.adj = []
        for i in range(v):
            self.adj.append([])

    def validateVertex(self, v):
        if (v < 0 or v >= self.vn):
            raise TypeError, 'invalid vertex number'
    
    def addEdge(self, v, w):
        '''add edge to graph'''
        self.validateVertex(v)
        self.validateVertex(w)
        self.en += 1
        self.adj[v].append(w)
        self.adj[w].append(v)  # comment out if not digraph
        
    def degree(self, v):
        self.validateVertex(v)
        return len(self.adj[v])
        
    def printConnections(self):
        print self.name, 'connections:',
        for i, a in enumerate(self.adj):
            for b in a:
                print i, '->', b, ',',
        print ''
        
class ConnComp:    # connected components, recursive dfs
    '''connected components of a gene graph'''
    # may exceed max recursion depth
    def __init__(self, gGraph):
        self.count = 0
        self.name = gGraph.name
        sz = gGraph.vn
        self.marked = [False] * sz
        self.id = [-1] * sz
        self.size = [0] * sz
        for v in range(sz):
            if (not self.marked[v]):
                self.dfs(gGraph, v)
                self.count += 1
    
    def dfs(self, gGraph, v):    # call stack too deep
        '''depth first search, recursive loop'''
        self.marked[v] = True
        self.id[v] = self.count
        self.size[self.count] += 1
        for w in gGraph.adj[v]:
            if (not self.marked[w]):
                self.dfs(gGraph, w)
    
    def printConn(self):
        print self.name, 'connected components count', self.count
        for i, a in enumerate(self.marked):
            if (a == False):
                print 'unmarked', i, ':', a
        print 'connComp id',
        for i, a in enumerate(self.id):
            print i, ':', a, ',',
        print ''

class ConnComp2:    # connected components, non-recursive dfs
    '''connected components of a gene graph'''
    def __init__(self, gGraph):
        self.count = 0
        self.name = gGraph.name
        sz = gGraph.vn
        self.marked = [False] * sz
        self.id = [-1] * sz
        self.size = [0] * sz
        self.lists = []
        self.listct = []
        for v in range(sz):
            if (not self.marked[v]):
                vct = 0
                vlist = [v]
                self.lists.append(vlist[0])
                while (vlist != []):
                    vlist = self.dfs2(gGraph, vlist)
                    if (vlist != []):
                        vct += 1
                self.listct.append(vct)
                self.count += 1

# how to do dfs() without recursion?  while loop, for loop?
# is this correct?
    
    def dfs2(self, gGraph, vlist):
        '''depth first search, non-recursive'''
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
        print self.name, 'connComp count', self.count
        for i, a in enumerate(self.marked):
            if (a == False):
                print 'unmarked', i, ':', a
        print 'connComp id',
        for i, a in enumerate(self.id):
            print i, ':', a, ',',
        print ''
        print self.name, 'connComp count', self.count
    
    def printSummary(self):
        print self.name, 'connected components count', self.count, 
        print '; start nodes', self.lists, '; comp sizes', self.listct
#       how to graph the connected components?

def processDnaFile(geneName, filename):
    '''process DNA file'''
    # nested methods
    
    rx = re.compile('[()\n]')  # split input file by regex
    
    def rmnull(x):
        return x <> ''
        
    def getsz(r1):
        '''get size'''
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
#            q0 = q2.pop(0)  # rotate q2[0] from start to end
#            q2.append(q0)
            for qq in zip(q1, q2):
                gene.addEdge(qq[1], qq[0])    
    
    print 'start reading', geneName
    f = open(filename,'r')
    s1 = f.readline()
    s2 = f.readline()
    r1 = filter(rmnull, rx.split(s1))
    r2 = filter(rmnull, rx.split(s2))
    sz1 = getsz(r1)
    sz2 = getsz(r2)
    print 'dna sizes:', sz1, sz2
    
    gene = GeneGraph(2 * sz2, geneName)
    addEdges(gene, r1)
    addEdges(gene, r2)
    return gene

def compGeneA():
    print 'geneA: one DNA strand, six linear blocks'
    geneA = GeneGraph(12, 'geneA')
    # add connections
    geneA.addEdge(1, 2)  # index+1 mod 12
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

def compGeneB():
    print 'geneB: 2 DNA strands length 4 & 2, construct as one'
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

def compGeneAB():
    print 'breakpoint graph of geneA, geneB'
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

# Above is a very manual process.  
# How to read and process a file?  processDnaFile()
#   e.g. 1st file line (+1 -3 -6 -5)(+2 -4)
#        2nd file line (+1 +2 +3 +4 +5 +6)

def compFileGeneC():
    geneC = processDnaFile('geneC','./data/data_gene1.txt')
    geneC.printConnections()
    compC = ConnComp(geneC)     # dfs with recursion
    compC.printConn()
    compC2 = ConnComp2(geneC)   # dfs without recursion
    compC2.printSummary()

def compFileGeneD():
    geneD = processDnaFile('geneD','./data/dataset_288_4.txt')
    # geneD.printConnections()  # may be very long
    # compD = ConnComp(geneD)   # max recursion depth exceeded
    compD = ConnComp2(geneD)    # dfs without recursion
    compD.printSummary()
    # compD.printConn()

def main():
    compGeneA()
    compGeneB()
    compGeneAB()
    compFileGeneC()
    compFileGeneD()

if __name__ == '__main__':
    main()

