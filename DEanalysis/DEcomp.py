import pandas as pd
import glob
from collections import defaultdict
import operator
import itertools
import argparse
import os
import time
import click

def loadResults2(contigToGene, fname):
    clust2genes = {}
    with open(fname, 'r') as f:
        for l in f:
            toks = l.rstrip().split()
            contig, clust = toks[0], toks[1]
            if clust not in clust2genes:
                clust2genes[clust] = set([])
            if contig in contigToGene:
                clust2genes[clust].add(contigToGene[contig])
    return clust2genes

def loadResults(contigToGene, fname):
    clust2gene = {}
    with open(fname, 'r') as f:
        line = f.readline()
        curclustID = line.split('\t')[1].strip('\n')
        contigIDs = [line.split('\t')[0]]
        for line in f:
            if (line.split('\t')[1].strip('\n') == curclustID):
                contigIDs.append(line.split('\t')[0])
            else:
                nomap = 0
                geneIDs = []
                for ID in contigIDs:
                    ID = ID.strip('\n')
                    try:
                        geneIDs.append((contigToGene.loc[ID, "cuffgene"]))
                    except KeyError:
                        nomap+=1
                if (nomap == len(contigIDs)):
                    clust2gene[curclustID] = 'None'
                else:
                    d = defaultdict(int)
                    for i in geneIDs:
                        d[i] += 1
                    result = max(d.iteritems(), key=lambda x: x[1])
                    clust2gene[curclustID] = result[0]
                curclustID = line.split('\t')[1].strip('\n')
                contigIDs = [line.split('\t')[0]]
        nomap = 0
        geneIDs = []
        for ID in contigIDs:
            ID = ID.strip('\n')
            try:
                geneIDs.append((contigToGene.loc[ID, "cuffgene"]))
            except KeyError:
                nomap+=1
        if (nomap == len(contigIDs)):
            clust2gene[curclustID] = 'None'
        else:
            d = defaultdict(int)
            for i in geneIDs:
                d[i] += 1
            result = max(d.iteritems(), key=lambda x: x[1])
            clust2gene[curclustID] = result[0]
    return clust2gene



@click.command()
@click.option("--method", default="sailfish", help="method to generate results for")
@click.option("--adir", help="analysis directory containing the clustering files and adjusted p-value files")
def genTPRate(method, adir):
    import os

    sigCutoff = 0.05
    path = adir
    print ("Data is in dir: " + path)

    usingSets = True
    contigToGeneFile = os.path.sep.join([path, "contig2cuffGene.tsv"])
    contigToGene = pd.read_table(contigToGeneFile, names=["contigs", "cuffgene"])
    contigToGene.set_index("contigs", inplace=True)
    contigToGene = contigToGene.to_dict()['cuffgene']

    if method == "sailfish":
        clust2gene = loadResults2(contigToGene, os.path.sep.join([path, "rapclust_clusters.flat"]))
    elif method == "corset":
        clust2gene = loadResults2(contigToGene, os.path.sep.join([path, "corset-clusters.txt"]))

    #print(clust2gene.keys())
    with open(os.path.sep.join([path,  method + "padj.txt"]), 'r') as f:
        data = pd.read_table(f, names = ['clust', 'pval'])
        clust2pval = data.set_index("clust").to_dict()['pval']

    sortedclusts = sorted(clust2pval, key=lambda key: clust2pval[key]) #sorted by p value

    ##
    # Find all genes that are DE with p-adj < 0.05 using the
    # true txp <-> gene mapping
    ##
    sigGenes = set([])
    with open(os.path.sep.join([path, "truthpadj.txt"]), 'r') as f:
        for l in f:
            toks = l.split('\t')
            if float(toks[1]) <= sigCutoff:
                sigGenes.add(toks[0])

    print ("Count of true positives for " + method)
    print ("Top ranked clusters" + "\t" + "Uniq True positives" + '\t' + "Uniq False positives")
    tp = set([]) # use a set so we don't count duplicates
    fp = set([])
    for i in range(len(sortedclusts)):
        if (i%500 == 0):
            print (str(i) + '\t' + str(len(tp)) + '\t' + str(len(fp)) + '\t' + str(clust2pval.get(sortedclusts[i])))
        #sigclustgene = clust2gene[sortedclusts[i]]
        #if sigclustgene in sigGenes:
        #    tp.add(sigclustgene)
        #else:
        #   fp.add(sigclustgene)
        #sigClustGenes = clust2gene[sortedclusts[i]]
        if usingSets:
            sigClustGenes = clust2gene[sortedclusts[i]]
            sigSet = sigGenes.intersection(sigClustGenes)
            if len(sigSet) > 0:
                if (len(sigSet.intersection(fp)) > 0):
                    print("contains some things already FP")
                tp |= (sigSet - fp)
            else:
                if (len(sigClustGenes.intersection(tp)) > 0):
                    print("contains some things already TP")
                fp |= (sigClustGenes - tp)
        if (clust2pval.get(sortedclusts[i]) > sigCutoff):
            break
    print (str(i) + '\t' + str(len(tp)) + '\t' + str(len(fp)) + '\t' + str(clust2pval.get(sortedclusts[i])))
    print ("precision " + str(float(len(tp))/(len(tp)+len(fp))))
    print ("recall " + str(float(len(tp))/len(sigGenes)))
if __name__ == "__main__":
    genTPRate()
