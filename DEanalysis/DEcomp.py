import pandas as pd
import glob
from collections import defaultdict
import operator
import itertools
import argparse
import os
import time
import click

def loadResults2(method, contigToGene, fname):
    clust2genes = {}
    with open(fname, 'r') as f:
        for l in f:
            toks = l.rstrip().split()
            if (method == "sailfish"):
                contig, clust = toks[1], toks[0]
            else:
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
        if (method == "sailfish"):
            curclustID = line.split('\t')[0].strip('\n')
            contigIDs = [line.split('\t')[1]].strip('\n')
        else:
            curclustID = line.split('\t')[1].strip('\n')
            contigIDs = [line.split('\t')[0]].strip('\n')
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
@click.option("--clustfile", help="path to the flat cluster file")
@click.option("--contig2gene", help="path to the flat contig2Cuffgene file")
@click.option("--adir", help="analysis directory containing the adjusted p-value files")
@click.option("--odir", default='.', help="directory where results will be written")
def genTPRate(method, clustfile, contig2gene, adir, odir):
    import os

    sigCutoff = 0.05
    path = adir
    print ("Data is in dir: " + path)

    usingSets = True
    contigToGeneFile = contig2gene
    contigToGene = pd.read_table(contigToGeneFile, names=["contigs", "cuffgene"])
    contigToGene.set_index("contigs", inplace=True)
    contigToGene = contigToGene.to_dict()['cuffgene']

    if method == "sailfish":
        clust2gene = loadResults2(method, contigToGene, clustfile)
    elif method == "corset":
        clust2gene = loadResults2(method, contigToGene, clustfile)
    elif method == "cdhit":
        clust2gene = loadResults2(method, contigToGene, clustfile)

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
    nonSigGenes = set([])
    with open(os.path.sep.join([path, "truthpadj.txt"]), 'r') as f:
        for l in f:
            toks = l.split('\t')
            if float(toks[1]) <= sigCutoff:
                sigGenes.add(toks[0])
            else:
                nonSigGenes.add(toks[0])

    # DESeq2 sometimes outputs p-values of NA,
    # but we drop these from the truthpadj.txt file.
    # Here, we consider any genes missing from truthadj.txt
    # as *not* differentially expressed.
    for k, v in contigToGene.iteritems():
        if v not in sigGenes:
            nonSigGenes.add(v)

    print ("Count of true positives for " + method)
    print ("Top ranked clusters" + "\t" + "Uniq True positives" + '\t' + "Uniq False positives")
    tp = set([]) # use a set so we don't count duplicates
    fp = set([])
    ofile = open(os.path.sep.join([odir, "{}_de_recovery.scored-label".format(method)]), 'w')
    for i in range(len(sortedclusts)):
        if (i%500 == 0):
            print (str(i) + '\t' + str(len(tp)) + '\t' + str(len(fp)) + '\t' + str(clust2pval.get(sortedclusts[i])))
        truePos = False
        falsePos = False
        newTP = set([])
        newFP = set([])
        if usingSets:
            sigClustGenes = clust2gene[sortedclusts[i]]
            sigSet = sigGenes.intersection(sigClustGenes)
            nonSigSet = nonSigGenes.intersection(sigClustGenes)
            assert(sigSet.union(nonSigSet) == sigClustGenes)
            if len(sigSet) > 0:
                if (len(sigSet.intersection(fp)) > 0):
                    print("contains some things already FP")
                newTP = (sigSet - fp) - tp
                tp |= newTP
                truePos = len(newTP) > 0
            if len(nonSigSet) > 0:
                if (len(nonSigSet.intersection(tp)) > 0):
                    print("contains some things already TP")
                newFP = (nonSigSet - tp) - fp
                fp |= newFP
                falsePos = len(newFP) > 0
        else: # Not using sets
            sigclustgene = clust2gene[sortedclusts[i]]
            if sigclustgene in sigGenes:
               tp.add(sigclustgene)
            else:
               fp.add(sigclustgene)
            sigClustGenes = clust2gene[sortedclusts[i]]

        # Only write output if the # of true positives or false positives
        # changed.  Unlabeled points are ignored
        if truePos or falsePos:
            # Every true positive is an example labled with 1
            for x in newTP:
                ofile.write("{}\t{}\n".format(-clust2pval[sortedclusts[i]], 1))
            # Every false positive is an example labled with 0
            for x in newFP:
                ofile.write("{}\t{}\n".format(-clust2pval[sortedclusts[i]], 0))
        if (clust2pval.get(sortedclusts[i]) > sigCutoff):
            break
    ofile.close()
    print (str(i) + '\t' + str(len(tp)) + '\t' + str(len(fp)) + '\t' + str(clust2pval.get(sortedclusts[i])))
    print ("precision " + str(float(len(tp))/(len(tp)+len(fp))))
    print ("recall " + str(float(len(tp))/len(sigGenes)))

if __name__ == "__main__":
    genTPRate()
