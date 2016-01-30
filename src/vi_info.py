import pandas as pd
from pyfasta import Fasta
from igraph import Clustering, compare_communities
import click
import sys

mclClustFile = "/home/rob/RapClust/DEanalysis/human_clusts/human_sum.clust"
#mclClustFile = sys.argv[1]
cdhit8File = "/mnt/scratch3/avi/clustering/human/trinity/cd_hit/0.8.clstr"
cdhit95File = "/mnt/scratch3/avi/clustering/human/trinity/cd_hit/0.95.clstr"
origClustFile = "/mnt/scratch3/avi/clustering/human/trinity/truth/contigs2genes.disambiguous.txt"
corsetClustFile = "/mnt/scratch3/avi/clustering/human/trinity/corset/corset-clusters.txt"

def read_original(filename):
    with open(filename) as f:
       data = pd.read_table(f, header=None, names=['contig', 'gene'])
       knownContigs = data['contig'].tolist()
       kcontig2Ind = {}
       for ind, contig in enumerate(knownContigs):
           kcontig2Ind[contig] = ind

       origClust = {k: list(v) for k,v in data.groupby("gene")["contig"]}
       origClustInd = [None] * len(knownContigs)
       clusterId = 0
       for _, val in origClust.items():
           for contig in val:
               origClustInd[kcontig2Ind[contig]] = clusterId
           clusterId += 1
    return (origClustInd, knownContigs, kcontig2Ind)

def read_corset(filename, knownContigs, kcontig2Ind):
    knownSet = set(knownContigs)
    with open(filename) as f:
        data = pd.read_table(f, header=None, names=['contig', 'cluster'])
        corsetClust = {k: list(v) for k,v in data.groupby("cluster")["contig"]}

        corsetClustInd = [None] * len(knownContigs)
        clusterId = 0
        for _, val in corsetClust.items():
            for contig in val:
                if contig in knownSet:
                    corsetClustInd[kcontig2Ind[contig]] = clusterId
            clusterId += 1
    return corsetClustInd

def read_mcl(filename, knownContigs, kcontig2Ind):
    knownSet = set(knownContigs)
    with open(filename) as f:
        mclClustInd = [None]*len(knownContigs)
        clusterId = 0
        for line in f:
            contigList = line.strip().split("\t")
            for contig in contigList:
                if contig in knownSet:
                    mclClustInd[kcontig2Ind[contig]] = clusterId
            clusterId += 1

    NoId = [ind for ind, x in enumerate(mclClustInd) if x == None]
    for ind in NoId:
        mclClustInd[ind] = clusterId
    return mclClustInd

def read_cdhit(filename, knownContigs, kcontig2Ind):
    knownSet = set(knownContigs)
    with open(filename) as f:
        cdhitClustInd = [None]*len(knownContigs)
        line = f.readline()
        clusterId = 0
        while(line):
            line = f.readline()
            while(line and line[0] != '>'):
                contig = line.strip().split("\t")[1].split(" ")[1].replace(">","").replace("...","")
                if contig in knownSet:
                    cdhitClustInd[kcontig2Ind[contig]] = clusterId
                line = f.readline()
            clusterId += 1
    NoId = [ind for ind, x in enumerate(cdhitClustInd) if x == None]
    for ind in NoId:
        cdhitClustInd[ind] = clusterId
    return cdhitClustInd

@click.command()
@click.option("--truth", help="true clustering file")
@click.option("--pred", type=(str, str), multiple=True, help="predicted clustering file and format {rap|corset|cdhit}")
def main(truth, pred):
    print("Reading true clustering")
    origClustInd, knownContigs, kcontig2Ind = read_original(truth)

    readDict = {"corset" : read_corset,
                "rap"    : read_mcl,
                "cdhit"  : read_cdhit}
    for p in pred:
        print("Reading {} clustering {}".format(p[1], p[0]))
        predInds = readDict[p[1]](p[0], knownContigs, kcontig2Ind)
        print ("truth vs {} : {}".format(p[1], compare_communities(Clustering(origClustInd), Clustering(predInds), method='vi')))

    #print "Reading corset"
    #corsetClustInd = read_corset(corsetClustFile, knownContigs, kcontig2Ind)
    #print "Reading MCL"
    #mclClustInd = read_mcl(mclClustFile, knownContigs, kcontig2Ind)
    #print "Reading cdhit8"
    #cdhit8ClustInd = read_cdhit(cdhit8File, knownContigs, kcontig2Ind)
    #print "Reading cdhit95"
    #cdhit95ClustInd = read_cdhit(cdhit95File, knownContigs, kcontig2Ind)
    #print "checking"
    #print ("Orig/corset {}".format(compare_communities(Clustering(origClustInd), Clustering(corsetClustInd), method='vi')))
    #print ("Orig/mcl {}".format(compare_communities(Clustering(origClustInd), Clustering(mclClustInd), method='vi')))
    #print ("mcl/corset {}".format(compare_communities(Clustering(mclClustInd), Clustering(corsetClustInd), method='vi')))
    #print ("mcl/cdhit8 {}".format(compare_communities(Clustering(mclClustInd), Clustering(cdhit8ClustInd), method='vi')))
    #print ("mcl/cdhit95 {}".format(compare_communities(Clustering(mclClustInd), Clustering(cdhit95ClustInd), method='vi')))
    #print ("orig/cdhit8 {}".format(compare_communities(Clustering(origClustInd), Clustering(cdhit8ClustInd), method='vi')))
    #print ("orig/cdhit95 {}".format(compare_communities(Clustering(origClustInd), Clustering(cdhit95ClustInd), method='vi')))


if __name__ == "__main__":
    main()


