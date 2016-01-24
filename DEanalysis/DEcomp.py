import pandas as pd
import glob
from collections import defaultdict
import operator
import itertools
import argparse
import os
import time

path = "/home/laraib/clust/DE_analysis/"
print ("Data is in dir: " + path)

parser = argparse.ArgumentParser()
parser.add_argument("method")
args = parser.parse_args()
method = args.method.lower()

with open("/home/laraib/clust/mappingData_avi/contig2cuffGene.txt", 'r') as f:
    data = pd.read_table(f, names = ['contig', 'cuffgene'])
    table = data.set_index("contig").to_dict()['cuffgene']

if (method == "sailfish"):
	clust2gene = {}
	clustID = 1

	with open("/mnt/scratch3/avi/clustering/geneModel/src/quant_human.clust", 'r') as f:
		for line in f:
			nomap = 0
			geneIDs = []
			contigIDs = line.split('\t')
			for ID in contigIDs:
				ID = ID.strip('\n')
				try:
					geneIDs.append((table[ID]))
				except KeyError:
					nomap+=1
			if (nomap == len(contigIDs)):
				clust2gene[("cluster"+str(clustID))] = 'None'
			else:
				d = defaultdict(int)
				for i in geneIDs:
					d[i] += 1
				result = max(d.iteritems(), key=lambda x: x[1])
				clust2gene[("cluster"+str(clustID))] = result[0]
			clustID += 1
elif (method == "corset"):
	clust2gene = {}
	with open("/mnt/scratch3/avi/clustering/data/corsetData/Human-Trinity/corset-clusters.txt", 'r') as f:
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
						geneIDs.append((table[ID]))
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
				geneIDs.append((table[ID]))
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

with open((path + method + "padj.txt"), 'r') as f:
    data = pd.read_table(f, names = ['clust', 'pval'])
    clust2pval = data.set_index("clust").to_dict()['pval']

sortedclusts = sorted(clust2pval, key=lambda key: clust2pval[key]) #sorted by p value

with open((path + "sigGenes_Human_Genome.txt"), 'r') as f:
    sigGenes = f.read().splitlines()

print ("Count of true positives for " + method)
print ("Top ranked clusters" + "\t" + "Uniq True positives")
tp = []
for i in range(0,35001):
	if (i%5000 == 0):
		print (str(i) + '\t' + str(len(set(tp))))
	sigclustgene = clust2gene[sortedclusts[i]]
	if sigclustgene in sigGenes:
		tp.append(sigclustgene)

print (str(i) + '\t' + str(len(set(tp))))


