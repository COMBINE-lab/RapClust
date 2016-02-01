import itertools
#import hirak
trinity_human = "/mnt/scratch3/avi/clustering/data/human/trinity/truth/contigs2genes.disambiguous.txt"
trinity_yeast = "/mnt/scratch3/avi/clustering/data/yeast/trinity/truth/contigs2genes.disambiguous.txt"
trinity_chicken = "/mnt/scratch3/avi/clustering/data/chicken/trinity/truth/contigs2genes.disambiguous.txt"
class Classification:
    TruePos, FalsePos, TrueNeg, FalseNeg = range(4)

def classType(true1, true2, pred1, pred2):
    if true1 == true2:
        if pred1 == pred2:
            return Classification.TruePos
        else: # truely the same, predicted different
            return Classification.FalseNeg
    else: # truly different
        if pred1 == pred2: #predicted same
            return Classification.FalsePos
        else:
            return Classification.TrueNeg


def accuracyExpressed(groundTruth_clust, tr_clust):
    #count true postive for each pair of transcripts O(N^2)
    tp, fp, tn, fn = 0, 0, 0, 0
    for tr_1, tr_2 in itertools.combinations(tr_clust.keys(), 2):
        if tr_1 not in groundTruth_clust or tr_2 not in groundTruth_clust:
            continue
        ct = classType(groundTruth_clust[tr_1], groundTruth_clust[tr_2], tr_clust[tr_1], tr_clust[tr_2])
        if ct == Classification.TruePos:
            tp += 1
        elif ct == Classification.TrueNeg:
            tn += 1
        elif ct == Classification.FalsePos:
            fp += 1
        elif ct == Classification.FalseNeg:
            fn += 1
    return tp, fp, tn, fn


def accuracyExpressedFast(groundTruth_clust, groundTruth_clust_inv,
                          tr_clust, tr_clust_inv):
    #num = len(set(tr_clust.keys()) & set(groundTruth_clust.keys()))
    num = len(set(groundTruth_clust.keys()))
    tp, fp, tn, fn = 0, 0, 0, 0
    for clustName, clustMems in tr_clust_inv.iteritems():
        for tr_1, tr_2 in itertools.combinations(clustMems,2):
            if tr_1 not in groundTruth_clust or tr_2 not in groundTruth_clust:
                continue
            if groundTruth_clust[tr_1] == groundTruth_clust[tr_2]:
                tp += 1
            else:
                fp += 1
    for clustName, clustMems in groundTruth_clust_inv.iteritems():
        for tr_1, tr_2 in itertools.combinations(clustMems,2):
            if tr_1 not in tr_clust or tr_2 not in tr_clust:
                continue
            if tr_clust[tr_1] != tr_clust[tr_2]:
                fn += 1
    nc2 = (num * (num-1)) / 2
    tn = nc2 - (fp + tp + fn)
    return tp, fp, tn, fn

def readCDHitClust(fn, filtDict=None):
    tr_clust = {}
    tr_clust_inv = {}
    fp = open(fn)
    cnum = None
    for l in fp:
        if l.startswith('>Cluster'):
            cnum = int(l.rstrip().split()[-1])
        else:
            e = l.split()[2].lstrip('>').rstrip('.')
            if not filtDict or e in filtDict:
                tr_clust[e] = cnum
    for k,v in tr_clust.iteritems():
        if v in tr_clust_inv:
            tr_clust_inv[v].append(k)
        else:
            tr_clust_inv[v] = [k]
    return tr_clust, tr_clust_inv

def readCorset(fn):
    ft = open(fn)
    clust_dict = {}
    clust_dict_inv = {}
    for line in ft:
        toks = line.rstrip().split()
        clust_dict[toks[0]] = toks[1]
    for k,v in clust_dict.iteritems():
        if v in clust_dict_inv:
            clust_dict_inv[v].append(k)
        else:
            clust_dict_inv[v] = [k]
    return clust_dict, clust_dict_inv

def readTrueLabels(fn):
    ft = open(fn)
    groundTruth_clust = {}
    groundTruth_clust_inv = {}
    gtClusterCount = {}
    for line in ft:
        tr_gn = line[:-1].split("\t")
        groundTruth_clust[tr_gn[0]] = tr_gn[1]
        if tr_gn[0] in gtClusterCount.keys():
            gtClusterCount[tr_gn[0]] += 1
        else:
            gtClusterCount[tr_gn[0]] = 1
    for k,v in groundTruth_clust.iteritems():
        if v in groundTruth_clust_inv:
            groundTruth_clust_inv[v].append(k)
        else:
            groundTruth_clust_inv[v] = [k]
    return groundTruth_clust, groundTruth_clust_inv
def readMCLClust(fn):
    fp = open(fn)
    tr_clust = {}
    tr_clust_inv = {}
    key = 1
    for cnum, line in enumerate(fp):
        same_cluster = line.rstrip().split('\t')
        for contig in same_cluster:
            tr_clust[contig] = cnum
    for k,v in tr_clust.iteritems():
        if v in tr_clust_inv:
            tr_clust_inv[v].append(k)
        else:
            tr_clust_inv[v] = [k]
    return tr_clust, tr_clust_inv

def readCommClust(fn):
    netdic = pickle.load(open(fn,'rb'))
    tr_clust = {}
    tr_clust_inv = {}
    for contig,cnum in netdic.iteritems():
        tr_clust[str(contig)] = cnum
    for k,v in tr_clust.iteritems():
        if v in tr_clust_inv:
            tr_clust_inv[v].append(k)
        else:
            tr_clust_inv[v] = [k]
    return tr_clust, tr_clust_inv

def measurePrecRecall(qsf,sp, ctype):
    if ctype == "mcl":
        tr_clust, tr_clust_inv = readMCLClust(qsf)
    elif ctype == "corset":
        tr_clust, tr_clust_inv = readCorset(qsf)

    if(sp == 'human'):
        ft = trinity_human
    elif (sp == 'yeast'):
        ft = trinity_yeast
    else:
        ft = trinity_chicken
    groundTruth_clust, ground_truth_clust_inv = readTrueLabels(ft)
    tp, fp, tn, fn = accuracyExpressedFast(groundTruth_clust, ground_truth_clust_inv, tr_clust, tr_clust_inv)
    print("tp : {}, fp : {}, tn : {}, fn : {}".format(tp, fp, tn, fn))
    print("prec: {}, recall: {}".format(tp / float(tp + fp), tp / float(tp + fn)))

import argparse

def main():
    parser = argparse.ArgumentParser(description="Give the SRR number")
    parser.add_argument('--clustfile',type = str, help="graph file")
    parser.add_argument('--sp',type = str, help="human or yeast or chicken")
    parser.add_argument('--ctype',type = str, default="mcl", help="human or yeast")
    args = parser.parse_args()
    measurePrecRecall(args.clustfile,args.sp, args.ctype)


if __name__ == "__main__":
    main()
