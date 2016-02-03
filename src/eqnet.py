from __future__ import print_function

def buildNetFile(sampdirs, netfile, cutoff, writecomponents=False):
    import itertools
    import pandas as pd
    import numpy as np
    import os

    sep = os.path.sep
    sffiles = [sep.join([sd, 'quant.sf']) for sd in sampdirs]

    quant = None
    for sffile in sffiles:
        if quant is None:
            quant = pd.read_table(sffile)
            quant.set_index('Name', inplace=True)
        else:
            quant2 = pd.read_table(sffile)
            quant2.set_index('Name', inplace=True)
            quant += quant2

    tnames = []
    weightDict = {}
    diagCounts = np.zeros(len(quant['TPM'].values))

    tot = 0
    eqfiles = [sep.join([sd, 'aux/eq_classes.txt']) for sd in sampdirs]

    firstSamp = True
    numSamp = 0
    eqClasses = {}
    for eqfile in eqfiles:
        with open(eqfile) as ifile:
            numSamp += 1
            numTran = int(ifile.readline().rstrip())
            numEq = int(ifile.readline().rstrip())
            print("file: {}; # tran = {}; # eq = {}".format(eqfile, numTran, numEq))
            if firstSamp:
                for i in xrange(numTran):
                    tnames.append(ifile.readline().rstrip())
            else:
                for i in xrange(numTran):
                    ifile.readline()

            for i in xrange(numEq):
                toks = map(int, ifile.readline().rstrip().split('\t'))
                nt = toks[0]
                tids = tuple(toks[1:-1])
                count = toks[-1]
                if tids in eqClasses:
                    eqClasses[tids] += count
                else:
                    eqClasses[tids] = count

            firstSamp = False

    tpm = quant.loc[tnames, 'TPM'].values / numSamp
    estCount = quant.loc[tnames, 'NumReads'].values
    efflens = quant.loc[tnames, 'EffectiveLength'].values
    epsilon =  np.finfo(float).eps
    for tids, count in eqClasses.iteritems():
        denom = sum([tpm[t] for t in tids])
        tot += count
        for t1, t2 in itertools.combinations(tids,2):
            #tpm1 = tpm[t1]
            #tpm2 = tpm[t2]
            #w = count * ((tpm1 + tpm2) / denom)
            if (t1, t2) in weightDict:
                weightDict[(t1, t2)] += count
            else:
                weightDict[(t1, t2)] = count
        for t in tids:
            #if (estCount[t] <= cutoff):
            #    continue
            #diagCounts[t] += count * (tpm[t] / denom)
            diagCounts[t] += count


    print("total reads = {}".format(tot))
    maxWeight = 0.0
    prior = 0.1
    edgesToRemove = []
    for k,v in weightDict.iteritems():
        c0, c1 = diagCounts[k[0]], diagCounts[k[1]]
        #w = (v + prior) / (min(c0, c1) + prior)
        if c0 + c1 > epsilon and c0 > cutoff and c1 > cutoff:
            w = v / min(c0, c1)
            weightDict[k] = w
            if w > maxWeight:
                maxWeight = w
        else:
            edgesToRemove.append(k)

    for e in edgesToRemove:
        del weightDict[e]

    tnamesFilt = []
    relabel = {}
    for i in xrange(len(estCount)):
        if (diagCounts[i] > cutoff):
            relabel[i] = len(tnamesFilt)
            tnamesFilt.append(tnames[i])
            weightDict[(i, i)] = 1.1

    import networkx as nx
    G = nx.Graph() if writecomponents else None
    with open(netfile, 'w') as ofile:
        writeEdgeList(weightDict, tnames, ofile, G)

    if G is not None:
        clustFile = netfile.split('.net')[0] + '.clust'
        print("Writing connected components as clusters to {}".format(clustFile))
        with open(clustFile, 'w') as ofile:
            cc = nx.connected_component_subgraphs(G)
            for c in cc:
                ofile.write('{}\n'.format('\t'.join(c.nodes())))

def writeEdgeList(weightDict, tnames, ofile, G):
    useGraph = G is not None
    for k,v in weightDict.iteritems():
        ofile.write("{}\t{}\t{}\n".format(tnames[k[0]], tnames[k[1]], v))
        if useGraph:
            G.add_edge(tnames[k[0]], tnames[k[1]])


def writePajek(weightDict, tnames, relabel, ofile):
    with open(netfile, 'w') as ofile:
        ofile.write("*Vertices\t{}\n".format(len(tnamesFilt)))
        for i, n in enumerate(tnamesFilt):
            ofile.write("{}\t\"{}\"\n".format(i, n))
        ofile.write("*Edges\n")
        print("There are {} edges\n".format(len(weightDict)))
        for k,v in weightDict.iteritems():
            ofile.write("{}\t{}\t{}\n".format(relabel[k[0]], relabel[k[1]], v))
            #ofile.write("{}\t{}\t{}\n".format(tnames[k[0]], tnames[k[1]], v))
            #if k[0] != k[1]:
            #    ofile.write("{}\t{}\t{}\n".format(tnames[k[1]], tnames[k[0]], v))

class EquivCollection(object):
    def __init__(self):
        self.tnames = []
        self.eqClasses = {}
        self.hasNames = False

    def setNames(self, names):
        self.tnames = names
        self.hasNames = True

    def add(self, tids, count):
        if tids in self.eqClasses:
            self.eqClasses[tids] += count
        else:
            self.eqClasses[tids] = count

def readEqClass(eqfile, eqCollection):
    with open(eqfile) as ifile:
        numTran = int(ifile.readline().rstrip())
        numEq = int(ifile.readline().rstrip())
        print("file: {}; # tran = {}; # eq = {}".format(eqfile, numTran, numEq))
        if not eqCollection.hasNames:
            tnames = []
            for i in xrange(numTran):
                tnames.append(ifile.readline().rstrip())
            eqCollection.setNames(tnames)
        else:
            for i in xrange(numTran):
                ifile.readline()

        for i in xrange(numEq):
            toks = map(int, ifile.readline().rstrip().split('\t'))
            nt = toks[0]
            tids = tuple(toks[1:-1])
            count = toks[-1]
            eqCollection.add(tids, count)

def getCountsFromEquiv(eqCollection):
    countDict = {}
    tn = eqCollection.tnames
    for tids, count in eqCollection.eqClasses.iteritems():
        for t in tids:
            if tn[t] in countDict:
                countDict[tn[t]] += count
            else:
                countDict[tn[t]] = count
    # ensure no division by 0
    for t in eqCollection.tnames:
        if t in countDict:
            countDict[t] += 1.0
        else:
            countDict[t] = 1.0
    return countDict

def flattenClusters(infile, outfile):
    with open(outfile, 'w') as ofile:
        with open(infile) as ifile:
            for i,l in enumerate(ifile):
                toks = l.rstrip().split()
                cname = "cluster{}".format(i)
                for t in toks:
                    ofile.write("{}\t{}\n".format(cname, t))

def filterGraph(expDict, netfile, ofile):
    import os
    import pandas as pd
    import math

    # Get just the set of condition names
    conditions = expDict.keys()
    print("conditions = {}".format(conditions))

    #for cond in conditions:
    #    sailfish[cond] = collections.defaultdict(float)
    #    for sample in samples:
    #        fileDir = cond+sample
    #        filePath = os.path.sep.join([samplesPath, fileDir, "quant.sf"])
    #        print(filePath)
    #        with open(filePath) as f:
    #            data = pd.read_table(f, header=0)
    #            for i in range(len(data)):
    #                name = data['Name'][i]
    #                numreads = data['NumReads'][i]
    #                sailfish[cond][name] += numreads
    #    # To avoid divide by 0 error
    #    for name in sailfish[cond]:
    #        sailfish[cond][name] += 1

    eqClasses = {}
    for cond in conditions:
        print(expDict[cond])
        for sampNum, sampPath in expDict[cond].iteritems():
            if cond not in eqClasses:
                eqClasses[cond] = EquivCollection()
            eqPath = os.path.sep.join([sampPath, "aux", "eq_classes.txt"])
            readEqClass(eqPath, eqClasses[cond])

    ambigCounts = {cond : getCountsFromEquiv(eqClasses[cond]) for cond in conditions}

    sailfish = {}
    for cond in conditions:
        sailfish[cond] = ambigCounts[cond]

    print ("Done Reading")
    count = 0
    numTrimmed = 0
    with open(netfile) as f, open(ofile, 'w') as ofile:
        data = pd.read_table(f, header=None)
        for i in range(len(data)):
            count += 1
            print("\r{} done".format(count), end="")
            #Alternative hypo
            x = data[0][i]
            y = data[1][i]
            non_null=0
            x_all=0
            y_all=0
            for cond in conditions:
                y_g = sailfish[cond][y]
                x_g = sailfish[cond][x]
                r = y_g / x_g
                non_null += (y_g * math.log(r*x_g)) - (r*x_g)
                non_null += (x_g * math.log(x_g)) - x_g
                x_all += x_g
                y_all += y_g
            #null hypothesis
            null = 0
            r_all = y_all / x_all
            for cond in conditions:
                y_g = sailfish[cond][y]
                x_g = sailfish[cond][x]
                mean_x = (x_g + y_g) / (1+r_all)
                null += (y_g * math.log(r_all * mean_x)) - (r_all * mean_x)
                null += (x_g * math.log(mean_x)) - mean_x
            D = 2*(non_null-null)
            if D <= 20:
                ofile.write("{}\t{}\t{}\n".format(x, y, data[2][i]))
            else:
                numTrimmed += 1
    print("\nTrimmed {} edges".format(numTrimmed))



