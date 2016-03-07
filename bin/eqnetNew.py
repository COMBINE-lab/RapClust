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
    import gzip
    import struct
    import numpy as np

    # Get just the set of condition names
    conditions = expDict.keys()
    print("conditions = {}".format(conditions))

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

    #NEW CODE BY AVI
    ipath = "/home/laraib/clust/sailfish/data/human/"

    import pandas as pd
    from scipy import stats
    #from rpy2.robjects.packages import importr
    #from rpy2.robjects.vectors import FloatVector

    #Rstats = importr('stats', on_conflict="warn")

    conditions = ['A', 'B']
    replicates = ['1', '2', '3']
    nboots = 50

    data = {}

    with open(ipath+"A1/aux/bootstrap/names.tsv") as names:
        contigIds = {k : v for v,k in enumerate(names.readline().strip().split("\t"))}

    n = len(contigIds)
    bsDict = {'A' : [], 'B' : []}
    s = struct.Struct('d'*n)
    for condition in conditions:
        for replicate in replicates:
            path = os.path.sep.join([ipath+condition+replicate,"aux","bootstrap"])
            bootstrapFile = os.path.sep.join([path, "bootstraps.gz"])

            with gzip.open(bootstrapFile, 'r') as boots:
                readCounts = [ np.array(s.unpack_from(boots.read(8*n))) for i in xrange(nboots)]
                bsDict[condition].append(np.vstack(readCounts).T)

    for k,a in bsDict.iteritems():
        for v in a:
            print("k = {}, size={}".format(k, v.shape))
    nsamp = 1 
    last = 0
    with open(netfile) as f, open(ofile, 'w') as ofile:
        netdata = pd.read_table(f, header=None)
        numTest = len(netdata)
        pvalsAll = []
        for i in range(len(netdata)):
            count += 1
            print("\r{} done".format(count), end="")

            xContig = netdata[0][i]
            yContig = netdata[1][i]

            xInd = contigIds[xContig]
            yInd = contigIds[yContig]

            fcx = []
            fcy = []
            xCondA = np.hstack([bsDict['A'][r][xInd,:] for r in xrange(3)])
            yCondA = np.hstack([bsDict['A'][r][yInd,:] for r in xrange(3)])
            xCondB = np.hstack([bsDict['B'][r][xInd,:] for r in xrange(3)])
            yCondB = np.hstack([bsDict['B'][r][yInd,:] for r in xrange(3)])
            for j in xrange(nsamp):
            #for bootIndA in range(nboots):
                np.random.shuffle(xCondA)
                np.random.shuffle(yCondA)
                np.random.shuffle(xCondB)
                np.random.shuffle(yCondB)
                condASamps = (xCondA + 1.0) / (yCondA + 1.0)
                condBSamps = (xCondB + 1.0) / (yCondB + 1.0)
                #for bootIndB in range(nboots):

                    #fcx.append((a1x+a2x+a3x)/(b1x+b2x+b3x))
                    #fcy.append((a1y+a2y+a3y)/(b1y+b2y+b3y))
                fcx += list(condASamps)
                fcy += list(condBSamps)
                    #fcy.append((a1y/b1y) + (a2y/b2y) + (a3y/b3y) + (a1y/b2y) + (a1y/b3y) + (a2y/b1y) + (a2y/b3y) + (a3y/b1y) + (a3y/b2y) )
            #(Dval,pval) = stats.ks_2samp(fcx, fcy)
            mux, stdx= stats.norm.fit(fcx)
            muy, stdy = stats.norm.fit(fcy)
            (stat, pval) = stats.ttest_ind(stats.norm.rvs(loc=mux, scale=stdx, size=50),
                                           stats.norm.rvs(loc=muy, scale=stdy, size=50), equal_var=False)

            diff = False
            if (muy >= mux + 2 * stdx) or  (mux >= muy  + 2 * stdy):
                diff = True

            #if not diff:
            #(stat, crit, pval) = stats.anderson_ksamp([fcx, fcy])
            if pval > (0.01 / numTest):
            #if pval > (0.05/len(netdata)):
            #if (Dval <= 0.8):
                ofile.write("{}\t{}\t{}\n".format(xContig, yContig, netdata[2][i]))
            else:
                numTrimmed += 1
                if numTrimmed - last >= 1000:
                    print("\r\r Trimmed {}% of edges".format(100.0 * numTrimmed / float(i+1), end=""))
                    print(mux, muy, pval)
                    last = numTrimmed

    print("\nTrimmed {} edges".format(numTrimmed))

#    conditionsAvi = ['A', 'B']
#    replicates = ['1', '2', '3']
#
#    data = {}
#
#    with open(ipath+"A1/aux/bootstrap/names.tsv") as names:
#        contigIds = names.readline().strip().split("\t")
#
#    for condition in conditionsAvi:
#        readCounts = collections.defaultdict(list)
#        for replicate in replicates:
#            path = ipath+condition+replicate+"/aux/bootstrap/"
#            bootstrapFile = path+"bootstraps.txt"
#
#            #shell("gunzip -c {bootstrapFile} > boots.txt")
#
#            with open(bootstrapFile) as boots:
#                for line in boots:
#                    counts = line.strip().split("\t")
#                    for ind, count in enumerate(counts):
#                        readCounts[contigIds[ind]].append(float(count)+1)
#        data[condition] = readCounts
#
#    #shell("rm boots.txt names.txt")
#    count = 0
#    with open(netfile) as f, open(ofile, 'w') as ofile:
#        netdata = pd.read_table(f, header=None)
#        pvalsAll = []
#        for i in range(len(netdata)):
#            count += 1
#            print("\r{} done".format(count), end="")
#
#            xContig = netdata[0][i]
#            yContig = netdata[1][i]
#
#            #xInd = contigIds.index(xContig)
#            #yInd = contigIds.index(yContig)
#
#            fcx = []
#            fcy = []
#            for _ in range(2500):
#                ind = randint(0,149)
#                ax = data['A'][xContig][ind]
#                ay = data['A'][yContig][ind]
#
#                ind = randint(0,149)
#                bx = data['B'][xContig][ind]
#                by = data['B'][yContig][ind]
#
#                if (bx == 0):
#                    fcx.append(0)
#                else:
#                    fcx.append( ax/bx )
#
#                if (by == 0):
#                    fcy.append(0)
#                else:
#                    fcy.append( ay/by )
#
#            np.sort(fcx)
#            np.sort(fcy)
#            (Dval,pval) = stats.ks_2samp(fcx, fcy)
#
#            #pvalsAll.append(pval)
#            #if pval > (0.05/len(netdata)):
#            #if (Dval <= 0.54):
#            ofile.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(xContig, yContig, netdata[2][i],D,Dval,pval))
#            #else:
#            #    numTrimmed += 1
#
#        #p_adjust = Rstats.p_adjust(FloatVector(pvalsAll), method = 'BH')
#        #for i in range(len(netdata)):
#        #    if p_adjust[i] > 0.05:
#        #        ofile.write("{}\t{}\t{}\n".format(xContig, yContig, netdata[2][i]))
#        #    else:
#        #        numTrimmed += 1
#
#            #if pvals > 0.05:
#           #if pval > (0.05/len(netdata)):
#            #    ofile.write("{}\t{}\t{}\n".format(xContig, yContig, netdata[2][i]))
#            #else:
#            #    numTrimmed += 1
#    print("\nTrimmed {} edges".format(numTrimmed))



