import click

@click.command()
@click.option("--sampdirs", help="equivalence class file")
@click.option("--netfile", help="output net file")
@click.option("--cutoff", default=10, help="filter contigs with fewer than this many reads")
def buildNetFromEq(sampdirs, netfile, cutoff):
    import itertools
    import pandas as pd
    import numpy as np
    import os

    sep = os.path.sep

    dirs = os.listdir(sampdirs)
    sampdirs = [sep.join([sampdirs, f]) for f in dirs]
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
         
    #quant.set_index('Name', inplace=True)

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
            print("# tran = {}\n# eq = {}".format(numTran, numEq))
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
            if (estCount[t1] <= cutoff or estCount[t2] <= cutoff):
                continue
            tpm1 = tpm[t1] 
            tpm2 = tpm[t2]
            assert(((tpm1 + tpm2) / denom) <= 1.0)
            w = count * ((tpm1 + tpm2) / denom)
            if (t1, t2) in weightDict:
                weightDict[(t1, t2)] += w 
            else:
                weightDict[(t1, t2)] = w 
        for t in tids:
            if (estCount[t] <= cutoff):
                continue
            diagCounts[t] += count * (tpm[t] / denom)

    print("total reads = {}".format(tot))
    maxWeight = 0.0
    prior = 0.1
    for k,v in weightDict.iteritems():
        c0, c1 = diagCounts[k[0]], diagCounts[k[1]]
        #w = (v + prior) / (min(c0, c1) + prior)
        if c0 + c1 > epsilon:
            w = (v) / ((c0 + c1))
            weightDict[k] = w 
            if w > maxWeight:
                maxWeight = w
        else:
            weightDict[k] = 0.0
    print("max weight was {}; rescaling\n".format(maxWeight))
    #for k,v in weightDict.iteritems():
    #    weightDict[k] /= maxWeight

    for i in xrange(len(estCount)):
        if (estCount[i] > cutoff):
            weightDict[(i, i)] = np.finfo(float).eps

    with open(netfile, 'w') as ofile:
        for k,v in weightDict.iteritems():
            ofile.write("{}\t{}\t{}\n".format(tnames[k[0]], tnames[k[1]], v))

if __name__ == "__main__":
    import sys
    buildNetFromEq()
