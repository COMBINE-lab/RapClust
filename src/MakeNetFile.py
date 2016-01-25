

def buildNetFromEq(eqfile, sfFile, netfile):
    import itertools
    import pandas as pd
    import numpy as np

    quant = pd.read_table(sfFile)
    quant['tname'] = quant['Name']
    quant.set_index('Name', inplace=True)
    tnames = []
    weightDict = {}
    diagCounts = np.zeros(len(quant['TPM'].values))
    tot = 0
    with open(eqfile) as ifile:
        numTran = int(ifile.readline().rstrip())
        numEq = int(ifile.readline().rstrip())
        print("# tran = {}\n# eq = {}".format(numTran, numEq))
        for i in xrange(numTran):
            tnames.append(ifile.readline().rstrip())

        tpm = quant.loc[tnames, 'TPM'].values
        estCount = quant.loc[tnames, 'NumReads'].values
        efflens = quant.loc[tnames, 'EffectiveLength'].values

        for i in xrange(numEq):
            toks = map(int, ifile.readline().rstrip().split('\t'))
            nt = toks[0]
            tids = toks[1:-1]
            count = toks[-1]
            denom = sum([tpm[t] for t in tids])
            tot += count
            for t1, t2 in itertools.combinations(tids,2):
                tpm1 = tpm[t1]
                tpm2 = tpm[t2]
                w = count * ((tpm1 + tpm2) / denom)
                if (t1, t2) in weightDict:
                    weightDict[(t1, t2)] += w
                else:
                    weightDict[(t1, t2)] = w
            for t in tids:
                w = count * (tpm[t]/ denom)
                diagCounts[t] += w

    print("total reads = {}".format(tot))
    for k,v in weightDict.iteritems():
        weightDict[k] = np.exp(-1 / ((v / (diagCounts[k[0]] + diagCounts[k[1]]) + 1.0))) # (0.5 * (diagCounts[k[0]] + diagCounts[k[1]]))

    for i in xrange(len(estCount)):
        weightDict[(i, i)] = 1.0

    with open(netfile, 'w') as ofile:
        for k,v in weightDict.iteritems():
            ofile.write("{}\t{}\t{}\n".format(tnames[k[0]], tnames[k[1]], v))

if __name__ == "__main__":
    import sys
    buildNetFromEq(sys.argv[1], sys.argv[2], sys.argv[3])
