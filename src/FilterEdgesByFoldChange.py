from __future__ import print_function
import pandas as pd
import collections
import math
import click
import os

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

@click.command()
@click.option("--sampdir", help="directory containing quantification results")
@click.option("--netfile", help="input network file")
@click.option("--ofile", help="output network file")
def filterGraph(sampdir, netfile, ofile):
    samplesPath = sampdir
    condSamplePairs = []
    for p in os.listdir(samplesPath):
        d = os.path.basename(os.path.normpath(p))
        # we assume that the sample directories 
        # have names of the form {condition}{replicate}
        # where condition is a single character.
        condSamplePairs.append((d[0], d[1:]))

    # Get just the set of condition names
    conditions = sorted(list(set(zip(*condSamplePairs)[0])))
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
    for cond,samp in condSamplePairs:
        if cond not in eqClasses:
            eqClasses[cond] = EquivCollection()
        filePath = os.path.sep.join([samplesPath, "{}{}".format(cond, samp), "aux", "eq_classes.txt"])
        readEqClass(filePath, eqClasses[cond])

    ambigCounts = {cond : getCountsFromEquiv(eqClasses[cond]) for cond in conditions}

    for cond in conditions:
        sailfish[cond] = ambigCounts[cond]

    print ("Done Reading")
    count = 0
    numTrimmed = 0
    with open(netfile) as f, open(ofile, 'w') as ofile:
        data = pd.read_table(f, header=None)
        for i in range(len(data)):
            count += 1
            print ("\r{} done".format(count), end= "")
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

if __name__ == "__main__":
    filterGraph()
