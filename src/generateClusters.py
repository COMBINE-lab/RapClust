import pandas as pd

from math import sqrt

def distance(a, b):
    return  sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2)

def point_line_distance(point, start, end):
    if (start == end):
        return distance(point, start)
    else:
        n = abs(
            (end[0] - start[0]) * (start[1] - point[1]) - (start[0] - point[0]) * (end[1] - start[1])
        )
        d = sqrt(
            (end[0] - start[0]) ** 2 + (end[1] - start[1]) ** 2
        )
        return n / d

def rdp(points, epsilon):
    """
    Reduces a series of points to a simplified version that loses detail, but
    maintains the general shape of the series.
    """
    dmax = 0.0
    index = 0
    for i in range(1, len(points) - 1):
        d = point_line_distance(points[i], points[0], points[-1])
        if d > dmax:
            index = i
            dmax = d
    return (index,dmax)
    if dmax >= epsilon:
        results = rdp(points[:index+1], epsilon)[:-1] + rdp(points[index:], epsilon)
    else:
        results = [points[0], points[-1]]
    return results

import os
from subprocess import call
def makeCluster(netfile,outfile):
    df = pd.read_csv(netfile,delimiter='\t', names = ['u','v','w'],header=None)
    df['w'] = df['w'].fillna(0)
    a = df['w'].values
    v,be = np.histogram(a,bins=100,density=True)
    (ind,maxx) = rdp(zip(range(0,len(v)),np.cumsum(v)),0.0)
    ind = ind/100.0
    cmd = "mcl " + netfile + " --abc -te 10 -o " + outfile + " -abc-tf 'gq("+str(ind)+")'"
    print cmd
    os.system(cmd)

import argparse

def main():
    parser = argparse.ArgumentParser(description="Give the SRR number")
    parser.add_argument('--netfile',type = str, help="input SRR path")
    parser.add_argument('--outfile',type=str, help="path to store temp files")
    args = parser.parse_args()
    makeCluster(args.netfile,args.outfile)


if __name__ == "__main__":
    main()
