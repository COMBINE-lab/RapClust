#Give input and ouput files as arguments
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('inputfile')
parser.add_argument('outputfile')
args = parser.parse_args()

infile = open(args.inputfile, 'r')
outfile = open(args.outputfile, 'w')

clusterID = 1

for line in iter(infile):
	for contigID in (line.split('\t')):
		outfile.write(contigID.strip('\n') + "\tcluster" + str(clusterID) + '\n')
	clusterID = clusterID + 1

infile.close()	
outfile.close()
