import matplotlib
matplotlib.use('PDF')
import click
import croc
import matplotlib.pyplot as plt
import seaborn as sns

@click.command()
@click.option("--dset", required=True, help="The dataset name used in the output pdf file")
@click.option("--label", type=(str, str), multiple=True, help="<filename> <method> that gives the scored-label list for a particular method")
def processCurves(dset, label):

    mdict = {}
    for fn, m in label:
        print("Processsing scored label file {} for method {}".format(fn, m))
        SD = None
        with open(fn) as scoreFile:
            SD = croc.ScoredData.read_from_file(scoreFile)
        curve = croc.ROC(SD.sweep_threshold())

        mdict[m] = curve
        print("AUC = {}".format(curve.area()))
        cfile = fn.split('.scored-label')[0] + '.curve'
        with open(cfile, 'w') as curveFile:
            curve.write_to_file(curveFile)

    sns.set_style("white")
    for m, c in mdict.iteritems():
        x, y = zip(*c)
        plt.plot(x, y, linewidth=5.0, label="{} (AUC = {:.3g})".format(m, c.area()))
    sns.despine()
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.legend(loc='lower right')
    plt.savefig("{}_dge_curves.pdf".format(dset))

if __name__ == "__main__":
    processCurves()
