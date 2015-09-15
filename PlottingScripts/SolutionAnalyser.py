"""
Plotting Tools for visualisation of FBA Scans by Thomas Pfau
The author does not take any responsibility whatsoever for any damage cause by the code.
The code is provided, as is, without any warranty.
This code can be used free of charge for any non commercial applications, as long as this notice is retained.
"""

from ScrumPy.Data import DataSets,Stats
from ScrumPy.Util import Sci

def GetReactionsFromTree(tree):
    res = []
    for r in tree.get_leaves():
        res.append(r.name)
    return res

def plot(x,ds,tree,count,exclude = []):
    """x - X-Axis, ds - dataset to plot, tree - the subtree to be plotted, count - #ofreactions per plot,
        exclude - reactions excluded from plot"""
    ds.RemoveAllFromPlot()
    type(tree)
    if not ( type(tree) == type([]) or type(tree) == type(set([]))):
        reacs = GetReactionsFromTree(tree)
    else:
        reacs = list(tree)
    pos = 0;
    items = -1
    Plot = True
    ds.SetPlotX(x)
    while pos < len(reacs):
        if items % count == 0:
            if not Plot:
                pass
            else:
                ds.Plot()
                Plot = False
                try:
                    temp = input("Next")
                except:
                    pass
                    
                try:
                    if temp == "X":
                        break
                except:
                    pass
                ds.RemoveAllFromPlot()
	if not reacs[pos] in exclude:
            ds.AddToPlot(reacs[pos])
            Plot = True
            if items < 0:
                items = 1
            else:
                items+=1
        pos+=1
    ds.Plot()
#    ds.RemoveAllFromPlot()


