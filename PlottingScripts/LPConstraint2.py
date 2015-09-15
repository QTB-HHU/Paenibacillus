from ScrumPy.Util import Set
import cplex
import ScrumPyQP 
reload(ScrumPyQP)
QPClass = ScrumPyQP.qp

from LP import ScrumPyLP
reload(ScrumPyLP)
LPClass = ScrumPyLP.lp


Inputs=[]
InOuts=[]
Outputs=[]
NoCost = []

def BuildLP(m, Input=Inputs, InOut=InOuts, Output=Outputs,NoCosts=NoCost,Test=[],target=None,Solver="CPLEX",inf = cplex.infinity):
    """ build lp form model m, satisfying constraints defined here """
    if Solver == "lpsolve":
        lp = LPClass(m)
    else:
        lp = QPClass(m,inf)
    
    for i in Input:
        lp.SetColBounds(i,0,None)
        
    for io in InOut:
        lp.SetColBounds(io,None,None)
        
    for o in Output:
        lp.SetColBounds(o,0,None)
    for n in NoCosts:
        lp.SetObjCoef(n,0)
    if target == None:
        lp.SetObjective(m.sm.cnames[:])
    else:
        lp.SetObjective(target)
    for n in NoCosts:
        lp.SetObjCoef(n,0)
    for t in Test:
        lp.SetFixedFlux({t:0})
    return lp


