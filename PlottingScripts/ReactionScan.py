"""
Implementations of FBA Scans by Thomas Pfau
The author does not take any responsibility whatsoever for any damage cause by the code.
The code is provided, as is, without any warranty.
This code can be used free of charge for any non commercial applications, as long as this notice is retained.
"""


from ScrumPy.Data import DataSets

reload(DataSets)
import numpy
def MaintenanceScan(m,lp,Biomass,Maintenance,lo,hi,steps,objectdirec="Min",addConst = {}, quad = False,objective = [], lpsolve = False,logscale = False):
    """ r = raction to alter , target = minimisation target (after which fluxmin will be performed)
       supply an lp with all necessary constraints set, along with the name of a reaction to scan
       The method will scan from steps values between lo and hi. if lo < hi scanning will be performed from hi to lo"""
    upper = max(lo,hi)
    lower = min(lo,hi)

    csteps = []
    if logscale:
        csteps = list(numpy.logspace(lo,hi,steps))
#    stepsize = int(float(upper-lower)/float(steps-1)*10000.0)
#    stepsize = float(stepsize) / 10000.0
    else:
        csteps= list(numpy.linspace(lo,hi,steps))
    lpobj = {}
    if len(objective) == 0:
        objective = m.sm.cnames[:]
    if lpsolve:
        for r in m.sm.cnames:
            lpobj[r] = 0
    scanlist = []
    lp.SetObjDirec(objectdirec)
    #(1) set the target for calculating first objective
    if lpsolve:
            lpobj[Biomass] = 1
            lp.SetObjective(lpobj)
            lpobj[Biomass] = 0
    else:
            lp.SetObjective([Biomass])        
            lp.ClearObjective()
            lp.SetObjective([Biomass])
    lp.UnboundFlux(Biomass)
    targetvals = {}
    #(1) set the target for calculating first objective
    for step in csteps:

        lp.SetFluxConstraint({Biomass : step, Maintenance : -1},0,0,"Maintenance")
        #calculate the optimal value for target
        lp.Solve(False)
#        bound = lower+stepsize*step
        if not lp.GetStatusMsg() == "Optimal":
            print "Error at step " +  str(step)
        temp = lp.GetPrimSol()
        temp["Maintenance"] = step
        targetvals[step] = lp.GetObjVal()
        lp.DelRow("Maintenance")
    if lpsolve:
            for val in objective:
                lpobj[val] = 1
            lp.SetObjective(lpobj)
            for val in objective:
                lpobj[val] = 0                        
    else:
            lp.ClearObjective()
            lp.SetObjective(objective,quad)
    lp.SetObjDirec("Min")

    for step in csteps:
#        bound = lower+stepsize*step
        lp.SetFluxConstraint({Biomass : step, Maintenance : -1},0,0,"Maintenance")
        lp.SetFixedFlux({Biomass : targetvals[step]})
        lp.Solve(False)
        if not lp.GetStatusMsg() == "Optimal":
            print "Error at step " +  str(step)
        temp = lp.GetPrimSol()
        temp["Maintenance"] = step
        scanlist.append(temp)
        lp.DelRow("Maintenance")
    ds = DataSets.DataSet(ItemNames=m.sm.cnames[:]+["Maintenance"])
    i = 0
    for sol in scanlist:
        ds.NewRow(name=str(i))
        for r in sol:
            ds[str(i),r] = sol[r]
        i+=1
    return ds,scanlist



def MaxMinScan(m,lp,r,target,lo,hi,steps,objectdirec="Max",addConst = {}, quad = False,objective = [], lpsolve = False):
    """ r = raction to alter , target = minimisation target (after which fluxmin will be performed)
       supply an lp with all necessary constraints set, along with the name of a reaction to scan
       The method will scan from steps values between lo and hi. if lo < hi scanning will be performed from hi to lo"""
    upper = max(lo,hi)
    lower = min(lo,hi)
    stepsize = int(float(upper-lower)/float(steps-1)*10000.0)
    stepsize = float(stepsize) / 10000.0
    lpobj = {}
    if len(objective) == 0:
        objective = m.sm.cnames[:]
    if lpsolve:
        for r in m.sm.cnames:
            lpobj[r] = 0
    scanlist = []
    optvals = []
    lp.SetObjDirec(objectdirec)
    lp.UnboundFlux(target)
    if lpsolve:
            lpobj[target] = 1
            lp.SetObjective(lpobj)
            lpobj[target] = 0
    else:
            lp.SetObjective([target])        
            lp.ClearObjective()
            lp.SetObjective([target])
    for step in range(steps):
        bound = lower+stepsize*step
        lp.SetFixedFlux({r : bound})
        #calculate the optimal value for target
        lp.Solve(False)
        if not lp.GetStatusMsg() == "Optimal":
            print "Error at step " +  str(step)
        optvals.append(lp.GetObjVal())
    if lpsolve:
            for val in objective:
                lpobj[val] = 1
            lp.SetObjective(lpobj)
            for val in objective:
                lpobj[val] = 0                        
    else:
            lp.ClearObjective()
            lp.SetObjective(objective,quad)
    lp.SetObjDirec("Min")

        
    for step in range(steps):
        #now set target to this optimal value
        bound = lower+stepsize*step
        lp.SetFixedFlux({r : bound, target : optvals[step]})
        lp.Solve(False)
        if not lp.GetStatusMsg() == "Optimal":
            print "Error at step " +  str(step)
        temp = lp.GetPrimSol()
        temp["Objective"] = lp.GetObjVal()
	temp['Maintenance'] = bound;
        scanlist.append(temp)

    ds = DataSets.DataSet(ItemNames=m.sm.cnames[:]+["Objective",'Maintenance'])
    i = 0
    for sol in scanlist:
        ds.NewRow(name=str(i))
        for r in sol:
            ds[str(i),r] = sol[r]
        i+=1
    return ds,scanlist

