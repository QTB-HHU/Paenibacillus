"""
Script to replicate the Plots (in gnuplot instead of the refined matlab plots) presented in the paper:
Model assisted genome annotation reveals Paenibacillus polymyxa as an efficient symbiont and a potential host for fuels and chemicals production.
Script created by Thomas Pfau
The author does not take any responsibility whatsoever for any damage cause by the code.
The code is provided, as is, without any warranty.
This code can be used free of charge for any non commercial applications, as long as this notice is retained.
"""

import ReactionScan
import SolutionAnalyser
import LPConstraint2
reload(ReactionScan)
reload(SolutionAnalyser)
import ScrumPy

def MaintenanceScan(m,lp,lpsolve,Count=200):
    ScanDS,sollist = ReactionScan.MaintenanceScan(m,lp,"Biomass","ATPase",0,20,Count,"Max",lpsolve=lpsolve)
    return ScanDS,sollist


def MaxMinScan(m,lp,lpsolve,Count=200):
    ScanDS,sollist = ReactionScan.MaxMinScan(m,lp,"R_TCE_CARBON-DIOXIDE","Biomass",1.5,2.1,Count,"Max",lpsolve=lpsolve)
    return ScanDS,sollist


def PlotRelevants(ds):
    SolutionAnalyser.plot("Maintenance",ds,["Biomass","FHL-RXN", "R_TCE_ACETONE", "R_TCE_ACET", "R_TCE_ETOH", "R_TCE_BUTANEDIOL", 'R_TCE_NITRITE', 'R_TEC_NITRATE', 'R_TEC_AMMONIA'],10)


def ScanAndPlot(m,lp,lpsolve=True,count=200):
    ScanDS,sollist = MaintenanceScan(m,lp,lpsolve,count)
    PlotRelevants(ScanDS)
    return ScanDS

def RedoxScan(m,lp,lpsolve,Count=200):
    ScanDS,sollist = ReactionScan.MaintenanceScan(m,lp,"Biomass","DEHOG",-1,1,Count,"Max",lpsolve=lpsolve)
    return ScanDS,sollist

def PlotDEHOG(ds):
    SolutionAnalyser.plot("Maintenance",ds,["Biomass","FHL-RXN", "R_TCE_ACET", "R_TCE_ETOH", "R_TCE_BUTANEDIOL"],10)


def PlotAmmonia(ds):
    SolutionAnalyser.plot("Maintenance",ds,["Biomass","FHL-RXN", "R_TCE_ACETONE", "R_TCE_ACET", "R_TCE_ETOH", "R_TCE_BUTANEDIOL", 'R_TEC_AMMONIA'],10)



def Figure4a(m,count=100):
    lp = LPConstraint2.BuildLP(m)
    setupAmmoniumLP(lp);
    lp.SetFluxConstraint({"Biomass" : 6.58, "ATPase" : -1},0,0,"ATPMaintenance")
    lp.UnboundFlux('DEHOG')
    ScanDS,sollist = RedoxScan(m,lp,False,count)
    PlotDEHOG(ScanDS)
    printRelevants(ScanDS,'RedoxScan.csv')
    return ScanDS


def Figure3(m,count=100):
	lp = LPConstraint2.BuildLP(m)
	setupAmmoniumLPNoFHL(lp);
	ds = ScanAndPlot(m,lp,False,200)
	PlotRelevants(ds)
	printRelevants(ds,'FinalAmmonium.csv')
	return ds

	
def Figure4b(m):
	lp = LPConstraint2.BuildLP(m)
	setupAmmoniumLP(lp);
	lp.SetFluxConstraint({'Biomass' : 0.63, 'DEHOG' : -1},0,0,"RedoxMaintenance")
	lp.UnboundFlux('DEHOG')
	ds = ScanAndPlot(m,lp,False,200)
	PlotRelevants(ds)
	printRelevants(ds,'FinalAmmoniumNoFHLDEHOG.csv')
	return ds


def Figure5(m):
	lp = LPConstraint2.BuildLP(m)
	setupNitrateLPTCAEC(lp);
	lp.SetFluxConstraint({'Biomass' : 0.63, 'DEHOG' : -1},0,0,"RedoxMaintenance")
	lp.UnboundFlux('DEHOG')
	ds = ScanAndPlot(m,lp,False,200)
	PlotRelevants(ds)
	printRelevants(ds,'FinalNitrateECTCADEHOG.csv')
	return ds


def SuppFigure6(m):
	lp = LPConstraint2.BuildLP(m)
	setupAmmoniumLP(lp);
#	lp.SetFluxConstraint({'Biomass' : 0.63, 'DEHOG' : -1},0,0,"RedoxMaintenance")
	lp.SetFixedFlux({'DEHOG' : 0})
	ds = ScanAndPlot(m,lp,False,200)
	PlotRelevants(ds)
	printRelevants(ds,'PREDEHOGFix.csv')
	return ds

def SuppFigure7(m,count=100):
    lp = LPConstraint2.BuildLP(m)
    setupAmmoniumLP(lp);
    lp.SetFixedFlux({'DEHOG' : 0});
    lp.SetFluxConstraint({"Biomass" : 6.58, "ATPase" : -1},0,0,"ATPMaintenance")
    ScanDS,sollist = MaxMinScan(m,lp,False,count)
    PlotRelevants(ScanDS)
    printRelevants(ScanDS,'CarbonScan.csv')
    return ScanDS


def printRelevants(ds,filename):
	rels = ["Maintenance","Biomass","FHL-RXN","R_TCE_ACETONE","R_TCE_ACET","R_TCE_ETOH","R_TCE_BUTANEDIOL",'R_TCE_NITRITE','R_TEC_NITRATE','R_TEC_AMMONIA']
	relnames = ["Maintenance","Biomass", "FHL", "Acetone production", "Acetate production", "Ethanol production", "Butanediol production", "Nitrite export", "Nitrate uptake", "Ammonia uptake/export"]
	f = open(filename,'w');
	for r in relnames:
		f.write(r + '\t')
	f.write('\n')
	for i in range(len(ds.GetCol(0))):
		for j in range(len(rels)):
			f.write(str(ds.GetCol(rels[j])[i]) + '\t')
		f.write('\n')
	f.close()


def setupNitrateLPTCAEC(lp):
#TCA open, Electron Transfer Chain usable with nitrate
	lp.SetFixedFlux({"R_TEC_GLC" : 1, "R_TEC_GLYCEROL" : 0, "R_TEC_XYLOSE" : 0, 
			 "R_TEC_CELLOBIOSE" : 0, "DEHOG" : 0,"2OXOGLUTARATEDEH-RXN":0 })
	lp.SetFluxBounds({"R_TCE_ACET": (0,None)})
	lp.SetFluxBounds({"PYRUVFORMLY-RXN":(None,0)})
	lp.SetFluxBounds({"PYRUFLAVREDUCT-RXN":(0,None)})
	lp.SetFluxBounds({"R_TCE_PROTONEFF":(0,0.122)})
	lp.SetFluxBounds({"R_TEC_AMMONIA":(-16/42.,0)})
	lp.SetFluxBounds({"R_TEC_NITRATE":(0,40.53/42.)})
	lp.SetFluxBounds({"R_TCE_NITRITE":(0,8.98/42.)})
	lp.SetColBounds("ACETOACETYL-COA-TRANSFER-RXN",-0.05,None)





def setupNitrateLP(lp):
#TCA interrupted, Electron Transfer Chain not usable
	lp.SetFixedFlux({"R_TEC_GLC" : 1, "R_TEC_GLYCEROL" : 0, "R_TEC_XYLOSE" : 0, 
			 "R_TEC_CELLOBIOSE" : 0, "DEHOG" : 0,"2OXOGLUTARATEDEH-RXN":0 })
	#See what could happen here...
#	lp.SetFixedFlux({"E-TransferChain": 0 })
	lp.SetFluxBounds({"R_TCE_ACET": (0,None)})
	lp.SetFluxBounds({"PYRUVFORMLY-RXN":(None,0)})
	lp.SetFluxBounds({"PYRUFLAVREDUCT-RXN":(0,None)})
	lp.SetFluxBounds({"R_TCE_PROTONEFF":(0,0.122)})
	lp.SetFluxBounds({"R_TEC_AMMONIA":(-16/42.,0)})
	lp.SetFluxBounds({"R_TEC_NITRATE":(0,40.53/42.)})
	lp.SetFluxBounds({"R_TCE_NITRITE":(0,8.98/42.)})
	lp.SetColBounds("SUCC-FUM-OXRED-RXN",None,0)
	lp.SetColBounds("2OXOGLUTARATEDEH-RXN",0,None)
	lp.SetColBounds("ACETOACETYL-COA-TRANSFER-RXN",-0.05,None)



def setupAmmoniumLPNoFHL(lp):
	lp.SetFixedFlux({"R_TEC_GLC" : 1, "R_TEC_GLYCEROL" : 0, "R_TEC_XYLOSE" : 0, 
			 "R_TEC_CELLOBIOSE" : 0, "DEHOG" : 0,"2OXOGLUTARATEDEH-RXN":0, 
			 "E-TransferChain": 0,  "FHL-RXN" : 0 })
	lp.SetFluxBounds({"R_TCE_ACET": (0,None)})
	lp.SetFluxBounds({"PYRUVFORMLY-RXN":(None,0)})
	lp.SetFluxBounds({"PYRUFLAVREDUCT-RXN":(0,None)})
	lp.SetFluxBounds({"R_TCE_PROTONEFF":(0,0.123)})
	lp.SetFluxBounds({"R_TEC_NITRATE":(0,0)})
	lp.SetFluxBounds({"R_TCE_NITRITE":(0,8.98/55)})
	lp.SetFluxBounds({"R_TEC_AMMONIA":(0,26.92/42)})
	lp.SetColBounds("SUCC-FUM-OXRED-RXN",None,0)
	lp.SetColBounds("2OXOGLUTARATEDEH-RXN",0,None)
	lp.SetColBounds("ACETOACETYL-COA-TRANSFER-RXN",-0.05,None)


def setupAmmoniumLP(lp):
	lp.SetFixedFlux({"R_TEC_GLC" : 1, "R_TEC_GLYCEROL" : 0, "R_TEC_XYLOSE" : 0, 
			 "R_TEC_CELLOBIOSE" : 0, "DEHOG" : 0,"2OXOGLUTARATEDEH-RXN":0, 
			 "E-TransferChain": 0,  "R_TEC_NITRATE" : 0 })
	lp.SetFluxBounds({"R_TCE_ACET": (0,None)})
	lp.SetFluxBounds({"PYRUVFORMLY-RXN":(None,0)})
	lp.SetFluxBounds({"PYRUFLAVREDUCT-RXN":(0,None)})
	lp.SetFluxBounds({"R_TCE_PROTONEFF":(0,0.123)})
	lp.SetFluxBounds({"R_TEC_NITRATE":(0,0)})
	lp.SetFluxBounds({"R_TCE_NITRITE":(0,8.98/42)})
	lp.SetFluxBounds({"R_TEC_AMMONIA":(0,26.92/42)})
	lp.SetColBounds("SUCC-FUM-OXRED-RXN",None,0)
	lp.SetColBounds("2OXOGLUTARATEDEH-RXN",0,None)
	lp.SetColBounds("ACETOACETYL-COA-TRANSFER-RXN",-0.05,None)


m = ScrumPy.Model("Paenibacillus.spy")

ds1 = Figure3(m)
ds2 = Figure4a(m)
ds3 = Figure4b(m)
ds4 = Figure5(m)
ds5 = SuppFigure6(m)
ds6 = SuppFigure7(m)
