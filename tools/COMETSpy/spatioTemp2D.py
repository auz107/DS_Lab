from __future__ import division
import sys,math
sys.path.append('../')
sys.path.append('../../')
from models import *
from metabolicModel import *
from fba import *
from DMMM import *
from tools.customError import *
from tools.pyomoSolverCreator import *

class spatioTemp2D(object):
    """
    A class for simulation of the spatio-temporal variations in microbial communities 

    Ali R. Zomorrodi - Daniel Segre Lab @ Boston University
    Last updated: July 23 2014
    """

    def __init__(self,layout2D,organisms,sharedMetabs,reactorType,timeRange,optSolverName = None,screenOutput = None, outputFile = None):
       
        """
        INPUTS:
        ------
                layout2D: Two dimensional layout (an instance of the class layout2D.See layout2D.py 
                          for more info)
               organisms: A dictionary whose keys are the nams of the species participating 
                          in the community and the values are the instances of class 
                          metabolicModel containing the metabolic models for that species.
                          The names chosen as keys in this dictionary coulud be different
                          from those in metabolicModel.microorganism_name
            sharedMetabs: A dictionary whose keys are the names of the shared metabolites
                          and the values are instances of class sharedMetab containing the 
                          information about that shared metabollite 
                          The names chosen as keys in this dictionary coulud be different
                          from those in sharedMetab.name
             reactorType: Type of reactor for simulations (batch or chemostat)
               timeRange: Range of the simulation time (in h). This is a list with 
                          either two or three elements.
                          Two-element format: [startTime,endTime]
                          Three-element format: [startTime,timeStep,endTime]. If the 
                          time step is not given a time step of 0.25 h (15 min o r 900 s) is used
           optSolverName: Name of the optimization solvers. This is a dictionary whose keys 
                          are 'LP', 'NLP', MILP, or 'MINLP' and the values are straings
                          containing the name of a solver for each class: 
                          Allowable choices for LP and MILP: cplex and gurobi
                          Allowable choices for NLP and MINLP: ** To be added ** 
            screenOutput: By default (None) writes  a summary including the 
                          solve  status, optimality status (if not optimal),
                          objective function value and the elapsed time on 
                          the screen if takes a value of 'silent' no resuults 
                          are written on the screen, in which case The user 
                          can instead specifiy  an output fiile using the 
                          option outputFile, or store them in the variable 
                          runOutput (see the 'run' method for details)
              outputFile: Optional input. It is a string containg the path to a
                          file and its name (e.g., 'results/fbaResults.txt'),
                          where the results should be written to.

        """

        # 2D layout
        self.layout2D = layout2D

        # List of the organisms
        self.organisms = organisms

        # Type of the reactor
        if reactorType.lower() not in ['batch','chemostat']:
           raise customError("**Error! Invalid reactorType (eligible choices are 'batch' and 'chemostat').")
        else:
            self.reactorType = reactorType

        # Start and end time of the simulation and the time step for growth and diffusion (in h) 
        self.t0 = timeRange[0]               # Start time
        if len(timeRange) == 2: 
            self.dtGrowth = 0.25             # Time step (for growth)
            self.tEndSim = timeRange[1]      # End of simulation time 
        else:
            self.dtGrowth = timeRange[1]     # Time step (for growth)
            self.tEndSim = timeRange[2]      # End of simulation time  

        # Time steps for diffusion process (in sec)
        self.dtDiff = 1 

        # Maximum time for diffusion simulaitons (in sec) 
        self.tmaxDiff = 30

        # Name of the optimization solver
        if optSolverName.keys() not in ['LP','MILP','NLP','MINLP']:
           raise customError("**Error! Invalid optSolverName (eligible choices are 'LP' and 'NLP','MILP','MINLP').")
        else:
           self.optSolverName = optSolverName

        # Output to the screen 
        if screenOutput != None and screenOutput.lower() != 'silent':
            raise customError("**Error! The only eligible value for screenOutput is 'silent'")
        else:
             self.screenOutput = screenOutput

        # Output file
        self.outputFile = outputFile

        self.orgNames = self.organisms.keys()
        self.sharedMetabNames = self.sharedMetabs.keys()

        # Extract the name of the shared metabolites each commmunity member exports and 
        # the names of correspondiing export reactions
        # self.orgExportRxnNames: A dictionary whose keys are the names of organisms and 
        #                      values are the names of the export reactions for the
        #                      shared metabolites it exports
        # self.orgSharedMetabExportRxn: A dictionary whose keys are the names of the 
        #                      community members and values are another dictionary
        #                      where the keys are the names of the shared metabolites
        #                      and values are the names of the corresponding export
        #                      reactions 
        self.orgExportRxnNames = dict((k,[]) for k in self.orgNames)
        self.orgSharedMetabExportRxn = dict((k,{}) for k in self.orgNames)
        for mname in self.sharedMetabNames:
            if self.sharedMetabs[mname].exportRxnNames != None:
                for orgName in self.sharedMetabs[mname].exportRxnNames.keys():
                    self.orgExportRxnNames[orgName].append(self.sharedMetabs[mname].exportRxnNames[orgName])
                    self.orgSharedMetabExportRxn[orgName][mname] = self.sharedMetabs[mname].exportRxnNames[orgName]


    def computeUptakeRates(self):
        """
        This function computes the uptake rates of the shared metabolites for a given 
        concentration 
        """

        for mname in self.sharedMetabNames:
            self.sharedMetabs[mname].evaluateUptakeRate()

            # Now assign the bounds on uptake reactions for each metabolic model
            for orgName in self.sharedMetabs[mname].uptakeRateVals.keys():
                uptakeRxnName = self.sharedMetabs[mname].uptakeRxnNames[orgName]
                self.organisms[orgName].flux_bounds[uptakeRxnName][0] = -self.sharedMetabs[mname].uptakeRateVals[orgName]

    def computeMuAndExport(self):
        """
        This function computes the specific growth rate of each species and the export
        rates of the shared metabolites using FBA for the beginning of the simulation time. 
        mu: A dictionary whose keys are the names of organisms and values are the specific 
            growth rates
        The export rates are stored in the instances of the class sharedMetab 
        """

        for orgName in self.orgNames:

            # Perform FBA
            modelFBA = fba(self.organisms[orgName],LPSolverName = 'gurobi',screenOutput = 'silent',reportLevel = self.orgExportRxnNames[orgName])
            [exitflag,mu,fbaExportRates] = modelFBA.run()

            if exitflag == 'globallyOptimal':
                self.mu[orgName] = mu
                for mname in self.orgSharedMetabExportRxn[orgName].keys():
                    expotRxnName = self.orgSharedMetabExportRxn[orgName][mname]
                    self.sharedMetabs[mname].exportRateVals[orgName] = fbaExportRates[expotRxnName]
            else:
                mu = self.deathPhase(orgName)
                for mname in self.orgSharedMetabExportRxn[orgName].keys():
                    expotRxnName = self.orgSharedMetabExportRxn[orgName][mname]
                    self.sharedMetabs[mname].exportRateVals[orgName] = 0

    def run(self):
        """
        Runs the spatio-temproal simulaitons

        OUTPUTS:
        """

        start_run = time.clock()

        # Compute the uptake rates first using the uptake kinetics and the initial metabolite
        # concentrations
        self.computeUptakeRates()

        # Compute the specific growth rate of each organism and export rates of
        # the limiting metabolites
        self.computeMuAndExport()

        done = 0

        # current time 
        self.t = 0

        while done == 0: 

            # Solve the diffusion equaiton in 2D for all extracellular metabolites
            diffModel = diffusion2D(layout2D,self.dtDiff,[0,self.dtDiff,self.tmaxDiff])
            diffModel.run()

            # Update the metabolite concentrations
            self.updateMetabConcs()

            # Update time in hours  
            self.t = self.t + self.dtmaxDiff/3600

            # Time range for DMMM simulations
            self.timeRangeDMMM = [self.t + self.tmaxDiff/3600,self.dtGrowth,self.t + self.tmaxDiff/3600 + self.dtGrowth]

            # Peform DMMM to get updated growth rates and export and uptate rates
            DMMM_current = DMMM(self.organisms,self.sharedMetabs,self.reactorType,timeRANGE = self.timeRangeDMMM,screenOutput = 'SILENT')
            DMMM_current.run()

            # Update time in hours  
            self.t = self.t + self.dtGrowth

            # Check termination condition
            if self.t >= self.tEndSim:
                done = 1



        # Elapsed to run
        elapsed_run = time.clock() - start_run


#--------- Sample implementation ------
if __name__ == "__main__":

