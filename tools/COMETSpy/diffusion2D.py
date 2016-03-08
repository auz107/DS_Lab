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

class diffusion2D(object):
    """
    A class for solving the 2-dimensional diffusion equation 

    Ali R. Zomorrodi - Daniel Segre Lab @ Boston University
    Last updated: July 23 2014
    """

    def __init__(self,layout2D,dtDiff,timeRange):
       
        """
        INPUTS:
        ------
                layout2D: Two dimensional layout (an instance of the class layout2D.See layout2D.py 
                          for more info)
               timeRange: Range of the simulation time (in h). This is a list with 
                          either two or three elements.
                          Two-element format: [startTime,endTime]
                          Three-element format: [startTime,timeStep,endTime]. If the 
                          time step is not given a time step of 1 second is used
            screenOutput: By default (on) writes  a summary including ??? 
                          on the screen.
                          if set to a value of 'off' nothing is written on the screen, 
                          in which case The user can instead specifiy 
                          an output fiile using the option outputFile, or store 
                          them in the variable runOutput (see the 'run' method for
                          details)
              outputFile: Optional input. It is a string containg the path to a
                          file and its name (e.g., 'results/fbaResults.txt'),
                          where the results should be written to.

        """

        # 2D layout
        self.layout2D = layout2D

        # Start and end time of the simulation and the time step for growth and diffusion (in h) 
        self.t0 = timeRange[0]               # Start time
        if len(timeRange) == 2: 
            self.dt = 1                      # Time step (for growth)
            self.tf = timeRange[1]           # End of simulation time 
        else:
            self.dt = timeRange[1]           # Time step (for growth)
            self.tf = timeRange[2]           # End of simulation time  

        # Output to the screen 
        if screenOutput != None and screenOutput.lower() != 'off':
            raise customError("**Error! The only eligible value for screenOutput is 'silent'")
        else:
             self.screenOutput = screenOutput

        # Output file
        self.outputFile = outputFile


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

    def updateMetabConcs(self):
        """
        Updates the concnetration of the shared metabolites 
        """

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

            # Solve the diffusion equaiton in 2D
            diffModel = diffusion2D(layout2D,self.dtDiff,[0,self.dtGrowth,self.tmaxDiff])
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

