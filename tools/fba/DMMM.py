from __future__ import division
import sys, time
sys.path.append('../../')
from datetime import timedelta  # To convert elapsed time to hh:mm:ss format
from tools.userError import userError 
from tools.core.compound import compound
from tools.core.reaction import reaction
from tools.core.organism import organism
from tools.core.model import model
from fba import fba
from death_rate import death_rate

class DMMM(object):
    """
    Performs dynamic multi-species metabolic modeling (DMMM) (PMID: 20668487) 

    Ali R. Zomorrodi - Daniel Segre Lab @ Boston University
    Last updated: 01-12-2015 (made sure dividing the time interval to smaller ones
                              to avoid negative compound concentrations works OK)
    """   

    def __init__(self,community_members,shared_compounds,reactor_type = 'batch',time_range = [0,1,10],carrying_capacity = {'cells_per_ml':float('inf'),'compounds_mM':float('inf')},serial_dilution_params = {'dilution_factor':None,'time_between_dilutions':None}, random_mortality_percentage = None, lag_phase_time = 0, store_all_fluxes = False, warnings = True, use_simplistic_time_update = False, stdout_msgs = True, **additional_args):

        """
        INPUTS:
        ------
          community_members: A list of objects of type model for each community member.
                             The object organism of the moell must have the following 
                             properties specified:
                             gDW_per_ml: Initial cell cencentration in gDW/ml in the form 
                                         of a dictionary where the key is the initial time
                                         point and the value is the initial cell concentraiton
                             The following optional properties can also be provided:
                             cells_per_ml: Initial cell cencentration in cells/ml in the form 
                                         of a dictionary where the key is the initial time
                                         point and the value is the initial cell concentraiton
                             gDW_per_cell: Gram of dry weight for one cell 
           shared_compounds: A list of objects of type compound containing the information
                             about shared compounds. The assigned compartment for these
                             compounds is 'shared', Each shared compound object must have the
                             following properties specified: 
                             reactant_reactions: Exchange reactions of community members taking
                             up this compound 
                             product_reactions: Echange reactions of community members producing 
                             this compound
                             reactions: reactant_reactions + product_reactions
                             concentration: Initial concentration in the form of a dictionary
                                            where the key is the initial time point and value
                                            is initial concentration in mM. Example: {0:0.1} 
               reactor_type: Type of reactor for simulations (batch, serial_dilution, life_cycle, chemostat)
                 time_range: Range of the simulation time. This is a list with 
                             either two or three elements.
                             Two-element format: [startTime,endTime]
                             Three-element format: [startTime,timeStep,endTime]. If the 
                             time step is not given a time step of one is used
          carrying_capacity: Carrying capacity for cells (in cells/ml) of the system 
                             (for logistic growth) and for compounds (in mM)
     serial_dilution_params: Ia needed if reactor_type is "serial_dilution". It is a dictionary
                             with two keys 'dilution_factor' (e.g., 1000) and time_between_dilutions in hour
     random_mortality_percentage: Percentage of biomass that dies out due to random mortality. If not value is provided
                             it is set to 0.1*dt such that for it is 0.1% for dt = 1, corresponding to mu_d = 0.001 1/h 
                             for dt = 1 (see self.random_mortality_percentage for details)
       store_all_fluxes: If True, dynamic changes in flux of all reactions are stored 
                             in the 'flux' field by setting store_opt_fluxes for fba to True.
                             The default for this input is False. Note that even when the value
                             of this parameter is True, dynamic changes in fluxes for reaction
                             with store_flux = True are stored. 
     use_simplistic_time_update: A parameter showing a simplistic time update methods should be used (True) or not. 
                             The simplistic method invovles setting all negative compound concentrations to zero if 
                             such a thing occurs. If this parameter is False, then the time interval is divided into 
                             smaller ones to avoid zero conetrations.
             lag_phase_time: The lag phase time (the time during which there is no growth). Default is zero
                   warnings: Can be True or False shwoing whether the warnings should be writtten to the 
                             screen or not
                stdout_msgs: Can be True or False shwoing whether any messages should be written to the
                             screen

        NOTE:
        -----
        Reaction fluxes and cell concentrations (organism.cells_per_ml) in the model as well as
        shared compound concentraitons (shared_cmp.concentration) should have already been
        stored in the form a dictionary. Keys of this dictionary are time points and 
        values are reactoin fluxes or concentraitons. This dictionary can be empty for
        fluxes at the begining but it must have one member for the other two 
        indicating the initial concentration of the shared compounds and initial cell 
        concentrations. For example, organsim.concentration = {0:10^5} or
        shared_cmp.concentraiton = [0:0.5]
        """
        self.community_members = community_members
        self.shared_compounds = shared_compounds

        if reactor_type not in ['batch','serial_dilution','life_cycle','chemostat']:
            raise ValueError('reactor_type must be batch, serial_dilution, life_cycle or chemostat.')
        else:
            self.reactor_type = reactor_type

        if store_all_fluxes not in [True,False]:
            raise userError('**ERROR! store_all_fluxes must be either True or False\n')
        else:
            self.store_all_fluxes = store_all_fluxes

        # Warnings and messages in the standard output
        if not isinstance(stdout_msgs,bool): 
            raise TypeError("stdout_msgs must be either True or False")
        else:
            self.stdout_msgs = stdout_msgs

        if not isinstance(warnings,bool): 
            raise TypeError("warnings must be either True or False")
        else:
            self.warnings = warnings

        if len(time_range) == 3:
            # Initial time
            self._t0 = time_range[0]

            # Original dt. This mught be adjusted during the simulations
            self._dt = time_range[1]

            # Final simulation time
            self._tf = time_range[2]   
        elif len(time_range) == 2:
            # Initial time
            self._t0 = time_range[0]

            # Original dt. This mught be adjusted during the simulations
            self._dt = 1

            # Final simulation time
            self._tf = time_range[1]   # Final simulation time
        else:
            userError("**ERROR! Invalid time_range (allowable size of the vectors are two or three)") 

        # Generate all time points
        decimal_digits_num_dt = len(str(self._dt - int(self._dt))[2:])
        self._time_points = [t/10 for t in range(self._t0,int(self._tf*(10**decimal_digits_num_dt) + self._dt*(10**decimal_digits_num_dt)),int(self._dt*(10**decimal_digits_num_dt)))]

        # Carrying capacity
        if not isinstance(carrying_capacity,dict):
            raise TypeError('carrying_capacity must be a dictionary')
        else:
            self.carrying_capacity = carrying_capacity

        if 'cells_per_ml' not in self.carrying_capacity.keys() and 'gDW_per_ml' not in self.carrying_capacity.keys() and 'compounds_mM' not in self.carrying_capacity.keys():
            raise ValueError('carrying_capacity must be a dictionary with keys cells_per_ml, gDW_per_ml, and/or compounds_mM')

        if 'cells_per_ml' not in self.carrying_capacity.keys():
            self.carrying_capacity['cells_per_ml'] = float('inf')
        elif type(self.carrying_capacity['cells_per_ml']) is not int and type(self.carrying_capacity['cells_per_ml']) is not float:
            raise ValueError("carrying_capacity['cells_per_ml'] must be float or integer.") 
        elif self.carrying_capacity['cells_per_ml'] != float('inf') and self.carrying_capacity['cells_per_ml'] <= 0:
            raise ValueError('Invalid carrying capacity for cells_per_ml(should be inf or a positive value)')

        if 'gDW_per_ml' not in self.carrying_capacity.keys() and 'cells_per_ml' not in self.carrying_capacity.keys():
            self.carrying_capacity['gDW_per_ml'] = float('inf')
        elif 'gDW_per_ml' not in self.carrying_capacity.keys() and 'cells_per_ml' in self.carrying_capacity.keys():
            # Average gDW/cell
            gDW_per_cell_ave = sum([member.organism.gDW_per_cell for member in self.community_members])/len(self.community_members)
            self.carrying_capacity['gDW_per_ml'] = self.carrying_capacity['cells_per_ml']*gDW_per_cell_ave
        elif type(self.carrying_capacity['gDW_per_ml']) is not int and type(self.carrying_capacity['gDW_per_ml']) is not float:
            raise ValueError("carrying_capacity['gDW_per_ml'] must be float or integer.") 
        elif self.carrying_capacity['gDW_per_ml'] != float('inf') and self.carrying_capacity['gDW_per_ml'] <= 0:
            raise ValueError('Invalid carrying capacity for gDW_per_ml(should be inf or a positive value)')

        if 'compounds_mM' not in self.carrying_capacity.keys():
            self.carrying_capacity['compounds_mM'] = float('inf')
        elif type(self.carrying_capacity['compounds_mM']) is not int and type(self.carrying_capacity['compounds_mM']) is not float:
            raise ValueError("carrying_capacity['compounds_mM'] must be float or integer.") 
        elif self.carrying_capacity['compounds_mM'] != float('inf') and self.carrying_capacity['compounds_mM'] <= 0:
            raise ValueError('Invalid carrying capacity for compounds_mM(should be inf or a positive value)')

        #-- Serial dilution parameters ---
        if not isinstance(serial_dilution_params,dict):
            raise TypeError('serial_dilution_params must a dictionary.')
        else:
            self.serial_dilution_params = serial_dilution_params

        # Percentage of biomass that dies out due to random mortality effects
        if random_mortality_percentage == None:
            self.random_mortality_percentage = 0.1*self._dt 
        else: 
            self.random_mortality_percentage = random_mortality_percentage
        
        # Convert random_mortality_percentage to a rate
        # Random mortality rate captured by random_mortality_percentage:
        # If this rate is shown by mu_d (mu_d < 0): we have dX/dt = mu_d*X
        # (X[t + dt] - X[t]) = mu_d*X[t]*dt
        # But, we have: X[t + dt] = (1 - random_mortality_percentage/100)*X[t]
        # Therefore, mu_d = -random_mortality_percentage/(100*dt)
        random_mortality_rate = -self.random_mortality_percentage/(100*self._dt)
        for member in self.community_members:
            if member.organism.random_mortality_rate == None:
                member.organism.random_mortality_rate = random_mortality_rate
            elif self.warnings:
                print '** WARNING (DMMM.py)! random_mortality_rate has already been assigned to organism ',member.id,' and the random_mortality_rate computed using random_mortality_percentage was not assigned to this member. Set member.mortality_rate to None in orderto assign random_mortality_rate computed using random_mortality_percentage'
 
        # Threshould beyond which a fractional delta_t is set to zero
        self._dt_frac_zero_thr = 1e-10

        # Lage phase time
        self.lag_phase_time = lag_phase_time
        self._lag_phase_time_orig = lag_phase_time

        #-- Check if an initial concentration for the shared compounds and cells --
        # have been provided. If the initial concentrations are provided but they 
        # are in the form of a single value, instead of a dictionary, they are assigned
        # to the first time point as a dictionary. 

        # The initial concentration can be given for a cell can be given in gDW/ml
        # or in cells_per_ml
        no_conc_members = []
        for member in self.community_members:
            if member.organism.cells_per_ml is None and member.organism.gDW_per_ml is None:
                no_conc_members.append(member)

            elif member.organism.cells_per_ml is not None and member.organism.gDW_per_ml is None:
                if type(member.organism.cells_per_ml) is not dict:
                    member.organism.cells_per_ml = {self._t0:member.organism.cells_per_ml}

                member.organism.gDW_per_ml = {self._t0:member.organism.cells_per_ml[self._t0]*member.organism.gDW_per_cell}

            elif member.organism.cells_per_ml is None and member.organism.gDW_per_ml is not None:
                if type(member.organism.gDW_per_ml) is not dict:
                    member.organism.gDW_per_ml = {self._t0:member.organism.gDW_per_ml}
                # If initial concentration is given in terms of gDW_per_ml, then there is no need
                # for DMMM to have gDW_per_cell
                if member.organism.gDW_per_cell is not None:
                    member.organism.cells_per_ml = {self._t0:member.organism.gDW_per_ml[self._t0]/member.organism.gDW_per_cell}

            elif member.organism.cells_per_ml is not None and member.organism.gDW_per_ml is not None:
                if type(member.organism.cells_per_ml) is not dict:
                    member.organism.cells_per_ml = {self._t0:member.organism.cells_per_ml}
                if type(member.organism.gDW_per_ml) is not dict:
                    member.organism.gDW_per_ml = {self._t0:member.organism.gDW_per_ml}

        if len(no_conc_members) > 0:
            print 'Initial biomass concnetration was not provided for the following community members:'
            for member in no_conc_members:
                print member.id
            raise userError()

        # Check initial concentration for shared compounds
        non_conc_cmp = []
        for shared_cmp in self.shared_compounds:
            if shared_cmp.concentration is None:
                non_conc_cmp.append(shared_meetab)
            # If the initial concentration has been provided but not in the form of
            # a dicitonary, convert it into a dicitonary and assigne it to the initial
            # time point specified by time_range[0]
            elif type(shared_cmp.concentration) is not dict:
                shared_cmp.concentration = {self._t0:shared_cmp.concentration}

        # Initialize mu for each community members
        for member in self.community_members:
            if member.organism.mu == None:
                member.organism.mu = {} 

        # Initialize the field "reactions" for all community members as a dictionary
        for member in self.community_members:
            for rxn in member.reactions:
                rxn.flux = {}

                # Set store_flux to False for all reactions 
                if self.store_all_fluxes:
                    rxn.store_flux = True 
                else: 
                    rxn.store_flux = False 

        # Set store_flux to True for biomass reacitons as well as for reactions in shared_cmp.reactions for shared compounds 
        for member in self.community_members:
            member.biomass_reaction.store_flux = True
        for shared_cmp in self.shared_compounds:
            for rxn in shared_cmp.reactions:
                rxn.store_flux = True
                
        # Check if each rxn has the field 'model' assigned
        for member in self.community_members:
            no_model_rxns =  [r for r in member.reactions if not hasattr(r,'model') or (hasattr(r,'model') and r.model == None)]
            if len(no_model_rxns) > 0:
                raise userError('No "model" has been assigned to reactions ' + str([r.id for r in no_model_rxns]) + ' related to community member ' + member.id)

        # Additoinal arguments. Additional arguments should be entered as normal but they are 
        # converted to a dictionary whose keys are the names of the arguments and values are 
        # the values of  those arguments
        argnames = additional_args.keys()
        argvals = additional_args.values()
        for argname in argnames:
           exec "self." + argname + " = " +"additional_args['" + argname + "']"

        # A parameter showing a simplistic time update methods should be used (True) or not. The simplistic method
        # invovles setting all negative compound concentrations to zero if such a thing occurs. If this parameter is
        # False, then the time interval is divided into smaller ones to avoid zero conetrations.
        if not isinstance(use_simplistic_time_update,bool):
            raise TypeError('use_simplistic_time_update must be either True or False')
        self.use_simplistic_time_update = use_simplistic_time_update

        # A parameter showing whether a time interval,which has already been divided into smaller ones should be 
        # divided again to avoid negative compounds concentrations (True) or not. Note that setting this parameter to
        # True may significanlty increase the runtime 
        self.redivide_time_interval = False

    def __setattr__(self,attr_name,attr_value):
        """
        Redefines funciton __setattr__
        INPUTS:
        -------
        attr_name: Attribute name
        attr_value: Attribute value
        """
        # lage phase time
        if attr_name.lower() == 'lag_phase_time' and not isinstance(attr_value,int) and not isinstance(attr_value,float):
            raise TypeError('lag_phase_time must be non-negative integer or float')
        elif attr_name.lower() == 'lag_phase_time' and attr_value < 0:
            raise ValueError('lah_phase_time must be a non-negative number') 

        self.__dict__[attr_name] = attr_value
         

    def uptake_rate_UB_calc(self):
        """ 
        Compute the upper bound on the uptake rates of the shared compounds 
        (LB on exchange fluxes) using kinetic expressions 
        """ 
        for shared_cmp in self.shared_compounds:
            for uptake_reaction in shared_cmp.reactant_reactions:         
                # Assume uptake reactions are exchange reactions with only one
                # reactant 
                uptake_reaction.reactants[0].concentration = shared_cmp.concentration[self._t] 
                uptake_reaction.kinetic_rate_calc(assignLB = True) 

    def update_fba_model(self):
        """ 
        This is an abstract method that can be defined by subclasses. This is useful
        for cases one needs to solve a different fba model (i.e., different objective
        function or different constriants) at each time point
        """ 
        pass 

    def mu_uptake_export_calc(self):
        """ 
        Compute the specific growth rate of community members and the uptake
        and export rates of shared compounds 
        """ 
        #-- Compute mu and uptake and export rate for any member who has not driven to extinsion --
        for member in [m for m in self.community_members if (m.organism.gDW_per_ml != None and m.organism.gDW_per_ml[self._t] > 0) or (m.organism.cells_per_ml != None and m.organism.cells_per_ml[self._t] > 0)]:

            # Lag phase
            if self._t < self.lag_phase_time:
                member.organism.mu[self._t] = member.organism.random_mortality_rate

                # Set the flux of exchange reactions for this member related to all shared compounds to zero
                for rxn in [r for c in self.shared_compounds for r in c.reactions if r in member.reactions]: 
                        rxn.flux[self._t] = 0

            # After lag phase
            else:
                if self._t == self._t0:
                    member.fba(build_new_optModel = True, reset_fluxes = True, store_opt_fluxes= True, flux_key = self._t, stdout_msgs = False)
                else:
                    member.fba(build_new_optModel = False, reset_fluxes = False, store_opt_fluxes = True, flux_key = self._t, stdout_msgs = False)
    
                if member.fba_model.solution['exit_flag'] == 'globallyOptimal':
                    member.organism.mu[self._t] = member.biomass_reaction.flux[self._t] + member.organism.random_mortality_rate
                    if self.stdout_msgs and self.organism.mu[self._t] < 0:
                        print '   Death rate for {} due to  having a growth rate less than random death rate ...'.format(member.id)
                else:
                    if self.stdout_msgs:
                        print '   Death rate for {} due to resource limitation...'.format(member.id)
                    # Find the limiting nutrients that this member takes up 
                    limiting_nutrients = [uptake_rxn.reactants[0] for shared_cmp in self.shared_compounds for uptake_rxn in shared_cmp.reactant_reactions if uptake_rxn in member.reactions]
                    if hasattr(member,'wildType_max_biomass'):
                        mu_death = death_rate(limiting_nutrients = limiting_nutrients,ms_biomassYield_calc = False, wildType_max_biomass = member.wildType_max_biomass, warnings = False, stdout_msgs = False) 
                    else:
                        mu_death = death_rate(limiting_nutrients = limiting_nutrients,ms_biomassYield_calc = False, wildType_max_biomass = None, warnings = False, stdout_msgs = False) 
                    member.organism.mu[self._t] = mu_death + member.organism.random_mortality_rate
    
                    # Set the flux of exchange reactions for this member related to all shared compounds to zero
                    for rxn in [r for c in self.shared_compounds for r in c.reactions if r in member.reactions]: 
                            rxn.flux[self._t] = 0

        #-- Set mu and uptake and export fluxes to zero for any member that has been extincted --
        for member in [m for m in self.community_members if (m.organism.gDW_per_ml != None and m.organism.gDW_per_ml[self._t] == 0) or (m.organism.cells_per_ml != None and m.organism.cells_per_ml[self._t] == 0)]:
            member.organism.mu[self._t] = 0 
            for rxn in [r for c in self.shared_compounds for r in c.reactions if r in member.reactions]: 
                    rxn.flux[self._t] = 0

    def update_concentrations_batch(self):
        """ 
        Updates cell concentrationns and the concentration of extracellular compounds
        for a batch reactor 
        """  
        #--- Update the cell concentrations ---
        # dX_i/dt = mu_i*X_i*(1 - sum(i,(1-rmp/100)X(i))/carrying_capacity) or 
        # (X_i(t+dt) - X_i(t))/dt = mu*X_i(t)*(1 - sum(i,(1-rmp/100)*X_i(t))/carrying_capacity)
        # If concentration is negative set it to zero
        members_gDW_per_ml_total = sum([member.organism.gDW_per_ml[self._t] for member in self.community_members])
        self._logisticG_factor_gDW_per_ml = 1 - members_gDW_per_ml_total/self.carrying_capacity['gDW_per_ml']
        if len([member for member in self.community_members if member.organism.cells_per_ml != None]) == len(self.community_members):
            members_cells_per_ml_total = sum([member.organism.cells_per_ml[self._t] for member in self.community_members])
            self._logisticG_factor_cells_per_ml = 1 - members_cells_per_ml_total/self.carrying_capacity['cells_per_ml']

        for member in self.community_members:
            # We always need gDW_per_ml to update compound concentrations but
            # providing cells_per_ml is optional
            member.organism.gDW_per_ml[self._t + self._dt] = max(member.organism.mu[self._t]*member.organism.gDW_per_ml[self._t]*self._logisticG_factor_gDW_per_ml*self._dt + member.organism.gDW_per_ml[self._t],0)

            if member.organism.cells_per_ml is not None:
                member.organism.cells_per_ml[self._t + self._dt] = max(member.organism.mu[self._t]*member.organism.cells_per_ml[self._t]*self._logisticG_factor_cells_per_ml*self._dt + member.organism.cells_per_ml[self._t],0)

        #--- Update shared compound concentrations ---
        # dC/dt = f*(1 - total_cmps_concs/carrying_capacity['compounds_mM']) 
        # where, f = sum(k,v_export_k*X_k) - sum(k,v_uptake_k*X_k)
        # (C(t+dt) - c(t))/dt = sum(k,v_export_k*X_k) - sum(k,v_uptake_k*X_k)
        # Total concentration of all shared compounds
        # Note: Since X_k are in gDW/ml they should be multipled by 1000 to be converted to gDW/l
        total_cmps_conc = sum([cmp.concentration[self._t] for cmp in self.shared_compounds])
        self._logisticG_factor_cmps = 1 - total_cmps_conc/self.carrying_capacity['compounds_mM']
        self._f = dict([(cmp,None) for cmp in self.shared_compounds])
        for shared_cmp in self.shared_compounds:
            f = sum([r.flux[self._t]*1000*r.model.organism.gDW_per_ml[self._t] for r in shared_cmp.reactions])
            self._f[shared_cmp] = f

            conc = f*self._logisticG_factor_cmps*self._dt + shared_cmp.concentration[self._t]

            if conc >= 0 or (conc < 0 and abs(conc) <= 1e-9):
                conc = max(conc,0)

            shared_cmp.concentration[self._t + self._dt] = conc


    def create_dt_frac(self):
        """
        This function creates the fractional time intervals to avoid getting 
        negative concentrations. This is done by divding the original time 
        interval into smaller ones, if any of the compound concentrations turns out 
        to be negative. This information is stored in dt_zero_conc. For example, the 
        original dt could be one, and dt's stored in dt_zero_conc can be 0.3 and 0.7. 
        Then we should do one simulation with a dt of 0.3, a second simulation with 
        dt of 0.7 - 0.3 = 0.4 and a last simulation with dt of 1 - 0.7 = 0.3
        """
        # A dictionary with keys and values as follows:
        #   Keys: compounds that get a negative concentration upon update
        # Values: The time interval needed for the concentration of compounds 
        #         to become zero
        dt_zero_conc = {}

        for shared_cmp in self._negative_conc_cmps: 
            # find the time interval needed to get a zero concentration
            # C(t + dt) = f*logisticG_factor*dt + C(t) = 0  ==> dt = (0 - C(t))/(f*logisticG_factor)
            dt_zero_conc_curr = (0 - shared_cmp.concentration[self._t])/(self._f[shared_cmp]*self._logisticG_factor_cmps)
    
            if dt_zero_conc_curr < 0:
                raise userError('\nNegative dt for compound ' + shared_cmp.id + '  (dt = ' + str(dt_zero_conc_curr) + '  conc(t) = ' + str(shared_cmp.concentration[self._t]) + '  f = ' + str(self._f[shared_cmp]) + '\n')
            elif dt_zero_conc_curr > self._dt:
                raise userError('The obtained modified dt for the shared compound ' + shared_cmp.id + ' (dt = ' + str(dt_zero_conc_curr) + ') is smaller than the original time interval (dt = ' + str(self._dt) + ')')
            else:
                dt_zero_conc[shared_cmp.id] = dt_zero_conc_curr

        # Add zero and the original dt to the list and sort it
        #**all_dts = sorted(dt_zero_conc.keys() + [0,self._dt])
        dt_zero_conc['start'] = 0
        if self._t in self._time_points:
            dt_zero_conc['end'] = self._dt_orig
        # Otherwise, find out between which two time points are now
        else:
            # The index of the smallest element of self._time_points, which is bigger than self._t 
            t_ind_timepoints = [i for i in range(1,len(self._time_points)) if self._time_points[i-1] < self._t and self._time_points[i] > self._t]
            if len(t_ind_timepoints) == 1:
                # Example: self._time_points[0,0.5,1,1.5,2] and self._t = 1.2, then 
                # dt_zero_conc['end'] = 1.5 - 1.2 = 0.3
                dt_zero_conc['end'] = self._time_points[t_ind_timepoints[0]] - self._t
            elif len(t_ind_timepoints) == 0:
                raise userError('t = ' + str(self._t) + ' is not between any two elements of self._time_points = ' + str(self._time_points))
            elif len(t_ind_timepoints) > 1:
                raise userError('Two indices found in self._time_points (' + str() + ') when trying to find out t = ' + str(self._t) + ' is between which two elements of self._time_points = ' + str(self._time_points))
 
        # Find out smaller dt's required for simulations (a fraction of the original dt)
        dt_zero_conc_sorted = sorted(dt_zero_conc.values())
        self._dt_fracs = [dt_zero_conc_sorted[i] - dt_zero_conc_sorted[i-1] for i in range(1,len(dt_zero_conc))]
        # Here we consider only the first element (smallest dt). This will resolve the negative concentraiton
        # only for the corresponding metabolites. The requires dt for other metabolites with negative 
        # conetrations should be updated after this udpate for dt (see funciton update_time)

        if self.stdout_msgs:
            print '\tDividing the original time interval into smaller ones to avoid getting negative concentrations ...'
            if False:
                print '\tdt_zero_conc = ',sorted([(shared_cmp,dt_zero_conc[shared_cmp]) for shared_cmp in dt_zero_conc.keys()],key=lambda x:x[1])
                print '\tdt_fracs = ',[([k for k in dt_zero_conc.keys() if dt_zero_conc[k] == dt_zero_conc_sorted[i]][0],dt_zero_conc_sorted[i] - dt_zero_conc_sorted[i-1]) for i in range(1,len(dt_zero_conc))]

    def update_time(self):
        """
        This function updates the time point
        """   
        # Shared compounds with negative concentrations
        self._negative_conc_cmps = [shared_cmp for shared_cmp in self.shared_compounds if shared_cmp.concentration[self._t + self._dt] < 0]

        # Update time, if there is no negative concentrations
        if len(self._negative_conc_cmps) == 0:
            self._t += self._dt
            # If too close to an entry in self._time_points set to that entry (to avoid numerical problems)
            if self._dt != self._dt_orig and self._t not in self._time_points and len([t for t in self._time_points if abs(t - self._t) <= 1e-6]) == 1:
                time_point = [t for t in self._time_points if abs(t - self._t) <= 1e-6][0] 
                for member in self.community_members:
                    member.organism.cells_per_ml[time_point] = member.organism.cells_per_ml[self._t] 
                    member.organism.gDW_per_ml[time_point] = member.organism.gDW_per_ml[self._t] 
                    del member.organism.cells_per_ml[self._t] 
                    del member.organism.gDW_per_ml[self._t] 

                for shared_cmp in self.shared_compounds:
                     shared_cmp.concentration[time_point] = shared_cmp.concentration[self._t] 
                     del shared_cmp.concentration[self._t] 
                self._t = time_point 
            if self.stdout_msgs:
                print 't = ',self._t

            #--- Sepcify the correct dt for the next time step, if needed ---
            if not hasattr(self,'_dt_fracs') or (hasattr(self,'_dt_fracs') and len(self._dt_fracs) == 0):
                self._dt = self._dt_orig

            # If the current fractional dt has already been examined remove it 
            # from dt_fracs
            elif len(self._dt_fracs) > 0 and self._dt in self._dt_fracs:
                del self._dt_fracs[self._dt_fracs.index(self._dt)]
                if len(self._dt_fracs) > 0:
                    self._dt = self._dt_fracs[0]
                else:
                    self._dt = self._dt_orig

            elif len(self._dt_fracs) > 0 and self._dt not in self._dt_fracs:
                 raise userError('Unknown dt (dt = ' + str(self._dt) + ') not in dt_fracs = ' + str(self._dt_fracs))
            
        # If there are negative concentrations
        else:
            # Create fractional dt's it hasn't been already created. Do not update time in this case
            if self._t in self._time_points: 
                # First remove all previously computed values for self._t + self._dt
                for member in self.community_members:
                    del member.organism.cells_per_ml[self._t + self._dt] 
                    del member.organism.gDW_per_ml[self._t + self._dt] 

                for shared_cmp in self.shared_compounds:
                     del shared_cmp.concentration[self._t + self._dt] 

                self.create_dt_frac()

                # Use dt_frac instead of dt_orig
                self._dt = self._dt_fracs[0]

            # If the current fractional dt has already been examined, we shouldn't get a negative concentraiton
            # in most cases. This, however, may sometimes happen, if there is an active produciton source for that
            # compound. Two strategies can be taken. The simplistic one setting all metabolites with a negative 
            # concetration to zero and move to the next dt. The other method is to divide the current time interal
            # into smaller ones. 
            elif len(self._dt_fracs) > 0 and self._dt in self._dt_fracs:
                if not self.redivide_time_interval:
                    for shared_cmp_ng in self._negative_conc_cmps: 
                        shared_cmp_ng.concentration[self._t + self._dt] = 0

                    # Update time
                    self._t += self._dt
                    # If too close to an entry in self._time_points set to that entry (to avoid numerical problems)
                    if self._dt != self._dt_orig and self._t not in self._time_points and len([t for t in self._time_points if abs(t - self._t) <= 1e-6]) == 1:
                        time_point = [t for t in self._time_points if abs(t - self._t) <= 1e-6][0] 
                        for member in self.community_members:
                            member.organism.cells_per_ml[time_point] = member.organism.cells_per_ml[self._t] 
                            member.organism.gDW_per_ml[time_point] = member.organism.gDW_per_ml[self._t] 
                            del member.organism.cells_per_ml[self._t] 
                            del member.organism.gDW_per_ml[self._t] 
        
                        for shared_cmp in self.shared_compounds:
                             shared_cmp.concentration[time_point] = shared_cmp.concentration[self._t] 
                             del shared_cmp.concentration[self._t] 
                        self._t = time_point 
                    if self.stdout_msgs:
                        print 't = ',self._t

                    # Then remove dt from dt_fracs
                    del self._dt_fracs[self._dt_fracs.index(self._dt)]

                    if len(self._dt_fracs) > 0:
                        self._dt = self._dt_fracs[0]
                    else:
                        self._dt = self._dt_orig
 
                # Divide the current time interval again to smaller ones and do not update time
                else:
                    for member in self.community_members:
                        del member.organism.cells_per_ml[self._t + self._dt] 
                        del member.organism.gDW_per_ml[self._t + self._dt] 
    
                    for shared_cmp in self.shared_compounds:
                         del shared_cmp.concentration[self._t + self._dt] 
    
                    self.create_dt_frac()
    
                    # Use dt_frac instead of dt_orig
                    self._dt = self._dt_fracs[0]

            else:
                raise userError('Unknown case happened at t = ' + str(self._t) + '! self._t in self._time_points = ' + str(self._t in self._time_points) + 'Check DMMM.py for more details')

    def dilute(self):
        """ 
        Performs dilution (used when reactor_type is serial_dilution 
        """ 
        #--- Cell concentrations are divided by the dilution factor ---
        for member in self.community_members:
            member.organism.gDW_per_ml[self._t] = member.organism.gDW_per_ml[self._t]/self.serial_dilution_params['dilution_factor']
            if member.organism.cells_per_ml is not None:
                member.organism.cells_per_ml[self._t] = int(member.organism.cells_per_ml[self._t]/self.serial_dilution_params['dilution_factor']) 

        #--- shared metabolite concentrations are set to those at the initial time point ---
        for shared_cmp in self.shared_compounds:
            shared_cmp.concentration[self._t] = shared_cmp.concentration[self._t0]

        # lag phase time
        self.lag_phase_time = self._t + self._lag_phase_time_orig

    def run_batch(self):
        """ 
        This function runs dynamic FBA ssimulations for a batch culture 
        """
        while self._t < self._tf:
            # Compute the upper bound on the uptake rates of 
            # the shared compounds using kinetic expressions
            self.uptake_rate_UB_calc()

            # Update the fba model for the initial point
            self.update_fba_model()

            # Compute the specific growth rate (mu) for each species using FBA 
            # as well as the export rates 
            self.mu_uptake_export_calc()

            # Update concentrations in the next time point 
            self.update_concentrations_batch()

            # Update time
            self.update_time()

    def run_serial_dilution(self):
        """ 
        This function runs dynamic FBA simulations for a batch culture with serial dilutions 
        """
        # Time points where dilutioning is performed. The final time should not be included. 
        # For example, if t0 = 0, tf = 118 h and 
        # serial_dilution_params['time_between_dilutions'] = 24 h, then 
        # dilution_times = [34, 58, 82, 106]. dilution_times may or may not include the last time point
        self._dilution_times = [self._t0 + self.serial_dilution_params['time_between_dilutions']*f for f in range(1,int((self._tf - self._t0)/self.serial_dilution_params['time_between_dilutions']) + 1)]
        if self._tf in self._dilution_times:
            del self._dilution_times[self._dilution_times.index(self._tf)]
         

        while self._t < self._tf:
            # Perform dilutioning
            # The second condition imposed in the following abs(t - self._t) <= 1e-6 is to take care of cases
            # where we make an adjustment to delta_t and self._t may not be exactly equal to the time points 
            # in self._dilution_times
            #if self._t in self._dilution_times or len([t for t in self._dilution_times if abs(t - self._t) <= 1e-6]) == 1:
            if self._t in self._dilution_times:
                if self.stdout_msgs:
                    print '\tDilutioning at time t = {} ...'.format(self._t)
                self.dilute()

            # Compute the upper bound on the uptake rates of 
            # the shared compounds using kinetic expressions
            self.uptake_rate_UB_calc()

            # Update the fba model for the initial point
            self.update_fba_model()

            # Compute the specific growth rate (mu) for each species using FBA 
            # as well as the export rates 
            self.mu_uptake_export_calc()

            # Update concentrations in the next time point 
            self.update_concentrations_batch()

            # Update time
            self.update_time()

    def run_life_cycle(self):
        """ 
        This function runs dynamic FBA ssimulations for microbial life cycle 
        """
        pass

    def run_chemostat(self):
        """ 
        This function runs FBA ssimulations for a chemostat 
        """
        pass

    def run(self):
        """ 
        This function runs the dynamic simulations. 
        """

        if self.stdout_msgs == 'on':
            print 'Start running DMMM ...'

        start_DMMM = time.clock()

        self._t = self._t0

        # Original delta_t 
        self._dt_orig = self._dt

        if self.reactor_type.lower() == 'batch':
            self.run_batch()
        elif self.reactor_type.lower() == 'serial_dilution':
            self.run_serial_dilution()
        elif self.reactor_type.lower() == 'life_cycle':
            self.run_life_cycle()
        elif self.reactor_type.lower() == 'chemostat':
            self.run_chemostat()

        # Compute the update rates and mu at the last time point in case they are 
        # needed to be reported 
        self.uptake_rate_UB_calc()
        self.update_fba_model()
        self.mu_uptake_export_calc()

        # Elapsed time to run DMMM
        elapsed_DMMM = str(timedelta(seconds = time.clock() - start_DMMM))

        if self.stdout_msgs:
           print '\nelapsed time (sec) for DMMM = ',elapsed_DMMM,'\n'

#---------- Test DMMM ---------------
if __name__ == "__main__":
    from copy import deepcopy
    from tools.io.read_sbml_model import read_sbml_model
    from set_specific_bounds import set_specific_bounds

    # Increse the recursion limit, otherwise deepcopy will complain
    sys.setrecursionlimit(10000)

    # Model path
    model_path = '/data/alizom/models/Escherichia_coli/iAF1260/'

    # Define the organism
    model_organism = organism(id = 'Ecoli', name = 'Escherichia coli',domain = 'Bacteria', genus = 'Escherichia', species = 'coli', strain = 'MG1655')

    iAF1260 = read_sbml_model(file_name = model_path + 'iAF1260_updated.xml', model_id = 'iAF1260',model_organism = model_organism,model_type = 'metabolic')
    
    #--- E. coli iAF1260 model ---
    print '\n--- Wild-type E.coli (iAF1260 model) ----'
    iAF1260.biomass_reaction = iAF1260.get_reactions({'Ec_biomass_iAF1260_core_59p81M':'id'}) 
    iAF1260.all_biomass_reactions = {'core':iAF1260.get_reactions({'Ec_biomass_iAF1260_core_59p81M':'id'}),'WT':iAF1260.get_reactions({'Ec_biomass_iAF1260_WT_59p81M':'id'})} 
    iAF1260.organism.gDW_per_cell = 2.8e-13
 
    # Assign a general Michaelis-Menten type uptake kinetics to all exchange reactions
    # Example: EX_glc(E): glc__D_e <==>    Vmax*C['glc__D_e']/(Km + C['glc__D_e']) 
    # Use a Vmax value of 10 mmole/gDW.h and a Km value of 10 micro-M 
    for reaction in [r for r in iAF1260.reactions if r.type.lower() == 'exchange']:
        # The id of compound participating in the exchange reaction
        metab_id = [m.id for m in reaction.compounds][0]
        reaction.kinetics = "10*C['" + metab_id + "']/(10 + C['" + metab_id + "'])"

    # Glucose uptake kinetics 
    exch_rxns = iAF1260.get_reactions({'EX_glc(e)':'id','EX_lys_L(e)':'id','EX_ile_L(e)':'id'})
    exch_rxns['EX_glc(e)'].kinetics = "10*C['glc__D_e']/(0.15 + C['glc__D_e'])"
    exch_rxns['EX_lys_L(e)'].kinetics = "0.1964*C['lys__L_e']/(5e-4 + C['lys__L_e']) + 0.3055*C['lys__L_e']/(1e-2 + C['lys__L_e'])"
    exch_rxns['EX_ile_L(e)'].kinetics = "0.0346*C['ile__L_e']/(1.22e-3 + C['ile__L_e'])"
   
    # Growth medium
    set_specific_bounds(model = iAF1260,file_name = model_path + 'iAF1260_minimal_glucose_aerobic.py',simulation_condition = 'minimal_glucose_aerobic')

    # Assign and objective function coefficients
    for rxn in iAF1260.reactions:
        rxn.objective_coefficient = 0
    iAF1260.biomass_reaction.objective_coefficient = 1

    # Perform FBA for the wild-type
    iAF1260.fba(assign_wildType_max_biomass = True)

    # Compute ms and biomass yield for glucose
    glc_D = iAF1260.get_compounds({'glc__D_e':'id'})
    glc_D.ms_calc() 
    glc_D.biomass_yield_calc() 
    print '\n  **glc-D ms = ',glc_D.ms
    print '  **glc-D biomass yield = ',glc_D.biomass_yield,'\n'
    # Check if the bounds and objective function were set back to original values
    # after calculating ms and biomass yield
    iAF1260.fba()

    #--- DMMM for the co-culture of lysA and ilvE mutants ---
    print '\n--- lysA mutant ----'
    # lysA mutant and related reaction whose flux should be set to zero 
    iAF1260_lysA = deepcopy(iAF1260)
    iAF1260_lysA.organism.id = 'lysA_Ecoli'
    iAF1260_lysA.organism.cells_per_ml = {0:7.5e6}
    iAF1260_lysA.organism.gDW_per_ml = {0:7.5e6*iAF1260_lysA.organism.gDW_per_cell}
    DAPDC = iAF1260_lysA.get_reactions({'DAPDC':'id'})
    DAPDC.flux_bounds = [0,0]
    iAF1260_lysA.fba(build_new_optModel = False)
    # Compute ms and biomass yield for lys_L
    lys_L = iAF1260_lysA.get_compounds({'lys__L_e':'id'})
    lys_L.ms_calc()    
    lys_L.biomass_yield_calc()    
    print '\n  **lys_L ms = ',lys_L.ms
    print '  **lys_L biomass yield = ',lys_L.biomass_yield,'\n'
    # Check if the bounds and objective function were set back to original values
    # after calculating ms and biomass yield
    iAF1260_lysA.fba()

    print '\n**Test cooperative behavior ....'
    for rxn in iAF1260_lysA.reactions:
        rxn.objective_coefficient = 0
    iAF1260_lysA.get_reactions({'EX_ile_L(e)':'id'}).objective_coefficient = 1
    iAF1260_lysA.fba(build_new_optModel = False)
    iAF1260_lysA.biomass_reaction.objective_coefficient = 1

    # ilvE mutant and related reaction whose flux should be set to zero 
    print '\n--- ilvE mutant ----'
    iAF1260_ilvE = deepcopy(iAF1260)
    iAF1260_ilvE.organism.id = 'ilvE_Ecoli'
    iAF1260_ilvE.organism.cells_per_ml = {0:7.5e6}
    iAF1260_ilvE.organism.gDW_per_ml = {0:7.5e6*iAF1260_ilvE.organism.gDW_per_cell}
    ILETA = iAF1260_ilvE.get_reactions({'ILETA':'id'})
    VALTA = iAF1260_ilvE.get_reactions({'VALTA':'id'})
    ILETA.flux_bounds = [0,0]
    VALTA.flux_bounds = [0,0]
    iAF1260_ilvE.fba(build_new_optModel = False)
    # Compute ms and biomass yield for ile_L
    ile_L = iAF1260_ilvE.get_compounds({'ile__L_e':'id'})
    ile_L.ms_calc()    
    ile_L.biomass_yield_calc()    
    print '\n  **ile_L ms = ',ile_L.ms
    print '  **ile_L biomass yield = ',ile_L.biomass_yield,'\n'
    # Check if the bounds and objective function were set back to original values
    # after calculating ms and biomass yield
    iAF1260_ilvE.fba()

    print '\n**Test cooperative behavior ....'
    for rxn in iAF1260_ilvE.reactions:
        rxn.objective_coefficient = 0
    iAF1260_ilvE.get_reactions({'EX_lys_L(e)':'id'}).objective_coefficient = 1
    iAF1260_ilvE.fba(build_new_optModel = False)
    iAF1260_ilvE.get_reactions({'EX_lys_L(e)':'id'}).objective_coefficient = 0
    iAF1260_ilvE.biomass_reaction.objective_coefficient = 1

    # --- Define the compounds available in the extracellular medium ---
    # These compounds include glucose, arg-L and ile_L
    # Glucose
    exch_rxns_lysA = iAF1260_lysA.get_reactions({'EX_glc(e)':'id','EX_lys_L(e)':'id','EX_ile_L(e)':'id'}) 
    exch_rxns_ilvE = iAF1260_ilvE.get_reactions({'EX_glc(e)':'id','EX_lys_L(e)':'id','EX_ile_L(e)':'id'}) 

    # Glucose as a shared compound (concentration in mM)
    glucose = compound(id = 'glc-D', name = 'D-Glucose', Kegg_id = 'C00031', reactant_reactions = [exch_rxns_lysA['EX_glc(e)'],exch_rxns_ilvE['EX_glc(e)']],reactions = [exch_rxns_lysA['EX_glc(e)'],exch_rxns_ilvE['EX_glc(e)']],concentration = {0:111.01})    
    # Lysine as a shared compound (concentration in mM)
    lys_L = compound (id = 'lys_L',name = 'L-Lysine', Kegg_id = 'C00047', reactant_reactions = [exch_rxns_lysA['EX_lys_L(e)']],product_reactions = [exch_rxns_ilvE['EX_lys_L(e)']], concentration = {0:1e-5}) 
   
    # Isoleucine as a shared compound (concentration in mM)
    ile_L = compound(id = 'ile_L', name = 'L-Isoleucine' , Kegg_id = 'C00407', reactant_reactions = [exch_rxns_ilvE['EX_ile_L(e)']], product_reactions = [exch_rxns_lysA['EX_ile_L(e)']], concentration = {0:1e-5})

    DMMM_test = DMMM(community_members = [iAF1260_lysA,iAF1260_ilvE],shared_compounds = [glucose,lys_L,ile_L],time_range = [0,0.5,10],reactor_type = 'batch',stdout_msgs = True)
    #DMMM_test = DMMM(community_members = [iAF1260_lysA,iAF1260_ilvE],shared_compounds = [glucose,lys_L,ile_L],time_range = [0,1.0,96],reactor_type = 'serial_dilution',serial_dilution_params = {'dilution_factor':1000,'time_between_dilutions':24},stdout_msgs = True)
    DMMM_test.run()
   
    # Print results in the output
    for t in sorted(glucose.concentration.keys()):
        print 't = ',t,
        for member in [iAF1260_lysA,iAF1260_ilvE]:
            print '  X(',member.organism.id,') = ',member.organism.cells_per_ml[t],

        for shared_cmp in [glucose,lys_L,ile_L]:
            print '   C(',shared_cmp.id,') = ',shared_cmp.concentration[t],

        print '\n'

