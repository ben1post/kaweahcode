from lmfit import minimize, Parameters # Parameter, report_fit
#


import numpy as np
import xsimlab as xs
from scipy.integrate import odeint
import matplotlib.pyplot as plt



@xs.process
class UniqueComponent:
    """This is for initializing a Nutrient (i.e. components that are unique (0D))"""
    label = xs.variable(intent='out', groups='label')
    init = xs.variable(intent='out', groups='init')
    
    def initialize(self):
        self.label = (self.component_label,)
        self.init = np.array(self.initVal)
        print('NutInit',self.init)
        
        
@xs.process
class Nutrient(UniqueComponent):
    component_label = xs.variable(default='N', groups='component-label')
    initVal = xs.variable(default=1)
    supplyrate = xs.variable(default=1)

    # Process variables
    supply = xs.variable(intent='out')
    
    def _supply(self):
        return self.supplyrate
    
    def initialize(self):
        super(Nutrient, self).initialize()
        self.supply = self._supply
        
        
@xs.process
class FlexSizeComponent:
    """This is for initializing both Phytoplankton and Zooplankton (i.e. components of variable number (1D))"""
    #flux = xs.variable(intent='out', groups='flux')
    label = xs.variable(intent='out', groups='label')
    init = xs.variable(intent='out', groups='init')
    size = xs.variable(intent='out', groups='size')
        
    def calculate_sizes(self):
        # must be implemented in subclasses
        raise NotImplementedError
        
    def initialize_alloparams(self):
        # must be implemented in subclasses
        raise NotImplementedError
        
    def initialize(self):
        self.label = tuple((self.component_label + str(i) for i in range(self.NP)))
        self.init = np.array([self.initVal / self.NP for i in range(self.NP)])
        
        self.size = self.calculate_sizes()

        self.initialize_alloparams()
        
        
        
@xs.process
class Phytoplankton(FlexSizeComponent):
    component_label = xs.variable(default='P', groups='component-label')
    # Phytoplankton params
    PMinEsd = xs.variable(default = 1,  description='P max growth', groups='parameter', static=True)
    PMaxEsd = xs.variable(default = 20, description='P max growth', groups='parameter', static=True)
    NP =      xs.variable(default = 40, description='number of P',  groups='parameter', static=True)
    
    # INPUT
    initVal = xs.variable()
    m = xs.variable(default=.1)
    
    # Rate params
    ks = xs.variable(intent='out', description='allometric half-saturation constant',  groups='parameter')
    mu0 = xs.variable(intent='out', description='allometric growth rate',  groups='parameter')
    
    # Process variables
    uptake = xs.variable(intent='out')
    mortality = xs.variable(intent='out')
    grazed = xs.variable(intent='out')
    
    def _uptake(self, N):
        return self.mu0 * N / (N + self.ks)
    
    def _mortality(self):
        """returns size-scaled mortality rate, calculated from mu0"""
        return self.m * self.mu0 
    
    def _grazed(self,FgrazP):
        return [sum(FgrazP[i,:]) for i in range(self.NP)]
    
    def calculate_sizes(self):
        numbers = np.array([i for i in range(self.NP)])
        sizes = (np.log(self.PMaxEsd) - np.log(self.PMinEsd))* numbers / (self.NP-1) + np.log(self.PMinEsd)
        return np.exp(sizes)
    
    def initialize_alloparams(self):
        # must be implemented in subclasses
        self.mu0 = 2.6 * (self.size) ** -0.45
        self.ks = (self.size) * .1
        
    def initialize(self):
        super(Phytoplankton, self).initialize()
        self.uptake = self._uptake
        self.mortality = self._mortality
        self.grazed = self._grazed

        
@xs.process
class Zooplankton(FlexSizeComponent):
    component_label = xs.variable(default='Z', groups='component-label')
    # SetupParams
    PMinEsd = xs.foreign(Phytoplankton, 'PMinEsd')
    PMaxEsd = xs.foreign(Phytoplankton, 'PMaxEsd')
    NP =      xs.foreign(Phytoplankton, 'NP')
    phytosize = xs.foreign(Phytoplankton, 'size')
   
    # Zooplankton params 
    initVal = xs.variable()
    zeta = xs.variable(default=1) #mortality rate quadratic
    deltaxprey = xs.variable(default=0.25) # log10 prey size tolerance
    KsZ = xs.variable(default=3) # grazing half saturation constant
    f_eg = xs.variable(default=.33) # egested food
    epsilon = xs.variable(default=.33) # assimilated food
    
    
    # Alloparams
    I0 = xs.variable(intent='out')
    xpreyopt = xs.variable(intent='out')
    
    phiP = xs.variable(intent='out')
    
    # Process variables
    grazingmatrix = xs.variable(intent='out')
    mortality = xs.variable(intent='out')
    ingestion = xs.variable(intent='out')
    excretion = xs.variable(intent='out')
    
    def _grazingmatrix(self,P,Z):
        PscaledAsFood = np.zeros((self.NP,self.NP))
        for j in range(self.NP):
            for i in range(self.NP):
                PscaledAsFood[i,j] = self.phiP[i,j] / self.KsZ * P[i]
        
        FgrazP = np.zeros((self.NP,self.NP))
        for j in range(self.NP):
            for i in range(self.NP):        
                FgrazP[i,j] = self.I0[j] * Z[j] * PscaledAsFood[i,j] / (1 + sum(PscaledAsFood[:,j]))
                
        return FgrazP
    
    def _ingestion(self,FgrazP):
        return [self.epsilon * sum(FgrazP[:,j]) for j in range(self.NP)]

    def _excretion(self,FgrazP):
        return [(1 - self.f_eg - self.epsilon) * sum(FgrazP[:,j]) for j in range(self.NP)]
    
    def _mortality(self,Z):
        return self.zeta * sum(Z)
        
    def calculate_sizes(self):
        zoosizes= 2.16 * self.phytosize ** 1.79
        return zoosizes
    
    def initialize_alloparams(self):
        # initializes allometric parameters as lists, based on sizes
        self.I0 = 26 * (self.size) ** -0.4
        self.xpreyopt = 0.65 * (self.size) ** .56 
    
    def init_phiP(self):
        """creates array of feeding preferences [P...P10] for each [Z]"""
        phiP= np.array([[np.exp(-((np.log10(xpreyi)-np.log10(xpreyoptj)) / self.deltaxprey)**2) 
               for xpreyi in self.phytosize] for xpreyoptj in self.xpreyopt])
        return phiP
    
    def initialize(self):        
        super(Zooplankton, self).initialize()
        self.grazingmatrix = self._grazingmatrix
        
        self.mortality = self._mortality
        self.ingestion = self._ingestion
        self.excretion = self._excretion
        
        self.phiP = self.init_phiP()

import itertools

def flatten(generatorlist):
    # returns 1D list from nested generator or multi-D list
    return list(itertools.chain.from_iterable(generatorlist))



@xs.process
class Chemostat:
    # foreign processes
    N_supply = xs.foreign(Nutrient, 'supply')
    
    P_uptake = xs.foreign(Phytoplankton, 'uptake')
    P_mortality = xs.foreign(Phytoplankton, 'mortality')
    
    grazingmatrix = xs.foreign(Zooplankton, 'grazingmatrix')
    Z_ingestion = xs.foreign(Zooplankton, 'ingestion')
    Z_excretion = xs.foreign(Zooplankton, 'excretion')
    Z_mortality = xs.foreign(Zooplankton, 'mortality')
    P_grazed = xs.foreign(Phytoplankton, 'grazed')
    
    # model construct labels
    complabels = xs.group('component-label')
    labels = xs.group('label')
    inits = xs.group('init')
    
    # output variables
    component = xs.variable(dims=('component'), intent='out')
    
    state = xs.variable(dims=('component'),intent='out')
    stateout = xs.variable(dims=('component'),intent='out')
    outflux = xs.variable(dims=('component'),intent='out')
    
    def initialize(self):
        self.time = 0

        # This is one of those things I am not sure about, it makes sense to use dicts to store model state (to allow calling different components by labels)
            # BUT this can not be used (afaik) in a nice way as output for plotting later on in the model output xarray.. do you have any idea how to do that more nicely?
        self.state = {label:val for label,val in zip(self.complabels,self.inits)}

        # in order to have output to plot, I also create an array that contains all output:
            # i have to use flatten function (defined above) to make group generators one dimensional (FlexComponents are Lists, so generators return List of Lists)
            # there should be a better solution to this.. but it relates to how xs.group() works in this model setup
        self.stateout = np.concatenate([self.state[label] for label in flatten(self.complabels)], axis=None)
        
        self.component = [i for i in flatten(self.labels)]
    
    def ode(self, state, time, labels, timestep):
        N, P, Z = state['N'], state['P'], state['Z']
        
        #P[P < 0.01] = 0.011
        #Z[Z < 0.01] = 0.011
        
        PGrazed = self.grazingmatrix(P,Z)

        dNdt = self.N_supply() - sum(self.P_uptake(N)*P) + sum(self.Z_excretion(PGrazed))
        
        dPdt = self.P_uptake(N)*P - self.P_mortality()*P - self.P_grazed(PGrazed)
        
        dZdt =  self.Z_ingestion(PGrazed) - self.Z_mortality(Z)*Z 
        
        # need to convert fluxes to correct time step by multiplying by dt        
        return {'N':dNdt*timestep, 'P':dPdt*timestep, 'Z':dZdt*timestep}

    @xs.runtime(args='step_delta')
    def run_step(self,dt):
        self.outflux = self.ode(self.state,self.time,self.labels, dt)
        # {label : self.funcs[label](self.state, self.time, self.labels) * dt for label in self.complabels}
        
        # keep track of timestep for forcing (+ other time dep. processes)
        self.time += dt
    
    def finalize_step(self):
        self.state = {label : self.state[label] + self.outflux[label] for label in self.complabels}
        
        self.stateout = np.concatenate(flatten([(self.state[label] for label in self.complabels)]), axis=None)


        
modmod = xs.Model({'env':Chemostat, 'N':Nutrient, 'P':Phytoplankton, 'Z':Zooplankton})


mod_in_x2 = xs.create_setup(
        model=modmod,
    clocks={
         'time': np.linspace(1,365*10,10*365*2)  # 20,1000)
     },
    input_vars={
        'P__initVal':1,
        'Z__initVal':1
    },
    output_vars={
        'env__component': None,
        'env__stateout': 'time'
    }
)

mod_in_x5 = xs.create_setup(
        model=modmod,
    clocks={
         'time': np.linspace(1,365*10,10*365*5)  # 20,1000)
     },
    input_vars={
        'P__initVal':1,
        'Z__initVal':1
    },
    output_vars={
        'env__component': None,
        'env__stateout': 'time'
    }
)

mod_in_x15 = xs.create_setup(
        model=modmod,
    clocks={
         'time': np.linspace(1,365*10,10*365*15)  # 20,1000)
     },
    input_vars={
        'P__initVal':1,
        'Z__initVal':1
    },
    output_vars={
        'env__component': None,
        'env__stateout': 'time'
    }
)

mod_in_x30 = xs.create_setup(
        model=modmod,
    clocks={
         'time': np.linspace(1,365*10,10*365*30)  # 20,1000)
     },
    input_vars={
        'P__initVal':1,
        'Z__initVal':1
    },
    output_vars={
        'env__component': None,
        'env__stateout': 'time'
    }
)

mod_in_x70 = xs.create_setup(
        model=modmod,
    clocks={
         'time': np.linspace(1,365*10,10*365*70)  # 20,1000)
     },
    input_vars={
        'P__initVal':1,
        'Z__initVal':1
    },
    output_vars={
        'env__component': None,
        'env__stateout': 'time'
    }
)


print('hehe')
mod_out_x2 = mod_in_x2.xsimlab.run(model=modmod)
mod_out_x2.to_netcdf('NEW_ASTroCAT_10yearoutput_x2res.nc')

print('huhe')
mod_out_x5 = mod_in_x5.xsimlab.run(model=modmod)
mod_out_x5.to_netcdf('NEW_ASTroCAT_10yearoutput_x5res.nc')


print('huhu')
mod_out_x15 = mod_in_x15.xsimlab.run(model=modmod)
mod_out_x15.to_netcdf('NEW_ASTroCAT_10yearoutput_x15res.nc')


print('huha')
mod_out_x30 = mod_in_x30.xsimlab.run(model=modmod)
mod_out_x30.to_netcdf('NEW_ASTroCAT_10yearoutput_x30res.nc')


print('haha')

mod_out_x70 = mod_in_x70.xsimlab.run(model=modmod)
mod_out_x70.to_netcdf('NEW_ASTroCAT_10yearoutput_x70res.nc')


print('hahahahahaha')