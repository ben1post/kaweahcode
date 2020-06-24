# 

#import pandas
import numpy as np
#import scipy.interpolate as intrp

import sys
sys.path.append('../../phydra_OSM/')

from phydra.aux import sliceparams, sliceoffparams, checkreplaceparam
from phydra.forcing import Forcing

from scipy.io import netcdf
import matplotlib

import numpy as np
import xsimlab as xs
from scipy.integrate import odeint
import matplotlib.pyplot as plt

import pickle 


@xs.process
class ExternalForcing:
    """Process supplies forcing"""
    #  -20, 47, 2.
    lat = xs.variable(default=47)
    lon = xs.variable(default=-20)
    rbb = xs.variable(default=2.)
    
    forcing = xs.variable(intent='out')
    
    def initialize(self):
        filehandler = open('forcing_47_minus20.obj', 'rb') 
        object_pickle = pickle.load(filehandler)
        self.forcing = object_pickle
        #Forcing('WOA2018', self.lat, self.lon, self.rbb)
        
@xs.process
class UniqueComponent:
    """This is for initializing a Nutrient (i.e. components that are unique (0D))"""
    label = xs.variable(intent='out', groups='label')
    init = xs.variable(intent='out', groups='init')
    
    def initialize(self):
        self.label = (self.component_label,) # labels are collected as tuples in initialization
        self.init = np.array(self.initVal)   # properties are collected as numpy arrays (might be better ways to do this)
        
        

@xs.process
class Nutrient(UniqueComponent):
    """There is only one Nutrient component in the current model setup, that is consumed by all Phytoplankton and get's replenished with a certain supply rate"""
    component_label = xs.variable(default='N', groups='component-label')
    
    initVal = xs.variable(default=1, description='initial Nutrient concentration')
    
    supplyrate = xs.variable(default=1, description='Nutrient supply rate', groups='parameter', static=True)

    # Process variables
    supply = xs.variable(intent='out')
    
    def _supply(self):
        return self.supplyrate
    
    def initialize(self):
        super(Nutrient, self).initialize() # for label & init
        self.supply = self._supply
        
        
@xs.process
class Detritus(UniqueComponent):
    """There is only one Detritus component in the current model setup, phytoplankton mortality and zooplankton excretion ends up in this pool that sinks and gets remineralized"""
    component_label = xs.variable(default='D', groups='component-label')
    
    initVal = xs.variable(default=1, description='initial Detritus concentration')
    
    remin_rate = xs.variable(default=.1, description='Detritus remineralization rate', groups='parameter', static=True)

    # Process variables
    remineralize = xs.variable(intent='out')
    
    def _remineralize(self, D):
        return self.remin_rate * D
    
    def initialize(self):
        super(Detritus, self).initialize() # for label & init
        self.remineralize = self._remineralize
        
        
@xs.process
class FlexSizeComponent:
    """This is for initializing both Phytoplankton and Zooplankton (i.e. components of variable number (1D))"""
    label = xs.variable(intent='out', groups='label')
    init = xs.variable(intent='out', groups='init')
    size = xs.variable(dims='flexComponent',intent='out', groups='size')
        
    def calculate_sizes(self):
        # must be implemented in subclasses
        raise NotImplementedError
        
    def initialize_alloparams(self):
        # must be implemented in subclasses
        raise NotImplementedError
        
    def initialize(self):
        # Note: i need some way to distinguish every individual component, even within FlexSizeComponents, so I am doing this:
        self.label = tuple((self.component_label + str(i) for i in range(self.NP))) 
        
        self.init = np.array([self.initVal / self.NP for i in range(self.NP)])
        
        self.size = self.calculate_sizes()

        self.initialize_alloparams()
        
@xs.process
class Phytoplankton(FlexSizeComponent):
    """The current model setup contains 40 Phytoplankton, ordered by Size (ESD = equivalent spherical diameter)"""
    component_label = xs.variable(default='P', groups='component-label')
    # Phytoplankton params
    PMinEsd = xs.variable(default = 1,  description='P size range minimum ESD', groups='parameter', static=True)
    PMaxEsd = xs.variable(default = 20, description='P size range maximum ESD', groups='parameter', static=True)
    NP =      xs.variable(default = 40, description='number of P - and Z',  groups='parameter', static=True)
    
    # INPUT
    initVal = xs.variable(description='initial total Phytoplankton concentration')
    m = xs.variable(default=.1, description='Phytoplankton mortality', groups='parameter', static=True)
    
    # Allometric params
    ks = xs.variable(intent='out', description='allometric half-saturation constant',  groups='alloparameter')
    mu0 = xs.variable(intent='out', description='allometric growth rate',  groups='alloparameter')
    
    # non-allometric params
    kw = xs.variable(default=0.04) # Light attenuation constant of water (m^-1)
    kc = xs.variable(default=0.03) # Light attenuation via phytoplankton pigment (m^-1)
    OptI = xs.variable(default=10)
    
    # Process variables
    uptake = xs.variable(intent='out')
    mortality = xs.variable(intent='out')
    grazed = xs.variable(intent='out')
    lightharvesting = xs.variable(intent='out')
    tempdepgrowth = xs.variable(intent='out')
    
    def _uptake(self, N):
        return self.mu0 * N / (N + self.ks)
    
    def _mortality(self):
        """returns size-scaled mortality rate, calculated from mu0"""
        return self.m * self.mu0 
    
    def _grazed(self,FgrazP):
        """returns grazed biomass (by all Z) per P"""
        return [sum(FgrazP[i,:]) for i in range(self.NP)]
    
    def _tempdepgrowth(self, Tmld):
        """temperature dependenc of growth rate"""
        return np.exp(0.063 * Tmld)
        
    def _lightharvesting(self, MLD, PAR, Psum):
        """Smith-type light modification of phytoplankton growth rate"""
        kPAR = self.kw + self.kc * Psum
        lighthrv = 1. / (kPAR * MLD) * \
                       (-np.exp(1. - PAR / self.OptI) - (
                           -np.exp((1. - (PAR * np.exp(-kPAR * MLD)) / self.OptI))))
        return lighthrv
    
    def calculate_sizes(self):
        """initializes array of sizes from ESD size range and number of P"""
        numbers = np.array([i for i in range(self.NP)])
        sizes = (np.log(self.PMaxEsd) - np.log(self.PMinEsd))* numbers / (self.NP-1) + np.log(self.PMinEsd)
        #npp = 3 #per sizeclass
        #npp_num = np.array([i for i in range(npp)])
        
        #pico_sizes = (np.log(1.999) - np.log(0.5))* npp_num / (npp-1) + np.log(0.5)
        
        #nano_sizes = (np.log(19.999) - np.log(2))* npp_num / (npp-1) + np.log(2)
        
        #micro_sizes = (np.log(49.999) - np.log(20))* npp_num / (npp-1) + np.log(0.5)
        #sizes = (self.PMaxEsd - self.PMinEsd)* numbers / (self.NP-1) + self.PMinEsd
        
        #all_sizes = np.concatenate([pico_sizes, nano_sizes, micro_sizes], axis=None)
        print(np.exp(sizes))
        return np.exp(sizes)
    
    def initialize_alloparams(self):
        """initializes allometric parameters based on array of sizes (ESD)
        allometric relationships are taken from meta-analyses of lab data"""
        self.mu0 = 2.6 * (self.size) ** -0.45
        self.ks = (self.size) * .1
        
    def initialize(self):
        """Phytoplankton initialization calls FlexSizeComponent initialize()
        and assigns process functions to be used in ODE"""
        super(Phytoplankton, self).initialize()
        self.uptake = self._uptake
        self.mortality = self._mortality
        self.grazed = self._grazed
        self.lightharvesting = self._lightharvesting
        self.tempdepgrowth = self._tempdepgrowth
        
@xs.process
class Zooplankton(FlexSizeComponent):
    """The current model setup contains 40 Zooplankton, ordered by Size (ESD = equivalent spherical diameter)"""
    component_label = xs.variable(default='Z', groups='component-label')
    # SetupParams
    PMinEsd = xs.foreign(Phytoplankton, 'PMinEsd')
    PMaxEsd = xs.foreign(Phytoplankton, 'PMaxEsd')
    NP =      xs.foreign(Phytoplankton, 'NP')
    phytosize = xs.foreign(Phytoplankton, 'size')
   
    # Zooplankton params 
    initVal =    xs.variable(description='initial total Zooplankton concentration')
    zeta =       xs.variable(default=.1,  description='Z quadratic mortality', groups='parameter', static=True)
    deltaxprey = xs.variable(default=0.25,  description='log10 prey size tolerance', groups='parameter', static=True)
    KsZ =        xs.variable(default=3,  description='grazing half saturation constant', groups='parameter', static=True) #3
    f_eg =       xs.variable(default=.33,  description='egested grazed biomass', groups='parameter', static=True)
    epsilon =    xs.variable(default=.33,  description='assimilated grazed biomass', groups='parameter', static=True)
    
    
    # Allometric params
    I0 = xs.variable(intent='out', description='Z allometric uptake rate',  groups='alloparameter')
    xpreyopt = xs.variable(intent='out', description='Z allometric optimal prey size',  groups='alloparameter')
    
    phiP = xs.variable(intent='out', description='grazing preference matrix')
    
    # Process variables
    grazingmatrix = xs.variable(intent='out')
    mortality = xs.variable(intent='out')
    ingestion = xs.variable(intent='out')
    excretion = xs.variable(intent='out')
    
    def _grazingmatrix(self,P,Z):
        """calculates biomass grazed as a matrix of P and Z"""
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
        """calculates biomass assimilated of P per Z"""
        return [self.epsilon * sum(FgrazP[:,j]) for j in range(self.NP)]

    def _excretion(self,FgrazP):
        """calculates biomass excreted back into N pool of P per Z"""
        return [(1 - self.f_eg - self.epsilon) * sum(FgrazP[:,j]) for j in range(self.NP)]
    
    def _mortality(self,Z):
        """returns quadratic mortality (multiplied by sum of all Z)"""
        return self.zeta * sum(Z)
        
    def calculate_sizes(self):
        """initializes array of sizes from ESD size range of Z based on size range of P
        so that each P has an optimized grazer"""
        zoosizes= 2.16 * self.phytosize ** 1.79
        return zoosizes
    
    def initialize_alloparams(self):
        """initializes allometric parameters based on array of sizes (ESD)"""
        self.I0 = 26 * (self.size) ** -0.4 #* .5
        self.xpreyopt = 0.65 * (self.size) ** .56 # which should equal = self.phytosize
    
    def init_phiP(self):
        """creates array of feeding preferences [P...P10] for each [Z] based on ESD and deltaxprey"""
        phiP= np.array([[np.exp(-((np.log10(xpreyi)-np.log10(xpreyoptj)) / self.deltaxprey)**2) 
               for xpreyi in self.phytosize] for xpreyoptj in self.xpreyopt])
        return phiP
    
    def initialize(self):         
        """Zooplankton initialization calls FlexSizeComponent initialize()
        and assigns process functions to be used in ODE"""       
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
class Slab:
    # added variables
    kappa = xs.variable(default=0.1)
    
    # foreign processes
    physx = xs.foreign(ExternalForcing, 'forcing') # new FORCING!
    
    N_supply = xs.foreign(Nutrient, 'supply')
    
    P_uptake = xs.foreign(Phytoplankton, 'uptake')
    P_mortality = xs.foreign(Phytoplankton, 'mortality')
    P_lightharvesting = xs.foreign(Phytoplankton, 'lightharvesting')
    P_tempdepgrowth = xs.foreign(Phytoplankton, 'tempdepgrowth')
    
    grazingmatrix = xs.foreign(Zooplankton, 'grazingmatrix')
    Z_ingestion = xs.foreign(Zooplankton, 'ingestion')
    Z_excretion = xs.foreign(Zooplankton, 'excretion')
    Z_mortality = xs.foreign(Zooplankton, 'mortality')
    P_grazed = xs.foreign(Phytoplankton, 'grazed')
    
    D_remin = xs.foreign(Detritus, 'remineralize')
    
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
        """model state is initialized"""
        self.time = 0

        # This is one of those things I am not sure about, it makes sense to use dicts to store model state (to allow calling different components by labels)
            # BUT this can not be used (afaik) in a nice way as output for plotting later on in the model output xarray.. do you have any idea how to do that more nicely?
        self.state = {label:val for label,val in zip(self.complabels,self.inits)}

        # in order to have output to plot, I also create an array that contains all output:
            # i have to use flatten function (defined above) to make group generators one dimensional (FlexComponents are Lists, so generators return List of Lists)
            # there should be a better solution to this.. but it relates to how xs.group() works in this model setup (i wrote more on this in "open questions" below)
        self.stateout = np.concatenate([self.state[label] for label in flatten(self.complabels)], axis=None)
        
        self.component = [i for i in flatten(self.labels)]
    
    def ode(self, state, time, labels, timestep):
        #slab addons:
        MLD = self.physx.MLD.return_interpvalattime(time)
        MLDderiv = self.physx.MLD.return_derivattime(time)
        N0 = self.physx.NOX.return_interpvalattime(time)
        PAR = self.physx.PAR.return_interpvalattime(time)
        Tmld = self.physx.SST.return_interpvalattime(time)
        
        MIX = (self.kappa + max(MLDderiv,0)) / MLD
        
        """the ode function defines the fluxes, based on functions passed from components
        state is split into the component parts and merged again, after calculating fluxes"""
        N, P, Z, D = state['N'], state['P'], state['Z'], state['D']
        
        # this was in Neil Bana's original code, however does not seem necessary here:
        #P[P < 0.01] = 0.01
        #Z[Z < 0.01] = 0.01
        
        PGrazed = self.grazingmatrix(P,Z)
        PGains = self.P_uptake(N)*self.P_lightharvesting(MLD,PAR,sum(P))*self.P_tempdepgrowth(Tmld)*P
        PMortality = self.P_mortality()*P
        
        ZIngested = self.Z_ingestion(PGrazed)
        
        DRemin = self.D_remin(D)
        
        dNdt = (N0-N) * MIX + sum(self.Z_excretion(PGrazed)) + DRemin - sum(PGains)
        
        # print('PAR', PAR) # 6 - 45 roughly
        dPdt = PGains - PMortality - self.P_grazed(PGrazed) - P * MIX
        
        dZdt =  ZIngested - self.Z_mortality(Z)*Z  - Z * MIX
        
        dDdt = sum(PMortality) + sum(ZIngested) - DRemin - D * MIX #* 5 # sinkrate v is 3, 0.33 digestion rate Z
        
        # need to convert fluxes to correct time step (rates per day) by multiplying by dt        
        return {'N':dNdt*timestep, 'P':dPdt*timestep, 'Z':dZdt*timestep, 'D':dDdt*timestep}

    @xs.runtime(args='step_delta')
    def run_step(self,dt):
        """fluxes are calculated per timestep"""
        self.outflux = self.ode(self.state,self.time,self.labels, dt)
        # keep track of timestep for forcing (+ other time dep. processes)
        self.time += dt
    
    def finalize_step(self):
        """new model state is calculated from fluxes"""
        self.state = {label : self.state[label] + self.outflux[label] for label in self.complabels}
        
        # state for computation, stateout for plotting later (this seems stupid)
        self.stateout = np.concatenate(flatten([(self.state[label] for label in self.complabels)]), axis=None)

ASTroCATclone = xs.Model({'physx':ExternalForcing, 'env':Slab, 'N':Nutrient, 'P':Phytoplankton, 'Z':Zooplankton, 'D':Detritus})










def new_plot(modeloutput, forcingverif, Pnum):
   
    data_coordfix = modeloutput.set_index(component = 'env__component', P__sizes='P__size', Z__sizes='Z__size')

    P_out = data_coordfix.sel(component=['P'+str(i) for i in range(Pnum)])
    Z_out = data_coordfix.sel(component=['Z'+str(i) for i in range(Pnum)])
    N_out = data_coordfix.sel(component='N')
    D_out = data_coordfix.sel(component='D')

    resolution = 10
    resmpl_step = 1


    P_pico = P_out.env__stateout[:,P_out.P__sizes.values<2]
    mask = ((P_out.P__sizes.values >= 2) & (P_out.P__sizes.values <= 20))
    P_nano = P_out.env__stateout[:,mask]
    mask2 = ((P_out.P__sizes.values >= 20) & (P_out.P__sizes.values <= 50))
    P_micro = P_out.env__stateout[:,mask2]
    ##########
    
    
    P_pico_sum = P_pico.sum(dim='component')
    P_nano_sum = P_nano.sum(dim='component')
    P_micro_sum = P_micro.sum(dim='component')
    
    Z_out_sum = Z_out.sum(dim='component')
    
    ##########
    P_pico_sum_ly = P_pico_sum[P_pico_sum.time > 365*5]
    P_nano_sum_ly = P_nano_sum[P_nano_sum.time > 365*5]
    P_micro_sum_ly = P_micro_sum[P_micro_sum.time > 365*5]

    N_out_ly = N_out.env__stateout[N_out.time > 365*5]
    D_out_ly = D_out.env__stateout[D_out.time > 365*5]
    Z_out_ly = Z_out_sum.env__stateout[Z_out.time > 365*5]

    P_pico_sum_ly_yday = P_pico_sum_ly.assign_coords(time=np.mod(P_pico_sum_ly.time, 365.))
    P_nano_sum_ly_yday = P_nano_sum_ly.assign_coords(time=np.mod(P_nano_sum_ly.time, 365.))
    P_micro_sum_ly_yday = P_micro_sum_ly.assign_coords(time=np.mod(P_micro_sum_ly.time, 365.))

    N_out_ly_yday = N_out_ly.assign_coords(time=np.mod(N_out_ly.time, 365.))
    D_out_ly_yday = D_out_ly.assign_coords(time=np.mod(D_out_ly.time, 365.))
    Z_out_ly_yday = Z_out_ly.assign_coords(time=np.mod(Z_out_ly.time, 365.))
    #################

    import datetime

    seconds = P_pico_sum_ly_yday.time.values * 60 * 60 * 24

    dt_array = [datetime.datetime.utcfromtimestamp(i) for i in seconds]
    #################

    P_pico_sum_ly_dt = P_pico_sum_ly_yday.assign_coords(time=dt_array).to_dataset()
    P_nano_sum_ly_dt = P_nano_sum_ly_yday.assign_coords(time=dt_array).to_dataset()
    P_micro_sum_ly_dt = P_micro_sum_ly_yday.assign_coords(time=dt_array).to_dataset()

    N_out_ly_dt = N_out_ly_yday.assign_coords(time=dt_array).to_dataset()
    D_out_ly_dt = D_out_ly_yday.assign_coords(time=dt_array).to_dataset()
    Z_out_ly_dt = Z_out_ly_yday.assign_coords(time=dt_array).to_dataset()

    P_pico_yagg = P_pico_sum_ly_dt.groupby('time.month').mean()
    P_pico_yagg_sd = P_pico_sum_ly_dt.groupby('time.month').std()

    P_nano_yagg = P_nano_sum_ly_dt.groupby('time.month').mean()
    P_nano_yagg_sd = P_nano_sum_ly_dt.groupby('time.month').std()

    P_micro_yagg = P_micro_sum_ly_dt.groupby('time.month').mean()
    P_micro_yagg_sd = P_micro_sum_ly_dt.groupby('time.month').std()

    N_out_yagg = N_out_ly_dt.groupby('time.month').mean()
    N_out_yagg_sd = N_out_ly_dt.groupby('time.month').std()

    D_out_yagg = D_out_ly_dt.groupby('time.month').mean()
    D_out_yagg_sd = D_out_ly_dt.groupby('time.month').std()
    
    Z_out_yagg = Z_out_ly_dt.groupby('time.month').mean()
    Z_out_yagg_sd = Z_out_ly_dt.groupby('time.month').std()

    # conversion
    CtoN = 106/16
    molarmassC = 12.0107
    CBMconvert = molarmassC * CtoN
    
    ##########
    dayspermonth = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    dpm_cumsum = np.cumsum(dayspermonth) - np.array(dayspermonth)/2 #- 15
    ##########
    
    
    plt.rcParams['figure.figsize'] = [8, 5]

    #c to N # 106/16
    labels = {'model_pico':'darkgreen','model_nano':'green','model_micro':'lightgreen', 'model_nuts':'blue','model_det':'brown','model_Z':'red','data':'grey'}
    markers = {'pico':'x', 'nano':'^', 'micro':'s', 'nuts':'s', 'det':'s'}

    #fig, ax = plt.subplots(1, 3, sharey='row')
    fig = plt.figure(constrained_layout=True)


    gs = fig.add_gridspec(nrows=2, ncols=3, hspace=0.1, wspace=0.1)#, height_ratios=[.1,1,1,.1,1])

    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[0, 1])#, sharey=ax0)
    ax2 = fig.add_subplot(gs[0, 2])#, sharey=ax0)
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])
    ax5 = fig.add_subplot(gs[1, 2])

    # Phyto MODEL
    ax0.errorbar(dpm_cumsum, P_pico_yagg.env__stateout, P_pico_yagg_sd.env__stateout, linestyle='None', c=labels['model_pico'])#, marker='^')
    ax0.scatter(dpm_cumsum, P_pico_yagg.env__stateout, label='model', c=labels['model_pico'], marker=markers['pico'])

    ax1.errorbar(dpm_cumsum, P_nano_yagg.env__stateout, P_nano_yagg_sd.env__stateout, linestyle='None', c=labels['model_nano'])#, marker='^')
    ax1.scatter(dpm_cumsum, P_nano_yagg.env__stateout, label='model', c=labels['model_nano'], marker=markers['nano'])

    ax2.errorbar(dpm_cumsum, P_micro_yagg.env__stateout, P_micro_yagg_sd.env__stateout, linestyle='None', c=labels['model_micro'])
    ax2.scatter(dpm_cumsum, P_micro_yagg.env__stateout, label='model', c=labels['model_micro'], marker=markers['micro'])

    # Phyto DATA

    ax0.errorbar(dpm_cumsum, np.array(forcingverif.verif.c_picop) /CBMconvert, np.array(forcingverif.verif.c_picop_sd) /CBMconvert, linestyle='None', capsize = 4,capthick = 1.5, linewidth=0.1, c=labels['data'])
    ax0.scatter(dpm_cumsum, np.array(forcingverif.verif.c_picop) /CBMconvert , label='satellite', c=labels['data'])

    ax1.errorbar(dpm_cumsum, np.array(forcingverif.verif.c_nanop) /CBMconvert, np.array(forcingverif.verif.c_nanop_sd) /CBMconvert, linestyle='None', capsize = 4,capthick = 1.5, linewidth=0.1, c=labels['data'])
    ax1.scatter(dpm_cumsum, np.array(forcingverif.verif.c_nanop) /CBMconvert , label='satellite', c=labels['data'])

    ax2.errorbar(dpm_cumsum, np.array(forcingverif.verif.c_microp) /CBMconvert, np.array(forcingverif.verif.c_microp_sd) /CBMconvert, linestyle='None', capsize = 4,capthick = 1.5, linewidth=0.1, c=labels['data'])
    ax2.scatter(dpm_cumsum, np.array(forcingverif.verif.c_microp) /CBMconvert , label='satellite', c=labels['data'])


    # NUTRIENTS
    ax3.errorbar(dpm_cumsum, N_out_yagg.env__stateout, N_out_yagg_sd.env__stateout, linestyle='None', c=labels['model_nuts'])#, marker='^')
    ax3.scatter(dpm_cumsum, N_out_yagg.env__stateout, label='model', c=labels['model_nuts'], marker=markers['nuts'])

    ax3.scatter(dpm_cumsum, np.array(forcingverif.verif.N) , label='WOA data', c=labels['data'])

    # DETITUS
    ax4.errorbar(dpm_cumsum, D_out_yagg.env__stateout, D_out_yagg_sd.env__stateout, linestyle='None', c=labels['model_det'])#, marker='^')
    ax4.scatter(dpm_cumsum, D_out_yagg.env__stateout, label='model', c=labels['model_det'], marker=markers['det'])

    # Zooplankton
    ax5.errorbar(dpm_cumsum, Z_out_yagg.env__stateout, Z_out_yagg_sd.env__stateout, linestyle='None', c=labels['model_Z'])#, marker='^')
    ax5.scatter(dpm_cumsum, Z_out_yagg.env__stateout, label='model', c=labels['model_Z'], marker=markers['det'])



    # LABELS
    ax0.set_title('Pico (0.5 - 2 µm)')
    ax1.set_title('Nano (2 - 20 µm)')
    ax2.set_title('Micro (20 - 50 µm')

    ax0.set_ylabel('Phytoplankton  [µM N]')
    ax1.set_xlabel('Day in year [days]')


    ax3.set_title('N')

    ax3.set_ylabel('Nitrogen [µM N]')
    ax3.set_xlabel('Day in year [days]')


    ax4.set_title('D')

    ax4.set_ylabel('Detritus [µM N]')
    ax4.set_xlabel('Day in year [days]')
    
    ax5.set_title('sum(Z)')

    ax5.set_ylabel('Zooplankton [µM N]')
    ax5.set_xlabel('Day in year [days]')

    #plt.plot(P_pico_yagg.month, P_pico_yagg.env__stateout, label=labels['pico'])
    #plt.plot(P_nano_yagg.month, P_nano_yagg.env__stateout, label=labels['nano'])
    #plt.plot(P_micro_yagg.month, P_micro_yagg.env__stateout, label=labels['micro'])

    #plt.scatter(P_nano_sum_ly_yday.time.values, P_nano_sum_ly_yday.values, label='nano')
    #plt.scatter(P_micro_sum_ly_yday.time.values, P_micro_sum_ly_yday.values, label='micro')
    xlim = (1, 365)
    from matplotlib.ticker import MaxNLocator
    for axe in [ax0,ax1,ax2,ax3, ax4, ax5]:
        axe.legend(fontsize='x-small')
        plt.setp(axe, xticks=[1,60,120,180,240,300,365], xlim=xlim)
        axe.grid(True, alpha=0.5)
        axe.get_yaxis().set_major_locator(MaxNLocator(nbins=4))
        axe.tick_params(top=False, right=True, direction="in")
        axe.set_ylim(bottom=0)

    #plt.setp(ax1.get_yticklabels(), visible=False)
    #plt.setp(ax2.get_yticklabels(), visible=False)

    # Defining custom 'xlim' and 'ylim' values.

    

    