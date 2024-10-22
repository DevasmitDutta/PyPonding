import sys
sys.path.append('/home/devasmit/Desktop/OpenSees/SRC/interpreter')
#import opensees as ops

import matplotlib.pyplot as plt
import openseespy.opensees as ops
import matplotlib.pyplot as plt

from math import pi,cos,cosh,ceil,log
from scipy.stats import norm
import pandas as pd
import scipy.stats
from scipy.stats import gumbel_r
import seaborn as sns
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import gumbel_r

sys.path.append('/../../../PyPonding/')
sys.path.append('/../../../PyPonding/PyPonding/')
sys.path.append('/../../../PyPonding/PyPonding/structures/')
from PyPonding.structures import wide_flange
from PyPonding.structures.wide_flange import wf
from libdenavit.section.database import wide_flange_database


from PyPonding.PondingLoadCell import PondingLoadCell2d
import numpy as np

class PondingLoadCell2d_OPS(PondingLoadCell2d):
    def __init__(self,id,nodeI,nodeJ,gamma,tw):
        self.id = id
        self.nodeI = nodeI
        self.nodeJ = nodeJ
        self.gamma = gamma
        self.tw = tw
        
        # Retreive Node Coordinates
        self.xI = ops.nodeCoord(self.nodeI,1)
        self.yI = ops.nodeCoord(self.nodeI,2)
        self.xJ = ops.nodeCoord(self.nodeJ,1)
        self.yJ = ops.nodeCoord(self.nodeJ,2)
        
    def update(self):
        # Code currently only updates y postion of nodes - @todo maybe update x position also
        # self.dxI = ops.nodeDisp(self.nodeI,1)
        self.dyI = ops.nodeDisp(self.nodeI,2)
        # self.dxJ = ops.nodeDisp(self.nodeJ,1)
        self.dyJ = ops.nodeDisp(self.nodeJ,2)

upper_limit = 0    
import numpy as np

def run_analysis(shape_name,L,slope,qD,label, IM, duration, window_size):

    Pf = fragility_assessment_copy(shape_name,L,slope,qD, label, IM, duration, window_size)
    return Pf

def fragility_assessment_copy(shape_name,L,slope,qD, label, IM, duration, window_size):

        # For creating a certain figure
        np.random.seed(15)
        
        inch = 1.0
        kip = 1.0
        minute = 1.0
        
        lb = kip/1000.0
        ft = 12.0*inch
        hr = 60*minute
        sec = minute/60.0

        ksi = kip/inch**2
        psf = lb/ft**2
        pcf = psf/ft
        
        g = 386.4*inch/sec**2
        gamma   = 62.4*pcf

        gal = 0.133681*ft**3
        cnt_fail = 0


        material_type = 'Elastic'
        material_type = 'ElasticPP'
        material_type = 'Hardening'
        
        Fy = 50.0*ksi # yield stress
        E = 29000.0*ksi # elastic modulus
        Hk = 29.0*ksi # kinematic hardening modulus
        Fr = 0.2*Fy # residual stress
        Fr = 0.001*Fy

        # Hourly rainfall rate (in/hr)
        locations = ['Denver','New York','New Orleans']
        rate = {}
        perc95 = {}
        df = pd.read_csv(os.path.join(os.getcwd(),'studies/reliability/LAX/final_precipitation_intensity_data.xlsx - Sheet1.csv'))
        intensity = df[duration][~np.isnan(df[duration])].values
        rate['Denver'] = np.mean(intensity)
        perc95['Denver'] = np.quantile(intensity,0.95)
 
        u_4 = norm.ppf(1 - 0.5)
        rate['New York'] = np.mean(intensity)
        perc95['New York'] = np.quantile(intensity,0.95)
        u_5 = norm.ppf(1 - 0.5)
        rate['New Orleans'] = np.mean(intensity)
        perc95['New Orleans'] = np.quantile(intensity,0.95)
        u_6 = norm.ppf(1 - 0.5)
        
        # "Drag coefficient" in scupper flow equation
        cd = 0.6
        
        # Scupper width
        ws = 6*inch
        ws = 12*inch

        # Scupper spacing
        Ss = 40*ft
        
        # Tributary width of beams
        tw = 10*ft

        # Number of analysis steps
        nsteps = 100
        nsteps = 250
        nsteps = 3000
        #nsteps = 1000
        
        # Number of Monte Carlo trials
        Ntrials = 200
        #Ntrials = 200
        # Ntrials = 99
        # Ntrials = 1
        # Ntrials = 500
        
        # From function input
        wfs = wide_flange_database[shape_name]
        d = wfs['d']*inch
        
        tweb = wfs['tw']*inch
        bf = wfs['bf']*inch
        tf = wfs['tf']*inch
        dw = d-2*tf
        
        wf_section = wf(d,tweb,bf,tf,Fy,E,Hk)
        wf_section.material_type = material_type

        A  = wf_section.A()
        Iz = wf_section.Iz()

        # From function input
        zi      = 0.0*inch
        zj      = slope*L

        # From function input
        wD = qD*tw # Dead load per length



        Atrib = L*tw
        wf_section.gamma    = gamma
        wf_section.TW    = tw
        wf_section.wD    = wD
        wf_section.L = L
        wf_section.zi = zi
        wf_section.zj = zj
        wf_section.frc = Fr

        # Nominal value for static head
        As = Ss*L
        dhnom = {}
        for city in locations:
            q = rate[city]*As
            # print('rate[city]',rate[city])
            dhnom[city] = (1.5*q/(cd*ws*(2*g)**0.5))**(2.0/3)
            # print('dhnom[city]',dhnom[city])
            # break
        
        methods = ['AISC Appendix 2','DAMP','Proposed for ASCE 7','Neglect Ponding']
        ds = {}
        zw_lim = {}
        # Set design method and compute zw_lim (ds+dh)
        for method in methods:
            zw_lim[method] = wf_section.maximum_permitted_zw(method)
            #ds[method] = zw_lim[method] - dh # Now computed in ReadMonteCarlo.py

        max_volume = (300*inch)*Atrib

        nsteps_vol = 30
        nele = 20
        vol_tol = max_volume/nsteps/100.0
        mid_node = int(nele/2)



                
        # set modelbuilder
        ops.wipe()
        ops.model('basic', '-ndm', 2, '-ndf', 3)
                
        # create nodes
        for i in range(nele+1):
            ops.node(i,L*i/(nele),zi+(zj-zi)*i/(nele))

        # set boundary condition
        ops.fix(   0, 1, 1, 0)
        ops.fix(nele, 0, 1, 0)

        # define coordinate transformation
        ops.geomTransf('Linear',1)

        # define cross section
        wf_section.define_fiber_section(1,1)
        Np = 4
        ops.beamIntegration('Lobatto', 1, 1, Np)
        #ops.beamIntegration('Legendre', 1, 1, 2)    

        EIplot = np.zeros((nele*Np,nsteps+1))

        # Time series for loads
        ops.timeSeries("Constant", 1)

        # Dead load
        ops.pattern('Plain',-1,1)

        
        ops.randomVariable(1,'lognormal','-mean',Fy,'-stdv',0.1*Fy); ops.parameter(1)
        ops.randomVariable(2,'lognormal','-mean', E,'-stdv',0.02*E); ops.parameter(2)
        ops.randomVariable(3,'normal','-mean',-wD,'-stdv',0.1*wD); ops.parameter(3)
        ops.randomVariable(4,'lognormal','-mean',Fr,'-stdv',0.15*Fr); ops.parameter(4)

        ops.parameter(8) # For the web
        Nregions = wf_section.num_regions # Number of fibers in each flange
        for i in range(Nregions):
            ops.parameter(9+i)
        
        ops.probabilityTransformation('Nataf')

        scale = {}
        loc = {}
        rateRVTag = 5
        for city in locations:
            #ops.randomVariable(rateRVTag,'lognormal','-mean',rate,'-stdv',0.2*rate)
            eulergamma = 0.57721566490153286061
            scale[city] = (rate[city]-perc95[city])/(log(-log(0.95))+eulergamma)
            loc[city] = rate[city] - scale[city]*eulergamma
            ops.randomVariable(rateRVTag,'type1LargestValue','-mean',rate[city],'-stdv',(scale[city]*6**0.5)/pi)

            ops.parameter(rateRVTag)
            ops.updateParameter(rateRVTag,rate[city])
            rateRVTag += 1
        
        
        # define elements
        for i in range(nele):

            ops.element("dispBeamColumn",i,i,i+1,1,1)
            ops.addToParameter(1,'element',i,'fy')
            ops.addToParameter(2,'element',i,'E')
            ops.eleLoad('-ele',i,'-type','beamUniform',-wD)
            ops.addToParameter(3,'loadPattern',-1,'elementLoad',i,'wy')
            ops.addToParameter(6,'element',i,'material',1,'F0') # Web
            for j in range(Nregions):
                ops.addToParameter(7+j,'element',i,'material',1+2*(j+1),'F0') # Flange fibers            

        legendLabel = {
            #9: 'gamma',
            9: 'rate',
            10: 'cd',
            1: 'Fy',
            2: 'E',
            3: 'd',
            4: 'tw',
            5: 'bf',
            6: 'tf',
            #8: 'Hkin',
            8: 'Fr',
            7: 'wD'
        }

        # ------------------------------
        # Start of analysis generation
        # ------------------------------

        # create SOE
        ops.system("SparseGeneral")
        
        # create DOF number
        ops.numberer("RCM")
        
        # create constraint handler
        ops.constraints("Plain")
        
        # create integrator
        ops.integrator("LoadControl", 0.0)
        
        ops.test('NormUnbalance',1.0e-6,20,0)
        
        # create algorithm
        ops.algorithm("Newton")
        
        ops.probabilityTransformation('Nataf')
        
        # create analysis object
        ops.analysis("Static")
        

        Nparam = len(ops.getParamTags())
        Nrv = len(ops.getRVTags())
        u = np.zeros(Nrv)

        duPlot = np.zeros((nsteps+1,Nparam))
        meanPlot = np.zeros((nsteps+1,2))

        # import matplotlib.pyplot as plt
        from matplotlib import rc
        # rc('text',usetex=True)
        # rc('font',family='serif')

        plt.figure()
        plt.subplot(2,1,1)
        # u_4 = intensity.quantile(0.95)      

        for j in range(Ntrials+1):
            print(j,label)
            
            ops.reset()

            # Transform random variables from standard normal to actual space
            # 1. Create random values in standard normal space
            jj = 0
            for rv in ops.getRVTags():
                if jj!=4 and jj!=5 and jj!=6:
                    u[jj] = norm.ppf(np.random.rand())   
                elif jj == 4:
                    u[jj] = u_4
                elif jj == 5:
                    u[jj] = u_5
                elif jj == 6:
                    u[jj] = u_6;          
                # if j == Ntrials: 
                #     u[jj] = 0 # mean realizations
                jj = jj+1

            # 2. Transform to real space
            x = ops.transformUtoX(*u)
            x[4] = x[5] = x[6] = IM
            # print('x[4] here %d' % x[4])
            # print('x for trial %d %s %s %s' % (j, x, u, [u_4]))

            # 3. Update parameters with random realizations
            jj = 0
            for rv in ops.getRVTags():
                if rv != 4:
                    ops.updateParameter(rv,x[jj])
                else:
                    Frj = x[jj]
                    Frt = -Frj*(bf*tf)/(bf*tf+tw*dw)
                    ops.updateParameter(8,Frt)
                    for j in range(Nregions):
                        Fri = Frj + (j+0.5)/Nregions*(Frt-Frj)
                        ops.updateParameter(9+j,Fri)
                jj = jj+1

            # Dead load analysis
            ops.analyze(1)

            # define ponding load cells
            # Inside the MC loop so we can make gamma random
            #
            #gamma = x[gammaRVTag-1]
            PondingLoadCells = dict()
            for i in range(nele):
                PondingLoadCells[i] = PondingLoadCell2d_OPS(id,i,i+1,gamma,tw)
        
            # ------------------------------
            # Finally perform the analysis
            # ------------------------------

            data_volume = np.zeros(nsteps+1)
            data_height = np.zeros(nsteps+1)
            end_step = nsteps
                    
            # Create dict of each node that can have ponding load applied and initilize load to zero
            EmptyPondingLoad = dict()
            for iCell in PondingLoadCells:
                if not PondingLoadCells[iCell].nodeI in EmptyPondingLoad:
                    EmptyPondingLoad[PondingLoadCells[iCell].nodeI] = 0.0    
                if not PondingLoadCells[iCell].nodeJ in EmptyPondingLoad:
                    EmptyPondingLoad[PondingLoadCells[iCell].nodeJ] = 0.0
            

            if j == Ntrials:
                meanHV = open(f'meanHV{label}.csv','w')
                #ops.sensitivityAlgorithm('-computeAtEachStep')
                
            # Perform analysis, ramping up volume      
            zw = 0.1
            CurrentPondingLoad = EmptyPondingLoad.copy()
            for iStep in range(nsteps):

                target_volume = (iStep+1)/nsteps*max_volume

                # Update ponding load cells
                for iCell in PondingLoadCells:
                    PondingLoadCells[iCell].update()

                # Estimate water height
                for i in range(nsteps_vol):
                    V = 0
                    dVdz = 0
                    for iCell in PondingLoadCells:
                        (iV,idVdz) = PondingLoadCells[iCell].get_volume(zw)
                        V += iV
                        dVdz += idVdz
                    zw = zw - (V-target_volume)/dVdz
                    if abs(target_volume-V) <= vol_tol:
                        break 

                # Compute load vector
                UpdatedPondingLoad = EmptyPondingLoad.copy()
                for iCell in PondingLoadCells:    
                    f = PondingLoadCells[iCell].get_load_vector(zw)
                    UpdatedPondingLoad[PondingLoadCells[iCell].nodeI] += f.item(0)
                    UpdatedPondingLoad[PondingLoadCells[iCell].nodeJ] += f.item(1)

                # Apply difference to model
                #
                ops.pattern("Plain", iStep, 1)
                for iNode in UpdatedPondingLoad:
                    fy = UpdatedPondingLoad[iNode] - CurrentPondingLoad[iNode]
                    ops.load(iNode, 0.0, fy, 0.0)
                CurrentPondingLoad = UpdatedPondingLoad

                # Run analysis
                ok = ops.analyze(1)
                if ok < 0:
                    print(target_volume,max_volume,zw)
                            
                # Store Data
                data_volume[iStep+1] = target_volume
                data_height[iStep+1] = zw
                if j == Ntrials:
                    meanPlot[iStep+1,0] = target_volume
                    meanPlot[iStep+1,1] = zw
                    meanHV.write(f'{target_volume},{zw}\n')
                    for rv in ops.getRVTags():
                        cov = ops.getStdv(rv)/ops.getMean(rv)
                        # Negative sign because disp is downward
                        #duPlot[iStep+1,rv-1] = -ops.sensNodeDisp(mid_node,2,rv)*abs(ops.getStdv(rv))


                # Stop analysis if water level too low or analysis failed
                if zw <= -4*inch or ok < 0:
                    end_step = iStep+1
                    break        

            lw = 0.1
            ls = 'grey'
            label = 'Safe Trial'
            method = 'DAMP'
            city = 'Denver'
            ds = zw_lim[method] - dhnom[city]
            print('x[4]',x[4])
            q = x[5-1]*As
            dh = (1.5*q/(cd*ws*(2*g)**0.5))**(2.0/3)
            print('dh',dh)
            # output.write(f'{max(data_height)},')
            if max(data_height) < ds + dh:
                lw = 1.5
                ls = 'k'
                label = 'Fail Trial'
                cnt_fail += 1
            if label=='Fail Trial':
                plt.plot(data_volume[:end_step]/(Ss*L)*25.4,data_height[:end_step]*25.4,ls,linewidth=lw)
                upper_limit = ds+dh
                # output.write(f'Fail ,')
            else:
                plt.plot(data_volume[:end_step]/(Ss*L)*25.4,data_height[:end_step]*25.4,ls,linewidth=lw)
                upper_limit = ds+dh
                # output.write(f'Safe ,')

            # output.write('\n')

            if j == Ntrials:
                meanHV.close()
                    
            # Remove patterns so there's not
            # duplicate tag errors in MC loop
            for iStep in range(nsteps):
                ops.remove('loadPattern',iStep)



        # output.close()

        # Remove random variables bc ops.wipe() doesn't do this yet
        ops.wipeReliability()

        print('zw_lim',zw_lim);
        print('dhnom',dhnom);
        print('ds',ds);
        print('dh',dh)
        print('ds + dhnom',ds+dhnom['Denver']);
        print('ds + dh',ds+dh);
        
        plt.plot([0,8*25.4],[upper_limit*25.4,upper_limit*25.4],'r--',linewidth=1.1,label='Design Limit')
        pf = cnt_fail/Ntrials
        plt.title('Pf = %f, for i = %s' % (pf, x[4]))
        plt.xlim(left=0,right=175)
        plt.ylim(bottom=0)
        plt.legend(ncol=2,frameon=False)
        plt.ylabel('Water Level (mm)')
        plt.xlabel('Normalized Water Volume, $V/SL$ (mm)')

        plt.savefig('MCtrials.pdf',bbox_inches='tight')
        plt.savefig(f'studies/reliability/LAX/hazard/Risk_assess_Frag_trials_{duration}_{window_size}/{IM}_1976.png', bbox_inches='tight')

        return pf   
        # plt.show()
