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


inch = 1.0
kip = 1.0
hr = 1.0

minute = hr/60.0
sec = minute/60.0

lb = kip/1000.0
ft = 12.0*inch

ksi = kip/inch**2
psf = lb/ft**2
pcf = psf/ft

g = 386.4*inch/sec**2
gamma   = 62.4*pcf

gal = 0.133681*ft**3


# Hourly rainfall rate (in/hr)
locations = ['Santa Barbara']
rate = {}
perc95 = {}


df_mean = pd.read_csv(os.path.join(os.getcwd(),f'studies/reliability/SB_Design_Specs/mean_intensity_conf_interval.csv'))
df_mean.set_index('Dur (hr)', inplace=True)

df_95_percentile = pd.read_csv(os.path.join(os.getcwd(),f'studies/reliability/SB_Design_Specs/95_percent_intensity_conf_interval.csv'))
df_95_percentile.set_index('Dur (hr)', inplace=True)

return_period = '100yr' # years
rate['Santa Barbara'] = df_mean.loc['x0.25hr',return_period]*inch/hr
print('Mean intensity for Santa Barbara',rate['Santa Barbara'],'in/hr for duration x0.25hr')

# "Drag coefficient" in scupper flow equation
cd = 0.6

# Scupper width
ws = 6*inch
ws = 12*inch

# Scupper spacing
Ss = 40*ft

# Tributary width of beams
tw = 10*ft

shape_name = 'W8X18'
L = 30.0*12.0 # in
slope = 0.25/12.0 # in/in
qD = 10.0*(1/1000)/(12**2) # kip
Fy = 50.0*ksi # psi
E = 29000.0*ksi # psi
Hk = 29.0*ksi # psi
Fr = 0.3*Fy
material_type = 'Hardening'  # 'Elastic' 'Hardening' 'Bilinear
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
    dhnom[city] = (1.5*q/(cd*ws*(2*g)**0.5))**(2.0/3)
    print('dhnom in',dhnom,'i in/hr',rate[city],'As in',As,'Ss in',Ss,'L in',L,'q in^3/hr',As*rate[city],'cd',cd,'ws in',ws,'g in/hr^2',g)

# methods = ['AISC Appendix 2','DAMP','Proposed for ASCE 7','Neglect Ponding']
ds = {}
zw_lim = {}
# Set design method and compute zw_lim (ds+dh)
# for method in methods:
#     zw_lim[method] = wf_section.maximum_permitted_zw(method)
#     #ds[method] = zw_lim[method] - dh # Now computed in ReadMonteCarlo.py


print('\n -------------Showing strength ratio useful for design purposes ----------\n')

wf_section_1 = wf_section
wf_section_1.zj = 0.0*L
zw_lim_no_slope = {}
zw_lim_no_slope['DAMP'] = wf_section_1.maximum_permitted_zw('DAMP')
print('zw_lim for no slope',zw_lim_no_slope,'\n')
design_z = 2 + (1.5*(rate['Santa Barbara']*As)/(cd*ws*(2*g)**0.5))**(2.0/3)#3.14
preliminary_design_step_SR_ratio_for_a_given_z = wf_section_1.preliminary_design_step_SR_ratio_for_a_given_z(design_z)
print('preliminary_step_without_slope_SR_ratio_for_a_given_z\n',preliminary_design_step_SR_ratio_for_a_given_z,'at z=%0.3f in' % design_z)

wf_section_2 = wf_section
wf_section_2.zj = 0.0*L
# zw_lim_2 = {}
# zw_lim_2['DAMP'] = wf_section_2.maximum_permitted_zw('DAMP')
design_z = 2 + (1.5*(rate['Santa Barbara']*As)/(cd*ws*(2*g)**0.5))**(2.0/3) #3.14
ponding_instability_step_without_slope_SR_ratio_for_a_given_z = wf_section_2.SR_ratio_for_a_given_z(design_z)
print('ponding_instability_step_without_slope_SR_ratio_for_a_given_z\n',ponding_instability_step_without_slope_SR_ratio_for_a_given_z,'at z=%0.3f in' % design_z)

wf_section_3 = wf_section
wf_section_3.zj = slope*L
zw_lim_with_slope = {}
zw_lim_with_slope['DAMP'] = wf_section_3.maximum_permitted_zw('DAMP')
print('zw_lim for slope',zw_lim_with_slope,'\n')
design_z = 2 + (1.5*(rate['Santa Barbara']*As)/(cd*ws*(2*g)**0.5))**(2.0/3) #3.14
preliminary_design_step_SR_ratio_for_a_given_z = wf_section_3.preliminary_design_step_SR_ratio_for_a_given_z(design_z)
print('preliminary_step_with_slope_SR_ratio_for_a_given_z\n',preliminary_design_step_SR_ratio_for_a_given_z,'at z=%0.3f in' % design_z)

wf_section_4 = wf_section
wf_section_4.zj = slope*L
# zw_lim_4 = {}
# zw_lim_4['DAMP'] = wf_section_4.maximum_permitted_zw('DAMP')
design_z = 2 + (1.5*(rate['Santa Barbara']*As)/(cd*ws*(2*g)**0.5))**(2.0/3) #3.14
ponding_instability_step_with_slope_SR_ratio_for_a_given_z = wf_section_4.SR_ratio_for_a_given_z(design_z)
print('ponding_instability_step_with_slope_SR_ratio_for_a_given_z',ponding_instability_step_with_slope_SR_ratio_for_a_given_z,'at z=%0.3f in' % design_z)
print('rate',rate['Santa Barbara'])

print('\n -------------End of strength ratio useful for design purposes ----------\n')
