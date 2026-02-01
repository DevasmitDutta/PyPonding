import openseespy.opensees as ops
from math import pi,cos,cosh,ceil,tan,tanh
import numpy as np
import matplotlib.pyplot as plt
from PyPonding.structures import wf#,wf_shapes
from libdenavit.section.database import wide_flange_database
import matplotlib as mpl
mpl.rcParams['text.usetex'] = False

# Define units
inch = 1.0
kip = 1.0
lb = kip/1000.0
ft = 12.0*inch
in_per_ft = inch/ft
ksi = kip/inch**2
psf = lb/ft**2
pcf = psf/ft
mm  = 1/25.4*inch
m   = 1000*mm
kN  = 1/4.448222*kip
kNm = kN*m

# Input parameters
Fy      = 50.0*ksi
E       = 29000.0*ksi
Hk      = 1.0e-3*E
TW      = 10*ft
shape_name = 'W14X22';
L1      = 30*ft
# L2      = 49.57*ft
# L3      = 54.85*ft
slope   = (0.02*12)*in_per_ft
qD      = 10.0*psf
gamma   = 62.4*pcf

# Lookup shape data
# shape_data = wf_shapes[shape_name]
shape_data = wide_flange_database[shape_name]
d  = shape_data['d']*inch
bf = shape_data['bf']*inch
tf = shape_data['tf']*inch
tw = shape_data['tw']*inch

# Create wide-flange object
wf_section = wf(d,tw,bf,tf,Fy,E,Hk)

# Set additional properties
wf_section.L        = L1
wf_section.gamma    = gamma
wf_section.TW       = TW
wf_section.zi       = 0.0*inch
wf_section.zj       = L1*slope
wf_section.wD       = qD*TW # Dead load per length

# Set OpenSees analysis options
wf_section.material_type = 'Elastic'
wf_section.num_steps = 3000
wf_section.nsteps_vol = 30
wf_section.geomTransfType = 'Linear'
wf_section.geomTransfTag = 1
# wf_section.percent_drop = 1.0

# Run OpenSees analyses 
# L1
wf_section.L  = L1
C1 = wf_section.C()
wf_section.max_volume = (300*inch)*L1*TW
wf_section.vol_tol = wf_section.max_volume/wf_section.num_steps/100.0; #wf_section.max_volume/wf_section.num_steps/10000.


wf_section.zj = slope*L1
(data_volume_L1s,data_height_L1s) = wf_section.perform_OpenSees_analysis();

wf_section.material_type = 'Hardening'
(data_volume_L1si,data_height_L1si) = wf_section.perform_OpenSees_analysis();

wf_section.frc = 0.3*Fy
(data_volume_L1sirs,data_height_L1sirs) = wf_section.perform_OpenSees_analysis();

line1, = plt.plot(   data_volume_L1s/(TW*L1)/mm,    data_height_L1s/mm, 'k-')
line2, = plt.plot(  data_volume_L1si/(TW*L1)/mm,   data_height_L1si/mm, 'k--')
line3, = plt.plot(data_volume_L1sirs/(TW*L1)/mm, data_height_L1sirs/mm, 'r-.')
plt.legend((line1, line2, line3), ('Elastic', 'Inelastic (no residual stress)', 'Inelastic'),frameon=False)
plt.xlabel('Normalized Water Volume (V/SL, mm)\n(a)')
plt.ylabel('Water Level (mm)')
plt.xlim(0,600)
plt.ylim(0,600)

# Save line data to CSV files with legends
np.savetxt('studies/reliability/line1.csv', np.column_stack((data_volume_L1s/(TW*L1)/mm, data_height_L1s/mm)), delimiter=',', header='Normalized Water Volume,Water Level,Legend', comments='', fmt='%f, %f, Elastic')
np.savetxt('studies/reliability/line2.csv', np.column_stack((data_volume_L1si/(TW*L1)/mm, data_height_L1si/mm)), delimiter=',', header='Normalized Water Volume,Water Level,Legend', comments='', fmt='%f, %f, Inelastic (no residual stress)')
np.savetxt('studies/reliability/line3.csv', np.column_stack((data_volume_L1sirs/(TW*L1)/mm, data_height_L1sirs/mm)), delimiter=',', header='Normalized Water Volume,Water Level,Legend', comments='', fmt='%f, %f, Inelastic')


# plt.savefig('Example_Plot_3.png',dpi=300)
plt.savefig('Example_Plot_3.pdf')
# print(wf_section.maximum_permitted_zw('DAMP'))

plt.show()
