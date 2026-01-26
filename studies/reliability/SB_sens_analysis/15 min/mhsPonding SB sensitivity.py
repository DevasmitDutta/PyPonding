
# run_analysis('W14X22',360.0,0.020833333333333332,6.944444444444444e-05,1976)
# print(360 == 30.0 * 12, 0.020833333333333332 == 0.25/12, 6.944444444444444e-05 == 10*((1/1000)/(12**2)), 1976 == 1976)

from matplotlib.pylab import gamma


exec(open('/home/devasmit/Desktop/PyPonding/studies/reliability/SB_sens_analysis/15 min/fragility_assessment SB sensitivity.py').read())
# '10x17'
for shape in ['W8X18']:
# [0.75, 0.5, 0.25] # [1.25, 1.5, 2]
  for beta in [1]:
      for gamm in [1]:
         for alpha in [1.25, 1.5, 2]:


            run_analysis_multiple_periods(shape,360.0,0.020833333333333332,6.944444444444444e-05,1976,beta, gamm, alpha)

# run_analysis_multiple_periods('W10X12',360.0,0.020833333333333332,6.944444444444444e-05,1976)