
# run_analysis('W14X22',360.0,0.020833333333333332,6.944444444444444e-05,1976)
# print(360 == 30.0 * 12, 0.020833333333333332 == 0.25/12, 6.944444444444444e-05 == 10*((1/1000)/(12**2)), 1976 == 1976)

exec(open('/home/devasmit/Desktop/PyPonding/studies/reliability/SB_Design_risk_analysis/15 min copy/fragility_assessment SB sensitivity.py').read())
# '10x17'
for shape in ['W14X22']:
# for shape in ['W10X39','W10X22','W8X21','W10X15']:


   run_analysis_multiple_periods(shape,360.0,0.020833333333333332,6.944444444444444e-05,1976)

# run_analysis_multiple_periods('W10X12',360.0,0.020833333333333332,6.944444444444444e-05,1976)