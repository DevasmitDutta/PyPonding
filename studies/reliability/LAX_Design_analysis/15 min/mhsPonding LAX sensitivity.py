
# run_analysis('W14X22',360.0,0.020833333333333332,6.944444444444444e-05,1976)
# print(360 == 30.0 * 12, 0.020833333333333332 == 0.25/12, 6.944444444444444e-05 == 10*((1/1000)/(12**2)), 1976 == 1976)

exec(open('/home/devasmit/Desktop/PyPonding/studies/reliability/LAX_Design_analysis/15 min/fragility_assessment LAX sensitivity.py').read())
# for shape in ['W12X35','W8X31']:
for shape in ['W8X24','W10X17']:
# for shape in ['W10X17']:
    run_analysis_multiple_periods(shape,360.0,0.020833333333333332,6.944444444444444e-05,1976)

# run_analysis_multiple_periods('W10X12',360.0,0.020833333333333332,6.944444444444444e-05,1976)