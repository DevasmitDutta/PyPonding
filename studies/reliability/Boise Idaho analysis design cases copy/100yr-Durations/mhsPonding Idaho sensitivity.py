
# run_analysis('W14X22',360.0,0.020833333333333332,6.944444444444444e-05,1976)
# print(360 == 30.0 * 12, 0.020833333333333332 == 0.25/12, 6.944444444444444e-05 == 10*((1/1000)/(12**2)), 1976 == 1976)

exec(open('/home/devasmit/Desktop/PyPonding/studies/reliability/Boise Idaho analysis design cases copy/100yr-Durations/fragility_assessment Idaho sensitivity.py').read())
# run_analysis_multiple_periods('W14X22',360.0,0.020833333333333332,6.944444444444444e-05,1976)
# shapes = ['W6X20','W10X12','W12X14']
# shapes = ['W10X12', 'W8X15','W6X25']
# shapes = ['W8X18']
shapes = ['W10X15']
# shapes = ['W8X21']

for shape in shapes:
    run_analysis_multiple_periods(shape,360.0,0.020833333333333332,6.944444444444444e-05,1976)