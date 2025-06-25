
# run_analysis('W14X22',360.0,0.020833333333333332,6.944444444444444e-05,1976)
# print(360 == 30.0 * 12, 0.020833333333333332 == 0.25/12, 6.944444444444444e-05 == 10*((1/1000)/(12**2)), 1976 == 1976)

exec(open('/home/devasmit/Desktop/PyPonding/studies/reliability/Boise Idaho analysis design cases copy/500yr-Durations/fragility_assessment MB sensitivity.py').read())
# shapes = ['W8X24','W10X17','W8X18']
# shapes = ['W10X17','W8X18','W10X15']
shapes = ['W8X18']
for shape in shapes:
    run_analysis_multiple_periods(shape,360.0,0.020833333333333332,6.944444444444444e-05,1976)
    # run_analysis_multiple_periods('W8X21',360.0,0.020833333333333332,6.944444444444444e-05,1976)
    # 360 == 30.0 * 12, 0.020833333333333332 == 0.25/12, 6.944444444444444e-05 == 10*((1/1000)/(12**2)), 1976 == 1976
