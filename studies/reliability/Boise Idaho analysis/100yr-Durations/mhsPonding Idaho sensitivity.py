
# run_analysis('W14X22',360.0,0.020833333333333332,6.944444444444444e-05,1976)
# print(360 == 30.0 * 12, 0.020833333333333332 == 0.25/12, 6.944444444444444e-05 == 10*((1/1000)/(12**2)), 1976 == 1976)

exec(open('/home/devasmit/Desktop/PyPonding/studies/reliability/Boise Idaho analysis/100yr-Durations/fragility_assessment Idaho sensitivity.py').read())
run_analysis_multiple_periods('W14X22',360.0,0.020833333333333332,6.944444444444444e-05,1976)