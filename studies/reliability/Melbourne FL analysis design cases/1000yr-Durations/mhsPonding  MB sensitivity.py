
# run_analysis('W14X22',360.0,0.020833333333333332,6.944444444444444e-05,1976)
# print(360 == 30.0 * 12, 0.020833333333333332 == 0.25/12, 6.944444444444444e-05 == 10*((1/1000)/(12**2)), 1976 == 1976)

exec(open('/home/devasmit/Desktop/PyPonding/studies/reliability/Melbourne FL analysis design cases/1000yr-Durations/fragility_assessment MB sensitivity.py').read())
# shapes = ['W10X15', 'W8X18', 'W6X25']
# shapes = ['W6X25']
shapes = ['W6X20']

# shapes = ['W8X18']
for shape in shapes:
    run_analysis_multiple_periods(shape,360.0,0.020833333333333332,6.944444444444444e-05,1976)
# run_analysis_multiple_periods('W8X18',360.0,0.020833333333333332,6.944444444444444e-05,1976)
# from libdenavit.section.database import wide_flange_database
# SR_ratios =[]
# shape_names = wide_flange_database.keys()
# cnt = 0
# for shape in shape_names:
#     cnt += 1
#     SR_Ratio = run_analysis_multiple_periods(shape,360.0,0.020833333333333332,6.944444444444444e-05,1976)
#     SR_ratios.append([SR_Ratio, shape])
#     print('finished shape', cnt, 'of', len(shape_names), 'SR_Ratio:', SR_Ratio, 'shape:', shape)

# print(SR_ratios)
