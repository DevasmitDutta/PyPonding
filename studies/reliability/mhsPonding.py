# exec(open('/home/devasmit/Desktop/PyPonding/studies/reliability/fragility_assessment copy.py').read())
# fragility_assessment_copy('W14X22',360.0,0.020833333333333332,6.944444444444444e-05,1976)

exec(open('/home/devasmit/Desktop/PyPonding/studies/reliability/fragility_assessment copy 2.py').read())
run_analysis('W14X22',360.0,0.020833333333333332,6.944444444444444e-05,1976)
# print(360 == 30.0 * 12, 0.020833333333333332 == 0.25/12, 6.944444444444444e-05 == 10*((1/1000)/(12**2)), 1976 == 1976)
# print(6.95e-5 * 2)
# run_analysis('W24X55',40.0*12,0.25/12,6.95e-5 * 2,1976)

# import os

# # Print the current working directory before executing the script
# print("Current working directory before exec:", os.getcwd())

# # Execute the script
# try:
#     exec(open('/home/devasmit/Desktop/PyPonding/studies/reliability/fragility_assessment copy.py').read())
# except Exception as e:
#     print("Error executing script:", e)

# # Print the current working directory after executing the script
# print("Current working directory after exec:", os.getcwd())

# # Call the function from the executed script
# try:
#     fragility_assessment_copy('W14X22', 360.0, 0.020833333333333332, 6.944444444444444e-05, 1976)
# except Exception as e:
#     print("Error calling fragility_assessment_copy:", e)