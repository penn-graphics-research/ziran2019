import subprocess

#Uncomment this block to run M comparisons!

M_array = [2.9, 9, 16, 34, 60] #these are the final values of M for this comparison (percent = 0.03)
for M in M_array:
    runCommand = './fracture -test 5 -M ' + str(M)
    subprocess.call([runCommand], shell=True)


#Uncomment this block to run sigmaC comparisons
"""
percent_array = [0.001, 0.009, 0.03, 0.079, 0.5] #these are the final values of sigma percent for this comparison (M = -1)
for percent in percent_array:
    runCommand = './fracture -test 5 -sigmaPercent ' + str(percent)
    subprocess.call([runCommand], shell=True)
"""