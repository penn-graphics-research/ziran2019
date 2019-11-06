import subprocess

Jp_array = ["1", "0.98"];
beta_array = ["4", "3", "2", "1"];
for Jp in Jp_array:
	for beta in beta_array:
		runCommand = './mpm -test 19 --3d -Jp ' + Jp + ' -beta ' + beta
		subprocess.call([runCommand], shell=True)
