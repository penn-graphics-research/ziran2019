import os
import subprocess

def render(inputFolder, outputFolder):
	# prepare data to render
	os.rename(inputFolder,'inputData')
	os.mkdir('frames')

	# render
	runCommand = '/opt/hfs16.5.571/bin/hbatch /home/minchen/Downloads/pumpkin_fanfu_testnig.hipnc -c \'render -V mantra2\' -c \'quit\''
	subprocess.call([runCommand], shell=True)

	# archive rendered data and input data
	os.rename('inputData', inputFolder)
	os.rename('frames', outputFolder)



inputFolder_array = [
	'pumpkin_smash_3d_0.94_1._0.39_2000.',
	'pumpkin_smash_3d_0.94_2._0.39_2000.',
	'pumpkin_smash_3d_0.94_3._0.39_2000.',
	'pumpkin_smash_3d_0.94_4._0.39_2000.',
	'pumpkin_smash_3d_0.96_1._0.39_2000.',
	'pumpkin_smash_3d_0.96_2._0.39_2000.',
	'pumpkin_smash_3d_0.96_3._0.39_2000.',
	'pumpkin_smash_3d_0.96_4._0.39_2000.'
];

for inputFolder in inputFolder_array:
	print('### rendering ' + inputFolder + " ...")
	render(inputFolder, 'frames_' + inputFolder)


