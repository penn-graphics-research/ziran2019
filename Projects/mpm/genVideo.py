import os
import subprocess

folder_array = [
	'pumpkin_smash_3d_0.94_1._0.39_2000.',
	'pumpkin_smash_3d_0.94_2._0.39_2000.',
	'pumpkin_smash_3d_0.94_3._0.39_2000.',
	'pumpkin_smash_3d_0.94_4._0.39_2000.',
	'pumpkin_smash_3d_0.96_1._0.39_2000.',
	'pumpkin_smash_3d_0.96_2._0.39_2000.',
	'pumpkin_smash_3d_0.96_3._0.39_2000.',
	'pumpkin_smash_3d_0.96_4._0.39_2000.',
	'pumpkin_smash_3d_0.98_1._0.39_2000.',
	'pumpkin_smash_3d_0.98_2._0.39_2000.',
	'pumpkin_smash_3d_0.98_3._0.39_2000.',
	'pumpkin_smash_3d_0.98_4._0.39_2000.',
	'pumpkin_smash_3d_1._1._0.39_2000.',
	'pumpkin_smash_3d_1._2._0.39_2000.',
	'pumpkin_smash_3d_1._3._0.39_2000.',
	'pumpkin_smash_3d_1._4._0.39_2000.'
];

for folder in folder_array:
	os.chdir('/home/minchen/CLionProjects/ziran/Projects/mpm/output/frames_' + folder)

	print('### ' + folder)

	print('turning transparent background to white...')
	runCommand = 'for i in `seq 0 80`; do convert $i.png -alpha remove -background white white$i.png; done'
	subprocess.call([runCommand], shell=True)

	print('generating video...')
	runCommand = 'ffmpeg -framerate 24 -i "white%d.png" -s 1280x720 -c:v libx264 -profile:v high -crf 10 -pix_fmt yuv420p -threads 11 ' + folder + '.mp4'
	subprocess.call([runCommand], shell=True)

