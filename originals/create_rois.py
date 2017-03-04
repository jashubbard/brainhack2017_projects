#!/usr/bin/env python

#Creates a spherical ROI from a set of coordinates in a text file
#


from __future__ import division
import argparse
import os
from subprocess import call
import glob

#parse command line arguments
parser = argparse.ArgumentParser(description="create spherical ROIs out of a list of coordinates in a text file")
parser.add_argument("infile",help="coordinate file. should be a text file with one x y z coordinate on each line, separated by whitespace")
parser.add_argument("radius",help="radius in mm")
parser.add_argument("--mm",help="coordinates are mm coordinates in MNI space instead of voxel coordinates",action="store_true")
args = parser.parse_args()

def mmToVox(mmcoords):
	#function to convert mm coordinates in the standard 2mm MNI atlas into voxel coordinates
	voxcoords = ['','','']
	voxcoords[0] = str(int(round(int(mmcoords[0])/2))*-1+45)
	voxcoords[1] = str(int(round(int(mmcoords[1])/2))+63)
	voxcoords[2] = str(int(round(int(mmcoords[2])/2))+36)
	return voxcoords

if os.path.exists(args.infile):
	cfile = open(args.infile)

	
	if not os.path.exists('./rois'):
		os.mkdir('rois')
	
	os.chdir('./rois')
	
	counter =1
	
	while 1:
		line = cfile.readline()
		if not line:
			break

		coords = line.split()

		if args.mm:
			#convert from mm to voxel coords
			voxcoords = mmToVox(coords)
		else:
			voxcoords = coords

		outfile = "roi-%s-%s-%s" % (coords[0],coords[1],coords[2])
		
		#only create if they don't exist yet-- good for long lists of ROIs
		if not os.path.exists(os.path.abspath('./'+outfile+".nii.gz")):
		  command = "fslmaths $FSLDIR/data/standard/MNI152_T1_2mm_brain_mask_dil1 -roi %s 1 %s 1 %s 1 0 1 tmp" % (voxcoords[0],voxcoords[1],voxcoords[2])
		  print command
		  call(command,shell=True)

		  command = "fslmaths tmp -kernel sphere %s -fmean tmp" % (args.radius)
		  print command
		  call(command,shell=True)

		  command = "fslmaths tmp -thr .00001 -bin %s" % outfile
		  print command
		  call(command,shell=True)
		
		  command = "fslmaths %s -mul %d %s" %(outfile,counter,outfile)
		  print command
		  call(command,shell=True)
		
	          command = "rm tmp.nii.gz"
		  call(command,shell=True)
		
		counter = counter + 1

	cfile.close()
	
	files = glob.glob('roi*')
	tmp = ["-add " + i for i in files[1:len(files)]]

	tmp = " ".join(tmp)

	command = "fslmaths " + files[0] + " " + tmp + " allrois"
	print command
	call(command,shell=True)
	
	command = "fslview $FSLDIR/data/standard/MNI152_T1_2mm_brain allrois"
	print command
	call(command,shell=True)


else:
	print "error: the file %s does not exist" % args.infile