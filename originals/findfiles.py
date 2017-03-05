#!/usr/bin/env python


from __future__ import division
import argparse
import os
from subprocess import call
import fnmatch



def findfiles(starting_dir,searchstring):
    fmatches = []
    dmatches = []

    
    for root, dirnames, filenames in os.walk(starting_dir):
      for filename in fnmatch.filter(filenames, searchstring):
          fmatches.append(os.path.join(root, filename))
      for dirname in fnmatch.filter(dirnames, searchstring):
          dmatches.append(os.path.join(root, dirname))
            
    return {'files': fmatches, 'dirs': dmatches}





if __name__ == '__main__':	
	#parse command line arguments
	parser = argparse.ArgumentParser(description="find files and subdirectories given a path and a search string")
	parser.add_argument("dir",help="directory to search within")
	parser.add_argument("search",help="regex-style search string to find files or directories",type=str)
	parser.add_argument("type",help="either 'files' (default) or 'dirs'",nargs='?', default='files')
	args = parser.parse_args()
	results = findfiles(args.dir,args.search)
	for p in results[args.type]: print(p)
	#print '\n\n'
	#print results['dirs']