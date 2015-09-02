#!/usr/bin/env python

#########################################################################
# BlueGene Mapfile Generation Script using RUBIK        Mon 01. Jul. 2013
#            Jens Henrik Goebbert, jens.henrik.goebbert(at)rwth-aachen.de
# Mapping 2D processor grid to 5D-torus
#########################################################################

import sys, getopt
from rubik import *

def main(argv):

  # get command line arguments
  bg_shape = '2x2x2x2x2'
  bg_nranks= 64
  mapfile  = 'mapfile.map'
  try:
    opts, args = getopt.getopt(argv,"s:n:m:",["shape=","nranks=","mapfile="])
  except getopt.GetoptError:
    print 'trans2map.py -s <bg_shape> -nr <node_ranks> -m <mapfile>'
    sys.exit(2)
  for opt, arg in opts:
    if opt in ("-s", "--shape"):
      bg_shape = arg
    elif opt in ("-n", "--nranks"):
      bg_nranks = arg
    elif opt in ("-m", "--mapfile"):
      mapfile = arg
  print 'bg_shape  is "', bg_shape, ' == ', bg_shape.split('x')
  print 'bg_nranks is "', int(bg_nranks)
  print 'map file  is "', mapfile
	
  # get vars from bg_shape
  shapeA = int(bg_shape.split('x')[0])
  shapeB = int(bg_shape.split('x')[1])
  shapeC = int(bg_shape.split('x')[2])
  shapeD = int(bg_shape.split('x')[3])
  shapeE = int(bg_shape.split('x')[4])
  shapeT = int(bg_nranks)	
  mpi_net = shapeA *shapeB *shapeC *shapeD *shapeE
  mpi_cnt = shapeA *shapeB *shapeC *shapeD *shapeE *shapeT
	
  # application topology, app(i*j) == torus(x*y*z)
  app = box([mpi_net, shapeT])
  app.tile( [1,shapeT]) 
	
  # processor topology
  torus = box([shapeA,shapeB,shapeC,shapeD,shapeE,shapeT])
  torus.tile( [1,1,1,1,1,shapeT])
  torus.map(app)
	
  f = open(mapfile, 'w')
  torus.write_map_file(f)
  f.close()

if __name__ == "__main__":
   main(sys.argv[1:])
