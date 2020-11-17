#!/usr/bin/python3

import numpy as np   # For array operations.
import sys

nx = 0 # number of points in Ox; if 0, then search dns.ini
ny = 0 # number of points in Oy; if 0, then search dns.ini
nz = 0 # number of points in Oz; if 0, then search dns.ini

sizeofmask = 6

# do not edit below this line

# getting grid size from dns.ini, if necessary
if ( nx == 0 ):
    for line in open('dns.ini'):
        if line.lower().replace(" ","").startswith("imax="):
            nx = int(line.split("=",1)[1])
            break

if ( ny == 0 ):
    for line in open('dns.ini'):
        if line.lower().replace(" ","").startswith("jmax="):
            ny = int(line.split("=",1)[1])
            break

if ( nz == 0 ):
    for line in open('dns.ini'):
        if line.lower().replace(" ","").startswith("kmax="):
            nz = int(line.split("=",1)[1])
            break

print("Grid size is {}x{}x{}.".format(nx,ny,nz))

def itnumber(filename):
    main = filename.split(".",1)[0]
    return int(main[len(main)-sizeofmask:len(main)])

if ( len(sys.argv) == 1 ):
    print("Usage: python $0 list-of-files.")
    quit()

filenames = sorted(sys.argv[1:],key=itnumber)

filetypes = []
for name in filenames:
    main = name.split(".",1)[0]
    type = main[:len(main)-sizeofmask]
    if not (any(type in s for s in filetypes)):
        filetypes.append(type)

filetimes = []
for name in filenames:
    main = name.split(".",1)[0]
    time = main[len(main)-sizeofmask:len(main)]
    if not (any(time in s for s in filetimes)):
        filetimes.append(time)

print("Processing %d times in fields %s..." % ( len(filetimes) , ', '.join(filetypes) ))

f = open('main'+'.times.xdmf', 'w')

# Definining entities depending on planesmode

f.write('''<?xml version="1.0" ?>

<!--
XDMF file to read time collections of plane data from TLAB into a data analysis and visualization application, like ParaView.

Add simply data items after the geometry block to read different files.

The structure of this file has been adapted from psOpen, from Jens Henrik Goebbert.
-->

<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" [
''')

data = ( nx,ny,nz, nx,ny,nz )
f.write('''
<!-- dimension of complete datasets -->
<!ENTITY DimsX   "%d">
<!ENTITY DimsY   "%d">
<!ENTITY DimsZ   "%d">

<!-- dimension of hyperslab to load -->
<!ENTITY HSDimsX  "%d">
<!ENTITY HSDimsY  "%d">
<!ENTITY HSDimsZ  "%d">

<!-- start of hyperslab in complete dataset -->
<!ENTITY HSDimsX_Start "0">
<!ENTITY HSDimsY_Start "0">
<!ENTITY HSDimsZ_Start "0">

<!-- stride of hyperslab in complete dataset -->
<!ENTITY HSStrideX "1">
<!ENTITY HSStrideY "1">
<!ENTITY HSStrideZ "1">

<!-- data precision (grid is always 8 bytes) -->
<!ENTITY Prec      "4">
''' % data )

data = (56+nx*8+8, 56+nx*8+8+ny*8+8)
f.write('''
<!-- offsets to grid blocks -->
<!ENTITY SeekGridX  "56">
<!ENTITY SeekGridY  "%d"> <!-- + DimX*8 + 8-->
<!ENTITY SeekGridZ  "%d"> <!-- + DimY*8 + 8-->

<!-- offsets to data -->
<!ENTITY SeekData   "0"> <!-- No header -->
<!-- <!ENTITY SeekData  "52"> --> <!-- Tlab header -->
<!-- <!ENTITY SeekData "244"> --> <!-- Ensight header -->
''' % data )

# code below is independent of filetype
f.write('''
]>

<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">
  <Domain Name="main">

    <!-- Hyperslab metadata referenced below -->
    <DataItem Name="HSMetaData" Dimensions="3 3" Format="XML">
      &HSDimsZ_Start; &HSDimsY_Start; &HSDimsX_Start;
      &HSStrideZ;     &HSStrideY;     &HSStrideX;
      &HSDimsZ;       &HSDimsY;       &HSDimsX;
    </DataItem>

    <!-- Defining common topology and common grid to all timeslices -->
    <Topology TopologyType="3DRectMesh" Dimensions="&HSDimsZ; &HSDimsY; &HSDimsX;">
    </Topology>

    <Geometry GeometryType="VXVYVZ">

      <DataItem Name="X" ItemType="HyperSlab" Dimensions="&HSDimsX;">
	<DataItem Dimensions="1 3" Format="XML">
          &HSDimsX_Start;
	  &HSStrideX;
	  &HSDimsX;
	</DataItem>
	<DataItem ItemType="Uniform" Format="Binary" Seek="&SeekGridX;" NumberType="Float" Precision="8" Endian="Little" Dimensions="&DimsX;">
	  grid
	</DataItem>
      </DataItem>

      <DataItem Name="Y" ItemType="HyperSlab" Dimensions="&HSDimsY;">
	<DataItem Dimensions="1 3" Format="XML">
	  &HSDimsY_Start;
	  &HSStrideY;
	  &HSDimsY;
	</DataItem>
	<DataItem ItemType="Uniform" Format="Binary" Seek="&SeekGridY;" NumberType="Float" Precision="8" Endian="Little" Dimensions="&DimsY;">
	  grid
	</DataItem>
      </DataItem>

      <DataItem Name="Z" ItemType="HyperSlab" Dimensions="&HSDimsZ;">
	<DataItem Dimensions="1 3" Format="XML">
	  &HSDimsZ_Start;
	  &HSStrideZ;
	  &HSDimsZ;
	</DataItem>
	<DataItem ItemType="Uniform" Format="Binary" Seek="&SeekGridZ;" NumberType="Float" Precision="8" Endian="Little" Dimensions="&DimsZ;">
	  grid
	</DataItem>
      </DataItem>

    </Geometry>

    <!-- Collection of timeslices -->
    <Grid GridType="Collection" CollectionType="Temporal">

      <Time TimeType="HyperSlab">
	<DataItem Format="XML" NumberType="Float" Dimensions="3"> <!-- start, stride, count-->
	  0.0 1.0 %d;
	</DataItem>
      </Time>
''' % (len(filetimes)) )

# Loop over timeslices
for time in filetimes:
    f.write('''
      <!-- Timeslice -->
      <Grid Name="It%d" GridType="Uniform">
	<Topology Reference="/Xdmf/Domain/Topology[1]"/>
	<Geometry Reference="/Xdmf/Domain/Geometry[1]"/>
    ''' % (int(time)) )

    for type in filetypes:
        if ( type in ['VelocityVector','VorticityVector'] ):
            f.write('''
        <Attribute AttributeType="Vector" Name="%s">
	  <DataItem ItemType="Function" Function="JOIN($0,$1,$2)" Dimensions="&HSDimsZ; &HSDimsY; &HSDimsX; 3">

	    <DataItem ItemType="HyperSlab" Dimensions="&HSDimsZ; &HSDimsY; &HSDimsX;">
	      <DataItem Reference="/Xdmf/Domain/DataItem[1]"/>
	      <DataItem ItemType="Uniform" Format="Binary" Seek="&SeekData;" NumberType="Float" Precision="4" Endian="Little" Dimensions="&DimsZ; &DimsY; &DimsX;">
		%s
	      </DataItem>
	    </DataItem>

	    <DataItem ItemType="HyperSlab" Dimensions="&HSDimsZ; &HSDimsY; &HSDimsX;">
	      <DataItem Reference="/Xdmf/Domain/DataItem[1]"/>
	      <DataItem ItemType="Uniform" Format="Binary" Seek="&SeekData;" NumberType="Float" Precision="4" Endian="Little" Dimensions="&DimsZ; &DimsY; &DimsX;">
		%s
	      </DataItem>
	    </DataItem>

	    <DataItem ItemType="HyperSlab" Dimensions="&HSDimsZ; &HSDimsY; &HSDimsX;">
	      <DataItem Reference="/Xdmf/Domain/DataItem[1]"/>
	      <DataItem ItemType="Uniform" Format="Binary" Seek="&SeekData;" NumberType="Float" Precision="4" Endian="Little" Dimensions="&DimsZ; &DimsY; &DimsX;">
		%s
	      </DataItem>
	    </DataItem>

	  </DataItem>
	</Attribute>
''' % (type, type+time+'.1',type+time+'.2',type+time+'.3') )
        elif ( type in ['StrainTensor','ReynoldsTensor'] ):
            f.write('''
        <Attribute AttributeType="Tensor6" Name="%s">
	  <DataItem ItemType="Function" Function="JOIN($0,$1,$2,$3,$4,$5)" Dimensions="&HSDimsZ; &HSDimsY; &HSDimsX; 6">

	    <DataItem ItemType="HyperSlab" Dimensions="&HSDimsZ; &HSDimsY; &HSDimsX;">
	      <DataItem Reference="/Xdmf/Domain/DataItem[1]"/>
	      <DataItem ItemType="Uniform" Format="Binary" Seek="&SeekData;" NumberType="Float" Precision="4" Endian="Little" Dimensions="&DimsZ; &DimsY; &DimsX;">
		%s
	      </DataItem>
	    </DataItem>

	    <DataItem ItemType="HyperSlab" Dimensions="&HSDimsZ; &HSDimsY; &HSDimsX;">
	      <DataItem Reference="/Xdmf/Domain/DataItem[1]"/>
	      <DataItem ItemType="Uniform" Format="Binary" Seek="&SeekData;" NumberType="Float" Precision="4" Endian="Little" Dimensions="&DimsZ; &DimsY; &DimsX;">
		%s
	      </DataItem>
	    </DataItem>

	    <DataItem ItemType="HyperSlab" Dimensions="&HSDimsZ; &HSDimsY; &HSDimsX;">
	      <DataItem Reference="/Xdmf/Domain/DataItem[1]"/>
	      <DataItem ItemType="Uniform" Format="Binary" Seek="&SeekData;" NumberType="Float" Precision="4" Endian="Little" Dimensions="&DimsZ; &DimsY; &DimsX;">
		%s
	      </DataItem>
	    </DataItem>

	    <DataItem ItemType="HyperSlab" Dimensions="&HSDimsZ; &HSDimsY; &HSDimsX;">
	      <DataItem Reference="/Xdmf/Domain/DataItem[1]"/>
	      <DataItem ItemType="Uniform" Format="Binary" Seek="&SeekData;" NumberType="Float" Precision="4" Endian="Little" Dimensions="&DimsZ; &DimsY; &DimsX;">
		%s
	      </DataItem>
	    </DataItem>

	    <DataItem ItemType="HyperSlab" Dimensions="&HSDimsZ; &HSDimsY; &HSDimsX;">
	      <DataItem Reference="/Xdmf/Domain/DataItem[1]"/>
	      <DataItem ItemType="Uniform" Format="Binary" Seek="&SeekData;" NumberType="Float" Precision="4" Endian="Little" Dimensions="&DimsZ; &DimsY; &DimsX;">
		%s
	      </DataItem>
	    </DataItem>

	    <DataItem ItemType="HyperSlab" Dimensions="&HSDimsZ; &HSDimsY; &HSDimsX;">
	      <DataItem Reference="/Xdmf/Domain/DataItem[1]"/>
	      <DataItem ItemType="Uniform" Format="Binary" Seek="&SeekData;" NumberType="Float" Precision="4" Endian="Little" Dimensions="&DimsZ; &DimsY; &DimsX;">
		%s
	      </DataItem>
	    </DataItem>

        </DataItem>
	</Attribute>
''' % (type, type+time+'.1',type+time+'.4',type+time+'.5', type+time+'.2',type+time+'.6',type+time+'.3') )
        else:
            f.write('''
	<Attribute Center="Node" Name="%s">
	  <DataItem ItemType="HyperSlab" Dimensions="&HSDimsZ; &HSDimsY; &HSDimsX;">
	    <DataItem Reference="/Xdmf/Domain/DataItem[1]"/>
	    <DataItem ItemType="Uniform" Format="Binary" Seek="&SeekData;" NumberType="Float" Precision="&Prec;" Endian="Little" Dimensions="&DimsZ; &DimsY; &DimsX;">
	      %s
	    </DataItem>
	  </DataItem>
	</Attribute>
''' % (type,type+time) )

    f.write('''
      </Grid>
''' )

# End the xmf file
f.write('''
    </Grid> <!-- End of time collection -->

  </Domain>
</Xdmf>
''')

f.close()
