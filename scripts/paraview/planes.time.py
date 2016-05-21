#!/usr/bin/python3

import numpy as np   # For array operations.
import sys

# edit as needed
npx = 1    # number of planesI
npy = 2    # number of planesJ
npz = 1    # number of planesK

nx = 128  # number of points in Ox
ny = 96   # number of points in Oy
nz = 128  # number of points in OZ

ns = 1    # number of scalar variables
nq = 3    # number of flow variables

# do not edit
def itnumber(filename):
    return int(filename.split(".",1)[1])

if ( len(sys.argv) == 1 ):
    print("Add filenames as arguments.")
    quit()

filetype  = sys.argv[1].split(".",1)[0]
print("Processing fiels %s..." % ( filetype ))

filenames = sorted(sys.argv[1:],key=itnumber)

f = open('TEST'+filetype+'.times.xdmf', 'w')

# Definining entities depending on planesmode

f.write('''<?xml version="1.0" ?>

<!--
XDMF file to read time collections of plane data from TLAB into a data analysis and visualization application, like ParaView.

Add simply data items after the geometry block to read different files.

The structure of this file has been adapted from psOpen, from Jens Henrik Goebbert.
-->

<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" [
''')

if   ( filetype == 'planesI'):
    data = ( 1,ny,nz, npx,ny,nz*(nq+ns) )
elif ( filetype == 'planesJ'):
    data = ( nx,1,nz, nx,npy,nz*(nq+ns) )
elif ( filetype == 'planesK'):
    data = ( nx,ny,1, nx,ny,npz*(nq+ns) )
f.write('''
<!-- dimension of complete datasets -->
<!ENTITY GridDimsX   "%d">
<!ENTITY GridDimsY   "%d">
<!ENTITY GridDimsZ   "%d">

<!ENTITY DimsX   "%d">
<!ENTITY DimsY   "%d">
<!ENTITY DimsZ   "%d">
''' % data )

if   ( filetype == 'planesI'):
    data = (  1,ny,nz, npx,1,nq+ns )
elif ( filetype == 'planesJ'):
    data = ( nx, 1,nz, 1,npy,nq+ns )
elif ( filetype == 'planesK'):
    data = ( nx,ny, 1, 1,1,npz*(nq+ns) )
f.write('''
<!-- dimension of hyperslab to load -->
<!ENTITY HSDimsX  "%d">
<!ENTITY HSDimsY  "%d">
<!ENTITY HSDimsZ  "%d">

<!-- start of hyperslab in complete dataset -->
<!ENTITY HSDimsX_Start "0">
<!ENTITY HSDimsY_Start "0">
<!ENTITY HSDimsZ_Start "0">

<!ENTITY HSGridDimsX_Start "0"> <!-- Set here the plane position -->
<!ENTITY HSGridDimsY_Start "0">
<!ENTITY HSGridDimsZ_Start "0">

<!-- stride of hyperslab in complete dataset -->
<!ENTITY HSStrideX "%d">
<!ENTITY HSStrideY "%d">
<!ENTITY HSStrideZ "%d">

<!ENTITY HSGridStrideX "1">
<!ENTITY HSGridStrideY "1">
<!ENTITY HSGridStrideZ "1">
''' % data )

if   ( filetype == 'planesI'):
    data = (56+nx*8+8, 56+nx*8+8+ny*8+8, ny*npx*4, ny*npx*4*2)
elif ( filetype == 'planesJ'):
    data = (56+nx*8+8, 56+nx*8+8+ny*8+8, nx*npy*4, nx*npy*4*2)
elif ( filetype == 'planesK'):
    data = (56+nx*8+8, 56+nx*8+8+ny*8+8, nx*ny *4, nx*ny *4*2)
f.write('''
<!-- offsets to grid blocks -->
<!ENTITY SeekGridX  "56"> 
<!ENTITY SeekGridY  "%d"> <!-- + DimX*8 + 8-->
<!ENTITY SeekGridZ  "%d"> <!-- + DimY*8 + 8-->

<!-- offsets to data -->
<!ENTITY SeekDataU  "0">
<!ENTITY SeekDataV  "%d"> <!-- + Dim*Prec-->
<!ENTITY SeekDataW  "%d"> <!-- + Dim*Prec-->
''' % data )

for i in range(ns):
    if   ( filetype == 'planesI'):
        data = (i+1,ny*npx*4*(i+3))
    elif ( filetype == 'planesJ'):
        data = (i+1,nx*npy*4*(i+3))
    elif ( filetype == 'planesK'):
        data = (i+1,nx*ny *4*(i+3))
    f.write('''<!ENTITY SeekDataS%d "%d"> <!-- + Dim*Prec-->
''' % data )

# code below is independent of filetype
f.write('''
]>

<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">
  <Domain Name="%s">
    
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
          &HSGridDimsX_Start;
	  &HSGridStrideX;
	  &HSDimsX;
	</DataItem>
	<DataItem ItemType="Uniform" Format="Binary" Seek="&SeekGridX;" NumberType="Float" Precision="8" Endian="Big" Dimensions="&GridDimsX;">
	  grid
	</DataItem>
      </DataItem>
      
      <DataItem Name="Y" ItemType="HyperSlab" Dimensions="&HSDimsY;">
	<DataItem Dimensions="1 3" Format="XML">
	  &HSGridDimsY_Start;
	  &HSGridStrideY;
	  &HSDimsY;
	</DataItem>
	<DataItem ItemType="Uniform" Format="Binary" Seek="&SeekGridY;" NumberType="Float" Precision="8" Endian="Big" Dimensions="&GridDimsY;">
	  grid
	</DataItem>
      </DataItem>
      
      <DataItem Name="Z" ItemType="HyperSlab" Dimensions="&HSDimsZ;">
	<DataItem Dimensions="1 3" Format="XML">
	  &HSGridDimsZ_Start;
	  &HSGridStrideZ;
	  &HSDimsZ;
	</DataItem>
	<DataItem ItemType="Uniform" Format="Binary" Seek="&SeekGridZ;" NumberType="Float" Precision="8" Endian="Big" Dimensions="&GridDimsZ;">
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
''' % (filetype,len(filenames)) )

# Loop over timeslices
for file in filenames:
    f.write('''
      <!-- Timeslice -->
      <Grid Name="It%d" GridType="Uniform">
	<Topology Reference="/Xdmf/Domain/Topology[1]"/>
	<Geometry Reference="/Xdmf/Domain/Geometry[1]"/>	
    ''' % (itnumber(file)) )
    
    for i in range(ns):
        f.write('''
	<Attribute Center="Node" Name="Scalar%d">
	  <DataItem ItemType="HyperSlab" Dimensions="&HSDimsZ; &HSDimsY; &HSDimsX;">
	    <DataItem Reference="/Xdmf/Domain/DataItem[1]"/>
	    <DataItem ItemType="Uniform" Format="Binary" Seek="&SeekDataS%d;" NumberType="Float" Precision="4" Endian="Big" Dimensions="&DimsZ; &DimsY; &DimsX;">
	      %s
	    </DataItem>
	  </DataItem>
	</Attribute>	
''' % (i+1,i+1,file) )

    f.write('''
        <Attribute AttributeType="Vector" Name="Velocity">
	  <DataItem ItemType="Function" Function="JOIN($0,$1,$2)" Dimensions="&HSDimsZ; &HSDimsY; &HSDimsX; 3">
	    
	    <DataItem ItemType="HyperSlab" Dimensions="&HSDimsZ; &HSDimsY; &HSDimsX;">
	      <DataItem Reference="/Xdmf/Domain/DataItem[1]"/>
	      <DataItem ItemType="Uniform" Format="Binary" Seek="&SeekDataU;" NumberType="Float" Precision="4" Endian="Big" Dimensions="&DimsZ; &DimsY; &DimsX;">
		%s
	      </DataItem>
	    </DataItem>

	    <DataItem ItemType="HyperSlab" Dimensions="&HSDimsZ; &HSDimsY; &HSDimsX;">
	      <DataItem Reference="/Xdmf/Domain/DataItem[1]"/>
	      <DataItem ItemType="Uniform" Format="Binary" Seek="&SeekDataV;" NumberType="Float" Precision="4" Endian="Big" Dimensions="&DimsZ; &DimsY; &DimsX;">
		%s
	      </DataItem>
	    </DataItem>
	    
	    <DataItem ItemType="HyperSlab" Dimensions="&HSDimsZ; &HSDimsY; &HSDimsX;">
	      <DataItem Reference="/Xdmf/Domain/DataItem[1]"/>
	      <DataItem ItemType="Uniform" Format="Binary" Seek="&SeekDataW;" NumberType="Float" Precision="4" Endian="Big" Dimensions="&DimsZ; &DimsY; &DimsX;">
		%s
	      </DataItem>
	    </DataItem>
	    
	  </DataItem>
	</Attribute>
	
      </Grid>
''' % (file,file,file) )

# End the xmf file
f.write('''
    </Grid> <!-- End of time collection -->

  </Domain>
</Xdmf>
''')

f.close()
