#!/usr/bin/python3
import numpy as np

nx    = 0   # number of points in Ox; if 0, then search tlab.ini
ny    = 0   # number of points in Oy; if 0, then search tlab.ini
nz    = 0   # number of points in OZ; if 0, then search tlab.ini
dtype = ''  # datatype of eps0.1 field; if empty, then search tlab.ini
fname = 'eps0.1'; path  = './'

# do not edit below this line

# getting grid size + data format from tlab.ini, if necessary
if ( nx == 0 ):
    for line in open(path + 'tlab.ini'):
        if line.lower().replace(" ","").startswith("imax="):
            nx = int(line.split("=",1)[1])
            break

if ( ny == 0 ):
    for line in open(path + 'tlab.ini'):
        if line.lower().replace(" ","").startswith("jmax="):
            ny = int(line.split("=",1)[1])
            break

if ( nz == 0 ):
    for line in open(path + 'tlab.ini'):
        if line.lower().replace(" ","").startswith("kmax="):
            nz = int(line.split("=",1)[1])
            break

print("Grid size is {}x{}x{}.".format(nx,ny,nz))

if ( len(dtype) == 0 ):
    for line in open(path + 'tlab.ini'):
        if line.lower().replace(" ","").startswith("datatypegeometry="):
            dtype = str(line.split("=",1)[1][:-1])
            break

print("Data type of eps-field: {}".format(dtype))

# convert eps bit field to int1 field 
if dtype == 'bit':
    def int2bit(out,data):
        bsize = data.size
        for i in range(bsize):
            ip = i * 8
            binary = f'{data[i]:08b}'
            if binary[0] == '-':
                binary = f'{data[i]+256:08b}'
            j = 0
            for k in range(-1,-9,-1):
                out[j+ip] = int(binary[k])
                j += 1
        return out
    # read header
    f = open(path + fname,'rb')
    f.seek(0,0)
    header = np.fromfile(f, '<i4', 5)
    f.close()
    print('Header size           :', header[0])
    print('Grid   size (nx*ny*nz):', header[1]*8,'x',header[2],'x',header[3])
    # data size 
    rsize = nx*ny*nz; bsize = rsize//8 
    # read eps field 
    f = open(path + fname,'rb')
    f.seek(header[0],0)
    data = np.fromfile(f, np.dtype('<i1'), bsize)
    f.close()
    # convert 
    eps = np.zeros(rsize)
    eps = int2bit(eps,data)
    # write converted eps field
    file_out = fname[:-1] + str(2)
    with open(file_out, "wb") as fout:
        print("Write out  file %s ..." % file_out)       
        header.astype('<i4').tofile(fout)
        eps.astype('<i1').tofile(fout)
        fout.close()
    
# types for .xdmf file
if dtype == 'real':
    prec      = 8
    seek_data = 52
    dtype     = 'Float'
elif dtype == 'int':
    prec      = 1
    seek_data = 20
    dtype     = 'Uchar'
elif dtype == 'bit':
    fname     = "eps0.2"
    prec      = 1
    seek_data = 20
    dtype     = 'Uchar'

print("Processing IBM geometry field... ", fname)

f = open(fname[:3]+'.xdmf', 'w')

# Definining entities depending on datatype

f.write('''<?xml version="1.0" ?>

<!--
XDMF file to read IBM geometry field 'eps0.1' from TLAB into a data analysis and visualization application, like ParaView.

The structure of this file has been adapted from psOpen, from Jens Henrik Goebbert.
-->

<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" [
''')

data = ( nx,ny,nz, nx,ny,nz, prec )
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
<!ENTITY Prec      "%d">
''' %data)

data = (56+nx*8+8, 56+nx*8+8+ny*8+8, seek_data)
f.write('''
<!-- offsets to grid blocks -->
<!ENTITY SeekGridX  "56"> 
<!ENTITY SeekGridY  "%d"> <!-- + DimX*8 + 8-->
<!ENTITY SeekGridZ  "%d"> <!-- + DimY*8 + 8-->

<!-- offsets to data -->
<!ENTITY SeekData   "%d"> <!-- Header -->
''' %data)

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
    
    <Grid Name="SingleFrame" GridType="Uniform">
    <Time TimeType="Single" Value=" 0.E+00"/>
    
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
''')

# file
data = (dtype,fname)
f.write('''
	<Attribute Center="Node" Name="eps">
	  <DataItem ItemType="HyperSlab" Dimensions="&HSDimsZ; &HSDimsY; &HSDimsX;">
	    <DataItem Reference="/Xdmf/Domain/DataItem[1]"/>
	    <DataItem ItemType="Uniform" Format="Binary" Seek="&SeekData;" NumberType="%s" Precision="&Prec;" Endian="Little" Dimensions="&DimsZ; &DimsY; &DimsX;">
	      %s
	    </DataItem>
	  </DataItem>
	</Attribute>	
''' %data)

f.write('''
      </Grid>
''' )

# end the xmf file
f.write('''
  </Domain>
</Xdmf>
''')

f.close()