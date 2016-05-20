import numpy   as np   # For array operations.

filenames = ["planesJ.5", "planesJ.10"]

nx = 1 # number of points in Ox
ny = 1 # number of points in Oy
nz = 1 # number of points in OZ
ns = 1 # number of scalar variables
nq = 3 # number of flow variables

f = open('TESTplanesJ.times.xdmf', 'w')

# Header for xml file
f.write('''<?xml version="1.0" ?>

<!--
XDMF file to read time collections of plane data from TLAB into a data analysis and visualization application, like ParaView.

Add simply data items after the geometry block to read different files.

The structure of this file has been adapted from psOpen, from Jens Henrik Goebbert.
-->

<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" [

<!-- number of timeslices -->
<!ENTITY ItMax "%d"> 

<!-- offsets to grid blocks -->
<!ENTITY SeekGridX   "56"> 
<!ENTITY SeekGridY "1088"> <!-- SeekGridX + DimX*8 + 8-->
<!ENTITY SeekGridZ "1864"> <!-- SeekGridY + DimY*8 + 8-->

<!-- offsets to data -->
<!ENTITY SeekDataU   "0">
<!ENTITY SeekDataV   "512">  <!-- SeekGridU  + DimX*Prec-->
<!ENTITY SeekDataW   "1024"> <!-- SeekGridV  + DimX*Prec-->
<!ENTITY SeekDataS1  "1536"> <!-- SeekGridW  + DimX*Prec-->
<!ENTITY SeekDataS2  "2048"> <!-- SeekGridS1 + DimX*Prec-->

<!-- dimension of complete datasets -->
<!ENTITY DimsX   "128">
<!ENTITY DimsY   "1">
<!ENTITY DimsZ   "512"><!-- #Vars*DimZ-->

<!ENTITY GridDimsX   "128">
<!ENTITY GridDimsY   "1">
<!ENTITY GridDimsZ   "128">

<!-- dimension of hyperslab to load -->
<!ENTITY HSDimsX  "128">
<!ENTITY HSDimsY  "1">
<!ENTITY HSDimsZ  "128">

<!-- start of hyperslab in complete dataset -->
<!ENTITY HSDimsX_Start "0">
<!ENTITY HSDimsY_Start "0">
<!ENTITY HSDimsZ_Start "0">

<!ENTITY HSGridDimsX_Start "0">
<!ENTITY HSGridDimsY_Start "11">
<!ENTITY HSGridDimsZ_Start "0">

<!-- stride of hyperslab in complete dataset -->
<!ENTITY HSStrideX "1">
<!ENTITY HSStrideY "1">
<!ENTITY HSStrideZ "4"><!-- #Vars-->

<!ENTITY HSGridStrideX "1">
<!ENTITY HSGridStrideY "1">
<!ENTITY HSGridStrideZ "1">

]>

<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">
  <Domain Name="PlanesJ">
    
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
	  0.0 1.0  &ItMax;
	</DataItem>
      </Time>
''' % (len(filenames)) )

# Loop over timeslices
for file in filenames:
    f.write('''
      <!-- Timeslice -->
      <Grid Name="It5" GridType="Uniform">
	<Topology Reference="/Xdmf/Domain/Topology[1]"/>
	<Geometry Reference="/Xdmf/Domain/Geometry[1]"/>	
''')
    
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
