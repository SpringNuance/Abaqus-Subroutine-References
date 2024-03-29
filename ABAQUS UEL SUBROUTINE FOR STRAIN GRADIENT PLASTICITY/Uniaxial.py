# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 2016 replay file
# Internal Version: 2015_09_24-22.31.09 126547
# Run by emimpa on Sun Nov 27 17:25:12 2016
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=304.055206298828, 
    height=237.870010375977)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from caeModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=ON)
o1 = session.openOdb(name='Uniaxial.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#: Model: D:/AbaqusTemp/Uniaxial.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       6
#: Number of Node Sets:          4
#: Number of Steps:              1
odb = session.odbs['Uniaxial.odb']
session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('RF', 
    NODAL, ((COMPONENT, 'RF2'), )), ), nodePick=(('PART-1-1', 41, (
    '[#1fffff #0:12 #22000000 #49249249 #92492 ]', )), ), )
odb = session.odbs['Uniaxial.odb']
session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('U', 
    NODAL, ((COMPONENT, 'U2'), )), ), nodePick=(('PART-1-1', 1, (
    '[#0:13 #10 ]', )), ), )
xy1 = session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 1']
xy2 = session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 2']
xy3 = session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 3']
xy4 = session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 4']
xy5 = session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 5']
xy6 = session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 6']
xy7 = session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 7']
xy8 = session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 8']
xy9 = session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 9']
xy10 = session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 10']
xy11 = session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 11']
xy12 = session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 12']
xy13 = session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 13']
xy14 = session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 14']
xy15 = session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 15']
xy16 = session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 16']
xy17 = session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 17']
xy18 = session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 18']
xy19 = session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 19']
xy20 = session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 20']
xy21 = session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 21']
xy22 = session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 442']
xy23 = session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 446']
xy24 = session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 449']
xy25 = session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 452']
xy26 = session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 455']
xy27 = session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 458']
xy28 = session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 461']
xy29 = session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 464']
xy30 = session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 467']
xy31 = session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 470']
xy32 = session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 473']
xy33 = session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 476']
xy34 = session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 479']
xy35 = session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 482']
xy36 = session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 485']
xy37 = session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 488']
xy38 = session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 491']
xy39 = session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 494']
xy40 = session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 497']
xy41 = session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 500']
xy42 = sum((xy1, xy2, xy3, xy4, xy5, xy6, xy7, xy8, xy9, xy10, xy11, xy12, 
    xy13, xy14, xy15, xy16, xy17, xy18, xy19, xy20, xy21, xy22, xy23, xy24, 
    xy25, xy26, xy27, xy28, xy29, xy30, xy31, xy32, xy33, xy34, xy35, xy36, 
    xy37, xy38, xy39, xy40, xy41))
xy42.setValues(
    sourceDescription='sum ( ( "RF:RF2 PI: PART-1-1 N: 1", "RF:RF2 PI: PART-1-1 N: 2", "RF:RF2 PI: PART-1-1 N: 3", "RF:RF2 PI: PART-1-1 N: 4", "RF:RF2 PI: PART-1-1 N: 5", "RF:RF2 PI: PART-1-1 N: 6", "RF:RF2 PI: PART-1-1 N: 7", "RF:RF2 PI: PART-1-1 N: 8", "RF:RF2 PI: PART-1-1 N: 9", "RF:RF2 PI: PART-1-1 N: 10", "RF:RF2 PI: PART-1-1 N: 11", "RF:RF2 PI: PART-1-1 N: 12", "RF:RF2 PI: PART-1-1 N: 13", "RF:RF2 PI: PART-1-1 N: 14", "RF:RF2 PI: PART-1-1 N: 15", "RF:RF2 PI: PART-1-1 N: 16", "RF:RF2 PI: PART-1-1 N: 17", "RF:RF2 PI: PART-1-1 N: 18", "RF:RF2 PI: PART-1-1 N: 19", "RF:RF2 PI: PART-1-1 N: 20", "RF:RF2 PI: PART-1-1 N: 21", "RF:RF2 PI: PART-1-1 N: 442", "RF:RF2 PI: PART-1-1 N: 446", "RF:RF2 PI: PART-1-1 N: 449", "RF:RF2 PI: PART-1-1 N: 452", "RF:RF2 PI: PART-1-1 N: 455", "RF:RF2 PI: PART-1-1 N: 458", "RF:RF2 PI: PART-1-1 N: 461", "RF:RF2 PI: PART-1-1 N: 464", "RF:RF2 PI: PART-1-1 N: 467", "RF:RF2 PI: PART-1-1 N: 470", "RF:RF2 PI: PART-1-1 N: 473", "RF:RF2 PI: PART-1-1 N: 476", "RF:RF2 PI: PART-1-1 N: 479", "RF:RF2 PI: PART-1-1 N: 482", "RF:RF2 PI: PART-1-1 N: 485", "RF:RF2 PI: PART-1-1 N: 488", "RF:RF2 PI: PART-1-1 N: 491", "RF:RF2 PI: PART-1-1 N: 494", "RF:RF2 PI: PART-1-1 N: 497", "RF:RF2 PI: PART-1-1 N: 500" ) )')
tmpName = xy42.name
session.xyDataObjects.changeKey(tmpName, 'RF')
xy1 = session.xyDataObjects['U:U2 PI: PART-1-1 N: 421']
xy2 = session.xyDataObjects['RF']
xy3 = combine(xy1, xy2/-0.001)
xy3.setValues(
    sourceDescription='combine ( "U:U2 PI: PART-1-1 N: 421","RF"/-0.001 )')
tmpName = xy3.name
session.xyDataObjects.changeKey(tmpName, 'Curve')
del session.xyDataObjects['RF']
del session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 1']
del session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 2']
del session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 3']
del session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 4']
del session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 5']
del session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 6']
del session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 7']
del session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 8']
del session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 9']
del session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 10']
del session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 11']
del session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 12']
del session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 13']
del session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 14']
del session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 15']
del session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 16']
del session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 17']
del session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 18']
del session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 19']
del session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 20']
del session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 21']
del session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 442']
del session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 446']
del session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 449']
del session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 452']
del session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 455']
del session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 458']
del session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 461']
del session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 464']
del session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 467']
del session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 470']
del session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 473']
del session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 476']
del session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 479']
del session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 482']
del session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 485']
del session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 488']
del session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 491']
del session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 494']
del session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 497']
del session.xyDataObjects['RF:RF2 PI: PART-1-1 N: 500']
del session.xyDataObjects['U:U2 PI: PART-1-1 N: 421']