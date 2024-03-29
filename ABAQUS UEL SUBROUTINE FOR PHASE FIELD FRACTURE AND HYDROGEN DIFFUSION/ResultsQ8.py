# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 2017 replay file
# Internal Version: 2016_09_27-22.54.59 126836
# Run by em642 on Wed Aug 08 17:18:39 2018
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=179.243743896484, 
    height=102.350006103516)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from caeModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=ON)
o1 = session.openOdb(name='Job-1.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#: Model: Job-1.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       4
#: Number of Node Sets:          9
#: Number of Steps:              1
session.viewports['Viewport: 1'].odbDisplay.basicOptions.setValues(
    pointElements=OFF)
session.viewports['Viewport: 1'].odbDisplay.commonOptions.setValues(
    visibleEdges=FEATURE, deformationScaling=UNIFORM, uniformScaleFactor=1)
odb = session.odbs['Job-1.odb']
xy0 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Reaction force: RF2 at Node 6 in NSET BOTTOMSIDE', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
xy1 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Reaction force: RF2 at Node 7 in NSET BOTTOMSIDE', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
xy2 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Reaction force: RF2 at Node 211 in NSET BOTTOMSIDE', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
xy3 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Reaction force: RF2 at Node 212 in NSET BOTTOMSIDE', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
xy4 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Reaction force: RF2 at Node 213 in NSET BOTTOMSIDE', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
xy5 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Reaction force: RF2 at Node 214 in NSET BOTTOMSIDE', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
xy6 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Reaction force: RF2 at Node 215 in NSET BOTTOMSIDE', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
xy7 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Reaction force: RF2 at Node 216 in NSET BOTTOMSIDE', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
xy8 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Reaction force: RF2 at Node 217 in NSET BOTTOMSIDE', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
xy9 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Reaction force: RF2 at Node 218 in NSET BOTTOMSIDE', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
xy10 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Reaction force: RF2 at Node 219 in NSET BOTTOMSIDE', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
xy11 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Reaction force: RF2 at Node 220 in NSET BOTTOMSIDE', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
xy12 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Reaction force: RF2 at Node 221 in NSET BOTTOMSIDE', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
xy13 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Reaction force: RF2 at Node 222 in NSET BOTTOMSIDE', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
xy14 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Reaction force: RF2 at Node 223 in NSET BOTTOMSIDE', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
xy15 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Reaction force: RF2 at Node 224 in NSET BOTTOMSIDE', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
xy16 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Reaction force: RF2 at Node 225 in NSET BOTTOMSIDE', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
xy17 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Reaction force: RF2 at Node 226 in NSET BOTTOMSIDE', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
xy18 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Reaction force: RF2 at Node 227 in NSET BOTTOMSIDE', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
xy19 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Reaction force: RF2 at Node 228 in NSET BOTTOMSIDE', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
xy20 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Reaction force: RF2 at Node 229 in NSET BOTTOMSIDE', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
xy21 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Reaction force: RF2 at Node 230 in NSET BOTTOMSIDE', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
xy22 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Reaction force: RF2 at Node 231 in NSET BOTTOMSIDE', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
xy23 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Reaction force: RF2 at Node 232 in NSET BOTTOMSIDE', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
xy24 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Reaction force: RF2 at Node 233 in NSET BOTTOMSIDE', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
xy25 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Reaction force: RF2 at Node 234 in NSET BOTTOMSIDE', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
xy26 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Reaction force: RF2 at Node 11153 in NSET BOTTOMSIDE', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
xy27 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Reaction force: RF2 at Node 11158 in NSET BOTTOMSIDE', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
xy28 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Reaction force: RF2 at Node 11160 in NSET BOTTOMSIDE', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
xy29 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Reaction force: RF2 at Node 11224 in NSET BOTTOMSIDE', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
xy30 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Reaction force: RF2 at Node 11228 in NSET BOTTOMSIDE', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
xy31 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Reaction force: RF2 at Node 11232 in NSET BOTTOMSIDE', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
xy32 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Reaction force: RF2 at Node 11234 in NSET BOTTOMSIDE', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
xy33 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Reaction force: RF2 at Node 11402 in NSET BOTTOMSIDE', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
xy34 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Reaction force: RF2 at Node 11404 in NSET BOTTOMSIDE', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
xy35 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Reaction force: RF2 at Node 11426 in NSET BOTTOMSIDE', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
xy36 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Reaction force: RF2 at Node 11450 in NSET BOTTOMSIDE', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
xy37 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Reaction force: RF2 at Node 11452 in NSET BOTTOMSIDE', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
xy38 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Reaction force: RF2 at Node 11487 in NSET BOTTOMSIDE', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
xy39 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Reaction force: RF2 at Node 11490 in NSET BOTTOMSIDE', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
xy40 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Reaction force: RF2 at Node 11503 in NSET BOTTOMSIDE', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
xy41 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Reaction force: RF2 at Node 11505 in NSET BOTTOMSIDE', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
xy42 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Reaction force: RF2 at Node 11518 in NSET BOTTOMSIDE', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
xy43 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Reaction force: RF2 at Node 11524 in NSET BOTTOMSIDE', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
xy44 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Reaction force: RF2 at Node 11594 in NSET BOTTOMSIDE', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
xy45 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Reaction force: RF2 at Node 11598 in NSET BOTTOMSIDE', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
xy46 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Reaction force: RF2 at Node 11633 in NSET BOTTOMSIDE', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
xy47 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Reaction force: RF2 at Node 11671 in NSET BOTTOMSIDE', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
xy48 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Reaction force: RF2 at Node 11675 in NSET BOTTOMSIDE', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
xy49 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Reaction force: RF2 at Node 11711 in NSET BOTTOMSIDE', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
xy50 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Reaction force: RF2 at Node 16400 in NSET BOTTOMSIDE', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
xy51 = sum((xy0, xy1, xy2, xy3, xy4, xy5, xy6, xy7, xy8, xy9, xy10, xy11, xy12, 
    xy13, xy14, xy15, xy16, xy17, xy18, xy19, xy20, xy21, xy22, xy23, xy24, 
    xy25, xy26, xy27, xy28, xy29, xy30, xy31, xy32, xy33, xy34, xy35, xy36, 
    xy37, xy38, xy39, xy40, xy41, xy42, xy43, xy44, xy45, xy46, xy47, xy48, 
    xy49, xy50, ), )
xy_result = session.XYData(name='RForces', objectToCopy=xy51, 
    sourceDescription='sum((RForces, RForces, RForces, RForces, RForces, RForces, RForces, RForces, RForces, RForces, RForces, RForces, RForces, RForces, RForces, RForces, RForces, RForces, RForces, RForces, RForces, RForces, RForces, RForces, RForces, RForces, RForces, RForces, RForces, RForces, RForces, RForces, RForces, RForces, RForces, RForces, RForces, RForces, RForces, RForces, RForces, RForces, RForces, RForces, RForces, RForces, RForces, RForces, RForces, RForces, RForces, ),)')
del session.xyDataObjects[xy0.name]
del session.xyDataObjects[xy1.name]
del session.xyDataObjects[xy2.name]
del session.xyDataObjects[xy3.name]
del session.xyDataObjects[xy4.name]
del session.xyDataObjects[xy5.name]
del session.xyDataObjects[xy6.name]
del session.xyDataObjects[xy7.name]
del session.xyDataObjects[xy8.name]
del session.xyDataObjects[xy9.name]
del session.xyDataObjects[xy10.name]
del session.xyDataObjects[xy11.name]
del session.xyDataObjects[xy12.name]
del session.xyDataObjects[xy13.name]
del session.xyDataObjects[xy14.name]
del session.xyDataObjects[xy15.name]
del session.xyDataObjects[xy16.name]
del session.xyDataObjects[xy17.name]
del session.xyDataObjects[xy18.name]
del session.xyDataObjects[xy19.name]
del session.xyDataObjects[xy20.name]
del session.xyDataObjects[xy21.name]
del session.xyDataObjects[xy22.name]
del session.xyDataObjects[xy23.name]
del session.xyDataObjects[xy24.name]
del session.xyDataObjects[xy25.name]
del session.xyDataObjects[xy26.name]
del session.xyDataObjects[xy27.name]
del session.xyDataObjects[xy28.name]
del session.xyDataObjects[xy29.name]
del session.xyDataObjects[xy30.name]
del session.xyDataObjects[xy31.name]
del session.xyDataObjects[xy32.name]
del session.xyDataObjects[xy33.name]
del session.xyDataObjects[xy34.name]
del session.xyDataObjects[xy35.name]
del session.xyDataObjects[xy36.name]
del session.xyDataObjects[xy37.name]
del session.xyDataObjects[xy38.name]
del session.xyDataObjects[xy39.name]
del session.xyDataObjects[xy40.name]
del session.xyDataObjects[xy41.name]
del session.xyDataObjects[xy42.name]
del session.xyDataObjects[xy43.name]
del session.xyDataObjects[xy44.name]
del session.xyDataObjects[xy45.name]
del session.xyDataObjects[xy46.name]
del session.xyDataObjects[xy47.name]
del session.xyDataObjects[xy48.name]
del session.xyDataObjects[xy49.name]
del session.xyDataObjects[xy50.name]
del session.xyDataObjects[xy51.name]
c1 = session.Curve(xyData=xy_result)
xyp = session.XYPlot('XYPlot-1')
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
chart.setValues(curvesToPlot=(c1, ), )
session.viewports['Viewport: 1'].setValues(displayedObject=xyp)
odb = session.odbs['Job-1.odb']
xy_result = session.XYDataFromHistory(name='U2', odb=odb, 
    outputVariableName='Spatial displacement: U2 at Node 4 in NSET TOPSIDE', 
    steps=('Step-1', ), __linkedVpName__='Viewport: 1')
c1 = session.Curve(xyData=xy_result)
xyp = session.xyPlots['XYPlot-1']
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
chart.setValues(curvesToPlot=(c1, ), )
xy1 = session.xyDataObjects['U2']
xy2 = session.xyDataObjects['RForces']
xy3 = combine(xy1, -1*xy2)
xyp = session.xyPlots['XYPlot-1']
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
c1 = session.Curve(xyData=xy3)
chart.setValues(curvesToPlot=(c1, ), )
xy1 = session.xyDataObjects['U2']
xy2 = session.xyDataObjects['RForces']
xy3 = combine(xy1, -1*xy2)
xy3.setValues(sourceDescription='combine ( "U2",-1*"RForces" )')
tmpName = xy3.name
session.xyDataObjects.changeKey(tmpName, 'FvsDisp')
xyp = session.xyPlots['XYPlot-1']
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
xy1 = session.xyDataObjects['FvsDisp']
c1 = session.Curve(xyData=xy1)
chart.setValues(curvesToPlot=(c1, ), )
odb = session.odbs['Job-1.odb']
session.viewports['Viewport: 1'].setValues(displayedObject=odb)
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_PHI', outputPosition=INTEGRATION_POINT, )
