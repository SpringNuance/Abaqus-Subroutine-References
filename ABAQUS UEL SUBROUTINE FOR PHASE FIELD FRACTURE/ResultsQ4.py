# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 2017 replay file
# Internal Version: 2016_09_27-22.54.59 126836
# Run by em642 on Tue Jul 31 15:37:28 2018
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
#: Warning: Permission was denied for "abaqus.rpy"; "abaqus.rpy.31" will be used for this session's replay file.
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=178.995819091797, 
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
#: Number of Element Sets:       10
#: Number of Node Sets:          8
#: Number of Steps:              1
session.viewports['Viewport: 1'].odbDisplay.basicOptions.setValues(
    pointElements=OFF)
session.viewports['Viewport: 1'].odbDisplay.commonOptions.setValues(
    visibleEdges=FEATURE, deformationScaling=UNIFORM, uniformScaleFactor=0)
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
xy26 = sum((xy0, xy1, xy2, xy3, xy4, xy5, xy6, xy7, xy8, xy9, xy10, xy11, xy12, 
    xy13, xy14, xy15, xy16, xy17, xy18, xy19, xy20, xy21, xy22, xy23, xy24, 
    xy25, ), )
session.XYData(name='RF', objectToCopy=xy26, 
    sourceDescription='sum((RF, RF, RF, RF, RF, RF, RF, RF, RF, RF, RF, RF, RF, RF, RF, RF, RF, RF, RF, RF, RF, RF, RF, RF, RF, RF, ),)')
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
odb = session.odbs['Job-1.odb']
session.XYDataFromHistory(name='U2 N: 4 NSET TOPSIDE-1', odb=odb, 
    outputVariableName='Spatial displacement: U2 at Node 4 in NSET TOPSIDE', 
    steps=('Step-1', ), __linkedVpName__='Viewport: 1')
xy1 = session.xyDataObjects['U2 N: 4 NSET TOPSIDE-1']
xy2 = session.xyDataObjects['RF']
xy3 = combine(xy1, -1*xy2)
xy3.setValues(sourceDescription='combine ( "U2 N: 4 NSET TOPSIDE-1",-1*"RF" )')
tmpName = xy3.name
session.xyDataObjects.changeKey(tmpName, 'ForceVSDisplacement')
xyp = session.XYPlot('XYPlot-1')
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
xy1 = session.xyDataObjects['ForceVSDisplacement']
c1 = session.Curve(xyData=xy1)
chart.setValues(curvesToPlot=(c1, ), )
session.viewports['Viewport: 1'].setValues(displayedObject=xyp)
odb = session.odbs['Job-1.odb']
session.viewports['Viewport: 1'].setValues(displayedObject=odb)
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_PHI', outputPosition=INTEGRATION_POINT, )
