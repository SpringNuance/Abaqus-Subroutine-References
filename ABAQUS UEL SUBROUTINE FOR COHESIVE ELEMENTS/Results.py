# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 2016 replay file
# Internal Version: 2015_09_24-22.31.09 126547
# Run by emimpa on Tue May 30 13:35:17 2017
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=279.234375, 
    height=131.760009765625)
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
#: Number of Node Sets:          2
#: Number of Steps:              1
odb = session.odbs['Job-1.odb']
xy0 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Spatial displacement: U2 at Node 1 in NSET _PICKEDSET6', 
    steps=('Step-1', ), suppressQuery=True)
xy1 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Stress components: S22 at Element 1 Int Point 1 in ELSET _PICKEDSET6', 
    steps=('Step-1', ), suppressQuery=True)
xy2 = combine(xy0, xy1, )
xy_result = session.XYData(name='XYData-1', objectToCopy=xy2, 
    sourceDescription='combine(Spatial displacement: U2 at Node 1 in NSET _PICKEDSET6, Stress components: S22 at Element 1 Int Point 1 in ELSET _PICKEDSET6, )')
del session.xyDataObjects[xy0.name]
del session.xyDataObjects[xy1.name]
del session.xyDataObjects[xy2.name]
c1 = session.Curve(xyData=xy_result)
xyp = session.XYPlot('XYPlot-1')
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
chart.setValues(curvesToPlot=(c1, ), )
session.viewports['Viewport: 1'].setValues(displayedObject=xyp)
