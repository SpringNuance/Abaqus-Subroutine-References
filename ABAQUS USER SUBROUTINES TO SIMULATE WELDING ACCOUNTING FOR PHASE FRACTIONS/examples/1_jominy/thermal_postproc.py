"""
abaqus python thermal_postproc.py filename.inp

Creates mutliple odb files where the phase fraction have the 
temperature label
"""

from odbAccess import *
import numpy
import shutil
import sys

filename = sys.argv[1]

# Open og and new databases
odb = openOdb(filename)

# Remove .odb suffix
filename = filename[:-4]

nodes = odb.rootAssembly.instances['PART-1-1'].nodes
labels = [node.label for node in nodes]
coordinates = [node.coordinates[:2] for node in nodes]
elements = odb.rootAssembly.instances['PART-1-1'].elements
connectivity_DC2D8 = [(i+1,) + element.connectivity for i, element in enumerate(elements) if element.type == 'DC2D8']
connectivity_DC2D6 = [(i+1,) + element.connectivity for i, element in enumerate(elements) if element.type == 'DC2D6']
connectivity_DC2D4 = [(i+1,) + element.connectivity for i, element in enumerate(elements) if element.type == 'DC2D4']
connectivity_DC2D3 = [(i+1,) + element.connectivity for i, element in enumerate(elements) if element.type == 'DC2D3']

steps = odb.steps

# Loop over phases
phases = ['f', 'p', 'b', 'm', 'a']

for iphase, phase in enumerate(phases):
    print "Phase", phase
    # Create prt files
    shutil.copyfile(filename+'.prt', filename+'_{:s}.prt'.format(phase))

    # Create new databases
    odb_phase = Odb(name='frac_{:s}'.format(phase),
        analysisTitle='derived data',
        description='{:s} phase fraction'.format(phase.upper()),
        path=filename+'_{:s}.odb'.format(phase))

    # Copy mesh from og odb
    part = odb_phase.Part('Part-1', TWO_D_PLANAR, DEFORMABLE_BODY)
    part.addNodes(labels, coordinates)

    if connectivity_DC2D8:
        part.addElements(connectivity_DC2D8, 'DC2D8')
    if connectivity_DC2D6:
        part.addElements(connectivity_DC2D6, 'DC2D6')
    if connectivity_DC2D4:
        part.addElements(connectivity_DC2D4, 'DC2D4')
    if connectivity_DC2D3:
        part.addElements(connectivity_DC2D3, 'DC2D3')

    inst = odb_phase.rootAssembly.Instance('PART-1-1', part)

    # Copy FV to NT11
    for stepkey in steps.keys():
        step = steps[stepkey]
        step_phase = odb_phase.Step(name=step.name, description=step.description, domain=step.domain,
            timePeriod=step.timePeriod, totalTime=step.totalTime)

        for iframe, frame in enumerate(step.frames):
            frame_phase = step_phase.Frame(incrementNumber=frame.incrementNumber, frameValue=frame.frameValue)
            field = frame_phase.FieldOutput(name='NT11', description='Nodal temperature', type=SCALAR)

            bulkdata = frame.fieldOutputs['FV{:d}'.format(iphase+1)].getSubset(position=NODAL).bulkDataBlocks[0]
            data = numpy.minimum(numpy.maximum(bulkdata.data, 0.0), 1.0)
            field.addData(position=NODAL, instance=inst, labels=bulkdata.nodeLabels, data=data)

            #datapoints = frame.fieldOutputs['FV{:d}'.format(iphase+1)].values
            # bulkdatablocks = frame.fieldOutputs['FV{:d}'.format(iphase+1)].bulkDataBlocks
           
            # for bulkdata in bulkdatablocks:
            #     if bulkdata.position == NODAL:
            #         data = numpy.minimum(numpy.maximum(bulkdata.data, 0.0), 1.0)
            #         field.addData(position=NODAL, instance=inst, labels=bulkdata.nodeLabels, data=data)

    odb_phase.save()
