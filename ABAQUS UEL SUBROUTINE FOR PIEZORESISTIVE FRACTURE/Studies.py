# Import working directory
import os


# Import modules and libraries
from abaqus import *
from abaqusConstants import *
from caeModules import *
import numpy as np
import sys
import fileinput
import math
from abaqus import getInput




def case_study(
    Vimp,
    Dimp,
    Inc,
    Eyoung,
    Nu,
    Gc,
    theta,
    cond,
    pi11,
    pi12,
    pi44,
    k_value,
    n_value,
    file_name,
    case_study_number,
    meshsize,
    meshsize2
):
    if case_study_number == 1:
        name_of_file = "Case1"
        Mdb()
        mdb.ModelFromInputFile(name=name_of_file, inputFileName=name_of_file + ".inp")

        mdb.JobFromInputFile(
            name=name_of_file,
            inputFileName=name_of_file + ".inp",
            type=ANALYSIS,
            atTime=None,
            waitMinutes=0,
            waitHours=0,
            queue=None,
            memory=90,
            memoryUnits=PERCENTAGE,
            getMemoryFromAnalysis=True,
            explicitPrecision=SINGLE,
            nodalOutputPrecision=SINGLE,
            userSubroutine="UEL_piezoresistive_phasefield.for",
            scratch="",
            resultsFormat=ODB,
            multiprocessingMode=DEFAULT,
            numCpus=8,
            numDomains=8,
            numGPUs=0,
        )

        jobname = mdb.jobs[name_of_file]
        jobname.submit()
        jobname.waitForCompletion()


    elif case_study_number == 2:
        return case_study_2(
            Vimp,
            Dimp,
            Inc,
            Eyoung,
            Nu,
            Gc,
            theta,
            cond,
            pi11,
            pi12,
            pi44,
            k_value,
            n_value,
            file_name,
            meshsize,
            meshsize2
        )
    elif case_study_number == 3:
        return case_study_3(
            Vimp,
            Dimp,
            Inc,
            Eyoung,
            Nu,
            Gc,
            theta,
            cond,
            pi11,
            pi12,
            pi44,
            k_value,
            n_value,
            file_name,
            meshsize,
            meshsize2
        )
    elif case_study_number == 4:
        return case_study_4(
            Vimp,
            Dimp,
            Inc,
            Eyoung,
            Nu,
            Gc,
            theta,
            cond,
            pi11,
            pi12,
            pi44,
            k_value,
            n_value,
            file_name,
            meshsize,
            meshsize2
        )
    elif case_study_number == 5:
        name_of_file = "Case5"
        Mdb()
        mdb.ModelFromInputFile(name=name_of_file, inputFileName=name_of_file + ".inp")

        mdb.JobFromInputFile(
            name=name_of_file,
            inputFileName=name_of_file + ".inp",
            type=ANALYSIS,
            atTime=None,
            waitMinutes=0,
            waitHours=0,
            queue=None,
            memory=90,
            memoryUnits=PERCENTAGE,
            getMemoryFromAnalysis=True,
            explicitPrecision=SINGLE,
            nodalOutputPrecision=SINGLE,
            userSubroutine="UEL_piezoresistive_phasefield.for",
            scratch="",
            resultsFormat=ODB,
            multiprocessingMode=DEFAULT,
            numCpus=8,
            numDomains=8,
            numGPUs=0,
        )

        jobname = mdb.jobs[name_of_file]
        jobname.submit()
        jobname.waitForCompletion()    


def case_study_2(
    Vimp,
    Dimp,
    Inc,
    Eyoung,
    Nu,
    Gc,
    theta,
    cond,
    pi11,
    pi12,
    pi44,
    k_value,
    n_value,
    file_name,
    meshsize,
    meshsize2
):
    name_of_file = file_name
    # Parameters
    ######################################################################
    lx = 0.1  # [m]
    ly = 0.2  # [m]
    lz = 0.005  # [m]
    delect = 0.01  # [m] Distance from the edge
    lelect = 0.08  # [m] Length of the electrode
    welect = 0.01  # [m] Width of the electrode
    lcrack = 0.02
    dcrack = 0.001
    theta = theta * (np.pi / 180)

    # Phase-field
    ######################################################################
    # Model
    modelSK = mdb.models["Model-1"]
    # Strip
    s = modelSK.ConstrainedSketch(name="__profile__", sheetSize=200.0)
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=STANDALONE)
    s.rectangle(point1=(0.0, 0.0), point2=(lx, ly))
    s.EllipseByCenterPerimeter(
        center=(lx / 2.0, ly / 2.0),
        axisPoint1=(
            lx / 2.0 + lcrack * np.cos(theta),
            ly / 2.0 + lcrack * np.sin(theta),
        ),
        axisPoint2=(
            lx / 2.0 + dcrack * np.cos(theta + np.pi / 2.0),
            ly / 2.0 + dcrack * np.sin(theta + np.pi / 2.0),
        ),
    )
    p = modelSK.Part(name="Part-1", dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p = modelSK.parts["Part-1"]
    p.BaseSolidExtrude(sketch=s, depth=lz)
    s.unsetPrimaryObject()
    p = modelSK.parts["Part-1"]

    # Cut the strip to allocate the electrodes
    p = modelSK.parts["Part-1"]
    datumid = p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=delect)
    datumid2 = p.DatumPlaneByPrincipalPlane(
        principalPlane=XZPLANE, offset=delect + welect
    )
    datumid3 = p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=ly - delect)
    datumid4 = p.DatumPlaneByPrincipalPlane(
        principalPlane=XZPLANE, offset=ly - delect - welect
    )
    datumid5 = p.DatumPlaneByPrincipalPlane(
        principalPlane=YZPLANE, offset=(lx - lelect) / 2.0
    )
    datumid6 = p.DatumPlaneByPrincipalPlane(
        principalPlane=YZPLANE, offset=lx - (lx - lelect) / 2.0
    )

    datumid7 = p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=ly * 0.4)
    datumid8 = p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=ly * 0.6)
    datumid9 = p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=ly * 0.7)
    datumid10 = p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=ly * 0.3)

    p = modelSK.parts["Part-1"]
    c = p.cells
    pickedCells = c
    d2 = p.datums
    p.PartitionCellByDatumPlane(datumPlane=p.datums[datumid.id], cells=pickedCells)

    p = modelSK.parts["Part-1"]
    c = p.cells
    pickedCells = c
    d2 = p.datums
    p.PartitionCellByDatumPlane(datumPlane=p.datums[datumid2.id], cells=pickedCells)

    p = modelSK.parts["Part-1"]
    c = p.cells
    pickedCells = c
    d2 = p.datums
    p.PartitionCellByDatumPlane(datumPlane=p.datums[datumid3.id], cells=pickedCells)

    p = modelSK.parts["Part-1"]
    c = p.cells
    pickedCells = c
    d2 = p.datums
    p.PartitionCellByDatumPlane(datumPlane=p.datums[datumid4.id], cells=pickedCells)

    p = modelSK.parts["Part-1"]
    c = p.cells
    pickedCells = c
    d2 = p.datums
    p.PartitionCellByDatumPlane(datumPlane=p.datums[datumid5.id], cells=pickedCells)

    p = modelSK.parts["Part-1"]
    c = p.cells
    pickedCells = c
    d2 = p.datums
    p.PartitionCellByDatumPlane(datumPlane=p.datums[datumid6.id], cells=pickedCells)

    p = modelSK.parts["Part-1"]
    c = p.cells
    pickedCells = c
    d2 = p.datums
    p.PartitionCellByDatumPlane(datumPlane=p.datums[datumid7.id], cells=pickedCells)

    p = modelSK.parts["Part-1"]
    c = p.cells
    pickedCells = c
    d2 = p.datums
    p.PartitionCellByDatumPlane(datumPlane=p.datums[datumid8.id], cells=pickedCells)

    p = modelSK.parts["Part-1"]
    c = p.cells
    pickedCells = c
    d2 = p.datums
    p.PartitionCellByDatumPlane(datumPlane=p.datums[datumid9.id], cells=pickedCells)

    p = modelSK.parts["Part-1"]
    c = p.cells
    pickedCells = c
    d2 = p.datums
    p.PartitionCellByDatumPlane(datumPlane=p.datums[datumid10.id], cells=pickedCells)

    # Master Node
    p = modelSK.Part(name="Part-2", dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p.ReferencePoint(point=(lx / 2, delect + welect / 2.0, -5.0 / 10000))

    p = modelSK.Part(name="Part-3", dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p.ReferencePoint(point=(lx / 2, ly - delect - welect / 2.0, -5.0 / 10000))

    # Material
    modelSK.Material(name="Material-1")
    mdb.models["Model-1"].materials["Material-1"].Depvar(n=18)
    mdb.models["Model-1"].materials["Material-1"].UserMaterial(
        mechanicalConstants=(0.0,)
    )
    modelSK.HomogeneousSolidSection(
        name="Section-1", material="Material-1", thickness=None
    )
    p = modelSK.parts["Part-1"]
    c = p.cells
    cells = c
    region = p.Set(cells=cells, name="Set-1")
    p = modelSK.parts["Part-1"]
    p.SectionAssignment(
        region=region,
        sectionName="Section-1",
        offset=0.0,
        offsetType=MIDDLE_SURFACE,
        offsetField="",
        thicknessAssignment=FROM_SECTION,
    )

    # Instance
    a1 = modelSK.rootAssembly
    p = modelSK.parts["Part-1"]
    a1.Instance(name="Part-1-1", part=p, dependent=OFF)
    p = modelSK.parts["Part-2"]
    a1.Instance(name="Part-2-1", part=p, dependent=OFF)
    p = modelSK.parts["Part-3"]
    a1.Instance(name="Part-3-1", part=p, dependent=OFF)

    a = modelSK.rootAssembly
    e1 = a.instances["Part-1-1"].edges
    edges1 = e1.getSequenceFromMask(
        mask=("[#0 #2d680 #800000 ]",),
    )
    a.Set(edges=edges1, name="Middle-1")
    e1 = a.instances["Part-1-1"].edges
    edges1 = e1.getSequenceFromMask(
        mask=("[#50a1a000 #a0000040 #104028a0 ]",),
    )
    a.Set(edges=edges1, name="Middle-2")

    # Mesh
    # ---- Strip
    elemType1 = mesh.ElemType(
        elemCode=C3D8,
        elemLibrary=STANDARD,
        secondOrderAccuracy=OFF,
        distortionControl=DEFAULT,
    )
    elemType2 = mesh.ElemType(elemCode=C3D6, elemLibrary=STANDARD)
    elemType3 = mesh.ElemType(elemCode=C3D4, elemLibrary=STANDARD)
    a = modelSK.rootAssembly
    c1 = a.instances["Part-1-1"].cells
    cells = c1
    pickedRegions = (cells,)
    a.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2, elemType3))

    e1 = a.instances["Part-1-1"].edges
    pickedEdges = e1.getSequenceFromMask(
        mask=("[#0 #2d680 #800000 ]",),
    )
    a.seedEdgeBySize(
        edges=pickedEdges,
        size=meshsize,
        deviationFactor=0.1,
        minSizeFactor=0.1,
        constraint=FINER,
    )
    e1 = a.instances["Part-1-1"].edges
    pickedEdges1 = e1.getSequenceFromMask(
        mask=("[#40a10000 #a0000000 #400800 ]",),
    )
    pickedEdges2 = e1.getSequenceFromMask(
        mask=("[#1000a000 #40 #100020a0 ]",),
    )
    a.seedEdgeByBias(
        biasMethod=SINGLE,
        end1Edges=pickedEdges1,
        end2Edges=pickedEdges2,
        minSize=meshsize,
        maxSize=meshsize,
        constraint=FINER,
    )

    a = modelSK.rootAssembly
    partInstances = (
        a.instances["Part-1-1"],
        a.instances["Part-2-1"],
        a.instances["Part-3-1"],
    )
    a.seedPartInstance(
        regions=partInstances, size=meshsize, deviationFactor=0.1, minSizeFactor=0.1
    )
    partInstances = (a.instances["Part-1-1"],)
    a.generateMesh(regions=partInstances)

    # Create Node Sets of electrodes
    a = modelSK.rootAssembly
    n1 = a.instances["Part-1-1"].nodes
    faces1 = n1.getByBoundingBox(-1000, -1000, -1000, 1000, 1000, 1000)
    a.Set(nodes=faces1, name="ALL")

    a = modelSK.rootAssembly
    n1 = a.instances["Part-1-1"].nodes
    faces1 = n1.getByBoundingBox(0, 0, 0, lx, ly, lz)
    a.Set(nodes=faces1, name="body")

    a = modelSK.rootAssembly
    n1 = a.instances["Part-1-1"].elements
    faces1 = n1.getByBoundingBox(-1000, -1000, -1000, 1000, 1000, 1000)
    a.Set(nodes=faces1, name="ALLEL")

    a = modelSK.rootAssembly
    n1 = a.instances["Part-1-1"].nodes
    faces1 = n1.getByBoundingBox(-100, -100, -100, 100, 0.0001, 100)
    a.Set(nodes=faces1, name="Base")

    a = modelSK.rootAssembly
    n1 = a.instances["Part-1-1"].nodes
    faces1 = n1.getByBoundingBox(-100, ly, -100, 100, ly + 0.0001, 100)
    a.Set(nodes=faces1, name="Top")

    a = modelSK.rootAssembly
    n1 = a.instances["Part-1-1"].nodes
    faces1 = n1.getByBoundingBox(
        delect - 0.0001,
        delect - 0.0001,
        0.004,
        lelect + delect + 0.0001,
        delect + welect + 0.0001,
        0.005,
    )
    a.Set(nodes=faces1, name="Electrode2")

    a = modelSK.rootAssembly
    n1 = a.instances["Part-1-1"].nodes
    faces1 = n1.getByBoundingBox(
        delect - 0.0001,
        ly - delect - welect - 0.0001,
        0.004,
        lelect + delect + 0.0001,
        ly - delect + 0.0001,
        0.005,
    )
    a.Set(nodes=faces1, name="Electrode1")

    # Simulation
    modelSK.StaticStep(
        name="Step-1",
        previous="Initial",
        timeIncrementationMethod=FIXED,
        initialInc=Inc,
        noStop=OFF,
        reformKernel=25,
        solutionTechnique=QUASI_NEWTON,
    )
    modelSK.steps["Step-1"].control.setValues(
        allowPropagation=OFF,
        resetDefaultValues=OFF,
        timeIncrementation=(
            25000.0,
            25000.0,
            25000.0,
            25000.0,
            25000.0,
            4.0,
            12.0,
            5.0,
            6.0,
            3.0,
            50.0,
        ),
        lineSearch=(5.0, 1.0, 0.0001, 0.25, 0.1),
    )
    # Amplitudes
    mdb.models["Model-1"].TabularAmplitude(
        data=((0.0, 0.9999999), (1.0, 1.0)),
        name="Const",
        smooth=SOLVER_DEFAULT,
        timeSpan=STEP,
    )
    mdb.models["Model-1"].TabularAmplitude(
        data=((0.0, 0.0), (1.0, 1.0)),
        name="Linear",
        smooth=SOLVER_DEFAULT,
        timeSpan=STEP,
    )

    # BCs
    # Electrical
    modelSK.ElectricPotentialBC(
        amplitude=UNSET,
        createStepName="Step-1",
        distributionType=UNIFORM,
        fieldName="Const",
        fixed=OFF,
        magnitude=0.0,
        name="BC-3",
        region=mdb.models["Model-1"].rootAssembly.sets["Electrode1"],
    )

    modelSK.ElectricPotentialBC(
        amplitude="Const",
        createStepName="Step-1",
        distributionType=UNIFORM,
        fieldName="",
        fixed=OFF,
        magnitude=Vimp,
        name="BC-4",
        region=mdb.models["Model-1"].rootAssembly.sets["Electrode2"],
    )

    # Mechanical
    a = modelSK.rootAssembly
    region = a.sets["Base"]
    modelSK.DisplacementBC(
        name="BC-1",
        createStepName="Initial",
        region=region,
        u1=0.0,
        u2=0.0,
        u3=0.0,
        ur1=UNSET,
        ur2=UNSET,
        ur3=UNSET,
        amplitude=UNSET,
        fixed=OFF,
        distributionType=UNIFORM,
        fieldName="",
        localCsys=None,
    )
    a = modelSK.rootAssembly
    region = a.sets["Top"]
    modelSK.DisplacementBC(
        name="BC-2",
        createStepName="Step-1",
        region=region,
        u1=0.0,
        u2=Dimp,
        u3=0.0,
        ur1=UNSET,
        ur2=UNSET,
        ur3=UNSET,
        amplitude="Linear",
        fixed=OFF,
        distributionType=UNIFORM,
        fieldName="",
        localCsys=None,
    )

    # Create job
    # Desired output information
    mdb.models["Model-1"].FieldOutputRequest(
        name="F-Output-1",
        createStepName="Step-1",
    )
    modelSK.fieldOutputRequests["F-Output-1"].setValues(
        variables=("S", "LE", "U", "RF", "CF", "E", "SDV")
    )
    # Submit job
    mdb.Job(
        name="PF",
        model="Model-1",
        description="",
        type=ANALYSIS,
        atTime=None,
        waitMinutes=0,
        waitHours=0,
        queue=None,
        memory=90,
        memoryUnits=PERCENTAGE,
        getMemoryFromAnalysis=True,
        explicitPrecision=SINGLE,
        nodalOutputPrecision=SINGLE,
        echoPrint=OFF,
        modelPrint=OFF,
        contactPrint=OFF,
        historyPrint=OFF,
        userSubroutine="",
        scratch="",
        multiprocessingMode=DEFAULT,
        numCpus=1,
        numGPUs=0,
    )
    mdb.jobs["PF"].writeInput(consistencyChecking=OFF)

    ## Prepare imput
    ######################################################################
    ######################################################################
    ######################################################################
    ######################################################################
    ######################################################################
    ######################################################################
    # --------------------------------------------------------------
    Beg = "*Heading"  # Beginning of the input file
    End = "*End Step"  # End of the input file.

    Eletype = float(2)

    if Eletype == 1.0:
        txtStart = "*Element, type=C3D10"  # Text before the element connectivity
    elif Eletype == 2.0:
        txtStart = "*Element, type=C3D8"  # Text before the element connectivity
    elif Eletype == 3.0:
        txtStart = "*Element, type=C3D20R"  # Text before the element connectivity
    elif Eletype == 4.0:
        txtStart = "*Element, type=CPE8R"  # Text before the element connectivity
    elif Eletype == 5.0:
        txtStart = "*Element, type=CPE4"  # Text before the element connectivity

    txtEnd = "*Nset, nset=Set-1, generate"  # Text after the element connectivity
    # --------------------------------------------------------------

    INP = "PF.inp"
    UEL1 = float(Eyoung)
    UEL2 = float(Nu)
    UEL3 = float(meshsize2)
    UEL4 = float(Gc)
    UEL5 = float(0)
    UEL6 = float(0)
    UEL7 = float(0)
    UEL8 = float(cond)
    UEL9 = float(pi11)
    UEL10 = float(pi12)
    UEL11 = float(pi44)
    UEL12 = float(k_value)
    UEL13 = float(n_value)

    UEL01 = str(UEL1)
    UEL02 = str(UEL2)
    UEL03 = str(UEL3)
    UEL04 = str(UEL4)
    UEL05 = str(UEL5)
    UEL06 = str(UEL6)
    UEL07 = str(UEL7)
    UEL08 = str(UEL8)
    UEL09 = str(UEL9)
    UEL010 = str(UEL10)
    UEL011 = "{:.7f}".format(UEL11)
    UEL012 = str(UEL12)
    UEL013 = str(UEL13)

    # Remove file Job-2.inp (if it exists)
    myfile = name_of_file + ".inp"
    if os.path.isfile(myfile):
        os.remove(myfile)

    # Generate the first file (F1)
    f = open(INP)
    content = f.read()
    text = content.split(Beg)[1].split(txtEnd)[0].strip()
    final_file = open("F1.txt", "w")
    final_file.write(text)
    final_file.close()
    f.close()

    # Change element type in F1
    fileToSearch = "F1.txt"
    tempFile = open(fileToSearch, "r+")
    for line in fileinput.input(fileToSearch):
        tempFile.write(line.replace(txtStart, "*ELEMENT, TYPE=U1, ELSET=SOLID"))

    tempFile.close()

    # Generate the second file (F2)
    f = open(INP)
    content = f.read()
    text = content.split(txtStart)[1].split(txtEnd)[0].strip()
    final_file = open("VisualElements.txt", "w")
    final_file.write(text)
    final_file.close()
    f.close()

    myDataFile = open("VisualElements.txt", "r+")
    myDataLines = myDataFile.readlines()
    myPointsList = [eval(dataLine) for dataLine in myDataLines]

    nelem = len(myPointsList)
    myDataFile.close()

    Order = math.floor(math.log10(nelem))
    num0 = 10 * 10 ** Order + 1

    DataFile = open("F2.txt", "a")
    if Eletype == 1.0:
        i = 0
        while i < nelem:
            num = i + num0
            linea_partida = myDataLines[i]
            nodes = tuple(map(str, linea_partida.split(",")))
            a = nodes[1]
            b = nodes[2]
            c = nodes[3]
            d = nodes[4]
            e = nodes[5]
            f = nodes[6]
            g = nodes[7]
            h = nodes[8]
            j = nodes[9]
            k = nodes[10]
            DataFile.writelines(
                "%d, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s"
                % (num, a, b, c, d, e, f, g, h, j, k)
            )
            i = i + 1
    elif Eletype == 2.0:
        i = 0
        while i < nelem:
            num = i + num0
            linea_partida = myDataLines[i]
            nodes = tuple(map(int, linea_partida.split(",")))
            a = nodes[1]
            b = nodes[2]
            c = nodes[3]
            d = nodes[4]
            e = nodes[5]
            f = nodes[6]
            g = nodes[7]
            h = nodes[8]
            if i < 1:
                DataFile.writelines(
                    "%d, %d, %d, %d, %d, %d, %d, %d, %d" % (num, a, b, c, d, e, f, g, h)
                )
            elif i > 0:
                DataFile.writelines(
                    "\n%d, %d, %d, %d, %d, %d, %d, %d, %d"
                    % (num, a, b, c, d, e, f, g, h)
                )
            i = i + 1
    elif Eletype == 3.0:
        i = 0
        while i < nelem:
            num = i / 2 + num0
            i1 = i + 1
            linea_partida = myDataLines[i]
            nodes = tuple(map(str, linea_partida.split(",")))
            a = nodes[1]
            b = nodes[2]
            c = nodes[3]
            d = nodes[4]
            e = nodes[5]
            f = nodes[6]
            g = nodes[7]
            h = nodes[8]
            j = nodes[9]
            k = nodes[10]
            a1 = nodes[11]
            b1 = nodes[12]
            c1 = nodes[13]
            d1 = nodes[14]
            e1 = nodes[15]
            linea_partida1 = myDataLines[i1]
            nodes1 = tuple(map(str, linea_partida1.split(",")))
            f1 = nodes1[0]
            g1 = nodes1[1]
            h1 = nodes1[2]
            j1 = nodes1[3]
            k1 = nodes1[4]
            DataFile.writelines(
                "%d, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, \n%s, %s, %s, %s, %s"
                % (
                    num,
                    a,
                    b,
                    c,
                    d,
                    e,
                    f,
                    g,
                    h,
                    j,
                    k,
                    a1,
                    b1,
                    c1,
                    d1,
                    e1,
                    f1,
                    g1,
                    h1,
                    j1,
                    k1,
                )
            )
            i = i + 2
    elif Eletype == 4.0:
        i = 0
        while i < nelem:
            num = i + num0
            linea_partida = myDataLines[i]
            nodes = tuple(map(int, linea_partida.split(",")))
            a = nodes[1]
            b = nodes[2]
            c = nodes[3]
            d = nodes[4]
            e = nodes[5]
            f = nodes[6]
            g = nodes[7]
            h = nodes[8]
            if i < 1:
                DataFile.writelines(
                    "%d, %d, %d, %d, %d, %d, %d, %d, %d" % (num, a, b, c, d, e, f, g, h)
                )
            elif i > 0:
                DataFile.writelines(
                    "\n%d, %d, %d, %d, %d, %d, %d, %d, %d"
                    % (num, a, b, c, d, e, f, g, h)
                )
            i = i + 1
    elif Eletype == 5.0:
        i = 0
        while i < nelem:
            num = i + num0
            linea_partida = myDataLines[i]
            nodes = tuple(map(int, linea_partida.split(",")))
            a = nodes[1]
            b = nodes[2]
            c = nodes[3]
            d = nodes[4]
            if i < 1:
                DataFile.writelines("%d, %d, %d, %d, %d" % (num, a, b, c, d))
            elif i > 0:
                DataFile.writelines("\n%d, %d, %d, %d, %d" % (num, a, b, c, d))
            i = i + 1

    DataFile.close()

    # Generate the second file (F3)
    f = open(INP)
    content = f.read()
    text = content.split(txtEnd)[1].split(End)[0].strip()
    final_file = open("F3.txt", "w")
    final_file.write(text)
    final_file.close()
    f.close()

    # Change DEPVAR in F3
    fileToSearch = "F3.txt"
    tempFile = open(fileToSearch, "r+")
    if Eletype < 4.0:
        for line in fileinput.input(fileToSearch):
            tempFile.write(
                line.replace(
                    "*User Material, constants=1",
                    "1, S11, S11\n2, S22, S22\n3, S33, S33\n4, S12, S12\n5, S13, S13\n6, S23, S23\n7, E11, E11\n8, E22, E22\n9, E33, E33\n10, E12, E12\n11, E13, E13\n12, E23, E23\n13, PHI, PHI\n14, H, H\n15, V, V\n16, JX, JX\n17, JY, JY\n18, JZ, JZ",
                )
            )
    elif Eletype > 3.0:
        for line in fileinput.input(fileToSearch):
            tempFile.write(
                line.replace(
                    "*User Material, constants=1",
                    "1, S11, S11\n2, S22, S22\n3, S33, S33\n4, S12, S12\n5, E11, E11\n6, E22, E22\n7, E33, E33\n8, E12, E12\n9, PHI, PHI\n10, H, H",
                )
            )

    tempFile.close()

    fileToSearch = "F3.txt"
    tempFile = open(fileToSearch, "r+")
    for line in fileinput.input(fileToSearch):
        tempFile.write(line.replace("0.,\n", "*User Material, constants=0\n"))
    tempFile.close()

    # Generate new input file
    DataF1 = open("F1.txt", "r+")
    F1 = DataF1.readlines()
    DataF1.close()

    DataF2 = open("F2.txt", "r+")
    F2 = DataF2.readlines()
    DataF2.close()

    DataF3 = open("F3.txt", "r+")
    F3 = DataF3.readlines()
    DataF3.close()

    DataFile = open(name_of_file + ".inp", "a")
    DataFile.writelines(Beg)
    DataFile.writelines("\n")
    DataFile.writelines(F1)
    DataFile.writelines("\n*UEL PROPERTY, ELSET=SOLID\n")
    DataFile.writelines(UEL01)
    DataFile.writelines(", ")
    DataFile.writelines(UEL02)
    DataFile.writelines(", ")
    DataFile.writelines(UEL03)
    DataFile.writelines(", ")
    DataFile.writelines(UEL04)
    DataFile.writelines(", ")
    DataFile.writelines(UEL05)
    DataFile.writelines(", ")
    DataFile.writelines(UEL06)
    DataFile.writelines(", ")
    DataFile.writelines(UEL07)
    DataFile.writelines(", ")
    DataFile.writelines(UEL08)
    DataFile.writelines(", ")
    DataFile.writelines("\n")
    DataFile.writelines(UEL09)
    DataFile.writelines(", ")
    DataFile.writelines(UEL010)
    DataFile.writelines(", ")
    DataFile.writelines(UEL011)
    DataFile.writelines(", ")
    DataFile.writelines(UEL012)
    DataFile.writelines(", ")
    DataFile.writelines(UEL013)

    if Eletype == 1.0:
        DataFile.writelines("\n*Element, type=C3D10, elset=Visualization\n")
    elif Eletype == 2.0:
        DataFile.writelines("\n*Element, type=C3D8, elset=Visualization\n")
    elif Eletype == 3.0:
        DataFile.writelines("\n*Element, type=C3D20R, elset=Visualization\n")
    elif Eletype == 4.0:
        DataFile.writelines("\n*Element, type=CPE8R, elset=Visualization\n")
    elif Eletype == 5.0:
        DataFile.writelines("\n*Element, type=CPE4, elset=Visualization\n")

    DataFile.writelines(F2)
    DataFile.writelines("\n")
    DataFile.writelines(txtEnd)
    DataFile.writelines("\n")
    DataFile.writelines(F3)
    DataFile.writelines("\n")
    DataFile.writelines(End)
    DataFile.writelines("\n")
    DataFile.close()

    # Change solid section
    fileToSearch = name_of_file + ".inp"
    tempFile = open(fileToSearch, "r+")
    if Eletype == 1.0:
        for line in fileinput.input(fileToSearch):
            tempFile.write(
                line.replace(
                    "*ELEMENT, TYPE=U1, ELSET=SOLID",
                    "*User element, nodes=10, type=U1, properties=13, coordinates=3, var=56\n1,2,3\n1,4\n1,9\n*ELEMENT, TYPE=U1, ELSET=SOLID",
                )
            )
    elif Eletype == 2.0:
        for line in fileinput.input(fileToSearch):
            tempFile.write(
                line.replace(
                    "*ELEMENT, TYPE=U1, ELSET=SOLID",
                    "*User element, nodes=8, type=U1, properties=13, coordinates=3, var=144\n1,2,3\n1,4\n1,9\n*ELEMENT, TYPE=U1, ELSET=SOLID",
                )
            )
    elif Eletype == 3.0:
        for line in fileinput.input(fileToSearch):
            tempFile.write(
                line.replace(
                    "*ELEMENT, TYPE=U1, ELSET=SOLID",
                    "*User element, nodes=20, type=U1, properties=13, coordinates=3, var=144\n1,2,3\n1,4\n1,9\n*ELEMENT, TYPE=U1, ELSET=SOLID",
                )
            )
    elif Eletype == 4.0:
        for line in fileinput.input(fileToSearch):
            tempFile.write(
                line.replace(
                    "*ELEMENT, TYPE=U1, ELSET=SOLID",
                    "*User element, nodes=8, type=U1, properties=13, coordinates=2, var=40\n1,2\n1,4\n*ELEMENT, TYPE=U1, ELSET=SOLID",
                )
            )
    elif Eletype == 5.0:
        for line in fileinput.input(fileToSearch):
            tempFile.write(
                line.replace(
                    "*ELEMENT, TYPE=U1, ELSET=SOLID",
                    "*User element, nodes=4, type=U1, properties=13, coordinates=2, var=40\n1,2\n1,4\n*ELEMENT, TYPE=U1, ELSET=SOLID",
                )
            )
    tempFile.close()

    fileToSearch = name_of_file + ".inp"
    tempFile = open(fileToSearch, "r+")
    for line in fileinput.input(fileToSearch):
        tempFile.write(
            line.replace(
                "*Solid Section, elset=Set-1, material=Material-1",
                "*Solid Section, elset=Visualization, material=Material-1",
            )
        )

    tempFile.close()

    # Eliminate existing files
    myfile = "F1.txt"
    if os.path.isfile(myfile):
        os.remove(myfile)

    myfile = "F2.txt"
    if os.path.isfile(myfile):
        os.remove(myfile)

    myfile = "F3.txt"
    if os.path.isfile(myfile):
        os.remove(myfile)

    myfile = "VisualElements.txt"
    if os.path.isfile(myfile):
        os.remove(myfile)

    myfile = "PF.inp"
    if os.path.isfile(myfile):
        os.remove(myfile)
    ##
    Mdb()
    mdb.ModelFromInputFile(name=name_of_file, inputFileName=name_of_file + ".inp")

    mdb.JobFromInputFile(
        name=name_of_file,
        inputFileName=name_of_file + ".inp",
        type=ANALYSIS,
        atTime=None,
        waitMinutes=0,
        waitHours=0,
        queue=None,
        memory=90,
        memoryUnits=PERCENTAGE,
        getMemoryFromAnalysis=True,
        explicitPrecision=SINGLE,
        nodalOutputPrecision=SINGLE,
        userSubroutine="UEL_piezoresistive_phasefield.for",
        scratch="",
        resultsFormat=ODB,
        multiprocessingMode=DEFAULT,
        numCpus=8,
        numDomains=8,
        numGPUs=0,
    )

    jobname = mdb.jobs[name_of_file]
    jobname.submit()
    jobname.waitForCompletion()


##################################
##################################
# Study case 2
def case_study_3(
    Vimp,
    Dimp,
    Inc,
    Eyoung,
    Nu,
    Gc,
    theta,
    cond,
    pi11,
    pi12,
    pi44,
    k_value,
    n_value,
    file_name,
    meshsize,
    meshsize2
):
    name_of_file = file_name

    # Parameters
    lx = 0.1
    ly = 0.2
    lz = 0.005
    delect = 0.01  # Distance from the edge
    lelect = 0.08  # Length of the electrode
    welect = 0.01  # Width of the electrode
    # Inc = 0.1 # Increment of steps
    holediam1 = 0.007
    sephor1 = 0.025
    spvert1 = 0.015
    holediam2 = 0.007
    sephor2 = 0.025
    spvert2 = 0.02
    lcrack = 0.02
    dcrack = 0.0002
    theta = theta * (np.pi / 180)

    # Phase-field
    ######################################################################
    ######################################################################
    ######################################################################
    ######################################################################
    ######################################################################
    ######################################################################
    # Model
    modelSK = mdb.models["Model-1"]
    # Strip
    s = modelSK.ConstrainedSketch(name="__profile__", sheetSize=200.0)
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=STANDALONE)
    s.rectangle(point1=(0.0, 0.0), point2=(lx, ly))
    s.EllipseByCenterPerimeter(
        center=(lx / 2.0, ly / 2.0),
        axisPoint1=(
            lx / 2.0 + lcrack * np.cos(theta),
            ly / 2.0 + lcrack * np.sin(theta),
        ),
        axisPoint2=(
            lx / 2.0 + dcrack * np.cos(theta + np.pi / 2.0),
            ly / 2.0 + dcrack * np.sin(theta + np.pi / 2.0),
        ),
    )
    s.CircleByCenterPerimeter(
        center=(lx / 2.0 - sephor2, ly / 2.0 + spvert2),
        point1=(lx / 2.0 - sephor2, ly / 2.0 + spvert2 + holediam2 / 2.0),
    )
    s.CircleByCenterPerimeter(
        center=(lx / 2.0 - sephor1, ly / 2.0 - spvert1),
        point1=(lx / 2.0 - sephor1, ly / 2.0 - spvert1 + holediam1 / 2.0),
    )
    s.CircleByCenterPerimeter(
        center=(lx / 2.0 + sephor2, ly / 2.0 + spvert2),
        point1=(lx / 2.0 + sephor2, ly / 2.0 + spvert2 + holediam2 / 2.0),
    )
    s.CircleByCenterPerimeter(
        center=(lx / 2.0 + sephor1, ly / 2.0 - spvert1),
        point1=(lx / 2.0 + sephor1, ly / 2.0 - spvert1 + holediam1 / 2.0),
    )
    p = modelSK.Part(name="Part-1", dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p = modelSK.parts["Part-1"]
    p.BaseSolidExtrude(sketch=s, depth=lz)
    s.unsetPrimaryObject()
    p = modelSK.parts["Part-1"]

    # Cut the strip to allocate the electrodes
    p = modelSK.parts["Part-1"]
    datumid = p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=delect)
    datumid2 = p.DatumPlaneByPrincipalPlane(
        principalPlane=XZPLANE, offset=delect + welect
    )
    datumid3 = p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=ly - delect)
    datumid4 = p.DatumPlaneByPrincipalPlane(
        principalPlane=XZPLANE, offset=ly - delect - welect
    )
    datumid5 = p.DatumPlaneByPrincipalPlane(
        principalPlane=YZPLANE, offset=(lx - lelect) / 2.0
    )
    datumid6 = p.DatumPlaneByPrincipalPlane(
        principalPlane=YZPLANE, offset=lx - (lx - lelect) / 2.0
    )

    datumid7 = p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=ly * 0.3)
    datumid8 = p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=ly * 0.7)
    datumid9 = p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=ly * 0.8)
    datumid10 = p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=ly * 0.2)

    p = modelSK.parts["Part-1"]
    c = p.cells
    pickedCells = c
    d2 = p.datums
    p.PartitionCellByDatumPlane(datumPlane=p.datums[datumid.id], cells=pickedCells)

    p = modelSK.parts["Part-1"]
    c = p.cells
    pickedCells = c
    d2 = p.datums
    p.PartitionCellByDatumPlane(datumPlane=p.datums[datumid2.id], cells=pickedCells)

    p = modelSK.parts["Part-1"]
    c = p.cells
    pickedCells = c
    d2 = p.datums
    p.PartitionCellByDatumPlane(datumPlane=p.datums[datumid3.id], cells=pickedCells)

    p = modelSK.parts["Part-1"]
    c = p.cells
    pickedCells = c
    d2 = p.datums
    p.PartitionCellByDatumPlane(datumPlane=p.datums[datumid4.id], cells=pickedCells)

    p = modelSK.parts["Part-1"]
    c = p.cells
    pickedCells = c
    d2 = p.datums
    p.PartitionCellByDatumPlane(datumPlane=p.datums[datumid5.id], cells=pickedCells)

    p = modelSK.parts["Part-1"]
    c = p.cells
    pickedCells = c
    d2 = p.datums
    p.PartitionCellByDatumPlane(datumPlane=p.datums[datumid6.id], cells=pickedCells)

    p = modelSK.parts["Part-1"]
    c = p.cells
    pickedCells = c
    d2 = p.datums
    p.PartitionCellByDatumPlane(datumPlane=p.datums[datumid7.id], cells=pickedCells)

    p = modelSK.parts["Part-1"]
    c = p.cells
    pickedCells = c
    d2 = p.datums
    p.PartitionCellByDatumPlane(datumPlane=p.datums[datumid8.id], cells=pickedCells)

    p = modelSK.parts["Part-1"]
    c = p.cells
    pickedCells = c
    d2 = p.datums
    p.PartitionCellByDatumPlane(datumPlane=p.datums[datumid9.id], cells=pickedCells)

    p = modelSK.parts["Part-1"]
    c = p.cells
    pickedCells = c
    d2 = p.datums
    p.PartitionCellByDatumPlane(datumPlane=p.datums[datumid10.id], cells=pickedCells)

    # Master Node
    p = modelSK.Part(name="Part-2", dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p.ReferencePoint(point=(lx / 2, delect + welect / 2.0, -5.0 / 10000))

    p = modelSK.Part(name="Part-3", dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p.ReferencePoint(point=(lx / 2, ly - delect - welect / 2.0, -5.0 / 10000))

    # Material
    modelSK.Material(name="Material-1")
    mdb.models["Model-1"].materials["Material-1"].Depvar(n=18)
    mdb.models["Model-1"].materials["Material-1"].UserMaterial(
        mechanicalConstants=(0.0,)
    )
    modelSK.HomogeneousSolidSection(
        name="Section-1", material="Material-1", thickness=None
    )
    p = modelSK.parts["Part-1"]
    c = p.cells
    cells = c
    region = p.Set(cells=cells, name="Set-1")
    p = modelSK.parts["Part-1"]
    p.SectionAssignment(
        region=region,
        sectionName="Section-1",
        offset=0.0,
        offsetType=MIDDLE_SURFACE,
        offsetField="",
        thicknessAssignment=FROM_SECTION,
    )

    # Instance
    a1 = modelSK.rootAssembly
    p = modelSK.parts["Part-1"]
    a1.Instance(name="Part-1-1", part=p, dependent=OFF)
    p = modelSK.parts["Part-2"]
    a1.Instance(name="Part-2-1", part=p, dependent=OFF)
    p = modelSK.parts["Part-3"]
    a1.Instance(name="Part-3-1", part=p, dependent=OFF)

    # Mesh
    elemType1 = mesh.ElemType(
        elemCode=C3D8,
        elemLibrary=STANDARD,
        secondOrderAccuracy=OFF,
        distortionControl=DEFAULT,
    )
    elemType2 = mesh.ElemType(elemCode=C3D6, elemLibrary=STANDARD)
    elemType3 = mesh.ElemType(elemCode=C3D4, elemLibrary=STANDARD)
    a = modelSK.rootAssembly
    c1 = a.instances["Part-1-1"].cells
    cells = c1
    pickedRegions = (cells,)
    a.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2, elemType3))

    partInstances = (
        a.instances["Part-1-1"],
        a.instances["Part-2-1"],
        a.instances["Part-3-1"],
    )
    a.seedPartInstance(
        regions=partInstances, size=meshsize, deviationFactor=0.1, minSizeFactor=0.1
    )
    e1 = a.instances["Part-1-1"].edges
    pickedEdges = e1.getSequenceFromMask(
        mask=("[#0:3 #e0000000 #3 ]",),
    )
    a.seedEdgeBySize(
        edges=pickedEdges,
        size=meshsize,
        deviationFactor=0.1,
        minSizeFactor=0.1,
        constraint=FINER,
    )

    partInstances = (
        a.instances["Part-1-1"],
        a.instances["Part-2-1"],
        a.instances["Part-3-1"],
    )
    a.seedPartInstance(
        regions=partInstances, size=meshsize, deviationFactor=0.1, minSizeFactor=0.1
    )
    partInstances = (a.instances["Part-1-1"],)
    a.generateMesh(regions=partInstances)

    # Create Node Sets of electrodes
    a = modelSK.rootAssembly
    n1 = a.instances["Part-1-1"].nodes
    faces1 = n1.getByBoundingBox(-1000, -1000, -1000, 1000, 1000, 1000)
    a.Set(nodes=faces1, name="ALL")

    a = modelSK.rootAssembly
    n1 = a.instances["Part-1-1"].nodes
    faces1 = n1.getByBoundingBox(0, 0, 0, lx, ly, lz)
    a.Set(nodes=faces1, name="body")

    a = modelSK.rootAssembly
    n1 = a.instances["Part-1-1"].elements
    faces1 = n1.getByBoundingBox(-1000, -1000, -1000, 1000, 1000, 1000)
    a.Set(nodes=faces1, name="ALLEL")

    a = modelSK.rootAssembly
    n1 = a.instances["Part-1-1"].nodes
    faces1 = n1.getByBoundingBox(-100, -100, -100, 100, 0.0001, 100)
    a.Set(nodes=faces1, name="Base")

    a = modelSK.rootAssembly
    n1 = a.instances["Part-1-1"].nodes
    faces1 = n1.getByBoundingBox(-100, ly, -100, 100, ly + 0.0001, 100)
    a.Set(nodes=faces1, name="Top")

    # a = modelSK.rootAssembly
    # n1 = a.instances['Part-1-1'].nodes
    # faces1 = n1.getByBoundingBox(-100,delect-0.0001,-100,100,delect+welect+0.0001,0.0001)
    # a.Set(nodes=faces1, name='Electrode1')

    # a = modelSK.rootAssembly
    # n1 = a.instances['Part-1-1'].nodes
    # faces1 = n1.getByBoundingBox(-100,ly-delect-welect-0.0001,-100,100,ly-delect+0.0001,0.0001)
    # a.Set(nodes=faces1, name='Electrode2')

    a = modelSK.rootAssembly
    n1 = a.instances["Part-1-1"].nodes
    faces1 = n1.getByBoundingBox(
        delect - 0.0001,
        delect - 0.0001,
        0.004,
        lelect + delect + 0.0001,
        delect + welect + 0.0001,
        0.005,
    )
    a.Set(nodes=faces1, name="Electrode2")

    a = modelSK.rootAssembly
    n1 = a.instances["Part-1-1"].nodes
    faces1 = n1.getByBoundingBox(
        delect - 0.0001,
        ly - delect - welect - 0.0001,
        0.004,
        lelect + delect + 0.0001,
        ly - delect + 0.0001,
        0.005,
    )
    a.Set(nodes=faces1, name="Electrode1")

    # Simulation
    modelSK.StaticStep(
        name="Step-1",
        previous="Initial",
        timeIncrementationMethod=FIXED,
        initialInc=Inc,
        noStop=OFF,
        reformKernel=25,
        solutionTechnique=QUASI_NEWTON,
    )
    modelSK.steps["Step-1"].control.setValues(
        allowPropagation=OFF,
        resetDefaultValues=OFF,
        timeIncrementation=(
            25000.0,
            25000.0,
            25000.0,
            25000.0,
            25000.0,
            4.0,
            12.0,
            5.0,
            6.0,
            3.0,
            50.0,
        ),
        lineSearch=(5.0, 1.0, 0.0001, 0.25, 0.1),
    )

    # Amplitudes
    mdb.models["Model-1"].TabularAmplitude(
        data=((0.0, 0.9999999), (1.0, 1.0)),
        name="Const",
        smooth=SOLVER_DEFAULT,
        timeSpan=STEP,
    )
    mdb.models["Model-1"].TabularAmplitude(
        data=((0.0, 0.0), (1.0, 1.0)),
        name="Linear",
        smooth=SOLVER_DEFAULT,
        timeSpan=STEP,
    )

    # BCs
    # Electrical
    modelSK.ElectricPotentialBC(
        amplitude=UNSET,
        createStepName="Step-1",
        distributionType=UNIFORM,
        fieldName="Const",
        fixed=OFF,
        magnitude=0.0,
        name="BC-3",
        region=mdb.models["Model-1"].rootAssembly.sets["Electrode1"],
    )

    modelSK.ElectricPotentialBC(
        amplitude="Const",
        createStepName="Step-1",
        distributionType=UNIFORM,
        fieldName="",
        fixed=OFF,
        magnitude=Vimp,
        name="BC-4",
        region=mdb.models["Model-1"].rootAssembly.sets["Electrode2"],
    )

    # Mechanical
    a = modelSK.rootAssembly
    region = a.sets["Base"]
    modelSK.DisplacementBC(
        name="BC-1",
        createStepName="Initial",
        region=region,
        u1=0.0,
        u2=0.0,
        u3=0.0,
        ur1=UNSET,
        ur2=UNSET,
        ur3=UNSET,
        amplitude=UNSET,
        fixed=OFF,
        distributionType=UNIFORM,
        fieldName="",
        localCsys=None,
    )
    a = modelSK.rootAssembly
    region = a.sets["Top"]
    modelSK.DisplacementBC(
        name="BC-2",
        createStepName="Step-1",
        region=region,
        u1=0.0,
        u2=Dimp,
        u3=0.0,
        ur1=UNSET,
        ur2=UNSET,
        ur3=UNSET,
        amplitude="Linear",
        fixed=OFF,
        distributionType=UNIFORM,
        fieldName="",
        localCsys=None,
    )

    # Create job
    # Desired output information
    mdb.models["Model-1"].FieldOutputRequest(
        name="F-Output-1",
        createStepName="Step-1",
    )
    modelSK.fieldOutputRequests["F-Output-1"].setValues(
        variables=("S", "LE", "U", "RF", "CF", "E", "SDV")
    )
    # Submit job
    mdb.Job(
        name="PF",
        model="Model-1",
        description="",
        type=ANALYSIS,
        atTime=None,
        waitMinutes=0,
        waitHours=0,
        queue=None,
        memory=90,
        memoryUnits=PERCENTAGE,
        getMemoryFromAnalysis=True,
        explicitPrecision=SINGLE,
        nodalOutputPrecision=SINGLE,
        echoPrint=OFF,
        modelPrint=OFF,
        contactPrint=OFF,
        historyPrint=OFF,
        userSubroutine="",
        scratch="",
        multiprocessingMode=DEFAULT,
        numCpus=1,
        numGPUs=0,
    )
    mdb.jobs["PF"].writeInput(consistencyChecking=OFF)

    ## Prepare imput
    ######################################################################
    ######################################################################
    ######################################################################
    ######################################################################
    ######################################################################
    ######################################################################

    # -*- coding: mbcs -*-
    # --------------------------------------------------------------
    Beg = "*Heading"  # Beginning of the input file
    End = "*End Step"  # End of the input file.

    Eletype = float(2)

    if Eletype == 1.0:
        txtStart = "*Element, type=C3D10"  # Text before the element connectivity
    elif Eletype == 2.0:
        txtStart = "*Element, type=C3D8"  # Text before the element connectivity
    elif Eletype == 3.0:
        txtStart = "*Element, type=C3D20R"  # Text before the element connectivity
    elif Eletype == 4.0:
        txtStart = "*Element, type=CPE8R"  # Text before the element connectivity
    elif Eletype == 5.0:
        txtStart = "*Element, type=CPE4"  # Text before the element connectivity

    txtEnd = "*Nset, nset=Set-1, generate"  # Text after the element connectivity
    # --------------------------------------------------------------

    INP = "PF.inp"
    UEL1 = float(Eyoung)
    UEL2 = float(Nu)
    UEL3 = float(meshsize2)
    UEL4 = float(Gc)
    UEL5 = float(0)
    UEL6 = float(0)
    UEL7 = float(0)
    UEL8 = float(cond)
    UEL9 = float(pi11)
    UEL10 = float(pi12)
    UEL11 = float(pi44)
    UEL12 = float(k_value)
    UEL13 = float(n_value)

    UEL01 = str(UEL1)
    UEL02 = str(UEL2)
    UEL03 = str(UEL3)
    UEL04 = str(UEL4)
    UEL05 = str(UEL5)
    UEL06 = str(UEL6)
    UEL07 = str(UEL7)
    UEL08 = str(UEL8)
    UEL09 = str(UEL9)
    UEL010 = str(UEL10)
    UEL011 = "{:.7f}".format(UEL11)
    UEL012 = str(UEL12)
    UEL013 = str(UEL13)

    # Remove file Job-2.inp (if it exists)
    myfile = name_of_file + ".inp"
    if os.path.isfile(myfile):
        os.remove(myfile)

    # Generate the first file (F1)
    f = open(INP)
    content = f.read()
    text = content.split(Beg)[1].split(txtEnd)[0].strip()
    final_file = open("F1.txt", "w")
    final_file.write(text)
    final_file.close()
    f.close()

    # Change element type in F1
    fileToSearch = "F1.txt"
    tempFile = open(fileToSearch, "r+")
    for line in fileinput.input(fileToSearch):
        tempFile.write(line.replace(txtStart, "*ELEMENT, TYPE=U1, ELSET=SOLID"))

    tempFile.close()

    # Generate the second file (F2)
    f = open(INP)
    content = f.read()
    text = content.split(txtStart)[1].split(txtEnd)[0].strip()
    final_file = open("VisualElements.txt", "w")
    final_file.write(text)
    final_file.close()
    f.close()

    myDataFile = open("VisualElements.txt", "r+")
    myDataLines = myDataFile.readlines()
    myPointsList = [eval(dataLine) for dataLine in myDataLines]

    nelem = len(myPointsList)
    myDataFile.close()

    Order = math.floor(math.log10(nelem))
    num0 = 10 * 10 ** Order + 1

    DataFile = open("F2.txt", "a")
    if Eletype == 1.0:
        i = 0
        while i < nelem:
            num = i + num0
            linea_partida = myDataLines[i]
            nodes = tuple(map(str, linea_partida.split(",")))
            a = nodes[1]
            b = nodes[2]
            c = nodes[3]
            d = nodes[4]
            e = nodes[5]
            f = nodes[6]
            g = nodes[7]
            h = nodes[8]
            j = nodes[9]
            k = nodes[10]
            DataFile.writelines(
                "%d, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s"
                % (num, a, b, c, d, e, f, g, h, j, k)
            )
            i = i + 1
    elif Eletype == 2.0:
        i = 0
        while i < nelem:
            num = i + num0
            linea_partida = myDataLines[i]
            nodes = tuple(map(int, linea_partida.split(",")))
            a = nodes[1]
            b = nodes[2]
            c = nodes[3]
            d = nodes[4]
            e = nodes[5]
            f = nodes[6]
            g = nodes[7]
            h = nodes[8]
            if i < 1:
                DataFile.writelines(
                    "%d, %d, %d, %d, %d, %d, %d, %d, %d" % (num, a, b, c, d, e, f, g, h)
                )
            elif i > 0:
                DataFile.writelines(
                    "\n%d, %d, %d, %d, %d, %d, %d, %d, %d"
                    % (num, a, b, c, d, e, f, g, h)
                )
            i = i + 1
    elif Eletype == 3.0:
        i = 0
        while i < nelem:
            num = i / 2 + num0
            i1 = i + 1
            linea_partida = myDataLines[i]
            nodes = tuple(map(str, linea_partida.split(",")))
            a = nodes[1]
            b = nodes[2]
            c = nodes[3]
            d = nodes[4]
            e = nodes[5]
            f = nodes[6]
            g = nodes[7]
            h = nodes[8]
            j = nodes[9]
            k = nodes[10]
            a1 = nodes[11]
            b1 = nodes[12]
            c1 = nodes[13]
            d1 = nodes[14]
            e1 = nodes[15]
            linea_partida1 = myDataLines[i1]
            nodes1 = tuple(map(str, linea_partida1.split(",")))
            f1 = nodes1[0]
            g1 = nodes1[1]
            h1 = nodes1[2]
            j1 = nodes1[3]
            k1 = nodes1[4]
            DataFile.writelines(
                "%d, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, \n%s, %s, %s, %s, %s"
                % (
                    num,
                    a,
                    b,
                    c,
                    d,
                    e,
                    f,
                    g,
                    h,
                    j,
                    k,
                    a1,
                    b1,
                    c1,
                    d1,
                    e1,
                    f1,
                    g1,
                    h1,
                    j1,
                    k1,
                )
            )
            i = i + 2
    elif Eletype == 4.0:
        i = 0
        while i < nelem:
            num = i + num0
            linea_partida = myDataLines[i]
            nodes = tuple(map(int, linea_partida.split(",")))
            a = nodes[1]
            b = nodes[2]
            c = nodes[3]
            d = nodes[4]
            e = nodes[5]
            f = nodes[6]
            g = nodes[7]
            h = nodes[8]
            if i < 1:
                DataFile.writelines(
                    "%d, %d, %d, %d, %d, %d, %d, %d, %d" % (num, a, b, c, d, e, f, g, h)
                )
            elif i > 0:
                DataFile.writelines(
                    "\n%d, %d, %d, %d, %d, %d, %d, %d, %d"
                    % (num, a, b, c, d, e, f, g, h)
                )
            i = i + 1
    elif Eletype == 5.0:
        i = 0
        while i < nelem:
            num = i + num0
            linea_partida = myDataLines[i]
            nodes = tuple(map(int, linea_partida.split(",")))
            a = nodes[1]
            b = nodes[2]
            c = nodes[3]
            d = nodes[4]
            if i < 1:
                DataFile.writelines("%d, %d, %d, %d, %d" % (num, a, b, c, d))
            elif i > 0:
                DataFile.writelines("\n%d, %d, %d, %d, %d" % (num, a, b, c, d))
            i = i + 1

    DataFile.close()

    # Generate the second file (F3)
    f = open(INP)
    content = f.read()
    text = content.split(txtEnd)[1].split(End)[0].strip()
    final_file = open("F3.txt", "w")
    final_file.write(text)
    final_file.close()
    f.close()

    # Change DEPVAR in F3
    fileToSearch = "F3.txt"
    tempFile = open(fileToSearch, "r+")
    if Eletype < 4.0:
        for line in fileinput.input(fileToSearch):
            tempFile.write(
                line.replace(
                    "*User Material, constants=1",
                    "1, S11, S11\n2, S22, S22\n3, S33, S33\n4, S12, S12\n5, S13, S13\n6, S23, S23\n7, E11, E11\n8, E22, E22\n9, E33, E33\n10, E12, E12\n11, E13, E13\n12, E23, E23\n13, PHI, PHI\n14, H, H\n15, V, V\n16, JX, JX\n17, JY, JY\n18, JZ, JZ",
                )
            )
    elif Eletype > 3.0:
        for line in fileinput.input(fileToSearch):
            tempFile.write(
                line.replace(
                    "*User Material, constants=1",
                    "1, S11, S11\n2, S22, S22\n3, S33, S33\n4, S12, S12\n5, E11, E11\n6, E22, E22\n7, E33, E33\n8, E12, E12\n9, PHI, PHI\n10, H, H",
                )
            )

    tempFile.close()

    fileToSearch = "F3.txt"
    tempFile = open(fileToSearch, "r+")
    for line in fileinput.input(fileToSearch):
        tempFile.write(line.replace("0.,\n", "*User Material, constants=0\n"))
    tempFile.close()

    # Generate new input file
    DataF1 = open("F1.txt", "r+")
    F1 = DataF1.readlines()
    DataF1.close()

    DataF2 = open("F2.txt", "r+")
    F2 = DataF2.readlines()
    DataF2.close()

    DataF3 = open("F3.txt", "r+")
    F3 = DataF3.readlines()
    DataF3.close()

    DataFile = open(name_of_file + ".inp", "a")
    DataFile.writelines(Beg)
    DataFile.writelines("\n")
    DataFile.writelines(F1)
    DataFile.writelines("\n*UEL PROPERTY, ELSET=SOLID\n")
    DataFile.writelines(UEL01)
    DataFile.writelines(", ")
    DataFile.writelines(UEL02)
    DataFile.writelines(", ")
    DataFile.writelines(UEL03)
    DataFile.writelines(", ")
    DataFile.writelines(UEL04)
    DataFile.writelines(", ")
    DataFile.writelines(UEL05)
    DataFile.writelines(", ")
    DataFile.writelines(UEL06)
    DataFile.writelines(", ")
    DataFile.writelines(UEL07)
    DataFile.writelines(", ")
    DataFile.writelines(UEL08)
    DataFile.writelines(", ")
    DataFile.writelines("\n")
    DataFile.writelines(UEL09)
    DataFile.writelines(", ")
    DataFile.writelines(UEL010)
    DataFile.writelines(", ")
    DataFile.writelines(UEL011)
    DataFile.writelines(", ")
    DataFile.writelines(UEL012)
    DataFile.writelines(", ")
    DataFile.writelines(UEL013)

    if Eletype == 1.0:
        DataFile.writelines("\n*Element, type=C3D10, elset=Visualization\n")
    elif Eletype == 2.0:
        DataFile.writelines("\n*Element, type=C3D8, elset=Visualization\n")
    elif Eletype == 3.0:
        DataFile.writelines("\n*Element, type=C3D20R, elset=Visualization\n")
    elif Eletype == 4.0:
        DataFile.writelines("\n*Element, type=CPE8R, elset=Visualization\n")
    elif Eletype == 5.0:
        DataFile.writelines("\n*Element, type=CPE4, elset=Visualization\n")

    DataFile.writelines(F2)
    DataFile.writelines("\n")
    DataFile.writelines(txtEnd)
    DataFile.writelines("\n")
    DataFile.writelines(F3)
    DataFile.writelines("\n")
    DataFile.writelines(End)
    DataFile.writelines("\n")
    DataFile.close()

    # Change solid section
    fileToSearch = name_of_file + ".inp"
    tempFile = open(fileToSearch, "r+")
    if Eletype == 1.0:
        for line in fileinput.input(fileToSearch):
            tempFile.write(
                line.replace(
                    "*ELEMENT, TYPE=U1, ELSET=SOLID",
                    "*User element, nodes=10, type=U1, properties=13, coordinates=3, var=56\n1,2,3\n1,4\n1,9\n*ELEMENT, TYPE=U1, ELSET=SOLID",
                )
            )
    elif Eletype == 2.0:
        for line in fileinput.input(fileToSearch):
            tempFile.write(
                line.replace(
                    "*ELEMENT, TYPE=U1, ELSET=SOLID",
                    "*User element, nodes=8, type=U1, properties=13, coordinates=3, var=144\n1,2,3\n1,4\n1,9\n*ELEMENT, TYPE=U1, ELSET=SOLID",
                )
            )
    elif Eletype == 3.0:
        for line in fileinput.input(fileToSearch):
            tempFile.write(
                line.replace(
                    "*ELEMENT, TYPE=U1, ELSET=SOLID",
                    "*User element, nodes=20, type=U1, properties=13, coordinates=3, var=144\n1,2,3\n1,4\n1,9\n*ELEMENT, TYPE=U1, ELSET=SOLID",
                )
            )
    elif Eletype == 4.0:
        for line in fileinput.input(fileToSearch):
            tempFile.write(
                line.replace(
                    "*ELEMENT, TYPE=U1, ELSET=SOLID",
                    "*User element, nodes=8, type=U1, properties=13, coordinates=2, var=40\n1,2\n1,4\n*ELEMENT, TYPE=U1, ELSET=SOLID",
                )
            )
    elif Eletype == 5.0:
        for line in fileinput.input(fileToSearch):
            tempFile.write(
                line.replace(
                    "*ELEMENT, TYPE=U1, ELSET=SOLID",
                    "*User element, nodes=4, type=U1, properties=13, coordinates=2, var=40\n1,2\n1,4\n*ELEMENT, TYPE=U1, ELSET=SOLID",
                )
            )
    tempFile.close()

    fileToSearch = name_of_file + ".inp"
    tempFile = open(fileToSearch, "r+")
    for line in fileinput.input(fileToSearch):
        tempFile.write(
            line.replace(
                "*Solid Section, elset=Set-1, material=Material-1",
                "*Solid Section, elset=Visualization, material=Material-1",
            )
        )

    tempFile.close()

    # Eliminate existing files
    myfile = "F1.txt"
    if os.path.isfile(myfile):
        os.remove(myfile)

    myfile = "F2.txt"
    if os.path.isfile(myfile):
        os.remove(myfile)

    myfile = "F3.txt"
    if os.path.isfile(myfile):
        os.remove(myfile)

    myfile = "VisualElements.txt"
    if os.path.isfile(myfile):
        os.remove(myfile)

    myfile = "PF.inp"
    if os.path.isfile(myfile):
        os.remove(myfile)
    ##
    Mdb()
    mdb.ModelFromInputFile(name=name_of_file, inputFileName=name_of_file + ".inp")

    mdb.JobFromInputFile(
        name=name_of_file,
        inputFileName=name_of_file + ".inp",
        type=ANALYSIS,
        atTime=None,
        waitMinutes=0,
        waitHours=0,
        queue=None,
        memory=90,
        memoryUnits=PERCENTAGE,
        getMemoryFromAnalysis=True,
        explicitPrecision=SINGLE,
        nodalOutputPrecision=SINGLE,
        userSubroutine="UEL_piezoresistive_phasefield.for",
        scratch="",
        resultsFormat=ODB,
        multiprocessingMode=DEFAULT,
        numCpus=8,
        numDomains=8,
        numGPUs=0,
    )

    jobname = mdb.jobs[name_of_file]
    jobname.submit()
    jobname.waitForCompletion()


##################################
##################################
# Study case 3
def case_study_4(
    Vimp,
    Dimp,
    Inc,
    Eyoung,
    Nu,
    Gc,
    theta,
    cond,
    pi11,
    pi12,
    pi44,
    k_value,
    n_value,
    file_name,
    meshsize,
    meshsize2
):
    name_of_file = file_name
    # Parameters
    lx = 0.1
    ly = 0.2
    lz = 0.005
    delect = 0.01  # Distance from the edge
    lelect = 0.08  # Length of the electrode
    welect = 0.01  # Width of the electrode
    holediam1 = 0.007
    sephor1 = 0.02
    spvert1 = 0.015
    holediam2 = 0.005
    sephor2 = 0.025
    spvert2 = 0.02
    lcrack = 0.02
    dcrack = 0.0002
    theta = theta * (np.pi / 180)

    # Phase-field
    ######################################################################
    ######################################################################
    ######################################################################
    ######################################################################
    ######################################################################
    ######################################################################

    # Model
    modelSK = mdb.models["Model-1"]
    # Strip
    s = modelSK.ConstrainedSketch(name="__profile__", sheetSize=200.0)
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=STANDALONE)
    s.rectangle(point1=(0.0, 0.0), point2=(lx, ly))
    p = modelSK.Part(name="Part-1", dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p = modelSK.parts["Part-1"]
    p.BaseSolidExtrude(sketch=s, depth=lz)
    s.unsetPrimaryObject()
    p = modelSK.parts["Part-1"]

    # Cut the strip to allocate the electrodes
    p = modelSK.parts["Part-1"]
    datumid = p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=delect)
    datumid2 = p.DatumPlaneByPrincipalPlane(
        principalPlane=XZPLANE, offset=delect + welect
    )
    datumid3 = p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=ly - delect)
    datumid4 = p.DatumPlaneByPrincipalPlane(
        principalPlane=XZPLANE, offset=ly - delect - welect
    )
    datumid5 = p.DatumPlaneByPrincipalPlane(
        principalPlane=YZPLANE, offset=(lx - lelect) / 2.0
    )
    datumid6 = p.DatumPlaneByPrincipalPlane(
        principalPlane=YZPLANE, offset=lx - (lx - lelect) / 2.0
    )

    # datumid7 = p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=ly*0.3)
    # datumid8 = p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=ly*0.7)
    # datumid9 = p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=ly*0.8)
    # datumid10 = p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=ly*0.2)

    p = modelSK.parts["Part-1"]
    c = p.cells
    pickedCells = c
    d2 = p.datums
    p.PartitionCellByDatumPlane(datumPlane=p.datums[datumid.id], cells=pickedCells)

    p = modelSK.parts["Part-1"]
    c = p.cells
    pickedCells = c
    d2 = p.datums
    p.PartitionCellByDatumPlane(datumPlane=p.datums[datumid2.id], cells=pickedCells)

    p = modelSK.parts["Part-1"]
    c = p.cells
    pickedCells = c
    d2 = p.datums
    p.PartitionCellByDatumPlane(datumPlane=p.datums[datumid3.id], cells=pickedCells)

    p = modelSK.parts["Part-1"]
    c = p.cells
    pickedCells = c
    d2 = p.datums
    p.PartitionCellByDatumPlane(datumPlane=p.datums[datumid4.id], cells=pickedCells)

    p = modelSK.parts["Part-1"]
    c = p.cells
    pickedCells = c
    d2 = p.datums
    p.PartitionCellByDatumPlane(datumPlane=p.datums[datumid5.id], cells=pickedCells)

    p = modelSK.parts["Part-1"]
    c = p.cells
    pickedCells = c
    d2 = p.datums
    p.PartitionCellByDatumPlane(datumPlane=p.datums[datumid6.id], cells=pickedCells)

    # p = modelSK.parts['Part-1']
    # c = p.cells
    # pickedCells = c
    # d2 = p.datums
    # p.PartitionCellByDatumPlane(datumPlane=p.datums[datumid7.id], cells=pickedCells)

    # p = modelSK.parts['Part-1']
    # c = p.cells
    # pickedCells = c
    # d2 = p.datums
    # p.PartitionCellByDatumPlane(datumPlane=p.datums[datumid8.id], cells=pickedCells)

    # p = modelSK.parts['Part-1']
    # c = p.cells
    # pickedCells = c
    # d2 = p.datums
    # p.PartitionCellByDatumPlane(datumPlane=p.datums[datumid9.id], cells=pickedCells)

    # p = modelSK.parts['Part-1']
    # c = p.cells
    # pickedCells = c
    # d2 = p.datums
    # p.PartitionCellByDatumPlane(datumPlane=p.datums[datumid10.id], cells=pickedCells)
    ##

    # Master Node
    p = modelSK.Part(name="Part-2", dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p.ReferencePoint(point=(lx / 2, delect + welect / 2.0, -5.0 / 10000))

    p = modelSK.Part(name="Part-3", dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p.ReferencePoint(point=(lx / 2, ly - delect - welect / 2.0, -5.0 / 10000))

    # Instance
    a1 = modelSK.rootAssembly
    p = modelSK.parts["Part-1"]
    a1.Instance(name="Part-1-1", part=p, dependent=OFF)
    p = modelSK.parts["Part-2"]
    a1.Instance(name="Part-2-1", part=p, dependent=OFF)
    p = modelSK.parts["Part-3"]
    a1.Instance(name="Part-3-1", part=p, dependent=OFF)

    # Create cracks
    initialvolume = lx * ly * lz
    crackdensity = 10 / 1000.0

    nparts = 0
    relvol = 0
    cont = 4
    part_volume = 0
    namepartorig = "Part-1"

    while relvol < crackdensity:
        mean1 = 0.2 / 100.0
        mean2 = 0.2 / 100.0
        std1 = mean1 * 12 / 100.0
        std2 = mean2 * 0.3 / 100.0

        namepart = "Part-" + str(cont)
        nameinst = "Part-" + str(cont) + "-1"
        namepartorig = "Part-" + str(cont + 1000)
        nameinstorig = "Part-" + str(cont + 1000) + "-1"
        namepartorigprev = "Part-" + str(cont + 1000 - 1)
        nameinstorigprev = "Part-" + str(cont + 1000 - 1) + "-1"
        semiaxis1 = float(np.random.randn(1) * std1 + mean1)
        s = modelSK.ConstrainedSketch(name="__profile__", sheetSize=200.0)
        # g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
        # s.setPrimaryObject(option=STANDALONE)
        # s.ConstructionLine(point1=(0.0, -100.0), point2=(0.0, 100.0))
        # s.FixedConstraint(entity=g[2])
        # s.EllipseByCenterPerimeter(center=(0.0, 0.0), axisPoint1=(0.0, semiaxis1),
        # axisPoint2=(semiaxis2, 0.0))
        # s.CoincidentConstraint(entity1=v[0], entity2=g[2], addUndoState=False)
        # s.autoTrimCurve(curve1=g[3], point1=(-semiaxis2, 0))
        # s.Line(point1=(0.0, semiaxis1), point2=(0.0, -semiaxis1))
        # s.VerticalConstraint(entity=g[6], addUndoState=False)
        # s.PerpendicularConstraint(entity1=g[5], entity2=g[6], addUndoState=False)
        p = modelSK.Part(name=namepart, dimensionality=THREE_D, type=DEFORMABLE_BODY)
        p = modelSK.parts[namepart]
        # p.BaseSolidRevolve(sketch=s, angle=360, flipRevolveDirection=OFF)
        s.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(semiaxis1, 0))
        p.BaseSolidExtrude(sketch=s, depth=lz)
        s.unsetPrimaryObject()
        datumid = p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=0.0)
        datumid2 = p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=0.0)
        datumid3 = p.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=0.0)
        # c = p.cells
        # pickedCells = c
        # d2 = p.datums
        # p.PartitionCellByDatumPlane(datumPlane=p.datums[datumid.id], cells=pickedCells)
        # p = modelSK.parts[namepart]
        # c = p.cells
        # pickedCells = c
        # d1 = p.datums
        # p.PartitionCellByDatumPlane(datumPlane=p.datums[datumid2.id], cells=pickedCells)
        # p = modelSK.parts[namepart]
        # c = p.cells
        # pickedCells = c
        # d2 = p.datums
        # p.PartitionCellByDatumPlane(datumPlane=p.datums[datumid3.id], cells=pickedCells)
        volellipsoid = (4.0 / 3.0) * np.pi * semiaxis2 * semiaxis2 * semiaxis1

        # Instance
        p = modelSK.parts[namepart]
        a1.Instance(name=nameinst, part=p, dependent=OFF)

        # Rotate and rotate crack
        angletorot = float(np.random.rand(1)) * 45
        angle1 = float(np.random.rand(1)) * 360
        angle2 = 2 * (float(np.random.rand(1)) - 0.5) * 180
        xpos = np.cos(angle1 * np.pi / 180.0)
        ypos = np.sin(angle1 * np.pi / 180.0)
        zpos = np.sin(angle2 * np.pi / 180.0)
        dx = lx / 6 + float(np.random.rand(1)) * (lx - 2 * lx / 6)  # *0+lx/2
        dy = ly / 6 + float(np.random.rand(1)) * (ly - 2 * ly / 6)  # *0+ly/2
        dz = lz / 6 + float(np.random.rand(1)) * (lz - 2 * lz / 6)  # *0+lz/2
        a = modelSK.rootAssembly
        # a.rotate(instanceList=(nameinst, ), axisPoint=(0.0, 0.0, 0.0),axisDirection=(xpos, ypos, zpos), angle=angletorot)
        # a = modelSK.rootAssembly
        a.translate(instanceList=(nameinst,), vector=(dx, dy, dz * 0))

        # CUT
        a = modelSK.rootAssembly
        if cont > 4:
            a.InstanceFromBooleanCut(
                name=namepartorig,
                instanceToBeCut=modelSK.rootAssembly.instances[nameinstorigprev],
                cuttingInstances=(a.instances[nameinst],),
                originalInstances=SUPPRESS,
            )
            # a.InstanceFromBooleanMerge(name=namepartorig, instances=(a.instances[nameinstorigprev],
            # a.instances[nameinst], ), keepIntersections=ON,
            # originalInstances=SUPPRESS, domain=GEOMETRY)
        else:
            a.InstanceFromBooleanCut(
                name=namepartorig,
                instanceToBeCut=modelSK.rootAssembly.instances["Part-1-1"],
                cuttingInstances=(a.instances[nameinst],),
                originalInstances=SUPPRESS,
            )
            # a.InstanceFromBooleanMerge(name=namepartorig, instances=(a.instances['Part-1-1'],
            # a.instances[nameinst], ), keepIntersections=ON,
            # originalInstances=SUPPRESS, domain=GEOMETRY)

        mask = modelSK.parts[namepartorig].cells.getMask()
        cellobj_sequence = modelSK.parts[namepartorig].cells.getSequenceFromMask(
            mask=mask
        )
        part_volume = modelSK.parts[namepartorig].getVolume(cells=cellobj_sequence)
        relvol = 1 - part_volume / initialvolume
        # mask=modelSK.parts[namepart].cells.getMask()
        # part_volume = part_volume+volellipsoid
        # relvol = part_volume/initialvolume
        cont = cont + 1
        print(relvol)
        print(namepart)
        # a = modelSK.rootAssembly
        # del mdb.models['Model-1'].parts[namepart]
        # namepart = namepart+'-1'
        # del a.features[namepart]
        nparts = nparts + 1

    cont = 4
    while cont < nparts + 3:
        nametodel = "Part-" + str(cont + 1000) + "-1"
        del a.features[nametodel]
        nametodel = "Part-" + str(cont) + "-1"
        del a.features[nametodel]
        nametodel = "Part-" + str(cont + 1000)
        del modelSK.parts[nametodel]
        nametodel = "Part-" + str(cont)
        del modelSK.parts[nametodel]
        cont = cont + 1
    if nparts > 1:
        nametodel = "Part-" + str(nparts + 3)
        del modelSK.parts[nametodel]

    # Suppress unnecesary assemblies
    a = modelSK.rootAssembly
    # a.features['Part-2-1'].suppress()
    # a = modelSK.rootAssembly
    # a.features['Part-3-1'].suppress()
    if nparts > 1:
        del a.features["Part-1-1"]
        a = mdb.models["Model-1"].rootAssembly
        a.regenerate()
        nametodel = "Part-" + str(nparts + 3) + "-1"
        del mdb.models["Model-1"].parts["Part-1"]
        del a.features[nametodel]

    # Material
    modelSK.Material(name="Material-1")
    mdb.models["Model-1"].materials["Material-1"].Depvar(n=18)
    mdb.models["Model-1"].materials["Material-1"].UserMaterial(
        mechanicalConstants=(0.0,)
    )
    modelSK.HomogeneousSolidSection(
        name="Section-1", material="Material-1", thickness=None
    )
    p = modelSK.parts[namepartorig]
    c = p.cells
    cells = c
    region = p.Set(cells=cells, name="Set-1")
    p = modelSK.parts[namepartorig]
    p.SectionAssignment(
        region=region,
        sectionName="Section-1",
        offset=0.0,
        offsetType=MIDDLE_SURFACE,
        offsetField="",
        thicknessAssignment=FROM_SECTION,
    )

    # Mesh
    p = modelSK.parts[namepartorig]
    c = p.cells
    pickedRegions = c
    # p.setMeshControls(regions=pickedRegions, elemShape=TET, technique=FREE)
    elemType1 = mesh.ElemType(
        elemCode=C3D8,
        elemLibrary=STANDARD,
        secondOrderAccuracy=OFF,
        distortionControl=DEFAULT,
    )
    elemType2 = mesh.ElemType(elemCode=C3D6, elemLibrary=STANDARD)
    elemType3 = mesh.ElemType(elemCode=C3D4, elemLibrary=STANDARD)
    p = modelSK.parts[namepartorig]
    c = p.cells
    cells = c
    pickedRegions = (cells,)
    p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2, elemType3))
    p = modelSK.parts[namepartorig]
    p.seedPart(size=meshsize, deviationFactor=0.1, minSizeFactor=0.01)
    p = modelSK.parts[namepartorig]
    e = p.edges
    pickedEdges = e.getByBoundingBox(-100, -100, -100, delect, 100, 100)
    p.seedEdgeBySize(
        edges=pickedEdges,
        size=meshsize,
        deviationFactor=0.1,
        minSizeFactor=0.001,
        constraint=FINER,
    )
    pickedEdges = e.getByBoundingBox(lx - delect, -100, -100, 100, 100, 100)
    p.seedEdgeBySize(
        edges=pickedEdges,
        size=meshsize,
        deviationFactor=0.1,
        minSizeFactor=0.001,
        constraint=FINER,
    )
    pickedEdges = e.getByBoundingBox(-100, -100, -100, 100, delect + welect, 100)
    p.seedEdgeBySize(
        edges=pickedEdges,
        size=meshsize,
        deviationFactor=0.1,
        minSizeFactor=0.001,
        constraint=FINER,
    )
    pickedEdges = e.getByBoundingBox(-100, ly - delect - welect, -100, 100, 100, 100)
    p.seedEdgeBySize(
        edges=pickedEdges,
        size=meshsize,
        deviationFactor=0.1,
        minSizeFactor=0.001,
        constraint=FINER,
    )
    p.generateMesh()
    a = mdb.models["Model-1"].rootAssembly
    a.regenerate()

    # Create Node Sets of electrodes
    a = modelSK.rootAssembly
    n1 = a.instances[nameinstorig].nodes
    faces1 = n1.getByBoundingBox(-100, -100, -100, 100, 0.0001, 100)
    a.Set(nodes=faces1, name="Base")

    a = modelSK.rootAssembly
    n1 = a.instances[nameinstorig].nodes
    faces1 = n1.getByBoundingBox(-100, ly, -100, 100, ly + 0.0001, 100)
    a.Set(nodes=faces1, name="Top")

    a = modelSK.rootAssembly
    n1 = a.instances[nameinstorig].nodes
    faces1 = n1.getByBoundingBox(
        delect - 0.0001,
        ly - delect - welect - 0.0001,
        0.004,
        lelect + delect + 0.0001,
        ly - delect + 0.0001,
        0.005,
    )
    a.Set(nodes=faces1, name="Electrode1")

    a = modelSK.rootAssembly
    n1 = a.instances[nameinstorig].nodes
    faces1 = n1.getByBoundingBox(
        delect - 0.0001,
        delect - 0.0001,
        0.004,
        lelect + delect + 0.0001,
        delect + welect + 0.0001,
        0.005,
    )
    a.Set(nodes=faces1, name="Electrode2")

    # Simulation
    # modelSK.StaticStep(name='Step-1', previous='Initial',
    #    timeIncrementationMethod=AUTOMATIC, initialInc=Inc, minInc=1e-05, maxInc=0.1, maxNumInc=50000, noStop=OFF,
    #    reformKernel=25, solutionTechnique=QUASI_NEWTON)

    modelSK.StaticStep(
        name="Step-1",
        previous="Initial",
        timeIncrementationMethod=FIXED,
        initialInc=Inc,
        noStop=OFF,
        reformKernel=25,
        solutionTechnique=QUASI_NEWTON,
    )

    modelSK.steps["Step-1"].control.setValues(
        allowPropagation=OFF,
        resetDefaultValues=OFF,
        timeIncrementation=(
            25000.0,
            25000.0,
            25000.0,
            25000.0,
            25000.0,
            4.0,
            12.0,
            5.0,
            6.0,
            3.0,
            50.0,
        ),
        lineSearch=(5.0, 1.0, 0.0001, 0.25, 0.1),
    )
    # Amplitudes
    mdb.models["Model-1"].TabularAmplitude(
        data=((0.0, 0.9999999), (1.0, 1.0)),
        name="Const",
        smooth=SOLVER_DEFAULT,
        timeSpan=STEP,
    )
    mdb.models["Model-1"].TabularAmplitude(
        data=((0.0, 0.0), (1.0, 1.0)),
        name="Linear",
        smooth=SOLVER_DEFAULT,
        timeSpan=STEP,
    )

    # BCs
    # Electrical
    modelSK.ElectricPotentialBC(
        amplitude=UNSET,
        createStepName="Step-1",
        distributionType=UNIFORM,
        fieldName="Const",
        fixed=OFF,
        magnitude=0.0,
        name="BC-3",
        region=mdb.models["Model-1"].rootAssembly.sets["Electrode1"],
    )

    modelSK.ElectricPotentialBC(
        amplitude="Const",
        createStepName="Step-1",
        distributionType=UNIFORM,
        fieldName="",
        fixed=OFF,
        magnitude=Vimp,
        name="BC-4",
        region=mdb.models["Model-1"].rootAssembly.sets["Electrode2"],
    )

    # Mechanical
    a = modelSK.rootAssembly
    region = a.sets["Base"]
    modelSK.DisplacementBC(
        name="BC-1",
        createStepName="Initial",
        region=region,
        u1=0.0,
        u2=0.0,
        u3=0.0,
        ur1=UNSET,
        ur2=UNSET,
        ur3=UNSET,
        amplitude=UNSET,
        fixed=OFF,
        distributionType=UNIFORM,
        fieldName="",
        localCsys=None,
    )
    a = modelSK.rootAssembly
    region = a.sets["Top"]
    modelSK.DisplacementBC(
        name="BC-2",
        createStepName="Step-1",
        region=region,
        u1=0.0,
        u2=Dimp,
        u3=0.0,
        ur1=UNSET,
        ur2=UNSET,
        ur3=UNSET,
        amplitude=UNSET,
        fixed=OFF,
        distributionType=UNIFORM,
        fieldName="",
        localCsys=None,
    )

    # Create job
    # Desired output information
    mdb.models["Model-1"].FieldOutputRequest(
        name="F-Output-1",
        createStepName="Step-1",
    )
    modelSK.fieldOutputRequests["F-Output-1"].setValues(
        variables=("S", "LE", "U", "RF", "CF", "E", "SDV")
    )
    # Submit job
    mdb.Job(
        name="PF",
        model="Model-1",
        description="",
        type=ANALYSIS,
        atTime=None,
        waitMinutes=0,
        waitHours=0,
        queue=None,
        memory=90,
        memoryUnits=PERCENTAGE,
        getMemoryFromAnalysis=True,
        explicitPrecision=SINGLE,
        nodalOutputPrecision=SINGLE,
        echoPrint=OFF,
        modelPrint=OFF,
        contactPrint=OFF,
        historyPrint=OFF,
        userSubroutine="",
        scratch="",
        multiprocessingMode=DEFAULT,
        numCpus=1,
        numGPUs=0,
    )
    mdb.jobs["PF"].writeInput(consistencyChecking=OFF)

    ## Prepare imput
    ######################################################################
    ######################################################################
    ######################################################################
    ######################################################################
    ######################################################################
    ######################################################################
    # --------------------------------------------------------------
    Beg = "*Heading"  # Beginning of the input file
    End = "*End Step"  # End of the input file.

    Eletype = float(2)

    if Eletype == 1.0:
        txtStart = "*Element, type=C3D10"  # Text before the element connectivity
    elif Eletype == 2.0:
        txtStart = "*Element, type=C3D8"  # Text before the element connectivity
    elif Eletype == 3.0:
        txtStart = "*Element, type=C3D20R"  # Text before the element connectivity
    elif Eletype == 4.0:
        txtStart = "*Element, type=CPE8R"  # Text before the element connectivity
    elif Eletype == 5.0:
        txtStart = "*Element, type=CPE4"  # Text before the element connectivity

    txtEnd = "*Nset, nset=Set-1, generate"  # Text after the element connectivity
    # --------------------------------------------------------------

    INP = "PF.inp"
    UEL1 = float(Eyoung)
    UEL2 = float(Nu)
    UEL3 = float(meshsize2)
    UEL4 = float(Gc)
    UEL5 = float(0)
    UEL6 = float(0)
    UEL7 = float(0)
    UEL8 = float(cond)
    UEL9 = float(pi11)
    UEL10 = float(pi12)
    UEL11 = float(pi44)
    UEL12 = float(k_value)
    UEL13 = float(n_value)

    UEL01 = str(UEL1)
    UEL02 = str(UEL2)
    UEL03 = str(UEL3)
    UEL04 = str(UEL4)
    UEL05 = str(UEL5)
    UEL06 = str(UEL6)
    UEL07 = str(UEL7)
    UEL08 = str(UEL8)
    UEL09 = str(UEL9)
    UEL010 = str(UEL10)
    UEL011 = "{:.7f}".format(UEL11)
    UEL012 = str(UEL12)
    UEL013 = str(UEL13)

    # Remove file .inp (if it exists)
    myfile = name_of_file + ".inp"
    if os.path.isfile(myfile):
        os.remove(myfile)

    # Generate the first file (F1)
    f = open(INP)
    content = f.read()
    text = content.split(Beg)[1].split(txtEnd)[0].strip()
    final_file = open("F1.txt", "w")
    final_file.write(text)
    final_file.close()
    f.close()

    # Change element type in F1
    fileToSearch = "F1.txt"
    tempFile = open(fileToSearch, "r+")
    for line in fileinput.input(fileToSearch):
        tempFile.write(line.replace(txtStart, "*ELEMENT, TYPE=U1, ELSET=SOLID"))

    tempFile.close()

    # Generate the second file (F2)
    f = open(INP)
    content = f.read()
    print(txtEnd)
    print(txtStart)

    text = content.split(txtStart)[1].split(txtEnd)[0].strip()

    final_file = open("VisualElements.txt", "w")
    final_file.write(text)
    final_file.close()
    f.close()

    myDataFile = open("VisualElements.txt", "r+")
    myDataLines = myDataFile.readlines()
    myPointsList = [eval(dataLine) for dataLine in myDataLines]

    nelem = len(myPointsList)
    myDataFile.close()

    Order = math.floor(math.log10(nelem))
    num0 = 10 * 10 ** Order + 1

    DataFile = open("F2.txt", "a")
    if Eletype == 1.0:
        i = 0
        while i < nelem:
            num = i + num0
            linea_partida = myDataLines[i]
            nodes = tuple(map(str, linea_partida.split(",")))
            a = nodes[1]
            b = nodes[2]
            c = nodes[3]
            d = nodes[4]
            e = nodes[5]
            f = nodes[6]
            g = nodes[7]
            h = nodes[8]
            j = nodes[9]
            k = nodes[10]
            DataFile.writelines(
                "%d, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s"
                % (num, a, b, c, d, e, f, g, h, j, k)
            )
            i = i + 1
    elif Eletype == 2.0:
        i = 0
        while i < nelem:
            num = i + num0
            linea_partida = myDataLines[i]
            nodes = tuple(map(int, linea_partida.split(",")))
            a = nodes[1]
            b = nodes[2]
            c = nodes[3]
            d = nodes[4]
            e = nodes[5]
            f = nodes[6]
            g = nodes[7]
            h = nodes[8]
            if i < 1:
                DataFile.writelines(
                    "%d, %d, %d, %d, %d, %d, %d, %d, %d" % (num, a, b, c, d, e, f, g, h)
                )
            elif i > 0:
                DataFile.writelines(
                    "\n%d, %d, %d, %d, %d, %d, %d, %d, %d"
                    % (num, a, b, c, d, e, f, g, h)
                )
            i = i + 1
    elif Eletype == 3.0:
        i = 0
        while i < nelem:
            num = i / 2 + num0
            i1 = i + 1
            linea_partida = myDataLines[i]
            nodes = tuple(map(str, linea_partida.split(",")))
            a = nodes[1]
            b = nodes[2]
            c = nodes[3]
            d = nodes[4]
            e = nodes[5]
            f = nodes[6]
            g = nodes[7]
            h = nodes[8]
            j = nodes[9]
            k = nodes[10]
            a1 = nodes[11]
            b1 = nodes[12]
            c1 = nodes[13]
            d1 = nodes[14]
            e1 = nodes[15]
            linea_partida1 = myDataLines[i1]
            nodes1 = tuple(map(str, linea_partida1.split(",")))
            f1 = nodes1[0]
            g1 = nodes1[1]
            h1 = nodes1[2]
            j1 = nodes1[3]
            k1 = nodes1[4]
            DataFile.writelines(
                "%d, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, \n%s, %s, %s, %s, %s"
                % (
                    num,
                    a,
                    b,
                    c,
                    d,
                    e,
                    f,
                    g,
                    h,
                    j,
                    k,
                    a1,
                    b1,
                    c1,
                    d1,
                    e1,
                    f1,
                    g1,
                    h1,
                    j1,
                    k1,
                )
            )
            i = i + 2
    elif Eletype == 4.0:
        i = 0
        while i < nelem:
            num = i + num0
            linea_partida = myDataLines[i]
            nodes = tuple(map(int, linea_partida.split(",")))
            a = nodes[1]
            b = nodes[2]
            c = nodes[3]
            d = nodes[4]
            e = nodes[5]
            f = nodes[6]
            g = nodes[7]
            h = nodes[8]
            if i < 1:
                DataFile.writelines(
                    "%d, %d, %d, %d, %d, %d, %d, %d, %d" % (num, a, b, c, d, e, f, g, h)
                )
            elif i > 0:
                DataFile.writelines(
                    "\n%d, %d, %d, %d, %d, %d, %d, %d, %d"
                    % (num, a, b, c, d, e, f, g, h)
                )
            i = i + 1
    elif Eletype == 5.0:
        i = 0
        while i < nelem:
            num = i + num0
            linea_partida = myDataLines[i]
            nodes = tuple(map(int, linea_partida.split(",")))
            a = nodes[1]
            b = nodes[2]
            c = nodes[3]
            d = nodes[4]
            if i < 1:
                DataFile.writelines("%d, %d, %d, %d, %d" % (num, a, b, c, d))
            elif i > 0:
                DataFile.writelines("\n%d, %d, %d, %d, %d" % (num, a, b, c, d))
            i = i + 1

    DataFile.close()

    # Generate the second file (F3)
    f = open(INP)
    content = f.read()
    text = content.split(txtEnd)[1].split(End)[0].strip()
    final_file = open("F3.txt", "w")
    final_file.write(text)
    final_file.close()
    f.close()

    # Change DEPVAR in F3
    fileToSearch = "F3.txt"
    tempFile = open(fileToSearch, "r+")
    if Eletype < 4.0:
        for line in fileinput.input(fileToSearch):
            tempFile.write(
                line.replace(
                    "*User Material, constants=1",
                    "1, S11, S11\n2, S22, S22\n3, S33, S33\n4, S12, S12\n5, S13, S13\n6, S23, S23\n7, E11, E11\n8, E22, E22\n9, E33, E33\n10, E12, E12\n11, E13, E13\n12, E23, E23\n13, PHI, PHI\n14, H, H\n15, V, V\n16, JX, JX\n17, JY, JY\n18, JZ, JZ",
                )
            )
    elif Eletype > 3.0:
        for line in fileinput.input(fileToSearch):
            tempFile.write(
                line.replace(
                    "*User Material, constants=1",
                    "1, S11, S11\n2, S22, S22\n3, S33, S33\n4, S12, S12\n5, E11, E11\n6, E22, E22\n7, E33, E33\n8, E12, E12\n9, PHI, PHI\n10, H, H",
                )
            )

    tempFile.close()

    fileToSearch = "F3.txt"
    tempFile = open(fileToSearch, "r+")
    for line in fileinput.input(fileToSearch):
        tempFile.write(line.replace("0.,\n", "*User Material, constants=0\n"))
    tempFile.close()

    # Generate new input file
    DataF1 = open("F1.txt", "r+")
    F1 = DataF1.readlines()
    DataF1.close()

    DataF2 = open("F2.txt", "r+")
    F2 = DataF2.readlines()
    DataF2.close()

    DataF3 = open("F3.txt", "r+")
    F3 = DataF3.readlines()
    DataF3.close()

    DataFile = open(name_of_file + ".inp", "a")
    DataFile.writelines(Beg)
    DataFile.writelines("\n")
    DataFile.writelines(F1)
    DataFile.writelines("\n*UEL PROPERTY, ELSET=SOLID\n")
    DataFile.writelines(UEL01)
    DataFile.writelines(", ")
    DataFile.writelines(UEL02)
    DataFile.writelines(", ")
    DataFile.writelines(UEL03)
    DataFile.writelines(", ")
    DataFile.writelines(UEL04)
    DataFile.writelines(", ")
    DataFile.writelines(UEL05)
    DataFile.writelines(", ")
    DataFile.writelines(UEL06)
    DataFile.writelines(", ")
    DataFile.writelines(UEL07)
    DataFile.writelines(", ")
    DataFile.writelines(UEL08)
    DataFile.writelines(", ")
    DataFile.writelines("\n")
    DataFile.writelines(UEL09)
    DataFile.writelines(", ")
    DataFile.writelines(UEL010)
    DataFile.writelines(", ")
    DataFile.writelines(UEL011)
    DataFile.writelines(", ")
    DataFile.writelines(UEL012)
    DataFile.writelines(", ")
    DataFile.writelines(UEL013)

    if Eletype == 1.0:
        DataFile.writelines("\n*Element, type=C3D10, elset=Visualization\n")
    elif Eletype == 2.0:
        DataFile.writelines("\n*Element, type=C3D8, elset=Visualization\n")
    elif Eletype == 3.0:
        DataFile.writelines("\n*Element, type=C3D20R, elset=Visualization\n")
    elif Eletype == 4.0:
        DataFile.writelines("\n*Element, type=CPE8R, elset=Visualization\n")
    elif Eletype == 5.0:
        DataFile.writelines("\n*Element, type=CPE4, elset=Visualization\n")

    DataFile.writelines(F2)
    DataFile.writelines("\n")
    DataFile.writelines(txtEnd)
    DataFile.writelines("\n")
    DataFile.writelines(F3)
    DataFile.writelines("\n")
    DataFile.writelines(End)
    DataFile.writelines("\n")
    DataFile.close()

    # Change solid section
    fileToSearch = name_of_file + ".inp"
    tempFile = open(fileToSearch, "r+")
    if Eletype == 1.0:
        for line in fileinput.input(fileToSearch):
            tempFile.write(
                line.replace(
                    "*ELEMENT, TYPE=U1, ELSET=SOLID",
                    "*User element, nodes=10, type=U1, properties=13, coordinates=3, var=56\n1,2,3\n1,4\n1,9\n*ELEMENT, TYPE=U1, ELSET=SOLID",
                )
            )
    elif Eletype == 2.0:
        for line in fileinput.input(fileToSearch):
            tempFile.write(
                line.replace(
                    "*ELEMENT, TYPE=U1, ELSET=SOLID",
                    "*User element, nodes=8, type=U1, properties=13, coordinates=3, var=144\n1,2,3\n1,4\n1,9\n*ELEMENT, TYPE=U1, ELSET=SOLID",
                )
            )
    elif Eletype == 3.0:
        for line in fileinput.input(fileToSearch):
            tempFile.write(
                line.replace(
                    "*ELEMENT, TYPE=U1, ELSET=SOLID",
                    "*User element, nodes=20, type=U1, properties=13, coordinates=3, var=144\n1,2,3\n1,4\n1,9\n*ELEMENT, TYPE=U1, ELSET=SOLID",
                )
            )
    elif Eletype == 4.0:
        for line in fileinput.input(fileToSearch):
            tempFile.write(
                line.replace(
                    "*ELEMENT, TYPE=U1, ELSET=SOLID",
                    "*User element, nodes=8, type=U1, properties=13, coordinates=2, var=40\n1,2\n1,4\n*ELEMENT, TYPE=U1, ELSET=SOLID",
                )
            )
    elif Eletype == 5.0:
        for line in fileinput.input(fileToSearch):
            tempFile.write(
                line.replace(
                    "*ELEMENT, TYPE=U1, ELSET=SOLID",
                    "*User element, nodes=4, type=U1, properties=13, coordinates=2, var=40\n1,2\n1,4\n*ELEMENT, TYPE=U1, ELSET=SOLID",
                )
            )
    tempFile.close()

    fileToSearch = name_of_file + ".inp"
    tempFile = open(fileToSearch, "r+")
    for line in fileinput.input(fileToSearch):
        tempFile.write(
            line.replace(
                "*Solid Section, elset=Set-1, material=Material-1",
                "*Solid Section, elset=Visualization, material=Material-1",
            )
        )

    tempFile.close()

    # Eliminate existing files
    myfile = "F1.txt"
    if os.path.isfile(myfile):
        os.remove(myfile)

    myfile = "F2.txt"
    if os.path.isfile(myfile):
        os.remove(myfile)

    myfile = "F3.txt"
    if os.path.isfile(myfile):
        os.remove(myfile)

    myfile = "VisualElements.txt"
    if os.path.isfile(myfile):
        os.remove(myfile)

    myfile = "PF.inp"
    if os.path.isfile(myfile):
        os.remove(myfile)
    ##
    Mdb()
    mdb.ModelFromInputFile(name=name_of_file, inputFileName=name_of_file + ".inp")

    mdb.JobFromInputFile(
        name=name_of_file,
        inputFileName=name_of_file + ".inp",
        type=ANALYSIS,
        atTime=None,
        waitMinutes=0,
        waitHours=0,
        queue=None,
        memory=90,
        memoryUnits=PERCENTAGE,
        getMemoryFromAnalysis=True,
        explicitPrecision=SINGLE,
        nodalOutputPrecision=SINGLE,
        userSubroutine="UEL_piezoresistive_phasefield.for",
        scratch="",
        resultsFormat=ODB,
        multiprocessingMode=DEFAULT,
        numCpus=8,
        numDomains=8,
        numGPUs=0,
    )

    jobname = mdb.jobs[name_of_file]
    jobname.submit()
    jobname.waitForCompletion()
