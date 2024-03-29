from abaqus import *
from abaqusConstants import *
def MatDef(MethodSol, ModelName, MatName, E, xnu, xl, Gc, kflagS, kflagC, kflagD, kflagH, ft):

    if (kflagS=='Monolithic'):
        kflagS1=0
    elif (kflagS=='Staggered'):
        kflagS1=1
        
    if (kflagC=='AT2'):
        kflagC1=0
    elif (kflagC=='AT1'):
        kflagC1=1
    elif (kflagC=='PF-CZM (linear)'):
        kflagC1=2
    elif (kflagC=='PF-CZM (exp)'):
        kflagC1=3
        
    if (kflagD=='Amor'):
        kflagD1=1
    elif (kflagD=='Miehe'):
        kflagD1=2
    elif (kflagD=='No split'):
        kflagD1=0
        
    if (kflagH=='Hybrid method'):
        kflagH1=1
    elif (kflagH=='Anistropic elasticty matrix'):
        kflagH1=2
        

    if (MethodSol=='UMAT'):
        mdb.models[ModelName].Material(name=MatName)
        mdb.models[ModelName].materials[MatName].Depvar(n=1)
        mdb.models[ModelName].materials[MatName].UserMaterial(
            mechanicalConstants=(E, xnu, xl, Gc, kflagS1, kflagC1, kflagD1, kflagH1, ft))
        mdb.models[ModelName].materials[MatName].Conductivity(table=((1.0, ), 
            ))
    elif (MethodSol=='UMAT+HETVAL'):
        mdb.models[ModelName].Material(name=MatName)
        mdb.models[ModelName].materials[MatName].Depvar(n=7)
        mdb.models[ModelName].materials[MatName].UserMaterial(
            mechanicalConstants=(E, xnu, xl, Gc, kflagS1, kflagC1, kflagD1, kflagH1, ft))
        mdb.models[ModelName].materials[MatName].Conductivity(table=((1.0, ), 
            ))
        mdb.models[ModelName].materials[MatName].HeatGeneration()
        
