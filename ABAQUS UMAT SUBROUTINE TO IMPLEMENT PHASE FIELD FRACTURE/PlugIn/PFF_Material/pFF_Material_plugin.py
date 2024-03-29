from abaqusGui import *
from abaqusConstants import ALL
import osutils, os


###########################################################################
# Class definition
###########################################################################

class PFF_Material_plugin(AFXForm):

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def __init__(self, owner):
        
        # Construct the base class.
        #
        AFXForm.__init__(self, owner)
        self.radioButtonGroups = {}

        self.cmd = AFXGuiCommand(mode=self, method='MatDef',
            objectName='plug', registerQuery=False)
        pickedDefault = ''
        self.ModelNameKw = AFXStringKeyword(self.cmd, 'ModelName', True, 'Model-1')
        self.MethodSolKw = AFXStringKeyword(self.cmd, 'MethodSol', True, 'UMAT')
        self.MatNameKw = AFXStringKeyword(self.cmd, 'MatName', True, 'Material-1')
        self.EKw = AFXFloatKeyword(self.cmd, 'E', True)
        self.xnuKw = AFXFloatKeyword(self.cmd, 'xnu', True)
        self.xlKw = AFXFloatKeyword(self.cmd, 'xl', True)
        self.GcKw = AFXFloatKeyword(self.cmd, 'Gc', True)
        self.ftKw = AFXFloatKeyword(self.cmd, 'ft', True, 0)
        self.kflagSKw = AFXStringKeyword(self.cmd, 'kflagS', True, 'Monolithic')
        self.kflagCKw = AFXStringKeyword(self.cmd, 'kflagC', True, 'AT2')
        self.kflagDKw = AFXStringKeyword(self.cmd, 'kflagD', True, 'No split')
        self.kflagHKw = AFXStringKeyword(self.cmd, 'kflagH', True, 'Hybrid method')

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def getFirstDialog(self):

        import pFF_MaterialDB
        return pFF_MaterialDB.PFF_MaterialDB(self)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def doCustomChecks(self):

        # Try to set the appropriate radio button on. If the user did
        # not specify any buttons to be on, do nothing.
        #
        for kw1,kw2,d in self.radioButtonGroups.values():
            try:
                value = d[ kw1.getValue() ]
                kw2.setValue(value)
            except:
                pass
        return True

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def okToCancel(self):

        # No need to close the dialog when a file operation (such
        # as New or Open) or model change is executed.
        #
        return False

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Register the plug-in
#
thisPath = os.path.abspath(__file__)
thisDir = os.path.dirname(thisPath)

toolset = getAFXApp().getAFXMainWindow().getPluginToolset()
toolset.registerGuiMenuButton(
    buttonText='PFF_Material', 
    object=PFF_Material_plugin(toolset),
    messageId=AFXMode.ID_ACTIVATE,
    icon=None,
    kernelInitString='import plug',
    applicableModules=ALL,
    version='N/A',
    author='N/A',
    description='N/A',
    helpUrl='N/A'
)
