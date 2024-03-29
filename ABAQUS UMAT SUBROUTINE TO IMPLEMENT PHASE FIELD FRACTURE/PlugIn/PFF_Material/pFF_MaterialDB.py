from abaqusConstants import *
from abaqusGui import *
from kernelAccess import mdb, session
import os

thisPath = os.path.abspath(__file__)
thisDir = os.path.dirname(thisPath)


###########################################################################
# Class definition
###########################################################################

class PFF_MaterialDB(AFXDataDialog):

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def __init__(self, form):

        # Construct the base class.
        #

        AFXDataDialog.__init__(self, form, 'Material Definition Plug-in for Phase Field Fracture',
            self.OK|self.APPLY|self.CANCEL, DIALOG_ACTIONS_SEPARATOR)
            

        okBtn = self.getActionButton(self.ID_CLICKED_OK)
        okBtn.setText('OK')
            

        applyBtn = self.getActionButton(self.ID_CLICKED_APPLY)
        applyBtn.setText('Apply')
            
        HFrame_3 = FXHorizontalFrame(p=self, opts=0, x=0, y=0, w=0, h=0,
            pl=0, pr=0, pt=0, pb=0)
        VFrame_2 = FXVerticalFrame(p=HFrame_3, opts=0, x=0, y=0, w=0, h=0,
            pl=0, pr=0, pt=0, pb=0)
        AFXTextField(p=VFrame_2, ncols=30, labelText='Model Name:', tgt=form.ModelNameKw, sel=0)
        ComboBox_1 = AFXComboBox(p=VFrame_2, ncols=0, nvis=1, text='UMAT or UMAT+HETVAL code:', tgt=form.MethodSolKw, sel=0)
        ComboBox_1.setMaxVisible(10)
        ComboBox_1.appendItem(text='UMAT')
        ComboBox_1.appendItem(text='UMAT+HETVAL')
        AFXTextField(p=VFrame_2, ncols=30, labelText='Material name:', tgt=form.MatNameKw, sel=0)
        l = FXLabel(p=VFrame_2, text='Warning: The UMAT version is applicable for Abaqus 2020 and later versions.', opts=JUSTIFY_LEFT)
        l.setFont( getAFXFont(FONT_BOLD) )
        AFXTextField(p=VFrame_2, ncols=12, labelText='Young s modulus (E):', tgt=form.EKw, sel=0)
        AFXTextField(p=VFrame_2, ncols=12, labelText='Poisson s ratio (xnu):', tgt=form.xnuKw, sel=0)
        AFXTextField(p=VFrame_2, ncols=12, labelText='Length scale (xl):', tgt=form.xlKw, sel=0)
        AFXTextField(p=VFrame_2, ncols=12, labelText='Toughness (Gc):', tgt=form.GcKw, sel=0)
        AFXTextField(p=VFrame_2, ncols=12, labelText='Tensile strength (ft):', tgt=form.ftKw, sel=0)
        ComboBox_2 = AFXComboBox(p=VFrame_2, ncols=0, nvis=1, text='Solution scheme (kflagS):', tgt=form.kflagSKw, sel=0)
        ComboBox_2.setMaxVisible(10)
        ComboBox_2.appendItem(text='Monolithic')
        ComboBox_2.appendItem(text='Staggered')
        ComboBox_3 = AFXComboBox(p=VFrame_2, ncols=0, nvis=1, text='Constitutive model (kflagC):', tgt=form.kflagCKw, sel=0)
        ComboBox_3.setMaxVisible(10)
        ComboBox_3.appendItem(text='AT1')
        ComboBox_3.appendItem(text='AT2')
        ComboBox_3.appendItem(text='PF-CZM (linear)')
        ComboBox_3.appendItem(text='PF-CZM (exp)')
        ComboBox_4 = AFXComboBox(p=VFrame_2, ncols=0, nvis=1, text='Strain energy split method (kflagD):', tgt=form.kflagDKw, sel=0)
        ComboBox_4.setMaxVisible(10)
        ComboBox_4.appendItem(text='No split')
        ComboBox_4.appendItem(text='Amor')
        ComboBox_4.appendItem(text='Miehe')
        ComboBox_5 = AFXComboBox(p=VFrame_2, ncols=0, nvis=1, text='Solution method (kflagH):', tgt=form.kflagHKw, sel=0)
        ComboBox_5.setMaxVisible(10)
        ComboBox_5.appendItem(text='Hybrid method')
        ComboBox_5.appendItem(text='Anistropic elasticty matrix')
        l = FXLabel(p=VFrame_2, text='Note: use anistropic elasticty matrix only with Mieh s split.', opts=JUSTIFY_LEFT)
        l.setFont( getAFXFont(FONT_BOLD) )
        VFrame_1 = FXVerticalFrame(p=HFrame_3, opts=0, x=0, y=0, w=0, h=0,
            pl=100, pr=0, pt=0, pb=0)
        fileName = os.path.join(thisDir, 'phi.png')
        icon = afxCreatePNGIcon(fileName)
        FXLabel(p=VFrame_1, text='', ic=icon)
        l = FXLabel(p=self, text='If using these codes for research or industrial purposes, please cite:  ', opts=JUSTIFY_LEFT)
        l.setFont( getAFXFont(FONT_BOLD) )
        l = FXLabel(p=self, text='Navidtehrani, Y.; Beteg\xf3n, C.; Mart\xednez-Pa\xf1eda, E. A simple and robust Abaqus implementation of the phase field fracture method, Applications in Engineering Science 6: 100050 (2021).', opts=JUSTIFY_LEFT)
        l = FXLabel(p=self, text='Navidtehrani, Y.; Beteg\xf3n, C.; Mart\xednez-Pa\xf1eda, E. A Unified Abaqus Implementation of the Phase Field Fracture Method Using Only a User Material Subroutine. Materials 14(8): 1913 (2021).', opts=JUSTIFY_LEFT)
