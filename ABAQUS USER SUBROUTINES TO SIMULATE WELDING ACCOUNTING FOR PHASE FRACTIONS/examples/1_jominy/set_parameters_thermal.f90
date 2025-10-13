! Material model names
param%matnames(1) = 'BASE'
param%matnames(2) = 'WELD'

! Chemical composition Base + Weld
param%C = 0.07
param%Si = 0.235
param%Mn = 1.26
param%P = 0.015
param%V = 0.002
param%Ti = 0.021
param%Cr = 0.015
param%Ni = 0.02 
param%Mo = 0.002
param%As = 0.002 
param%Al = 0.024
param%Gsize_min = 10.0
param%Gsize_max = 100.0

! Change C content in weld
! note that the other alloying elements of the weld
! are similar to the base, as assigned above
param%C(2) = 0.14

! Transformation temperatures of the Base metal.
! Since we dont specify the temperatures of the Weld
! metal, they are calctulated by empirical equations.
param%Ae3(1) = 809.1d0
param%Ae1(1) = 711.5d0
param%Bs(1) = 618.0d0
param%Ms(1) = 408.4d0

! Absolute path to prop file
! Note that below is a linux path, a Windows path uses backward slash
param%propfile = '~/simulations/weld_model/examples/1_jominy/properties_thermal.inp'
