%% Function to extract data from COMSOL

function [desRate] = funCOMSOL(DuCu,time)
import com.comsol.model.*
import com.comsol.model.util.*
DuCu(1) % Hydrogen diffusivity *1e-11, D [m^2/s] (called D_iso in the paper) in the current curve fitting iteration step
DuCu(2) % Diffusible hydrogen concentration, C_0 [mol/m^3] (called C_0iso in the paper) in the current curve fitting iteration step
model = mphopen('isothermalTDS.mph');
model.param.set('D',strcat(string(1e-11*DuCu(1)),' [m^2/s]'));
model.param.set('C_0',strcat(string(DuCu(2)),' [mol/m^3]'));
model.study('std1').run;
desRate = mphinterp(model,'-d(intop1(c),t)/(L/2)','coord',0,'t',time);
mphsave(model,'isothermalTDS.mph');
end


