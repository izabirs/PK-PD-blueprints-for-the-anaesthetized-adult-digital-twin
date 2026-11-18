function m = pkEleveldRemifentanil(patient)
%PKELEVELDREMIFENTANIL  Remifentanil PK/PD, allometric covariate model.
%
%   Source: Eleveld DJ, Proost JH, Vereecke H, Absalom AR, Olofsen E,
%   Vuyk J, Struys MMRF. An allometric model of remifentanil
%   pharmacokinetics and pharmacodynamics.
%   Anesthesiology 2017;126:1005-1018.  doi:10.1097/ALN.0000000000001634
%
%   Volumes in L, clearances in L/min, rate constants in 1/min.

age = patient.age;  wgt = patient.weight;  bmi = patient.bmi;
sex = double(upper(char(patient.sex)) == 'M');

AGE_ref = 35;  WGT_ref = 70;  HGT_ref = 1.7;  BMI_ref = WGT_ref/HGT_ref^2;
faging = @(x) exp(x*(age - AGE_ref));
fsig   = @(x,c50,g) x.^g ./ (c50.^g + x.^g);

SIZE = alSallamiFFM(sex, wgt, age, bmi) / alSallamiFFM(1, WGT_ref, AGE_ref, BMI_ref);
th   = [2.88, -0.00554, -0.00327, -0.0315, 0.470, -0.0260];

KMAT    = fsig(wgt,     th(1), 2);
KMATref = fsig(WGT_ref, th(1), 2);
if sex
    KSEX = 1;
else
    KSEX = 1 + th(5)*fsig(age,12,6)*(1 - fsig(age,45,6));
end

m.V1 = 5.81  * SIZE * faging(th(2));
m.V2 = 8.882 * SIZE * faging(th(3));
m.V3 = 5.03  * SIZE * faging(th(4)) * exp(th(6)*(wgt - WGT_ref));
m.CL = 2.58  * SIZE^0.75 * (KMAT/KMATref) * KSEX * faging(th(3));
m.Q2 = 1.72  * (m.V2/8.882)^0.75 * faging(th(2)) * KSEX;
m.Q3 = 0.124 * (m.V3/5.03)^0.75  * faging(th(2));
m.ke0 = 1.09 * faging(-0.0289);

m.omega2 = struct('V1',0.104,'V2',0.115,'V3',0.810, ...
                  'CL',0.0197,'Q2',0.0547,'Q3',0.285,'ke0',1.26);
m.source = 'Eleveld 2017, Anesthesiology 126:1005-1018';
end
