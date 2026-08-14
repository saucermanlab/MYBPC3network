function dydt=ODEfun(t,y,params) 
% Assign names for parameters 
[rpar,tau,ymax,speciesNames]=params{:}; 
Mass = 1; 
Eccentricity = 2; 
HCM_Mutation = 3; 
DCM_Mutation = 4; 
BAdrenergic_R = 5; 
MYBPC3_Mutation = 6; 
Ca_Sensitivity = 7; 
Calcium = 8; 
Active_F = 9; 
Stiffness_T = 10; 
CaMK = 11; 
CaN = 12; 
NFAT = 13; 
GBG = 14; 
Ras = 15; 
PI3K_AKT = 16; 
ERK12 = 17; 
Pressure_Overload = 18; 
Volume_Overload = 19; 
PLB = 20; 
PKA = 21; 
SERCA = 22; 
LCC = 23; 
RyR = 24; 
SR_LOAD = 25; 
Ca_Decay_T = 26; 
Ca_Amplitude = 27; 
Ca_Diastolic = 28; 
PKD = 29; 
PKG = 30; 
GF_R = 31; 
PKC = 32; 
Gq_R = 33; 
NO_BNP_CNP = 34; 
Gq = 35; 
FHL1 = 36; 
mTOR = 37; 
MEF2 = 38; 
Ca_Buffer = 39; 
RBM20 = 40; 
cMYBP_C = 41; 
dydt = zeros(41,1); 
dydt(Mass) = (OR(act(y(NFAT),rpar(:,38)),OR(act(y(ERK12),rpar(:,39)),OR(act(y(mTOR),rpar(:,75)),act(y(MEF2),rpar(:,76)))))*ymax(Mass) - y(Mass))/tau(Mass); 
dydt(Eccentricity) = (inhib(y(ERK12),rpar(:,77))*ymax(Eccentricity) - y(Eccentricity))/tau(Eccentricity); 
dydt(HCM_Mutation) = (rpar(1,1)*ymax(HCM_Mutation) - y(HCM_Mutation))/tau(HCM_Mutation); 
dydt(DCM_Mutation) = (rpar(1,2)*ymax(DCM_Mutation) - y(DCM_Mutation))/tau(DCM_Mutation); 
dydt(BAdrenergic_R) = (rpar(1,5)*ymax(BAdrenergic_R) - y(BAdrenergic_R))/tau(BAdrenergic_R); 
dydt(MYBPC3_Mutation) = (rpar(1,8)*ymax(MYBPC3_Mutation) - y(MYBPC3_Mutation))/tau(MYBPC3_Mutation); 
dydt(Ca_Sensitivity) = (OR(inhib(y(DCM_Mutation),rpar(:,14)),OR(act(y(HCM_Mutation),rpar(:,16)),OR(act(y(Stiffness_T),rpar(:,21)),inhib(y(cMYBP_C),rpar(:,90)))))*ymax(Ca_Sensitivity) - y(Ca_Sensitivity))/tau(Ca_Sensitivity); 
dydt(Calcium) = (OR(act(y(Ca_Diastolic),rpar(:,58)),OR(act(y(Ca_Amplitude),rpar(:,59)),OR(act(y(Ca_Decay_T),rpar(:,60)),act(y(Ca_Buffer),rpar(:,61)))))*ymax(Calcium) - y(Calcium))/tau(Calcium); 
dydt(Active_F) = (AND(rpar(:,17),act(y(Ca_Sensitivity),rpar(:,17)),act(y(Calcium),rpar(:,17)))*ymax(Active_F) - y(Active_F))/tau(Active_F); 
dydt(Stiffness_T) = (OR(act(y(Calcium),rpar(:,24)),OR(act(y(RBM20),rpar(:,79)),OR(act(y(PKC),rpar(:,80)),OR(inhib(y(PKA),rpar(:,81)),OR(inhib(y(PKD),rpar(:,82)),OR(inhib(y(CaMK),rpar(:,83)),OR(inhib(y(PKG),rpar(:,84)),AND(rpar(:,85),inhib(y(ERK12),rpar(:,85)),act(y(FHL1),rpar(:,85))))))))))*ymax(Stiffness_T) - y(Stiffness_T))/tau(Stiffness_T); 
dydt(CaMK) = (act(y(Calcium),rpar(:,22))*ymax(CaMK) - y(CaMK))/tau(CaMK); 
dydt(CaN) = (act(y(Calcium),rpar(:,23))*ymax(CaN) - y(CaN))/tau(CaN); 
dydt(NFAT) = (OR(AND(rpar(:,26),act(y(Stiffness_T),rpar(:,26)),act(y(CaN),rpar(:,26))),act(y(PI3K_AKT),rpar(:,35)))*ymax(NFAT) - y(NFAT))/tau(NFAT); 
dydt(GBG) = (OR(act(y(BAdrenergic_R),rpar(:,28)),act(y(Gq),rpar(:,30)))*ymax(GBG) - y(GBG))/tau(GBG); 
dydt(Ras) = (OR(AND(rpar(:,31),act(y(CaMK),rpar(:,31)),act(y(GBG),rpar(:,31)),inhib(y(Volume_Overload),rpar(:,31))),act(y(GF_R),rpar(:,36)))*ymax(Ras) - y(Ras))/tau(Ras); 
dydt(PI3K_AKT) = (OR(act(y(DCM_Mutation),rpar(:,15)),OR(act(y(GBG),rpar(:,32)),OR(act(y(GF_R),rpar(:,33)),act(y(PKD),rpar(:,64)))))*ymax(PI3K_AKT) - y(PI3K_AKT))/tau(PI3K_AKT); 
dydt(ERK12) = (AND(rpar(:,37),act(y(Ras),rpar(:,37)),act(y(FHL1),rpar(:,37)))*ymax(ERK12) - y(ERK12))/tau(ERK12); 
dydt(Pressure_Overload) = (rpar(1,3)*ymax(Pressure_Overload) - y(Pressure_Overload))/tau(Pressure_Overload); 
dydt(Volume_Overload) = (rpar(1,4)*ymax(Volume_Overload) - y(Volume_Overload))/tau(Volume_Overload); 
dydt(PLB) = (OR(inhib(y(CaMK),rpar(:,41)),OR(inhib(y(PKA),rpar(:,42)),inhib(y(PKG),rpar(:,43))))*ymax(PLB) - y(PLB))/tau(PLB); 
dydt(PKA) = (act(y(BAdrenergic_R),rpar(:,27))*ymax(PKA) - y(PKA))/tau(PKA); 
dydt(SERCA) = (inhib(y(PLB),rpar(:,44))*ymax(SERCA) - y(SERCA))/tau(SERCA); 
dydt(LCC) = (OR(act(y(CaMK),rpar(:,25)),act(y(PKA),rpar(:,45)))*ymax(LCC) - y(LCC))/tau(LCC); 
dydt(RyR) = (OR(act(y(Gq),rpar(:,29)),OR(act(y(GF_R),rpar(:,34)),act(y(CaMK),rpar(:,40))))*ymax(RyR) - y(RyR))/tau(RyR); 
dydt(SR_LOAD) = (OR(inhib(y(RyR),rpar(:,46)),act(y(SERCA),rpar(:,50)))*ymax(SR_LOAD) - y(SR_LOAD))/tau(SR_LOAD); 
dydt(Ca_Decay_T) = (OR(act(y(RyR),rpar(:,47)),OR(inhib(y(SERCA),rpar(:,51)),act(y(Ca_Sensitivity),rpar(:,55))))*ymax(Ca_Decay_T) - y(Ca_Decay_T))/tau(Ca_Decay_T); 
dydt(Ca_Amplitude) = (OR(act(y(RyR),rpar(:,49)),OR(act(y(LCC),rpar(:,53)),OR(act(y(SR_LOAD),rpar(:,54)),inhib(y(Ca_Buffer),rpar(:,57)))))*ymax(Ca_Amplitude) - y(Ca_Amplitude))/tau(Ca_Amplitude); 
dydt(Ca_Diastolic) = (OR(act(y(RyR),rpar(:,48)),OR(inhib(y(SERCA),rpar(:,52)),act(y(PKA),rpar(:,69))))*ymax(Ca_Diastolic) - y(Ca_Diastolic))/tau(Ca_Diastolic); 
dydt(PKD) = (act(y(PKC),rpar(:,63))*ymax(PKD) - y(PKD))/tau(PKD); 
dydt(PKG) = (act(y(NO_BNP_CNP),rpar(:,66))*ymax(PKG) - y(PKG))/tau(PKG); 
dydt(GF_R) = (OR(rpar(1,6),OR(act(y(Pressure_Overload),rpar(:,9)),act(y(Active_F),rpar(:,18))))*ymax(GF_R) - y(GF_R))/tau(GF_R); 
dydt(PKC) = (OR(AND(rpar(:,62),act(y(Calcium),rpar(:,62)),act(y(GF_R),rpar(:,62))),AND(rpar(:,70),act(y(Calcium),rpar(:,70)),act(y(Gq),rpar(:,70))))*ymax(PKC) - y(PKC))/tau(PKC); 
dydt(Gq_R) = (OR(rpar(1,7),OR(act(y(Pressure_Overload),rpar(:,10)),OR(act(y(Volume_Overload),rpar(:,12)),act(y(Active_F),rpar(:,19)))))*ymax(Gq_R) - y(Gq_R))/tau(Gq_R); 
dydt(NO_BNP_CNP) = (OR(act(y(Pressure_Overload),rpar(:,11)),OR(act(y(Volume_Overload),rpar(:,13)),act(y(Active_F),rpar(:,20))))*ymax(NO_BNP_CNP) - y(NO_BNP_CNP))/tau(NO_BNP_CNP); 
dydt(Gq) = (act(y(Gq_R),rpar(:,71))*ymax(Gq) - y(Gq))/tau(Gq); 
dydt(FHL1) = (OR(act(y(Stiffness_T),rpar(:,67)),OR(act(y(ERK12),rpar(:,68)),act(y(Gq),rpar(:,72))))*ymax(FHL1) - y(FHL1))/tau(FHL1); 
dydt(mTOR) = (OR(act(y(PI3K_AKT),rpar(:,73)),act(y(ERK12),rpar(:,74)))*ymax(mTOR) - y(mTOR))/tau(mTOR); 
dydt(MEF2) = (AND(rpar(:,65),act(y(CaMK),rpar(:,65)),act(y(PKD),rpar(:,65)))*ymax(MEF2) - y(MEF2))/tau(MEF2); 
dydt(Ca_Buffer) = (act(y(Ca_Sensitivity),rpar(:,56))*ymax(Ca_Buffer) - y(Ca_Buffer))/tau(Ca_Buffer); 
dydt(RBM20) = (act(y(mTOR),rpar(:,78))*ymax(RBM20) - y(RBM20))/tau(RBM20); 
dydt(cMYBP_C) = (OR(AND(rpar(:,86),inhib(y(MYBPC3_Mutation),rpar(:,86)),act(y(PKD),rpar(:,86))),OR(AND(rpar(:,87),inhib(y(MYBPC3_Mutation),rpar(:,87)),act(y(PKA),rpar(:,87))),OR(AND(rpar(:,88),inhib(y(MYBPC3_Mutation),rpar(:,88)),act(y(PKC),rpar(:,88))),AND(rpar(:,89),inhib(y(MYBPC3_Mutation),rpar(:,89)),act(y(CaMK),rpar(:,89))))))*ymax(cMYBP_C) - y(cMYBP_C))/tau(cMYBP_C); 

% utility functions 
function fact = act(x,rpar) 
% hill activation function with parameters w (weight), n (Hill coeff), EC50 
    w = rpar(1); 
    n = rpar(2); 
    EC50 = rpar(3); 
    beta = (EC50.^n - 1)./(2*EC50.^n - 1); 
    K = (beta - 1).^(1./n); 
    fact = w.*(beta.*x.^n)./(K.^n + x.^n); 
    if fact>w,                 % cap fact(x)<= 1 
        fact = w; 
    end
 
function finhib = inhib(x,rpar) 
% inverse hill function with parameters w (weight), n (Hill coeff), EC50 
    finhib = 1 - act(x,rpar);
 
function z = OR(x,y) 
% OR logic gate 
    z = 1-(1-x)*(1-y);
 
function z = AND(rpar,varargin) 
% AND logic gate, multiplying all of the reactants together 
    w = rpar(1); 
        v = real(cell2mat(varargin)); 
        h = (sum(v)/(nargin-1))^(nargin-2); 
    if w == 0 
        z = 0; 
    elseif h == 0 
        z = 0; 
    else 
        z = max(prod(v)/h,min(v)); 
    end 
