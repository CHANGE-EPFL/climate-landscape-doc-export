%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction BIOGEOCHEMISTRY_DYNAMIC   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[LEAK_nt]= Biogeo_mass_leak(B,Lk,Vtotal,BiogeoPar)  
%%%%%%%%%%%%%% Leakage of Nitrogen Compound 
%%%% Proportional to leakage from the nitrogen zone with fraction of
%%%% dissolved ammonium and nitrate 
%%% For simplicity leakage only at the base of soil column 
%%%
VSUM=Vtotal(:);
Lk=Lk(:);


%%%% Porporato et al., 2003 
aNH4 =BiogeoPar.aNH4;  
aNO3 = BiogeoPar.aNO3;  
aP = BiogeoPar.aP;
aK = BiogeoPar.aK; 
aDON= BiogeoPar.aDON;
aDOP= BiogeoPar.aDOP;
aDOC = BiogeoPar.aDOC; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% we are calculating the mass
% here we may have a unit problem?? this is
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% actually [gN/ gsoil], should *(ZBIOG*rsd)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% to get g/m2
LEAK_NH4 = zeros(size(VSUM));
LEAK_NO3 = zeros(size(VSUM));
LEAK_P = zeros(size(VSUM));
LEAK_K = zeros(size(VSUM));
LEAK_DON = zeros(size(VSUM));
LEAK_DOP = zeros(size(VSUM));
LEAK_DOC = zeros(size(VSUM));

valid = isfinite(VSUM) & VSUM > 0;
if any(valid)
LEAK_NH4(valid) = max(0,min(B(valid,1),aNH4*B(valid,1).*(Lk(valid))./(VSUM(valid)))); %% [gN/m2]
LEAK_NO3(valid) = max(0,min(B(valid,2),aNO3*B(valid,2).*(Lk(valid))./(VSUM(valid)))); %% [gN/m2]
LEAK_P(valid) = max(0,min(B(valid,3),aP*B(valid,3).*(Lk(valid))./(VSUM(valid)))); %% [gP/m2]
LEAK_K(valid) = max(0,min(B(valid,4),aK*B(valid,4).*(Lk(valid))./(VSUM(valid)))); %% [gK/m2]
LEAK_DON(valid) = max(0,min(B(valid,5),aDON*B(valid,5).*(Lk(valid))./(VSUM(valid)))); %% [gN/m2]
LEAK_DOP(valid) = max(0,min(B(valid,6),aDOP*B(valid,6).*(Lk(valid))./(VSUM(valid)))); %% [gP/m2]
LEAK_DOC(valid) = max(0,min(B(valid,7)+B(valid,8),aDOC*(B(valid,7)+B(valid,8)).*(Lk(valid))./(VSUM(valid)))); %% [gC/m2]
end
LEAK_nt=[LEAK_NH4,LEAK_NO3,LEAK_P,LEAK_K,LEAK_DON,LEAK_DOP,LEAK_DOC];












return 
