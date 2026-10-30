%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  FRACTURED_ROCK             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[TR_rock,TR_ch,Q_channel,FROCK,Qflow_rock]= FRACTURED_ROCK_BG(cellsize,N_nt,TR_rock,TR_ch,LEAK_nt,Q_channeltm1,FROCKtm1,FROCKtm11,SNn,dth,m_cell,n_cell,num_cell,Kres_Rock)
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Reference: 
%%%%%% OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% INPUT
%%% P [g/m2]: 8 dissolved P pool to rout, the last two are all DOC. 
%%% LEAK_nt [g/m2]: nutrient leakage to bedrock
%%% FROCK [mm] Volume in fractured rock vertical 
%%% Q_channel [mm] Volume in the channels 
%%% SNn [0/1] Stream network identifier 
%%% dth [h] time step 
%%% FROCKtm11 [mm] Volume in fractured rock vertical before the hydrology
%%% unit, to calculate nutrient concentration
%Kres_Rock [h] Aquifer constant 
%%%%% 

Q_channeltm1=reshape(Q_channeltm1,num_cell,1);
Qrock = dth*FROCKtm1/Kres_Rock; %%% [mm]   
FROCK = FROCKtm1 - Qrock; %%[mm]
Qrock_TOT = sum(Qrock); 
Q_channel = Q_channeltm1 + SNn*(Qrock_TOT/sum(SNn)); %%[mm]
Qflow_rock = Qrock_TOT;
Q_channel=reshape(Q_channel,m_cell,n_cell);

% %%%% update TR_rock, due to leakage
TR_rock= (TR_rock.*sum(FROCKtm11) +sum(LEAK_nt)*cellsize^2) / (sum(FROCKtm1)+1e-5); %%% TR_rock for next time step
TR_rock(isnan(TR_rock))=0;

%%%%%%%%%%%%%%% Sequence is inverse, should firstly have lk in, then add to
%%%%%%%%%%%%%%% river
%%%% update TR_ch, due to baseflow input, use TR_rock from last time step
TR_ch=(Q_channeltm1.*reshape(TR_ch,m_cell*n_cell,N_nt)+SNn*(Qrock_TOT/sum(SNn)).*TR_rock) ./ (reshape(Q_channel,num_cell,1)+1e-5);
TR_ch=reshape(TR_ch,m_cell,n_cell,N_nt);
TR_ch(isnan(TR_ch))=0;

% %%%% no need to update P, the P input is into the channel, not soil.

end 