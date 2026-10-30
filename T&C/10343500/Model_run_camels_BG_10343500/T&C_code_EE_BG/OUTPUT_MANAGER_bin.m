%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% OUTPUT WRITING   %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if t==2
    Qlat_in_tgtm1 = 0;
    q_runon_tgtm1 = 0;
    Q_channel_tgtm1 = 0;
    QE_cum = 0;
    H_cum = 0;
    Rn_cum = 0; 
    LAI_H_cum = LAI_Htm1;
    LAI_L_cum  =LAI_Ltm1;
    NPP_L_cum = 0; 
    NPP_H_cum = 0; 
    RA_L_cum = 0; 
    RA_H_cum = 0; 
    T_L_cum = 0; 
    T_H_cum = 0; 
    G_cum = 0; 
    GPP_L_cum = 0; 
    GPP_H_cum =0;
    SND_cum = 0; 
    

else
    V_tgtm1 = V_tg;
    Vice_tgtm1 = Vice_tg;
    SWE_tgtm1 = SWE_tg;
    In_tgtm1 =In_tg;
    ICE_tgtm1 =  ICE_tg; %%
    WAT_tgtm1 =  WAT_tg; %%
    FROCK_tgtm1 = FROCK_tg; %%
    %%%%---
    Qlat_in_tgtm1 = Qlat_in_tg;
    q_runon_tgtm1 = q_runon_tg;
    Q_channel_tgtm1 = Q_channel_tg;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% SPATIAL AVERAGE over the watershed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 


%%%%%
Pr_tg = mean(Pr_S(Kinde)); %%% [mm]
Ta_tg = mean(Ta_S(Kinde)); %%% [°C]
Ws_tg = mean(Ws_S(Kinde)); %%% [m/s]
Ds_tg = mean(Ds_S(Kinde)); %%% [°C]
ea_tg = mean(ea_S(Kinde)); %%% [Pa]
N_tg = mean(N_S(Kinde)); %%% [-]
Pre_tg = mean(Pre_S(Kinde)); %%% [mbar]
Tdew_tg = mean(Tdew_S(Kinde)); %%% [°C]
Rsw_space = SAB1_S(:)+ SAB2_S(:) + SAD1_S(:)+ SAD2_S(:);
Rsw_tg = mean(Rsw_space(Kinde)); %%% [W/m^2]
PAR_space = PARB_S(:) + PARD_S(:); %%% %% [W/m^2]
PAR_tg= mean(PAR_space(Kinde));%% [W/m^2]
Ca_tg = mean(Ca_S(Kinde)); %% [ppm]
%%%%%
Ts_tg  =  mean(Ts(Kinde));%% [°C]
Tdamp_tg  = mean(Tdamp(Kinde));%% [°C]
Csno_tg = mean(Csno(Kinde)); %%[]
Cice_tg = mean(Cice(Kinde)); %%[]
Csnow_tg = mean(Csnow(Kinde)); %%[]
Cicew_tg = mean(Cicew(Kinde)); %%[]
Pr_sno_tg  =  mean(Pr_sno(Kinde));%%[mm]
Pr_liq_tg = mean(Pr_liq(Kinde)); %%[mm]
%%%%
Rn_tg=mean(Rn(Kinde));%% [W/m^2]
H_tg=mean(H(Kinde));%% [W/m^2]
G_tg=mean(G(Kinde));%% [W/m^2]
Gfin_tg=mean(Gfin(Kinde)); % [W/m^2]
QE_tg=mean(QE(Kinde));%% [W/m^2]
Qv_tg = mean(Qv(Kinde)); %% [W/m^2]
Qfm_tg =  mean(Qfm(Kinde));% [W/m^2]
%%%%
SWE_tg = mean(SWE(Kinde)); %% [mm]
SND_tg = mean(SND(Kinde)); %% [m]
WR_SP_tg = mean(WR_SP(Kinde)); %%[mm]
U_SWE_tg = mean(U_SWE(Kinde)); %%[mm]
NIn_SWE_tg  =  mean(NIn_SWE(Kinde));%% [mm]
dw_SNO_tg=mean(dw_SNO(Kinde)); %%[]
ros_tg = mean(ros(Kinde)); %% [kg/m^3]
In_SWE_tg =  mean(In_SWE(Kinde)); %%[mm]
SP_wc_tg =  mean(SP_wc(Kinde));%%[mm]
%%%%
ICE_tg = mean(ICE(Kinde)); %% [mm]
ICE_D_tg = mean(ICE_D(Kinde)); %% [m]
WR_IP_tg = mean(WR_IP(Kinde)); %%[mm]
IP_wc_tg =  mean(IP_wc(Kinde));%%[mm]
NIce_tg  =  mean(NIce(Kinde));%% [mm]
%%%
Imelt_tg = mean(Imelt(Kinde)); %% [mm]
Smelt_tg = mean(Smelt(Kinde)); %% [mm]
Tice_tg =  mean(Tice(Kinde)); %% [C]
%%%%%%%%%
T_H_space = sum(T_H,2);
T_L_space = sum(T_L,2);
T_H_tg = mean(T_H_space(Kinde)); %%%[mm/h]
T_L_tg = mean(T_L_space(Kinde)); %%%[mm/h]
T_tg= T_L_tg +T_H_tg  ; %%% [mm/h]
EIn_H_space = sum(EIn_H,2);
EIn_L_space = sum(EIn_L,2);
EIn_L_tg  = mean(EIn_L_space(Kinde)); %%% [mm/h]
EIn_H_tg =  mean(EIn_H_space(Kinde)); %%% [mm/h]
EIn_tg =EIn_H_tg+ EIn_L_tg; %%% [mm/h]
EG_tg =mean(EG(Kinde)); %%% [mm/h]
ESN_tg =mean(ESN(Kinde)+ESN_In(Kinde)); %%% [mm/h]
EWAT_tg = mean(EWAT(Kinde)); %%[mm/h]
EICE_tg = mean(EICE(Kinde)); %%[mm/h]
Dr_H_space = sum(Dr_H,2);
Dr_L_space = sum(Dr_L,2);
Dr_L_tg  = mean(Dr_L_space(Kinde)); %%% [mm]
Dr_H_tg =  mean(Dr_H_space(Kinde)); %%% [mm]
%%%%%%%%%%%%%%%%
EIn_urb_tg = mean(EIn_urb(Kinde)); %%[mm/h]
EIn_rock_tg = mean(EIn_rock(Kinde)); %%[mm/h]
SE_rock_tg = mean(SE_rock(Kinde)); %%[]
SE_urb_tg = mean(SE_urb(Kinde)); %%[]
%%%
In_H_space = sum(In_H,2);
In_L_space = sum(In_L,2);
In_urb_tg =  mean(In_urb(Kinde)); %%[mm]
In_rock_tg =  mean(In_rock(Kinde)); %%[mm]
In_tg = mean(In_H_space(Kinde) +  In_L_space(Kinde) +  SP_wc(Kinde) + In_SWE(Kinde) + In_urb(Kinde) + In_rock(Kinde)+ IP_wc(Kinde) ); %%% [mm]
Inveg_tg = mean(In_H_space(Kinde) +  In_L_space(Kinde));
%%%
WAT_tg= mean(WAT(Kinde)) ;%%%[mm]
FROCK_tg = mean(FROCK(Kinde)); %%%[mm]
%%%
f_tg = mean(f(Kinde)*dth); %%[mm]
WIS_tg = mean(WIS(Kinde)); %%[mm]
Rd_tg =mean(Rd(Kinde)); %% [mm]
Rh_tg =mean(Rh(Kinde)); %%[mm]
Lk_tg =mean(Lk(Kinde)*dth); %% [mm]
Lk_wat_tg =mean(Lk_wat(Kinde)*dth); %% [mm]
Lk_rock_tg =mean(Lk_rock(Kinde)*dth); %% [mm]
%%%
OF_tg = mean(OF(Kinde)); %%[]
OS_tg = mean(OS(Kinde)); %%[]
ZWT_tg= mean(ZWT(Kinde)); %% [mm]
Fract_sat_tg = sum(ZWT(Kinde)==0)/num_cell;%% [-]
%%%%%%%%
Qi_in_space = sum(Qi_in*dth,2);
Qi_out_space = sum(Qi_out*dth,2);
Qlat_in_tg=  mean(Qi_in_space(Kinde)); %% [mm]
Qlat_out_tg =  mean(Qi_out_space(Kinde)); %% [mm]
q_runon_tg = mean(q_runon(Kinde)*dth); %% [mm]
Q_channel_tg =  mean(Q_channel(Kinde)); %% [mm]
%%%%%%%%%
V_space =  sum(V,2); %%%
Vice_space =  sum(Vice,2); %%%
O_space =  V_space./Zs_OUT + Ohy_OUT ; %%[-]
V_tg = mean(Asur(Kinde).*V_space(Kinde)); %% [mm]
Vice_tg = mean(Asur(Kinde).*Vice_space(Kinde)); %% [mm]
O_tg = mean(O_space(Kinde)); %% [-]
%%%%%%
ra_tg  =  mean(ra(Kinde));%% [s/m]
r_soil_tg  =  mean(r_soil(Kinde));%% [s/m]
alp_soil_tg = mean(alp_soil(Kinde)); %%[-]
%%%%
Tdp_space =  mean(Tdp,2); %%%
Tdp_tg = mean(Tdp_space(Kinde)); %% % [°C]
%%%%
er_tg =  mean(er(Kinde)); %%  % [kg/s m^2]
%%%
TsVEG_tg  =  mean(TsVEG(Kinde));%% [°C]
DQ_S_tg = mean(DQ_S(Kinde)); %%[]
DT_S_tg = mean(DT_S(Kinde)); %%[]
dQ_S_tg = mean(dQ_S(Kinde)); %%[]
CK1_tg=  mean(CK1(Kinde));%% [mm]
%CK2_tg =  mean(CK2(Kinde));%% [mm]
%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% SIMPLE MASS CHECK CONTROL
if sum(SNn)>=1
    CKt = (V_tgtm1 - V_tg) + (Vice_tgtm1 - Vice_tg) + Pr_tg - EG_tg*dth - T_tg*dth - EIn_tg*dth ...
        - ESN_tg*dth - EIn_urb_tg*dth - EWAT_tg*dth - EIn_rock_tg*dth - EICE_tg*dth ...
        + (SWE_tgtm1 -SWE_tg) + (In_tgtm1 -In_tg) ...
        +  (ICE_tgtm1 -ICE_tg) +  (WAT_tgtm1 -WAT_tg)  +  (FROCK_tgtm1 -FROCK_tg) ...
        + Qlat_in_tgtm1 + q_runon_tgtm1 + Q_channel_tgtm1  ...
        - Qlat_in_tg - Q_exit - Qsub_exit - q_runon_tg - Q_channel_tg -Swe_exit;
else
    CKt = (V_tgtm1 - V_tg) + (Vice_tgtm1 - Vice_tg) + Pr_tg - EG_tg*dth - T_tg*dth - EIn_tg*dth ...
        - ESN_tg*dth - EIn_urb_tg*dth - EWAT_tg*dth - EIn_rock_tg*dth - EICE_tg*dth ...
        + (SWE_tgtm1 -SWE_tg) + (In_tgtm1 -In_tg) ...
        +  (ICE_tgtm1 -ICE_tg) +  (WAT_tgtm1 -WAT_tg)  +  (FROCK_tgtm1 -FROCK_tg) ...
        + Qlat_in_tgtm1 + q_runon_tgtm1   ...
        - Qlat_in_tg - Q_exit - Qsub_exit - q_runon_tg  -Swe_exit;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ckk1=  Pr_tg  + q_runon_tgtm1 - EIn_tg*dth - WIS_tg + (In_tgtm1 -In_tg); %%Not correct
%Ckk2 = WIS_tg - f_tg - Rh_tg;
%Ckk3 = Rh_tg - Q_exit - q_runon_tg;
%Ckk4 = f_tg + (V_tgtm1 - V_tg) - EG_tg*dth - T_tg*dth - Lk_tg ...
%  - Qlat_in_tg -Rd_tg + Qlat_in_tgtm1 - Qsub_exit ;
%Ckk5 = Pr_sno_tg*dth + Pr_liq_tg*dth*Csno_tg ...
%    - ESN_tg*dth + (SWE_tgtm1 -SWE_tg) ...
%    + (In_tgtm1 -In_tg) -sum(Dr_H_tg)*Csno_tg ; %% Not correct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if t==2
%     tit{1}=strcat([Out.OutDir, '/OUTPUT_',TITLE_SAVE,'_AVG.dat']);
%     fid(1)=fopen(tit{1},'a');
% end


SIM_AVG(t,	1	)=	Pr_tg	;
SIM_AVG(t,	2	)=	Ta_tg	;
SIM_AVG(t,	3	)=	Ws_tg	;
SIM_AVG(t,	4	)=	Ds_tg	;
SIM_AVG(t,	5	)=	ea_tg	;
SIM_AVG(t,	6	)=	N_tg	;
SIM_AVG(t,	7	)=	Pre_tg	;
SIM_AVG(t,	8	)=	Tdew_tg	;
SIM_AVG(t,	9	)=	Rsw_tg	;
SIM_AVG(t,	10	)=	PAR_tg	;
SIM_AVG(t,	11	)=	Ca_tg	;
SIM_AVG(t,	12	)=	Ts_tg	;
SIM_AVG(t,	13	)=	Tdamp_tg	;
SIM_AVG(t,	14	)=	Csno_tg	;
SIM_AVG(t,	15	)=	Cice_tg	;
SIM_AVG(t,	16	)=	Csnow_tg	;
SIM_AVG(t,	17	)=	Cicew_tg	;
SIM_AVG(t,	18	)=	Pr_sno_tg	;
SIM_AVG(t,	19	)=	Pr_liq_tg	;
SIM_AVG(t,	20	)=	Rn_tg	;
SIM_AVG(t,	21	)=	H_tg	;
SIM_AVG(t,	22	)=	G_tg	;
SIM_AVG(t,	23	)=	Gfin_tg	;
SIM_AVG(t,	24	)=	QE_tg	;
SIM_AVG(t,	25	)=	Qv_tg	;
SIM_AVG(t,	26	)=	Qfm_tg	;
SIM_AVG(t,	27	)=	SWE_tg	;
SIM_AVG(t,	28	)=	SND_tg	;
SIM_AVG(t,	29	)=	WR_SP_tg	;
SIM_AVG(t,	30	)=	dw_SNO_tg	;
SIM_AVG(t,	31	)=	ros_tg	;
SIM_AVG(t,	32	)=	In_SWE_tg	;
SIM_AVG(t,	33	)=	SP_wc_tg	;
SIM_AVG(t,	34	)=	ICE_tg	;
SIM_AVG(t,	35	)=	ICE_D_tg	;
SIM_AVG(t,	36	)=	WR_IP_tg	;
SIM_AVG(t,	37	)=	NIce_tg	;
SIM_AVG(t,	38	)=	IP_wc_tg	;
SIM_AVG(t,	39	)=	T_H_tg	;
SIM_AVG(t,	40	)=	T_L_tg 	;
SIM_AVG(t,	41	)=	EIn_H_tg	;
SIM_AVG(t,	42	)=	EIn_L_tg	;
SIM_AVG(t,	43	)=	EG_tg	;
SIM_AVG(t,	44	)=	ESN_tg	;
SIM_AVG(t,	45	)=	EWAT_tg	;
SIM_AVG(t,	46	)=	EICE_tg	;
SIM_AVG(t,	47	)=	Dr_L_tg	;
SIM_AVG(t,	48	)=	Dr_H_tg	;
SIM_AVG(t,	49	)=	EIn_urb_tg	;
SIM_AVG(t,	50	)=	EIn_rock_tg	;
SIM_AVG(t,	51	)=	In_tg	;
SIM_AVG(t,	52	)=	Inveg_tg	;
SIM_AVG(t,	53	)=	WAT_tg	;
SIM_AVG(t,	54	)=	FROCK_tg	;
SIM_AVG(t,	55	)=	f_tg	;
SIM_AVG(t,	56	)=	WIS_tg	;
SIM_AVG(t,	57	)=	Rd_tg	;
SIM_AVG(t,	58	)=	Rh_tg	;
SIM_AVG(t,	59	)=	Lk_tg	;
SIM_AVG(t,	60	)=	Lk_wat_tg	;
SIM_AVG(t,	61	)=	Lk_rock_tg	;
SIM_AVG(t,	62	)=	OF_tg	;
SIM_AVG(t,	63	)=	OS_tg	;
SIM_AVG(t,	64	)=	ZWT_tg	;
SIM_AVG(t,	65	)=	Fract_sat_tg	;
SIM_AVG(t,	66	)=	Qlat_in_tg	;
SIM_AVG(t,	67	)=	Qlat_out_tg	;
SIM_AVG(t,	68	)=	q_runon_tg	;
SIM_AVG(t,	69	)=	Q_channel_tg	;
SIM_AVG(t,	70	)=	V_tg	;
SIM_AVG(t,	71	)=	O_tg	;
SIM_AVG(t,	72	)=	Q_exit	;
SIM_AVG(t,	73	)=	Qsub_exit	;
SIM_AVG(t,	74	)=	Swe_exit	;
SIM_AVG(t,	75	)=	r_soil_tg	;
SIM_AVG(t,	76	)=	alp_soil_tg	;
SIM_AVG(t,	77	)=	ra_tg 	;
SIM_AVG(t,	78	)=	Tdp_tg 	;
SIM_AVG(t,	79	)=	er_tg	;
SIM_AVG(t,	80	)=	TsVEG_tg	;
SIM_AVG(t,	81	)=	DQ_S_tg	;
SIM_AVG(t,	82	)=	DT_S_tg	;
SIM_AVG(t,	83	)=	dQ_S_tg	;
SIM_AVG(t,	84	)=	Imelt_tg	;
SIM_AVG(t,	85	)=	Smelt_tg	;
SIM_AVG(t,	86	)=	Tice_tg	;
SIM_AVG(t,	87	)=	Vice_tg	;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% SPATIAL STANDARD DEVIATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
std_Pr_tg = std(Pr_S(Kinde)); %%% [mm]
std_Ta_tg = std(Ta_S(Kinde)); %%% [°C]
std_Ws_tg = std(Ws_S(Kinde)); %%% [m/s]
std_Ds_tg = std(Ds_S(Kinde)); %%% [°C]
std_ea_tg = std(ea_S(Kinde)); %%% [Pa]
std_N_tg = std(N_S(Kinde)); %%% [-]
std_Pre_tg = std(Pre_S(Kinde)); %%% [mbar]
std_Tdew_tg = std(Tdew_S(Kinde)); %%% [°C]
std_Rsw_tg = std(Rsw_space(Kinde)); %%% [W/m^2]
std_PAR_tg= std(PAR_space(Kinde));%% [W/m^2]
std_Ca_tg= std(Ca_S(Kinde));%% [W/m^2]
std_Ts_tg  =  std(Ts(Kinde));%% [°C]
std_Tdamp_tg  = std(Tdamp(Kinde));%% [°C]
std_Pr_sno_tg  =  std(Pr_sno(Kinde));%% [mm]
std_Pr_liq_tg = std(Pr_liq(Kinde)); %%[mm]
%%%%
%%%%
std_Rn_tg=std(Rn(Kinde));%% [W/m^2]
std_H_tg=std(H(Kinde));%% [W/m^2]
std_G_tg=std(G(Kinde));%% [W/m^2]
std_Gfin_tg=std(Gfin(Kinde));%% [W/m^2]
std_QE_tg=std(QE(Kinde));%% [W/m^2]
std_Qv_tg = std(Qv(Kinde)); %% [W/m^2]
std_Qfm_tg =  std(Qfm(Kinde));% [W/m^2]
%%%%
std_SWE_tg = std(SWE(Kinde)); %% [mm]
std_SND_tg = std(SND(Kinde)); %% [m]
std_WR_SP_tg = std(WR_SP(Kinde)); %%[]
std_U_SWE_tg = std(U_SWE(Kinde)); %%[]
std_NIn_SWE_tg  =  std(NIn_SWE(Kinde));%% [m/s]
std_dw_SNO_tg=std(dw_SNO(Kinde)); %%[]
std_ros_tg = std(ros(Kinde)); %% [kg/m^3]
std_In_SWE_tg =  std(In_SWE(Kinde)); %%[mm]
std_SP_wc_tg =  std(SP_wc(Kinde));%%[mm]
%%%%
std_ICE_tg = std(ICE(Kinde)); %% [mm]
std_ICE_D_tg = std(ICE_D(Kinde)); %% [m]
std_WR_IP_tg = std(WR_IP(Kinde)); %%[]
std_NIce_tg  =  std(NIce(Kinde));%%
std_IP_wc_tg =  std(IP_wc(Kinde));%%[mm]
%%%
std_T_H_tg = std(T_H_space(Kinde)); %%%[mm/h]
std_T_L_tg = std(T_L_space(Kinde)); %%%[mm/h]
std_EIn_L_tg  = std(EIn_L_space(Kinde)); %%% [mm/h]
std_EIn_H_tg =  std(EIn_H_space(Kinde)); %%% [mm/h]
std_EG_tg =std(EG(Kinde)); %%% [mm/h]
std_ESN_tg =std(ESN(Kinde)+ESN_In(Kinde)); %%% [mm/h]
std_EWAT_tg = std(EWAT(Kinde)); %%[mm/h]
std_EICE_tg = std(EICE(Kinde)); %%[mm/h]
%%%%%%%%%%%%%%%%
std_EIn_urb_tg = std(EIn_urb(Kinde)); %%[mm/h]
std_EIn_rock_tg = std(EIn_rock(Kinde)); %%[mm/h]
std_In_tg = std(In_H_space(Kinde) +  In_L_space(Kinde) +  SP_wc(Kinde) + In_SWE(Kinde) +  In_urb(Kinde) + In_rock(Kinde)+ IP_wc(Kinde) ); %%% [mm]
std_Inveg_tg = std(In_H_space(Kinde) +  In_L_space(Kinde));
%%%
std_WAT_tg = std(WAT(Kinde)); %% [mm]
std_FROCK_tg = std(FROCK(Kinde)); %% [mm]
%%%
std_f_tg = std(f(Kinde)*dth); %%[mm]
std_WIS_tg = std(WIS(Kinde)); %%[mm]
std_Rd_tg =std(Rd(Kinde)); %% [mm]
std_Rh_tg =std(Rh(Kinde)); %%[mm]
std_Lk_wat_tg =std(Lk_wat(Kinde)*dth); %% [mm]
std_Lk_rock_tg =std(Lk_rock(Kinde)*dth); %% [mm]
std_Lk_tg =std(Lk(Kinde)*dth); %% [mm]
std_OF_tg = std(OF(Kinde)); %%[]
std_OS_tg = std(OS(Kinde)); %%[]
std_ZWT_tg= std(ZWT(Kinde)); %% [mm]
%%%%%%%%
std_Qlat_in_tg=  std(Qi_in_space(Kinde)); %% [mm]
std_Qlat_out_tg =  std(Qi_out_space(Kinde)); %% [mm]
std_q_runon_tg = std(q_runon(Kinde)*dth); %% [mm]
std_Q_channel_tg = std(Q_channel(Kinde)); %% [mm]
%%%%%%%%%
std_O_tg =std(O_space(Kinde)); %%
std_V_tg =std(V_space(Kinde)); %%
std_Tdp_tg = std(Tdp_space(Kinde));
std_TsVEG_tg  = std(TsVEG(Kinde));%% [°C]
%%%%%%%%
std_er_tg =  std(er(Kinde)); %%  % [kg/s m^2]
std_r_soil_tg  =  std(r_soil(Kinde));%% [°C]
std_alp_soil_tg = std(alp_soil(Kinde)); %%[]
std_ra_tg  =  std(ra(Kinde));%% [m/s]
%%%%%%%%%



SIM_STD(t,	1	)= std_Pr_tg	;
SIM_STD(t,	2	)= std_Ta_tg	;
SIM_STD(t,	3	)= std_Ws_tg	;
SIM_STD(t,	4	)= std_Ds_tg	;
SIM_STD(t,	5	)= std_ea_tg	;
SIM_STD(t,	6	)= std_N_tg	;
SIM_STD(t,	7	)= std_Pre_tg	;
SIM_STD(t,	8	)= std_Tdew_tg	;
SIM_STD(t,	9	)= std_Rsw_tg	;
SIM_STD(t,	10	)= std_PAR_tg	;
SIM_STD(t,	11	)= std_Ca_tg	;
SIM_STD(t,	12	)= std_Ts_tg	;
SIM_STD(t,	13	)= std_Tdamp_tg	;
SIM_STD(t,	14	)= std_Pr_sno_tg	;
SIM_STD(t,	15	)= std_Pr_liq_tg	;
SIM_STD(t,	16	)= std_Rn_tg	;
SIM_STD(t,	17	)= std_H_tg	;
SIM_STD(t,	18	)= std_G_tg	;
SIM_STD(t,	19	)= std_Gfin_tg	;
SIM_STD(t,	20	)= std_QE_tg	;
SIM_STD(t,	21	)= std_Qv_tg	;
SIM_STD(t,	22	)= std_Qfm_tg	;
SIM_STD(t,	23	)= std_SWE_tg	;
SIM_STD(t,	24	)= std_SND_tg	;
SIM_STD(t,	25	)= std_WR_SP_tg	;
SIM_STD(t,	26	)= std_dw_SNO_tg	;
SIM_STD(t,	27	)= std_ros_tg	;
SIM_STD(t,	28	)= std_In_SWE_tg	;
SIM_STD(t,	29	)= std_SP_wc_tg	;
SIM_STD(t,	30	)= std_ICE_tg	;
SIM_STD(t,	31	)= std_ICE_D_tg	;
SIM_STD(t,	32	)= std_WR_IP_tg	;
SIM_STD(t,	33	)= std_NIce_tg	;
SIM_STD(t,	34	)= std_IP_wc_tg	;
SIM_STD(t,	35	)= std_T_H_tg	;
SIM_STD(t,	36	)= std_T_L_tg 	;
SIM_STD(t,	37	)= std_EIn_H_tg	;
SIM_STD(t,	38	)= std_EIn_L_tg	;
SIM_STD(t,	39	)= std_EG_tg	;
SIM_STD(t,	40	)= std_ESN_tg	;
SIM_STD(t,	41	)= std_EWAT_tg	;
SIM_STD(t,	42	)= std_EICE_tg	;
SIM_STD(t,	43	)= std_EIn_urb_tg	;
SIM_STD(t,	44	)= std_EIn_rock_tg	;
SIM_STD(t,	45	)= std_In_tg	;
SIM_STD(t,	46	)= std_Inveg_tg	;
SIM_STD(t,	47	)= std_WAT_tg	;
SIM_STD(t,	48	)= std_FROCK_tg	;
SIM_STD(t,	49	)= std_f_tg	;
SIM_STD(t,	50	)= std_WIS_tg	;
SIM_STD(t,	51	)= std_Rd_tg	;
SIM_STD(t,	52	)= std_Rh_tg	;
SIM_STD(t,	53	)= std_Lk_tg	;
SIM_STD(t,	54	)= std_Lk_wat_tg	;
SIM_STD(t,	55	)= std_Lk_rock_tg	;
SIM_STD(t,	56	)= std_OF_tg	;
SIM_STD(t,	57	)= std_OS_tg	;
SIM_STD(t,	58	)= std_ZWT_tg	;
SIM_STD(t,	59	)= std_Qlat_in_tg	;
SIM_STD(t,	60	)= std_Qlat_out_tg	;
SIM_STD(t,	61	)= std_q_runon_tg	;
SIM_STD(t,	62	)= std_Q_channel_tg	;
SIM_STD(t,	63	)= std_V_tg	;
SIM_STD(t,	64	)= std_O_tg	;
SIM_STD(t,	65	)= std_r_soil_tg	;
SIM_STD(t,	66	)= std_alp_soil_tg	;
SIM_STD(t,	67	)= std_ra_tg 	;
SIM_STD(t,	68	)= std_Tdp_tg 	;
SIM_STD(t,	69	)= std_er_tg	;
SIM_STD(t,	70	)= std_TsVEG_tg	;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


OPT_VAR_SM_ON=0;


if OPT_VAR_SM_ON == 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% SPATIAL VARIABILITY OF SOIL MOISTURE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% V, f in n-coordinate
    %%% O_tg O_tg_prec dth O_space O_space_prec  Qi_in_prec  std_O_tg
    f_ver_tg = mean(Asur(Kinde).*f(Kinde)*dth); %%[mm]
    %%%
    T_H_tg_n = mean(T_H_space(Kinde)./Asur(Kinde)); %%%[mm/h]
    T_L_tg_n = mean(T_L_space(Kinde)./Asur(Kinde)); %%%[mm/h]
    EG_tg_n =mean(EG(Kinde)./Asur(Kinde)); %%% [mm/h]
    Rd_tg_n =mean(Rd(Kinde)./Asur(Kinde)); %% [mm]
    Lk_tg_n =mean(Lk(Kinde)./Asur(Kinde)*dth); %% [mm]
    Qlat_out_tg_n =  mean(Qi_out_space(Kinde)./Asur(Kinde)); %% [mm]
    %%%%%%%%%%%
    V_space_prec = sum(V_prec,2);
    O_space_prec =  V_space_prec./Zs_OUT + Ohy_OUT ; %%[-]
    %V_tg_prec = mean(Asur(Kinde).*V_space_prec(Kinde)); %% [mm]
    O_tg_prec = mean(O_space_prec(Kinde)); %% [-]
    Qi_in_space_prec = sum(Qi_in_prec*dth,2);
    Qlat_in_tg_prec =  mean(Qi_in_space_prec(Kinde)); %% [mm]
    %%%
    Qlat_in_tg_prec_n=  mean(Qi_in_space_prec(Kinde)./Asur(Kinde)); %% [mm]
    %%%%%
    Qi_in_prec = Qi_in ; %%%
    V_prec = V ;
    %V_space =  sum(V,2); %%%
    %O_space =  V_space./Zs_OUT + Ohy_OUT ; %%[-]
    %V_tg = mean(Asur(Kinde).*V_space(Kinde)); %% [mm]
    %O_tg = mean(O_space(Kinde)); %% [-]
    %%%%%%%%%% Computation in n-coordinate
    Te1 = (O_space-O_space_prec)/dth - (O_tg - O_tg_prec)/dth; %%[1/h]
    Te2 = 2*(O_space-O_tg); %% [-]
    dOspacep2dt = mean(Te2(Kinde).*Te1(Kinde));
    %dOspacedt2 =  mean( (1./Zs_OUT(Kinde)).*(Te2(Kinde).*(Asur.*f(Kinde)-f_ver_tg/dth) - Te2(Kinde).*(EG(Kinde)-EG_tg) - Te2(Kinde).*(T_H_space(Kinde)+T_L_space(Kinde)-T_H_tg-T_L_tg) + ...
    %     -Te2(Kinde).*(Lk(Kinde)-Lk_tg/dth) + Te2(Kinde).*(Qi_in_space_prec(Kinde)/dth-Qlat_in_tg_prec/dth) - Te2(Kinde).*(Qi_out_space(Kinde)/dth-Qlat_out_tg/dth) +...
    %    -Te2(Kinde).*(Rd(Kinde)/dth-Rd_tg/dth)) );
    %%%%%%%%%
    %Cov_f  = mean( (1./Zs_OUT(Kinde)).*( Te2(Kinde).*(f(Kinde)-f_tg/dth) ));
    %Cov_Eg =  mean( (1./Zs_OUT(Kinde)).*( Te2(Kinde).*(EG(Kinde)-EG_tg) ));
    %Cov_T =  mean( (1./Zs_OUT(Kinde)).*( Te2(Kinde).*(T_H_space(Kinde)+T_L_space(Kinde)-T_H_tg-T_L_tg) ));
    %Cov_Lk = mean( (1./Zs_OUT(Kinde)).*( Te2(Kinde).*(Lk(Kinde)-Lk_tg/dth) ));
    %Cov_Qin = mean( (1./Zs_OUT(Kinde)).*( Te2(Kinde).*(Qi_in_space_prec(Kinde)/dth-Qlat_in_tg_prec/dth) ));
    %Cov_Qout = mean( (1./Zs_OUT(Kinde)).*( Te2(Kinde).*(Qi_out_space(Kinde)/dth-Qlat_out_tg/dth) ));
    %Cov_Rd = mean( (1./Zs_OUT(Kinde)).*( Te2(Kinde).*(Rd(Kinde)/dth-Rd_tg/dth) ));
    %%%%%%%
    Cov_f  = mean( (1./Zs_OUT(Kinde)).*( Te2(Kinde).*(f(Kinde)-f_tg/dth) ));
    Cov_Eg =  mean( (1./Zs_OUT(Kinde)).*( Te2(Kinde).*(EG(Kinde)-EG_tg_n) ));
    Cov_T =  mean( (1./Zs_OUT(Kinde)).*( Te2(Kinde).*(T_H_space(Kinde)./Asur+T_L_space(Kinde)./Asur-T_H_tg_n-T_L_tg_n) ));
    Cov_Lk = mean( (1./Zs_OUT(Kinde)).*( Te2(Kinde).*(Lk(Kinde)./Asur-Lk_tg_n/dth) ));
    Cov_Qin = mean( (1./Zs_OUT(Kinde)).*( Te2(Kinde).*(Qi_in_space_prec(Kinde)/dth./Asur-Qlat_in_tg_prec_n/dth) ));
    Cov_Qout = mean( (1./Zs_OUT(Kinde)).*( Te2(Kinde).*(Qi_out_space(Kinde)/dth./Asur-Qlat_out_tg_n/dth) ));
    Cov_Rd = mean( (1./Zs_OUT(Kinde)).*( Te2(Kinde).*(Rd(Kinde)/dth./Asur-Rd_tg_n/dth) ));
    %%%%%%
    dOspacep2dt2 = Cov_f - Cov_Eg - Cov_T - Cov_Lk + Cov_Qin - Cov_Qout - Cov_Rd;
    A_var = Cov_f - Cov_Eg  - Cov_Lk + Cov_Qin - Cov_Qout - Cov_Rd;
    B_var = Cov_T;
    %%%%%%%%%%%%%%%% %% Budget in v - Coordinate
    dOspacedt =  (V_tg-V_tgtm1) ;
    dOspacedt2 = f_ver_tg  - EG_tg*dth - T_tg*dth - Lk_tg +  Qlat_in_tgtm1 - Qlat_out_tg - Rd_tg ;
    Ckk4 = (V_tgtm1 - V_tg) + f_ver_tg  - EG_tg*dth - T_tg*dth - Lk_tg +  Qlat_in_tgtm1 + (- Qlat_in_tg  +  Qlat_in_tg) - Qlat_out_tg - Rd_tg ;
    %Ckk44 = f_tg + (V_tgtm1 - V_tg) - EG_tg*dth - T_tg*dth - Lk_tg - Qlat_in_tg -Rd_tg + Qlat_in_tgtm1 - Qsub_exit ;
    A_mu =  mean(1./Zs_OUT(Kinde)).*(f_ver_tg  - EG_tg*dth - Lk_tg +  Qlat_in_tgtm1 - Qlat_out_tg - Rd_tg);
    B_mu =  mean(1./Zs_OUT(Kinde)).*(T_tg*dth);
    %%%%%%%%%%%%%
    OCv_tg = std_O_tg/O_tg;
    dCv_1 =  (A_var/(2*O_tg*std_O_tg));
    dCv_2 =  A_mu*OCv_tg/O_tg;
    dCv_3 =  B_mu*OCv_tg/O_tg;
    dCv_4 =  (B_var/(2*O_tg*std_O_tg));
    dCvdt = dCv_1 - dCv_2 + dCv_3 - dCv_4;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    deltaTSO = (O_space - O_tg)./O_tg ; %%% [-] Mean Relative difference
    deltaTSO2 = (deltaTSO).^2 ; %%% [-] Square of Mean Relative difference
    %%%%Ohy_OUT Osat_OUT
    SEspace = (O_space-Ohy_OUT)./(Osat_OUT-Ohy_OUT);
    SE_tg = mean( SEspace(Kinde) ); %% Effective saturation [-]
    std_SE_tg =std(SEspace(Kinde)); %%
    deltaTSS = (SEspace - SE_tg)./SE_tg ; %%% [-] Mean Relative difference
    deltaTSS2 = (deltaTSS).^2 ; %%% [-] Square of Mean Relative difference
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [WL_Mdep1]=Evaporation_layers(Zs,Mdep1); %%% Weighiting layers
    OMdep1 = (O*WL_Mdep1');%% Weighting layers 1
    OMdep1_tg = mean(OMdep1(Kinde));
    std_OMdep1_tg =std(OMdep1(Kinde)); %%
    [WL_Mdep2]=Evaporation_layers(Zs,Mdep2); %%% Weighiting layers
    OMdep2 = (O*WL_Mdep2');%% Weighting layers 2
    OMdep2_tg = mean(OMdep2(Kinde));
    std_OMdep2_tg =std(OMdep2(Kinde)); %%
    [WL_Mdep3]=Evaporation_layers(Zs,Mdep3); %%% Weighiting layers
    OMdep3 = (O*WL_Mdep3');%% Weighting layers 3
    OMdep3_tg = mean(OMdep3(Kinde));
    std_OMdep3_tg =std(OMdep3(Kinde)); %%
    %%%%%%%%%%%%%%%%

    SIM_VAR_SM(t,	1	)=	dOspacep2dt	;
    SIM_VAR_SM(t,	2	)=	Cov_f	;
    SIM_VAR_SM(t,	3	)=	Cov_Eg	;
    SIM_VAR_SM(t,	4	)=	Cov_T	;
    SIM_VAR_SM(t,	5	)=	Cov_Lk	;
    SIM_VAR_SM(t,	6	)=	Cov_Qin	;
    SIM_VAR_SM(t,	7	)=	Cov_Qout	;
    SIM_VAR_SM(t,	8	)=	Cov_Rd	;
    SIM_VAR_SM(t,	9	)=	dOspacep2dt2	;
    SIM_VAR_SM(t,	10	)=	A_var	;
    SIM_VAR_SM(t,	11	)=	B_var	;
    SIM_VAR_SM(t,	12	)=	dOspacedt	;
    SIM_VAR_SM(t,	13	)=	dOspacedt2	;
    SIM_VAR_SM(t,	14	)=	Ckk4	;
    SIM_VAR_SM(t,	15	)=	A_mu	;
    SIM_VAR_SM(t,	16	)=	B_mu	;
    SIM_VAR_SM(t,	17	)=	dCv_1	;
    SIM_VAR_SM(t,	18	)=	dCv_2	;
    SIM_VAR_SM(t,	19	)=	dCv_3	;
    SIM_VAR_SM(t,	20	)=	dCv_4	;
    SIM_VAR_SM(t,	21	)=	dCvdt	;
    SIM_VAR_SM(t,	22	)=	SE_tg	;
    SIM_VAR_SM(t,	23	)=	std_SE_tg	;
    SIM_VAR_SM(t,	24	)=	OMdep1_tg	;
    SIM_VAR_SM(t,	25	)=	std_OMdep1_tg	;
    SIM_VAR_SM(t,	26	)=	OMdep2_tg	;
    SIM_VAR_SM(t,	27	)=	std_OMdep2_tg	;
    SIM_VAR_SM(t,	28	)=	OMdep3_tg	;
    SIM_VAR_SM(t,	29	)=	std_OMdep3_tg	;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% SPATIAL AVERAGE over the Vegetation type
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ijki=1:cc_max

    %%%%%%
    OH_tg(1,ijki) =  mean(OH(:,ijki)); %%[]
    std_OH_tg(1,ijki) = std(OH(:,ijki)); %%[]
    OL_tg(1,ijki) =  mean(OL(:,ijki)); %%[]
    std_OL_tg(1,ijki) = std(OL(:,ijki)); %%[]
    An_H_tg(1,ijki) =  mean(An_H(:,ijki)); %%[]
    std_An_H_tg(1,ijki) = std(An_H(:,ijki)); %%[]
    An_L_tg(1,ijki) =  mean(An_L(:,ijki)); %%[]
    std_An_L_tg(1,ijki) = std(An_L(:,ijki)); %%[]
    Tdp_H_tg(1,ijki) =  mean(Tdp_H(:,ijki)); %%[]
    std_Tdp_H_tg(1,ijki) = std(Tdp_H(:,ijki)); %%[]
    Tdp_L_tg(1,ijki) =  mean(Tdp_L(:,ijki)); %%[]
    std_Tdp_L_tg(1,ijki) = std(Tdp_L(:,ijki)); %%[]
    Rdark_H_tg(1,ijki) =  mean(Rdark_H(:,ijki)); %%[]
    std_Rdark_H_tg(1,ijki) = std(Rdark_H(:,ijki)); %%[]
    Rdark_L_tg(1,ijki) =  mean(Rdark_L(:,ijki)); %%[]
    std_Rdark_L_tg(1,ijki) = std(Rdark_L(:,ijki)); %%[]
    LAI_H_tg(1,ijki) =  mean(LAI_H(:,ijki)); %%[]
    std_LAI_H_tg(1,ijki) = std(LAI_H(:,ijki)); %%[]
    LAI_L_tg(1,ijki) =  mean(LAI_L(:,ijki)); %%[]
    std_LAI_L_tg(1,ijki) = std(LAI_L(:,ijki)); %%[]
    NPP_H_tg(1,ijki) =  mean(NPP_H(:,ijki)); %%[]
    std_NPP_H_tg(1,ijki) = std(NPP_H(:,ijki)); %%[]
    NPP_L_tg(1,ijki) =  mean(NPP_L(:,ijki)); %%[]
    std_NPP_L_tg(1,ijki) = std(NPP_L(:,ijki)); %%[]
    ANPP_H_tg(1,ijki) =  mean(ANPP_H(:,ijki)); %%[]
    std_ANPP_H_tg(1,ijki) = std(ANPP_H(:,ijki)); %%[]
    ANPP_L_tg(1,ijki) =  mean(ANPP_L(:,ijki)); %%[]
    std_ANPP_L_tg(1,ijki) = std(ANPP_L(:,ijki)); %%[]
    RA_H_tg(1,ijki) =  mean(RA_H(:,ijki)); %%[]
    std_RA_H_tg(1,ijki) = std(RA_H(:,ijki)); %%[]
    RA_L_tg(1,ijki) =  mean(RA_L(:,ijki)); %%[]
    std_RA_L_tg(1,ijki) = std(RA_L(:,ijki)); %%[]
    Rg_H_tg(1,ijki) =  mean(Rg_H(:,ijki)); %%[]
    std_Rg_H_tg(1,ijki) = std(Rg_H(:,ijki)); %%[]
    Rg_L_tg(1,ijki) =  mean(Rg_L(:,ijki)); %%[]
    std_Rg_L_tg(1,ijki) = std(Rg_L(:,ijki)); %%[]
    LAIdead_H_tg(1,ijki) =  mean(LAIdead_H(:,ijki)); %%[]
    std_LAIdead_H_tg(1,ijki) = std(LAIdead_H(:,ijki)); %%[]
    LAIdead_L_tg(1,ijki) =  mean(LAIdead_L(:,ijki)); %%[]
    std_LAIdead_L_tg(1,ijki) = std(LAIdead_L(:,ijki)); %%[]
    hc_H_tg(1,ijki) =  mean(hc_H(:,ijki)); %%[]
    std_hc_H_tg(1,ijki) = std(hc_H(:,ijki)); %%[]
    hc_L_tg(1,ijki) =  mean(hc_L(:,ijki)); %%[]
    std_hc_L_tg(1,ijki) = std(hc_L(:,ijki)); %%[]
    AgeL_H_tg(1,ijki) =  mean(AgeL_H(:,ijki)); %%[]
    std_AgeL_H_tg(1,ijki) = std(AgeL_H(:,ijki)); %%[]
    AgeL_L_tg(1,ijki) =  mean(AgeL_L(:,ijki)); %%[]
    std_AgeL_L_tg(1,ijki) = std(AgeL_L(:,ijki)); %%[]
    SAI_H_tg(1,ijki) =  mean(SAI_H(:,ijki)); %%[]
    std_SAI_H_tg(1,ijki) = std(SAI_H(:,ijki)); %%[]
    SAI_L_tg(1,ijki) =  mean(SAI_L(:,ijki)); %%[]
    std_SAI_L_tg(1,ijki) = std(SAI_L(:,ijki)); %%[]
    PHE_S_H_tg(1,ijki) =  mean(PHE_S_H(:,ijki)); %%[]
    std_PHE_S_H_tg(1,ijki) = std(PHE_S_H(:,ijki)); %%[]
    PHE_S_L_tg(1,ijki) =  mean(PHE_S_L(:,ijki)); %%[]
    std_PHE_S_L_tg(1,ijki) = std(PHE_S_L(:,ijki)); %%[]


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ijki=1:cc_max

    PFT_AVG(t,	1	,ijki)=	OH_tg(1,	ijki)	;
    PFT_AVG(t,	2	,ijki)=	OL_tg(1,	ijki)	;
    PFT_AVG(t,	3	,ijki)=	An_H_tg(1,	ijki)	;
    PFT_AVG(t,	4	,ijki)=	An_L_tg(1,	ijki)	;
    PFT_AVG(t,	5	,ijki)=	Rdark_H_tg(1,	ijki)	;
    PFT_AVG(t,	6	,ijki)=	Rdark_L_tg(1,	ijki)	;
    PFT_AVG(t,	7	,ijki)=	Tdp_H_tg(1,	ijki)	;
    PFT_AVG(t,	8	,ijki)=	Tdp_L_tg(1,	ijki)	;
    PFT_AVG(t,	9	,ijki)=	LAI_H_tg(1,	ijki)	;
    PFT_AVG(t,	10	,ijki)=	LAI_L_tg(1,	ijki)	;
    PFT_AVG(t,	11	,ijki)=	NPP_H_tg(1,	ijki)	;
    PFT_AVG(t,	12	,ijki)=	NPP_L_tg(1,	ijki)	;
    PFT_AVG(t,	13	,ijki)=	ANPP_H_tg(1,	ijki)	;
    PFT_AVG(t,	14	,ijki)=	ANPP_L_tg(1,	ijki)	;
    PFT_AVG(t,	15	,ijki)=	RA_H_tg(1,	ijki)	;
    PFT_AVG(t,	16	,ijki)=	RA_L_tg(1,	ijki)	;
    PFT_AVG(t,	17	,ijki)=	Rg_H_tg(1,	ijki)	;
    PFT_AVG(t,	18	,ijki)=	Rg_L_tg(1,	ijki)	;
    PFT_AVG(t,	19	,ijki)=	PHE_S_H_tg(1,	ijki)	;
    PFT_AVG(t,	20	,ijki)=	PHE_S_L_tg(1,	ijki)	;
    PFT_AVG(t,	21	,ijki)=	AgeL_H_tg(1,	ijki)	;
    PFT_AVG(t,	22	,ijki)=	AgeL_L_tg(1,	ijki)	;
    PFT_AVG(t,	23	,ijki)=	SAI_H_tg(1,	ijki)	;
    PFT_AVG(t,	24	,ijki)=	SAI_L_tg(1,	ijki)	;
    PFT_AVG(t,	25	,ijki)=	LAIdead_H_tg(1,	ijki)	;
    PFT_AVG(t,	26	,ijki)=	LAIdead_L_tg(1,	ijki)	;
    PFT_AVG(t,	27	,ijki)=	hc_H_tg(1,	ijki)	;
    PFT_AVG(t,	28	,ijki)=	hc_L_tg(1,	ijki)	;

    PFT_STD(t,	1	,ijki)=	std_OH_tg(1,	ijki)	;
    PFT_STD(t,	2	,ijki)=	std_OL_tg(1,	ijki)	;
    PFT_STD(t,	3	,ijki)=	std_An_H_tg(1,	ijki)	;
    PFT_STD(t,	4	,ijki)=	std_An_L_tg(1,	ijki)	;
    PFT_STD(t,	5	,ijki)=	std_Rdark_H_tg(1,	ijki)	;
    PFT_STD(t,	6	,ijki)=	std_Rdark_L_tg(1,	ijki)	;
    PFT_STD(t,	7	,ijki)=	std_Tdp_H_tg(1,	ijki)	;
    PFT_STD(t,	8	,ijki)=	std_Tdp_L_tg(1,	ijki)	;
    PFT_STD(t,	9	,ijki)=	std_LAI_H_tg(1,	ijki)	;
    PFT_STD(t,	10	,ijki)=	std_LAI_L_tg(1,	ijki)	;
    PFT_STD(t,	11	,ijki)=	std_NPP_H_tg(1,	ijki)	;
    PFT_STD(t,	12	,ijki)=	std_NPP_L_tg(1,	ijki)	;
    PFT_STD(t,	13	,ijki)=	std_ANPP_H_tg(1,	ijki)	;
    PFT_STD(t,	14	,ijki)=	std_ANPP_L_tg(1,	ijki)	;
    PFT_STD(t,	15	,ijki)=	std_RA_H_tg(1,	ijki)	;
    PFT_STD(t,	16	,ijki)=	std_RA_L_tg(1,	ijki)	;
    PFT_STD(t,	17	,ijki)=	std_Rg_H_tg(1,	ijki)	;
    PFT_STD(t,	18	,ijki)=	std_Rg_L_tg(1,	ijki)	;
    PFT_STD(t,	19	,ijki)=	std_PHE_S_H_tg(1,	ijki)	;
    PFT_STD(t,	20	,ijki)=	std_PHE_S_L_tg(1,	ijki)	;
    PFT_STD(t,	21	,ijki)=	std_AgeL_H_tg(1,	ijki)	;
    PFT_STD(t,	22	,ijki)=	std_AgeL_L_tg(1,	ijki)	;
    PFT_STD(t,	23	,ijki)=	std_SAI_H_tg(1,	ijki)	;
    PFT_STD(t,	24	,ijki)=	std_SAI_L_tg(1,	ijki)	;
    PFT_STD(t,	25	,ijki)=	std_LAIdead_H_tg(1,	ijki)	;
    PFT_STD(t,	26	,ijki)=	std_LAIdead_L_tg(1,	ijki)	;
    PFT_STD(t,	27	,ijki)=	std_hc_H_tg(1,	ijki)	;
    PFT_STD(t,	28	,ijki)=	std_hc_L_tg(1,	ijki)	;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% TEMPORAL AVERAGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if t==2
    toutp=2;
else
    toutp=toutp+1;
end
%%%
if toutp==2
    %%%%%
    Pr_spatial=Pr_S;
    Ta_spatial=Ta_S;
    Ws_spatial=Ws_S;
    Ds_spatial=Ds_S;
    ea_spatial=ea_S;
    N_spatial=N_S;
    Pre_spatial=Pre_S;
    Tdew_spatial=Tdew_S;
    Rsw_spatial=Rsw_space;
    PAR_spatial=PAR_space;
    Ca_spatial = Ca_S;
    %%%
    Ts_spatial=Ts;
    Tdamp_spatial=Tdamp;
    Csno_spatial=Csno;
    Cice_spatial=Cice;
    Csnow_spatial=Csnow;
    Cicew_spatial=Cicew;
    Pr_sno_spatial=Pr_sno;
    Pr_liq_spatial=Pr_liq;
    %%%%
    Rn_spatial=Rn;
    H_spatial=H;
    G_spatial=G;
    Gfin_spatial=Gfin;
    QE_spatial=QE;
    Qv_spatial=Qv;
    Qfm_spatial=Qfm;
    %%%%
    SWE_spatial = SWE;
    SND_spatial = SND ;
    WR_SP_spatial = WR_SP;
    U_SWE_spatial = U_SWE;
    dw_SNO_spatial = dw_SNO;
    ros_spatial = ros;
    In_SWE_spatial=  In_SWE;
    SP_wc_spatial =  SP_wc ;
    %SWE_avalanched_spatial = SWE_avalanched;
    %%%%
    ICE_spatial = ICE;
    ICE_D_spatial = ICE_D ;
    WR_IP_spatial = WR_IP;
    NICe_spatial = NIce;
    IP_wc_spatial =  IP_wc ;
    %%%%
    Imelt_spatial = Imelt;
    Smelt_spatial = Smelt;
    Tice_spatial = Tice;
    T_H_spatial=T_H_space;
    T_L_spatial=T_L_space;
    EIn_H_spatial=EIn_H_space;
    EIn_L_spatial=EIn_L_space;
    EG_spatial = EG;
    ESN_spatial = ESN + ESN_In;
    EWAT_spatial = EWAT ;
    EICE_spatial = EICE ;
    Dr_H_spatial=Dr_H_space;
    Dr_L_spatial=Dr_L_space;
    EIn_rock_spatial = EIn_rock;
    %%%%
    In_spatial =  In_H_space + In_L_space + SP_wc + In_SWE  + In_urb + In_rock + IP_wc;
    Inveg_spatial = In_H_space + In_L_space;
    In_rock_spatial = In_rock;
    WAT_spatial = WAT;
    FROCK_spatial = FROCK ;
    %%%%%
    f_spatial = f;
    WIS_spatial = WIS;
    Rd_spatial =Rd;
    Rh_spatial =Rh;
    Lk_spatial =Lk;
    Lk_rock_spatial =Lk_rock;
    Lk_wat_spatial =Lk_wat;
    OF_spatial =OF;
    OS_spatial = OS;
    ZWT_spatial=ZWT;
    %%%%%%%%
    Qlat_in_spatial=  Qi_in_space;
    Qlat_out_spatial =  Qi_out_space;
    q_runon_spatial = q_runon;
    Q_channel_spatial = Q_channel;
    %%%%%%%%%
    O_spatial=O_space;
    V_spatial=V_space;
    Vice_spatial=Vice_space;
    Tdp_spatial = Tdp_space;
    SAT_spatial = (ZWT==0);
    TsVEG_spatial = TsVEG;
    %%%%%%%%%%%%%%%%
    er_spatial =  er ;
    DQ_S_spatial =DQ_S;
    DT_S_spatial = DT_S;
    dQ_S_spatial =dQ_S;
    CK1_spatial =  CK1;
    %CK2_spatial = CK2;
    %%%%
    OH_spatial = OH;
    OL_spatial = OL;
    An_H_spatial = An_H;
    An_L_spatial=  An_L ;
    Rdark_H_spatial =Rdark_H;
    Rdark_L_spatial = Rdark_L;
    Tdp_H_spatial =Tdp_H;
    Tdp_L_spatial = Tdp_L;
    Ci_sunH_spatial =Ci_sunH;
    Ci_sunL_spatial = Ci_sunL;
    Ci_shdH_spatial =Ci_shdH;
    Ci_shdL_spatial = Ci_shdL;
    LAI_H_spatial = LAI_H;
    LAI_L_spatial = LAI_L;
    GPP_H_spatial = NPP_H+RA_H;
    GPP_L_spatial = NPP_L +RA_L;
    NPP_H_spatial = NPP_H;
    NPP_L_spatial = NPP_L;
    ANPP_H_spatial = ANPP_H;
    ANPP_L_spatial = ANPP_L;
    B_H_spatial = B_H;
    B_L_spatial = B_L;
    RA_H_spatial =RA_H;
    RA_L_spatial = RA_L;
    Rg_H_spatial =Rg_H;
    Rg_L_spatial = Rg_L;
    Rms_H_spatial =Rms_H;
    Rms_L_spatial = Rms_L;
    Rmr_H_spatial =Rmr_H;
    Rmr_L_spatial = Rmr_L;
    Rmc_H_spatial =Rmc_H;
    Rmc_L_spatial = Rmc_L;
    SAI_H_spatial = SAI_H;
    SAI_L_spatial = SAI_L;
    hc_H_spatial = hc_H;
    hc_L_spatial = hc_L;
    LAIdead_H_spatial = LAIdead_H;
    LAIdead_L_spatial = LAIdead_L;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    %%%%%
    Pr_spatial = Pr_spatial + Pr_S;
    Ta_spatial = ((toutp-2)*Ta_spatial + Ta_S)/(toutp-1) ;
    Ws_spatial = ((toutp-2)*Ws_spatial + Ws_S)/(toutp-1) ;
    Ds_spatial = ((toutp-2)*Ds_spatial + Ds_S)/(toutp-1) ;
    ea_spatial = ((toutp-2)*ea_spatial + ea_S)/(toutp-1) ;
    N_spatial = ((toutp-2)*N_spatial + N_S)/(toutp-1) ;
    Pre_spatial = ((toutp-2)*Pre_spatial + Pre_S)/(toutp-1) ;
    Tdew_spatial = ((toutp-2)*Tdew_spatial + Tdew_S)/(toutp-1) ;
    Rsw_spatial = ((toutp-2)*Rsw_spatial(:) + Rsw_space(:))/(toutp-1) ;
    PAR_spatial = ((toutp-2)*PAR_spatial + PAR_space(:))/(toutp-1) ;
    Ca_spatial = ((toutp-2)*Ca_spatial + Ca_S)/(toutp-1) ;
    %%%
    Ts_spatial = ((toutp-2)*Ts_spatial + Ts)/(toutp-1) ;
    Tdamp_spatial = ((toutp-2)*Tdamp_spatial + Tdamp)/(toutp-1) ;
    Csno_spatial = Csno_spatial + Csno ;
    Csnow_spatial = Csnow_spatial + Csnow ;
    Cice_spatial = Cice_spatial + Cice ;
    Cicew_spatial = Cicew_spatial + Cicew ;
    Pr_sno_spatial = Pr_sno_spatial + Pr_sno;
    Pr_liq_spatial = Pr_liq_spatial + Pr_liq;
    %%%
    Rn_spatial = ((toutp-2)*Rn_spatial + Rn)/(toutp-1) ;
    H_spatial = ((toutp-2)*H_spatial + H)/(toutp-1) ;
    G_spatial = ((toutp-2)*G_spatial + G)/(toutp-1) ;
    Gfin_spatial = ((toutp-2)*Gfin_spatial + Gfin)/(toutp-1) ;
    QE_spatial = ((toutp-2)*QE_spatial + QE)/(toutp-1) ;
    Qv_spatial = ((toutp-2)*Qv_spatial + Qv)/(toutp-1) ;
    Qfm_spatial = ((toutp-2)*Qfm_spatial + Qfm)/(toutp-1) ;
    %%%%
    SWE_spatial = ((toutp-2)*SWE_spatial + SWE)/(toutp-1) ;
    SND_spatial = ((toutp-2)*SND_spatial + SND)/(toutp-1) ;
    WR_SP_spatial = ((toutp-2)*WR_SP_spatial + WR_SP)/(toutp-1) ;
    U_SWE_spatial = ((toutp-2)*U_SWE_spatial + U_SWE)/(toutp-1) ;
    dw_SNO_spatial = ((toutp-2)*dw_SNO_spatial + dw_SNO)/(toutp-1) ;
    ros_spatial = ((toutp-2)*ros_spatial + ros)/(toutp-1) ;
    In_SWE_spatial= ((toutp-2)*In_SWE_spatial + In_SWE)/(toutp-1) ;
    SP_wc_spatial = ((toutp-2)*SP_wc_spatial + SP_wc)/(toutp-1) ;
    % SWE_avalanched_spatial = ((toutp-2)*SWE_avalanched_spatial + SWE_avalanched)/(toutp-1) ;
    %%%%
    ICE_spatial = ((toutp-2)*ICE_spatial + ICE)/(toutp-1) ;
    ICE_D_spatial= ((toutp-2)*ICE_D_spatial + ICE_D)/(toutp-1) ;
    WR_IP_spatial = ((toutp-2)*WR_IP_spatial + WR_IP)/(toutp-1) ;
    NICe_spatial = ((toutp-2)*NICe_spatial + NIce)/(toutp-1) ;
    IP_wc_spatial = ((toutp-2)*IP_wc_spatial + IP_wc)/(toutp-1) ;
    Imelt_spatial = ((toutp-2)*Imelt_spatial + Imelt)/(toutp-1) ;
    Smelt_spatial = ((toutp-2)*Smelt_spatial + Smelt)/(toutp-1) ;
    Tice_spatial = ((toutp-2)*Tice_spatial + Tice)/(toutp-1) ;
    %%%%
    T_H_spatial = ((toutp-2)*T_H_spatial + T_H_space)/(toutp-1) ;
    T_L_spatial = ((toutp-2)*T_L_spatial + T_L_space)/(toutp-1) ;
    EIn_H_spatial = ((toutp-2)*EIn_H_spatial + EIn_H_space)/(toutp-1) ;
    EIn_L_spatial = ((toutp-2)*EIn_L_spatial + EIn_L_space)/(toutp-1) ;
    EG_spatial = ((toutp-2)*EG_spatial + EG)/(toutp-1) ;
    ESN_spatial = ((toutp-2)*ESN_spatial + ESN + ESN_In)/(toutp-1) ;
    EWAT_spatial = ((toutp-2)*EWAT_spatial + EWAT)/(toutp-1) ;
    EICE_spatial = ((toutp-2)*EICE_spatial + EICE)/(toutp-1) ;
    Dr_H_spatial = ((toutp-2)*Dr_H_spatial + Dr_H_space)/(toutp-1) ;
    Dr_L_spatial = ((toutp-2)*Dr_L_spatial + Dr_L_space)/(toutp-1) ;
    EIn_rock_spatial = ((toutp-2)*EIn_rock_spatial + EIn_rock)/(toutp-1) ;
    %%%%%
    In_spatial = ((toutp-2)*In_spatial + In_H_space + In_L_space + SP_wc + In_SWE + In_urb + In_rock + IP_wc)/(toutp-1) ;
    Inveg_spatial = ((toutp-2)*Inveg_spatial + In_H_space + In_L_space)/(toutp-1) ;
    In_rock_spatial =  ((toutp-2)*In_rock_spatial + In_rock)/(toutp-1) ;
    FROCK_spatial = ((toutp-2)*FROCK_spatial + FROCK)/(toutp-1) ;
    WAT_spatial = ((toutp-2)*WAT_spatial + WAT)/(toutp-1) ;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    f_spatial = ((toutp-2)*f_spatial +  f)/(toutp-1) ;
    WIS_spatial = ((toutp-2)*WIS_spatial +  WIS)/(toutp-1) ;
    Rd_spatial =((toutp-2)*Rd_spatial +  Rd)/(toutp-1) ;
    Rh_spatial =((toutp-2)*Rh_spatial +  Rh)/(toutp-1) ;
    Lk_spatial =((toutp-2)*Lk_spatial +  Lk)/(toutp-1) ;
    Lk_rock_spatial =((toutp-2)*Lk_rock_spatial +  Lk_rock)/(toutp-1) ;
    Lk_wat_spatial =((toutp-2)*Lk_wat_spatial +  Lk_wat)/(toutp-1) ;
    OF_spatial =((toutp-2)*OF_spatial +  OF)/(toutp-1) ;
    OS_spatial = ((toutp-2)*OS_spatial +  OS)/(toutp-1) ;
    ZWT_spatial=((toutp-2)*ZWT_spatial +  ZWT)/(toutp-1) ;
    %%%%%%%%
    Qlat_in_spatial=  ((toutp-2)*Qlat_in_spatial +  Qi_in_space)/(toutp-1) ;
    Qlat_out_spatial =  ((toutp-2)*Qlat_out_spatial +  Qi_out_space)/(toutp-1) ;
    q_runon_spatial = ((toutp-2)*q_runon_spatial +  q_runon)/(toutp-1) ;
    Q_channel_spatial = ((toutp-2)*Q_channel_spatial +  Q_channel)/(toutp-1) ;
    %%%%%%%%%
    V_spatial= ((toutp-2)*V_spatial + V_space)/(toutp-1) ;
    Vice_spatial=((toutp-2)*Vice_spatial + Vice_space)/(toutp-1) ;
    O_spatial = ((toutp-2)*O_spatial + O_space)/(toutp-1) ;
    SAT_spatial = SAT_spatial + (ZWT==0) ;
    Tdp_spatial = ((toutp-2)*Tdp_spatial + Tdp_space)/(toutp-1) ;
    TsVEG_spatial = ((toutp-2)*TsVEG_spatial + TsVEG)/(toutp-1) ;
    %%%%%%%%%
    er_spatial = ((toutp-2)*er_spatial + er)/(toutp-1) ;
    DQ_S_spatial = DQ_S_spatial + DQ_S;
    DT_S_spatial = DT_S_spatial + DT_S;
    dQ_S_spatial = dQ_S_spatial + dQ_S;
    CK1_spatial = CK1_spatial + CK1;
    %CK2_spatial = CK2_spatial + CK2;
    %%%%%
    OH_spatial = ((toutp-2)*OH_spatial + OH)/(toutp-1) ;
    OL_spatial = ((toutp-2)*OL_spatial + OL)/(toutp-1) ;
    An_H_spatial = ((toutp-2)*An_H_spatial + An_H)/(toutp-1) ;
    An_L_spatial = ((toutp-2)*An_L_spatial + An_L)/(toutp-1) ;
    Rdark_H_spatial = ((toutp-2)*Rdark_H_spatial + Rdark_H)/(toutp-1) ;
    Rdark_L_spatial = ((toutp-2)*Rdark_L_spatial + Rdark_L)/(toutp-1) ;
    Ci_sunH_spatial = ((toutp-2)*Ci_sunH_spatial + Ci_sunH)/(toutp-1) ;
    Ci_sunL_spatial = ((toutp-2)*Ci_sunL_spatial + Ci_sunL)/(toutp-1) ;
    Ci_shdH_spatial = ((toutp-2)*Ci_shdH_spatial + Ci_shdH)/(toutp-1) ;
    Ci_shdL_spatial = ((toutp-2)*Ci_shdL_spatial + Ci_shdL)/(toutp-1) ;
    LAI_H_spatial = ((toutp-2)*LAI_H_spatial + LAI_H)/(toutp-1) ;
    LAI_L_spatial = ((toutp-2)*LAI_L_spatial + LAI_L)/(toutp-1) ;
    GPP_H_spatial = ((toutp-2)*GPP_H_spatial + NPP_H+RA_H)/(toutp-1) ;
    GPP_L_spatial = ((toutp-2)*GPP_L_spatial + NPP_L+RA_L)/(toutp-1) ;
    NPP_H_spatial = ((toutp-2)*NPP_H_spatial + NPP_H)/(toutp-1) ;
    NPP_L_spatial = ((toutp-2)*NPP_L_spatial + NPP_L)/(toutp-1) ;
    ANPP_H_spatial = ((toutp-2)*ANPP_H_spatial + ANPP_H)/(toutp-1) ;
    ANPP_L_spatial = ((toutp-2)*ANPP_L_spatial + ANPP_L)/(toutp-1) ;
    B_H_spatial = ((toutp-2)*B_H_spatial + B_H)/(toutp-1) ;
    B_L_spatial = ((toutp-2)*B_L_spatial + B_L)/(toutp-1) ;
    RA_H_spatial = ((toutp-2)*RA_H_spatial + RA_H)/(toutp-1) ;
    RA_L_spatial = ((toutp-2)*RA_L_spatial + RA_L)/(toutp-1) ;
    Rg_H_spatial = ((toutp-2)*Rg_H_spatial + Rg_H)/(toutp-1) ;
    Rg_L_spatial = ((toutp-2)*Rg_L_spatial + Rg_L)/(toutp-1) ;
    Rms_H_spatial = ((toutp-2)*Rms_H_spatial + Rms_H)/(toutp-1) ;
    Rms_L_spatial = ((toutp-2)*Rms_L_spatial + Rms_L)/(toutp-1) ;
    Rmr_H_spatial = ((toutp-2)*Rmr_H_spatial + Rmr_H)/(toutp-1) ;
    Rmr_L_spatial = ((toutp-2)*Rmr_L_spatial + Rmr_L)/(toutp-1) ;
    Rmc_H_spatial = ((toutp-2)*Rmc_H_spatial + Rmc_H)/(toutp-1) ;
    Rmc_L_spatial = ((toutp-2)*Rmc_L_spatial + Rmc_L)/(toutp-1) ;
    SAI_H_spatial = ((toutp-2)*SAI_H_spatial +  SAI_H)/(toutp-1) ;
    SAI_L_spatial = ((toutp-2)*SAI_L_spatial +  SAI_L)/(toutp-1) ;
    hc_H_spatial = ((toutp-2)*hc_H_spatial +  hc_H)/(toutp-1) ;
    hc_L_spatial = ((toutp-2)*hc_L_spatial +  hc_L)/(toutp-1) ;
    LAIdead_H_spatial = ((toutp-2)*LAIdead_H_spatial +  LAIdead_H)/(toutp-1) ;
    LAIdead_L_spatial = ((toutp-2)*LAIdead_L_spatial +  LAIdead_L)/(toutp-1) ;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if  length(intersect(t,tstore))==1 ||  t==N_time_step
    Title_save = strcat([Out.OutDir, '/OUTPUT_',TITLE_SAVE,'_SPATIAL_',num2str(t)]);
    save(Title_save);
    %%%%%%%%%%
    toutp = 1;
    %%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% TRACKED PIXELS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ipo=1:npoint
    ij = sub2ind(size(DTM),Yout(ipo),Xout(ipo));
    TrackPoint(t,ipo,	1	)=(	Pr_S(ij))	;
    TrackPoint(t,ipo,	2	)=(	Ta_S(ij))	;
    TrackPoint(t,ipo,	3	)=(	Ws_S(ij))	;
    TrackPoint(t,ipo,	4	)=(	Ds_S(ij))	;
    TrackPoint(t,ipo,	5	)=(	ea_S(ij))	;
    TrackPoint(t,ipo,	6	)=(	N_S(ij))	;
    TrackPoint(t,ipo,	7	)=(	Pre_S(ij))	;
    TrackPoint(t,ipo,	8	)=(	Tdew_S(ij))	;
    TrackPoint(t,ipo,	9	)=(	Rsw_space(ij))	;
    TrackPoint(t,ipo,	10	)=(	PAR_space(ij))	;
    TrackPoint(t,ipo,	11	)=(	Ca_S(ij))	;
    TrackPoint(t,ipo,	12	)=(	Ts(ij))	;
    TrackPoint(t,ipo,	13	)=(	Tdamp(ij))	;
    TrackPoint(t,ipo,	14	)=(	Csno(ij))	;
    TrackPoint(t,ipo,	15	)=(	Cice(ij))	;
    TrackPoint(t,ipo,	16	)=(	Csnow(ij))	;
    TrackPoint(t,ipo,	17	)=(	Cicew(ij))	;
    TrackPoint(t,ipo,	18	)=(	Pr_sno(ij))	;
    TrackPoint(t,ipo,	19	)=(	Pr_liq(ij))	;
    TrackPoint(t,ipo,	20	)=(	Rn(ij))	;
    TrackPoint(t,ipo,	21	)=(	H(ij))	;
    TrackPoint(t,ipo,	22	)=(	G(ij))	;
    TrackPoint(t,ipo,	23	)=(	Gfin(ij))	;
    TrackPoint(t,ipo,	24	)=(	QE(ij))	;
    TrackPoint(t,ipo,	25	)=(	Qv(ij))	;
    TrackPoint(t,ipo,	26	)=(	Qfm(ij))	;
    TrackPoint(t,ipo,	27	)=(	SWE(ij))	;
    TrackPoint(t,ipo,	28	)=(	SND(ij))	;
    TrackPoint(t,ipo,	29	)=(	WR_SP(ij))	;
    TrackPoint(t,ipo,	30	)=(	U_SWE(ij))	;
    TrackPoint(t,ipo,	31	)=(	NIn_SWE(ij))	;
    TrackPoint(t,ipo,	32	)=(	dw_SNO(ij))	;
    TrackPoint(t,ipo,	33	)=(	ros(ij))	;
    TrackPoint(t,ipo,	34	)=(	In_SWE(ij))	;
    TrackPoint(t,ipo,	35	)=(	SP_wc(ij))	;
    TrackPoint(t,ipo,	36	)=(	ICE(ij))	;
    TrackPoint(t,ipo,	37	)=(	ICE_D(ij))	;
    TrackPoint(t,ipo,	38	)=(	WR_IP(ij))	;
    TrackPoint(t,ipo,	39	)=(	IP_wc(ij))	;
    TrackPoint(t,ipo,	40	)=(	NIce(ij))	;
    TrackPoint(t,ipo,	41	)=(	EG(ij))	;
    TrackPoint(t,ipo,	42	)=(	ESN(ij))	;
    TrackPoint(t,ipo,	43	)=(	ESN_In(ij))	;
    TrackPoint(t,ipo,	44	)=(	EWAT(ij))	;
    TrackPoint(t,ipo,	45	)=(	EICE(ij))	;
    TrackPoint(t,ipo,	46	)=(	EIn_urb(ij))	;
    TrackPoint(t,ipo,	47	)=(	EIn_rock(ij))	;
    TrackPoint(t,ipo,	48	)=(	SE_rock(ij))	;
    TrackPoint(t,ipo,	49	)=(	SE_urb(ij))	;
    TrackPoint(t,ipo,	50	)=(	In_urb(ij))	;
    TrackPoint(t,ipo,	51	)=(	In_rock(ij))	;
    TrackPoint(t,ipo,	52	)=(	WAT(ij))	;
    TrackPoint(t,ipo,	53	)=(	FROCK(ij))	;
    TrackPoint(t,ipo,	54	)=(	f(ij))	;
    TrackPoint(t,ipo,	55	)=(	WIS(ij))	;
    TrackPoint(t,ipo,	56	)=(	Rd(ij))	;
    TrackPoint(t,ipo,	57	)=(	Rh(ij))	;
    TrackPoint(t,ipo,	58	)=(	Lk(ij))	;
    TrackPoint(t,ipo,	59	)=(	Lk_wat(ij))	;
    TrackPoint(t,ipo,	60	)=(	Lk_rock(ij))	;
    TrackPoint(t,ipo,	61	)=(	OF(ij))	;
    TrackPoint(t,ipo,	62	)=(	OS(ij))	;
    TrackPoint(t,ipo,	63	)=(	ZWT(ij))	;
    TrackPoint(t,ipo,	64	)=(	q_runon(ij))	;
    TrackPoint(t,ipo,	65	)=(	Q_channel(ij))	;
    TrackPoint(t,ipo,	66	)=(	r_soil(ij))	;
    TrackPoint(t,ipo,	67	)=(	alp_soil(ij))	;
    TrackPoint(t,ipo,	68	)=(	ra(ij))	;
    TrackPoint(t,ipo,	69	)=(	er(ij))	;
    TrackPoint(t,ipo,	70	)=(	DQ_S(ij))	;
    TrackPoint(t,ipo,	71	)=(	DT_S(ij))	;
    TrackPoint(t,ipo,	72	)=(	dQ_S(ij))	;
    TrackPoint(t,ipo,	73	)=(	QpointH(ipo))	;
    TrackPoint(t,ipo,	74	)=(	QpointC(ipo))	;
    TrackPoint(t,ipo,	75	)=(	UpointH(ipo))	;
    TrackPoint(t,ipo,	76	)=(	UpointC(ipo))	;
    TrackPoint(t,ipo,	77	)=(	CK1(ij))	;


    for ijki = 1:cc_max
        TrackPointVeg(t,ipo,ijki,	1	)=(	T_H(ij	,	ijki))	;
        TrackPointVeg(t,ipo,ijki,	2	)=(	T_L(ij	,	ijki))	;
        TrackPointVeg(t,ipo,ijki,	3	)=(	EIn_H(ij	,	ijki))	;
        TrackPointVeg(t,ipo,ijki,	4	)=(	EIn_L(ij	,	ijki))	;
        TrackPointVeg(t,ipo,ijki,	5	)=(	Dr_H(ij	,	ijki))	;
        TrackPointVeg(t,ipo,ijki,	6	)=(	Dr_L(ij	,	ijki))	;
        TrackPointVeg(t,ipo,ijki,	7	)=(	In_H(ij	,	ijki))	;
        TrackPointVeg(t,ipo,ijki,	8	)=(	In_L(ij	,	ijki))	;
        TrackPointVeg(t,ipo,ijki,	9	)=(	Tdp_H(ij	,	ijki))	;
        TrackPointVeg(t,ipo,ijki,	10	)=(	Tdp_L(ij	,	ijki))	;
        TrackPointVeg(t,ipo,ijki,	11	)=(	OH(ij	,	ijki))	;
        TrackPointVeg(t,ipo,ijki,	12	)=(	OL(ij	,	ijki))	;
        TrackPointVeg(t,ipo,ijki,	13	)=(	rb_H(ij	,	ijki))	;
        TrackPointVeg(t,ipo,ijki,	14	)=(	rb_L(ij	,	ijki))	;
        TrackPointVeg(t,ipo,ijki,	15	)=(	rs_sunH(ij	,	ijki))	;
        TrackPointVeg(t,ipo,ijki,	16	)=(	rs_sunL(ij	,	ijki))	;
        TrackPointVeg(t,ipo,ijki,	17	)=(	rs_shdH(ij	,	ijki))	;
        TrackPointVeg(t,ipo,ijki,	18	)=(	rs_shdL(ij	,	ijki))	;
        TrackPointVeg(t,ipo,ijki,	19	)=(	rap_H(ij	,	ijki))	;
        TrackPointVeg(t,ipo,ijki,	20	)=(	rap_L(ij	,	ijki))	;
        TrackPointVeg(t,ipo,ijki,	21	)=(	An_H(ij	,	ijki))	;
        TrackPointVeg(t,ipo,ijki,	22	)=(	An_L(ij	,	ijki))	;
        TrackPointVeg(t,ipo,ijki,	23	)=(	Rdark_H(ij	,	ijki))	;
        TrackPointVeg(t,ipo,ijki,	24	)=(	Rdark_L(ij	,	ijki))	;
        TrackPointVeg(t,ipo,ijki,	25	)=(	Ci_sunH(ij	,	ijki))	;
        TrackPointVeg(t,ipo,ijki,	26	)=(	Ci_sunL(ij	,	ijki))	;
        TrackPointVeg(t,ipo,ijki,	27	)=(	Ci_shdH(ij	,	ijki))	;
        TrackPointVeg(t,ipo,ijki,	28	)=(	Ci_shdL(ij	,	ijki))	;
        TrackPointVeg(t,ipo,ijki,	29	)=(	LAI_H(ij	,	ijki))	;
        TrackPointVeg(t,ipo,ijki,	30	)=(	LAI_L(ij	,	ijki))	;
        TrackPointVeg(t,ipo,ijki,	31	)=(	NPP_H(ij	,	ijki))	;
        TrackPointVeg(t,ipo,ijki,	32	)=(	NPP_L(ij	,	ijki))	;
        TrackPointVeg(t,ipo,ijki,	33	)=(	ANPP_H(ij	,	ijki))	;
        TrackPointVeg(t,ipo,ijki,	34	)=(	ANPP_L(ij	,	ijki))	;
        TrackPointVeg(t,ipo,ijki,	35	)=(	RA_H(ij	,	ijki))	;
        TrackPointVeg(t,ipo,ijki,	36	)=(	RA_L(ij	,	ijki))	;
        TrackPointVeg(t,ipo,ijki,	37	)=(	Rg_H(ij	,	ijki))	;
        TrackPointVeg(t,ipo,ijki,	38	)=(	Rg_L(ij	,	ijki))	;
        TrackPointVeg(t,ipo,ijki,	39	)=(	Rms_H(ij	,	ijki))	;
        TrackPointVeg(t,ipo,ijki,	40	)=(	Rms_L(ij	,	ijki))	;
        TrackPointVeg(t,ipo,ijki,	41	)=(	Rmr_H(ij	,	ijki))	;
        TrackPointVeg(t,ipo,ijki,	42	)=(	Rmr_L(ij	,	ijki))	;
        TrackPointVeg(t,ipo,ijki,	43	)=(	Rmc_H(ij	,	ijki))	;
        TrackPointVeg(t,ipo,ijki,	44	)=(	Rmc_L(ij	,	ijki))	;
        TrackPointVeg(t,ipo,ijki,	45	)=(	PHE_S_H(ij	,	ijki))	;
        TrackPointVeg(t,ipo,ijki,	46	)=(	PHE_S_L(ij	,	ijki))	;
        TrackPointVeg(t,ipo,ijki,	47	)=(	dflo_H(ij	,	ijki))	;
        TrackPointVeg(t,ipo,ijki,	48	)=(	dflo_L(ij	,	ijki))	;
        TrackPointVeg(t,ipo,ijki,	49	)=(	AgeL_H(ij	,	ijki))	;
        TrackPointVeg(t,ipo,ijki,	50	)=(	AgeL_L(ij	,	ijki))	;
        TrackPointVeg(t,ipo,ijki,	51	)=(	SAI_H(ij	,	ijki))	;
        TrackPointVeg(t,ipo,ijki,	52	)=(	SAI_L(ij	,	ijki))	;
        TrackPointVeg(t,ipo,ijki,	53	)=(	LAIdead_H(ij	,	ijki))	;
        TrackPointVeg(t,ipo,ijki,	54	)=(	LAIdead_L(ij	,	ijki))	;
        TrackPointVeg(t,ipo,ijki,	55	)=(	hc_H(ij	,	ijki))	;
        TrackPointVeg(t,ipo,ijki,	56	)=(	hc_L(ij	,	ijki))	;
        TrackPointVeg(t,ipo,ijki,	57	)=(	e_rel_L(ij	,	ijki))	;
        TrackPointVeg(t,ipo,ijki,	58	)=(	e_rel_H(ij	,	ijki))	;
        TrackPointVeg(t,ipo,ijki,	59	)=(	B_H(ij	,	ijki	,1));
        TrackPointVeg(t,ipo,ijki,	60	)=(	B_H(ij	,	ijki	,2));
        TrackPointVeg(t,ipo,ijki,	61	)=(	B_H(ij	,	ijki	,3));
        TrackPointVeg(t,ipo,ijki,	62	)=(	B_H(ij	,	ijki	,4));
        TrackPointVeg(t,ipo,ijki,	63	)=(	B_H(ij	,	ijki	,5));
        TrackPointVeg(t,ipo,ijki,	64	)=(	B_H(ij	,	ijki	,6));
        TrackPointVeg(t,ipo,ijki,	65	)=(	B_H(ij	,	ijki	,7));
        TrackPointVeg(t,ipo,ijki,	66	)=(	B_L(ij	,	ijki	,1));
        TrackPointVeg(t,ipo,ijki,	67	)=(	B_L(ij	,	ijki	,2));
        TrackPointVeg(t,ipo,ijki,	68	)=(	B_L(ij	,	ijki	,3));
        TrackPointVeg(t,ipo,ijki,	69	)=(	B_L(ij	,	ijki	,4));
        TrackPointVeg(t,ipo,ijki,	70	)=(	B_L(ij	,	ijki	,5));
        TrackPointVeg(t,ipo,ijki,	71	)=(	B_L(ij	,	ijki	,6));
        TrackPointVeg(t,ipo,ijki,	72	)=(	B_L(ij	,	ijki	,7));

    end



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


% update memory items

QE_cum = QE_cum + QE;
Rn_cum = Rn_cum + Rn;
H_cum = H_cum + H;
LAI_H_cum = LAI_Htm1+LAI_H;
LAI_L_cum  =LAI_Ltm1+ LAI_L;
NPP_L_cum = NPP_L_cum +NPP_L; 
NPP_H_cum = NPP_L_cum +NPP_L; 
RA_L_cum = RA_L_cum + RA_L; 
RA_H_cum = RA_H_cum + RA_H; 
T_L_cum = T_L_cum + T_L; 
T_H_cum = T_H_cum + T_H;
G_cum = G_cum + G; 
GPP_L_cum = GPP_L_cum+RA_L+NPP_L; 
GPP_H_cum = GPP_H_cum+RA_H+NPP_H;
SND_cum = SND_cum +SND; 




if OPT_VAR_SM_ON == 1
    save([Out.OutDir, '/OUTPUT_stats'],'TrackPointVeg','TrackPoint','PFT_STD','PFT_AVG','SIM_AVG','SIM_STD','SIM_VAR_SM','QE_cum','Rn_cum','H_cum','LAI_H_cum', ...
        'LAI_L_cum','NPP_L_cum','NPP_H_cum','RA_H_cum','RA_L_cum','T_L_cum','T_H_cum','G_cum','GPP_L_cum','GPP_H_cum','SND_cum')
else
    save([Out.OutDir, '/OUTPUT_stats'],'TrackPointVeg','TrackPoint','PFT_STD','PFT_AVG','SIM_AVG','SIM_STD','QE_cum','Rn_cum','H_cum','LAI_H_cum', ...
        'LAI_L_cum','NPP_L_cum','NPP_H_cum','RA_H_cum','RA_L_cum','T_L_cum','T_H_cum','G_cum','GPP_L_cum','GPP_H_cum','SND_cum')
end











