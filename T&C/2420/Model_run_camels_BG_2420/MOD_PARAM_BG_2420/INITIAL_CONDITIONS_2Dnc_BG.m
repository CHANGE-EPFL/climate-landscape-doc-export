%% Initial Conditions

%vi = 0.3*[10 30 50 400 1500]; %irrelevant
% vi=[8.4480 12.7442 21.4297 43.5578 44.4718 91.7058 96.2888];
FROCKtm1= 3.*MASKn;
Zs = [ 0 20 50 100 200 300 500 Zsmax];
dz= diff(Zs); %%%% [mm]  Thickness of the Layers
for ii=1:ms_max
    O_t(:,ii,1)=  reshape(Ofc,m_cell*n_cell,1);
end
V_t(:,:,1)= (O_t(:,:,1) -reshape(Ohy,m_cell*n_cell,1)).*dz;
Vtm1=V_t(:,:,1);
% 
% vi = [0.2 0.2 0.2 0.2 0.3 0.3 0.8].*[10 30 50 100 100 200 1000];
% 
% Vtm1=zeros(num_cell,ms_max);
% for jk=1:ms_max
%     Vtm1(:,jk)=vi(jk)*MASKn;
% end
%%%%%%%%
Vicetm1=zeros(num_cell,ms_max);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% HYDROLOGY PART %%%%%
Tstm1 =1*MASKn;
Tdptm1 = 1*ones(num_cell,ms_max);
Tdpsnowtm1 = -2*zeros(num_cell,ms_max);
Tdamptm1 = Tstm1;
Ts_undertm1 = 5*ones(num_cell,ms_max);
TdpI_Htm1= Tdamptm1*ones(1,cc_max);
TdpI_Ltm1= Tdamptm1*ones(1,cc_max);
%%%
SWEtm1 = 0;%0.7*(DTM-4800);
SWEtm1(SWEtm1<0)=0;
SWEtm1 =0.*MASKn;
SNDtm1 =0.*MASKn;
snow_albedotm1 = 0.8*ones(num_cell,4) ;
e_snotm1 = 0.97*MASKn;
t_slstm1= 0*MASKn;
rostm1=0*MASKn;
SP_wctm1=0.*MASKn;
tau_snotm1 = 0*MASKn;
In_Htm1 = 0*ones(num_cell,cc_max);
In_Ltm1 = 0*ones(num_cell,cc_max);
In_Littertm1 = 0*MASKn;
In_SWEtm1 = 0*MASKn;
In_urbtm1 = 0*MASKn;
In_rocktm1 = 0*MASKn;
EKtm1 = 0*MASKn;
WATtm1= 0*MASKn;
ICEtm1 = 0;%0.5*(DTM-4800);
ICEtm1(ICEtm1<0)=0;
ICEtm1= 0.*MASKn;
IP_wctm1= 0.*MASKn;
ICE_Dtm1= 0.*MASKn;
Cicewtm1= 0*MASKn;
% FROCKtm1= 0*MASKn;
Ticetm1= 0*MASKn;
Oicetm1 =zeros(num_cell,ms_max);
Tdebtm1 = 0*ones(num_cell,md_max);
Ws_undertm1 = 1*MASKn;
%%%
Ccrown_t_tm1 = ones(num_cell,cc_max);
for iii = 1:size(Geo.PFT_def, 3)
    temp_f = Geo.PFT_def(:,:,iii);
    Ccrown_t_tm1(:,iii) = temp_f(:);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tstm0 =Tstm1;
Slo_head = reshape(Slo_top,num_cell,1)*ones(1,ms_max); %%%
q_runon =  0*MASKn;
Q_channel = 10*ones(m_cell,n_cell);
Qi_in = zeros(num_cell,ms_max);
Qi_out = zeros(num_cell,ms_max);
Qi_out_Rout= zeros(m_cell,n_cell,ms_max);
Qi_in_Rout= zeros(m_cell,n_cell,ms_max);
Q_exit=0;
Qsub_exit=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%% HIGH VEGETATION
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
Citm1_sunH = 370*ones(num_cell,cc_max);
Citm1_shdH = 370*ones(num_cell,cc_max);
LAI_Htm1 = zeros(num_cell,cc_max);
B_Htm1 =  zeros(num_cell,cc_max,8);
PHE_S_Htm1=1*ones(num_cell,cc_max);
NPP_Htm1=zeros(num_cell,cc_max);
NPPI_Htm1=zeros(num_cell,cc_max);
Bfac_weekHtm1 = ones(num_cell,cc_max);
dflo_Htm1= zeros(num_cell,cc_max);
AgeL_Htm1= zeros(num_cell,cc_max);
AgeDL_Htm1= zeros(num_cell,cc_max);
e_rel_Htm1=ones(num_cell,cc_max);
e_relN_Htm1=ones(num_cell,cc_max);
SAI_Htm1=zeros(num_cell,cc_max)+0.1;
hc_Htm1=zeros(num_cell,cc_max);
Tden_Htm1 =zeros(num_cell,cc_max);
AgePl_Htm1 =zeros(num_cell,cc_max);
Nreserve_Htm1= 1000*ones(num_cell,cc_max);
Preserve_Htm1= 1000*ones(num_cell,cc_max);
Kreserve_Htm1= 1000*ones(num_cell,cc_max);
NupI_Htm1 = zeros(num_cell,cc_max,3);
Bfac_weekLtm1 = ones(num_cell,cc_max);
FNC_Htm1=1*ones(num_cell,cc_max);
Vx_Htm1= zeros(num_cell,cc_max);
Vl_Htm1= zeros(num_cell,cc_max);
Psi_x_Htm1= zeros(num_cell,cc_max);
Psi_l_Htm1 = zeros(num_cell,cc_max);

LAI_Htm1(:,1)=1;
% introduce automated initlial carbon pools

% B_H(1,:,:)= [304 461 507 309 35 10 76 10 ].*ones(n_PFTs,1);
% for a = 1:n_PFTs
%     if aSE_H(a)==2
%         B_H(1,a,:) = [66 0 165 118 10 0 2 0 ];
%     end
% end
cc_idx = find(IIpar==1);
for i = 1:cc_max
    if (cc_idx(i) ==4) ||(cc_idx(i) ==5)
        B_Htm1(:,i,:) = [66 0 165 118 10 0 2 0 ].*ones(num_cell,1,1);
    else
        B_Htm1(:,i,:)= [304 461 507 309 35 10 76 10 ].*ones(num_cell,1,1);
    end
end
%
%  % B_Htm1 = zeros(:,cc_max,)


PHE_S_Htm1=1*ones(num_cell,cc_max);
dflo_Htm1= zeros(num_cell,cc_max);
AgeL_Htm1= zeros(num_cell,cc_max);
AgeL_Htm1(:,1)=0;
SAI_Htm1=zeros(num_cell,cc_max);
SAI_Htm1(: ,1)=0.1;
hc_Htm1=zeros(num_cell,cc_max);
hc_Htm1(:,1)=20;
e_rel_Htm1=ones(num_cell,cc_max);
%%%%%%


%%%% LOW VEGETATION
Citm1_sunL = 370*ones(num_cell,cc_max);
Citm1_shdL = 370*ones(num_cell,cc_max);
LAI_Ltm1= zeros(num_cell,cc_max);
B_Ltm1= zeros(num_cell,cc_max,8);
PHE_S_Ltm1=1*ones(num_cell,cc_max);
NPP_Ltm1=zeros(num_cell,cc_max);
NPPI_Ltm1=zeros(num_cell,cc_max);
dflo_Ltm1= zeros(num_cell,cc_max);
AgeL_Ltm1= zeros(num_cell,cc_max);
AgeDL_Ltm1= zeros(num_cell,cc_max);
SAI_Ltm1=zeros(num_cell,cc_max);
hc_Ltm1=zeros(num_cell,cc_max);
Tden_Ltm1= zeros(num_cell,cc_max);
AgePl_Ltm1= zeros(num_cell,cc_max);
e_rel_Ltm1=ones(num_cell,cc_max);
e_relN_Ltm1=ones(num_cell,cc_max);
Nreserve_Ltm1= 1000*ones(num_cell,cc_max);
Preserve_Ltm1= 1000*ones(num_cell,cc_max);
Kreserve_Ltm1= 1000*ones(num_cell,cc_max);
NupI_Ltm1 = zeros(num_cell,cc_max,3);
FNC_Ltm1=1*ones(num_cell,cc_max);
Vx_Ltm1= zeros(num_cell,cc_max);
Vl_Ltm1= zeros(num_cell,cc_max);
Psi_x_Ltm1 = zeros(num_cell,cc_max);
Psi_l_Ltm1 = zeros(num_cell,cc_max);

%%%%%%
NuLit_Htm1= zeros(num_cell,cc_max,3);
NuLit_Ltm1= zeros(num_cell,cc_max,3);
NBLeaf_Htm1 =zeros(num_cell,cc_max);
NBLeaf_Ltm1 =zeros(num_cell,cc_max);
PARI_Htm1 =zeros(num_cell,cc_max,3);
NBLI_Htm1  =zeros(num_cell,cc_max);
PARI_Ltm1 =zeros(num_cell,cc_max,3);
NBLI_Ltm1  =zeros(num_cell,cc_max);
%%%
TBio_H = zeros(num_cell,cc_max);
TBio_L = zeros(num_cell,cc_max);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BG added %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RexmyI=zeros(num_cell,3);
LEAK_NH4=zeros(num_cell,1);
LEAK_NO3=zeros(num_cell,1);
LEAK_P=zeros(num_cell,1);
LEAK_K=zeros(num_cell,1);
LEAK_DOC=zeros(num_cell,1);
LEAK_DON=zeros(num_cell,1);
LEAK_DOP=zeros(num_cell,1);
R_litter=zeros(num_cell,1);
R_microbe=zeros(num_cell,1);
R_litter_sur=zeros(num_cell,1);
R_ew=zeros(num_cell,1);
VOL=zeros(num_cell,1);
N2flx=zeros(num_cell,1);
Min_N=zeros(num_cell,1);
Min_P=zeros(num_cell,1);
R_bacteria=zeros(num_cell,1);
RmycAM=zeros(num_cell,1);
RmycEM=zeros(num_cell,1);
Prod_B=zeros(num_cell,1);
Prod_F=zeros(num_cell,1);
BfixN=zeros(num_cell,1);
LitFirEmi=zeros(num_cell,2);
NEE=zeros(num_cell,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if OPT_SoilBiogeochemistry == 1
    load('Rule_based_cluster_map')
    ksv=veg_class_map;

    ksv_list=unique(veg_class_map(~isnan(veg_class_map)));
    % P=zeros(num_cell,55);
    % Ptm1=zeros(num_cell,55);
    % PHs=5;
    % Upl=0.008;
    % N_nt=7;
    TR_ice=zeros(m_cell*n_cell,N_nt);
    TR_sur=zeros(m_cell,n_cell,N_nt);
    TR_ch=zeros(m_cell,n_cell,N_nt);%nutrient channal concentration [g/mm]
    TR_rock=zeros(1,N_nt);
    for ii=1:num_cell
        for kk=1:length(ksv_list)
            if MASKn(ii)== 1 && ksv(ii) == ksv_list(kk)

                load([Param.basin_name '_Plant_' num2str(ksv(ii)) '_2D_initial.mat'])
                B_Htm1(ii,:,:)=B_H_avg;
                Nreserve_Htm1(ii,:)= Nreserve_avg;
                Preserve_Htm1(ii,:)= Preserve_avg; %
                Kreserve_Htm1(ii,:)= Kreserve_avg; %%%% cc_max = 1; so its (ii,1,1:3)
                NupI_Htm1(ii,:,1:3)= NupI_avg;
                RexmyItm1(ii,1:3)= RexmyI_avg;
                NavlItm1(ii,1:3)=NavlI_avg;
                Ptm1(ii,:)=P_avg;
            end
        end
    end
    clear ii Nreserve_avg Preserve_avg Kreserve_avg NupI_avg RexmyI_avg NavlI_avg P_avg
    if isfield(Param, 'spin') && Param.spin == 1
        clear Ptm1
        disp("RF model spin up applied")
        load(['RF_Ptm1_' Param.basin_name '_revised.mat']) %%%load Ptm1
    else
        disp('RF model spin up disabled')
    end
    Zs = [ 0 20 50 100 200 300 500 Zsmax];
    dz= diff(Zs); %%%% [mm]  Thickness of the Layers
    for ii=1:ms_max
        O_t(:,ii,1)=  reshape(Ofc,m_cell*n_cell,1);
    end
    V_t(:,:,1)= (O_t(:,:,1) -reshape(Ohy,m_cell*n_cell,1)).*dz;
    Vtm1=V_t(:,:,1);
    clear ii
    clear Ofc Ohy
end