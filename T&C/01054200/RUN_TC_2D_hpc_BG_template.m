
function [] = RUN_TC_2D_hpc_BG_template(Meteo, Geo, Param, Out, Date)
%% set start end year
start_year=year(Date(1));
end_year=year(Date(end));
basin_name=Param.basin_name;
time_index=[num2str(start_year) num2str(end_year) 'yr'];
%% use other rainfall product
fn = [pwd filesep 'meteo' filesep basin_name filesep basin_name '.nc'];

% GMT_OFFSET = '-06:00';   % <-- change this to your local GMT offset
GMT_OFFSET = sprintf('%+03d:00', Geo.DeltaGMT);
% ---- 1) Read time & precipitation ----
t_hours = ncread(fn,'DateTime');   % hours since 1980-01-01 00:00:00 (UTC)
Pr      = ncread(fn,'Rainf');      % precipitation variable

% ---- 2) Build UTC datetime ----
t0_utc = datetime(1980,1,1,0,0,0,'TimeZone','UTC');
t_utc  = t0_utc + hours(double(t_hours));  % UTC timeline

% ---- 3) Create "local" timeline using fixed GMT offset (no DST) ----
t_local = t_utc;
t_local.TimeZone = GMT_OFFSET;

% ---- 4) Define target period in LOCAL time (fixed GMT offset) ----
t_start_local = datetime(start_year,1,1,0,0,0,'TimeZone',GMT_OFFSET);
t_end_local   = datetime(end_year,12,31,23,0,0,'TimeZone',GMT_OFFSET);

% ---- 5) Find indices based on LOCAL time ----
idx = (t_local >= t_start_local) & (t_local <= t_end_local);

% ---- 6) Subset outputs ----
t_sel_local = t_local(idx);        % local timestamps (fixed GMT offset)
t_sel_utc   = t_utc(idx);          % corresponding UTC timestamps (optional)
Pr_sel      = double(Pr(idx));     % precipitation values aligned with t_sel_local
clear Pr



%% Unpack
% Meteo - to change with netcdf
x1 = Meteo.x1;
x2 = Meteo.x2;
mkdir(Meteo.FolderLoc)

Date_dv = datevec(Date);
%%%%%% simplification for 2D, but 1D can be different
t_bef = 0;
t_aft = 1;
 DeltaGMT = Geo.DeltaGMT; 

Date_time = datetime(Date', 'ConvertFrom','datenum');
if DeltaGMT<0
    Date_time = Date_time - hours(abs(DeltaGMT));
    Date = datenum(Date_time);

    [YE,MO,DA,HO,MI,SE] = datevec(Date);
    Datam(:,1) = YE; Datam(:,2)= MO; Datam(:,3)= DA; Datam(:,4)= HO;
    clear YE MO DA HO MI SE

    x1 =abs(DeltaGMT)+1;  
    x2_max=length(Date)-abs(DeltaGMT);
else
    Date_time = Date_time + hours(abs(DeltaGMT));
    Date = datenum(Date_time);
    [YE,MO,DA,HO,MI,SE] = datevec(Date);
    Datam(:,1) = YE; Datam(:,2)= MO; Datam(:,3)= DA; Datam(:,4)= HO;
    clear YE MO DA HO MI SE

    x1 =24-abs(DeltaGMT)+1;  
    x2_max=length(Date)-abs(DeltaGMT);

end
%%% if the mannually lenth is not long enough
if x2_max<x2
    x2=x2_max;
end
Date=Date(x1:x2);Datam = Datam(x1:x2,:); 

try
    %ZATM = ncread([Meteo.FolderLoc, '\Ta_' time_index '.nc'], 'zatm');
    ZATM = ncread([Meteo.FolderLoc, filesep 'Ta_' time_index '.nc'], 'zatm'); %% mannualy?

catch
    ZATM = 30*ones(size(Geo.DTM));
end

N_time_step= x2 - x1 + 1; 
Nd_time_step = ceil(N_time_step/24)+1;
tstore= [40,8766:8766:N_time_step];

% Geo
DTM = Geo.DTM;

num_cell = numel(DTM);

%% Load the first month of ERA5 - Data
% 
% WorldClim = py.ModuleERAmonthly142225.wc_get(Geo.shapefile_path, Meteo.Meteo_sc);
% Load_EE_ERA5(Date_dv(1,1), Date_dv(1,2), Geo.shapefile_path, Meteo.Meteo_sc, WorldClim,...
%     Meteo.FolderLoc, t_bef, t_aft, DeltaGMT )

ix_rd = x1+1; %netcdf Reading index
%%%%%% to speed up
blockSize = 192; %%% load 8 days

% optional: save file names once
fnTa   = [Meteo.FolderLoc, filesep, 'Ta_' time_index '.nc'];
fnU    = [Meteo.FolderLoc, filesep, 'U_' time_index '.nc'];
fnN    = [Meteo.FolderLoc, filesep, 'N.nc'];
fnPr   = [Meteo.FolderLoc, filesep, 'Pr_' time_index '.nc'];
fnPre  = [Meteo.FolderLoc, filesep, 'Pre_' time_index '.nc'];
fnWs   = [Meteo.FolderLoc, filesep, 'Ws_' time_index '.nc'];
fnSAB1 = [Meteo.FolderLoc, filesep, 'SAB1.nc'];
fnSAB2 = [Meteo.FolderLoc, filesep, 'SAB2.nc'];
fnSAD1 = [Meteo.FolderLoc, filesep, 'SAD1.nc'];
fnSAD2 = [Meteo.FolderLoc, filesep, 'SAD2.nc'];
fnPARB = [Meteo.FolderLoc, filesep, 'PARB.nc'];
fnPARD = [Meteo.FolderLoc, filesep, 'PARD.nc'];
fnCa   = [Meteo.FolderLoc, filesep, 'Ca.nc'];

% block cache
blk_start_ix = -1;
blk_end_ix   = -1;
% py.ModuleERAmonthly142225.ERA5Land_mon_retrieve_save(Date_dv(1,1), Date_dv(1,2), Geo.shapefile_path, Meteo.Meteo_sc, WorldClim, Meteo.FolderLoc )

% info = ncinfo([Meteo.FolderLoc, '\Ta_' time_index '.nc']);
info = ncinfo([Meteo.FolderLoc, filesep 'Ta_' time_index '.nc']);
projM = projcrs(int32(str2num(info.Attributes.Value)));

% Meteo.x = ncread([Meteo.FolderLoc, '\Ta_' time_index '.nc'], 'x');
% Meteo.y = ncread([Meteo.FolderLoc, '\Ta_' time_index '.nc'], 'y');
Meteo.x = ncread([Meteo.FolderLoc, filesep 'Ta_' time_index '.nc'], 'x');
Meteo.y = ncread([Meteo.FolderLoc, filesep 'Ta_' time_index '.nc'], 'y');
[x_, y_] = meshgrid(Meteo.x, Meteo.y);
% Reproject to Geo System. Should be ideantical, but as a fail safe
[latM, lonM] = projinv(projM, x_, y_);
[xm, ym] = projfwd(Geo.proj, latM, lonM);

Meteo.x = xm;
Meteo.y = ym;

SPRINGn = 0;
% xllcorner = Geo.xllcorner;
% yllcorner = Geo.yllcorner;
PSAN = Geo.PSAN;
PCLA = Geo.PCLA;
PORG = Geo.PORG;
%KBOT = Geo.KBOT;
KROCK = Geo.KROCK;
PFTdef = Geo.PFT_def;
proj = Geo.proj;
leaf_sncs_T =  Geo.leaf_sncs_T;
leaf_onset_T =  Geo.leaf_onset_T;
leaf_onset_day =Geo.leaf_onset_day;
leaf_sncs_day = Geo.leaf_sncs_day;
eco_idx =double(Geo.eco_idx);
lon = Geo.lon;
lat = Geo.lat;
IIpar = Geo.IIpar; 

UALB = Geo.UALB;
BuildH = Geo.Bh;
Ue_sur = Geo.Ue_sur;

%ksv= Geo.ksv; % Land use type
cellsize = Geo.sc;
Xout = Geo.Xout; % geo outlet index, use Out.Xout instead
Yout = Geo.Yout;
% Xout = Out.Xout;
% Yout = Out.Yout;


Slo_top = Geo.Slo_top;
Aspect = Geo.Aspect;
T_flow = Geo.T;
%WC = Geo.WC;
SN = Geo.SN;
NMAN_C = Geo.NMAN_C;
NMAN_H = Geo.NMAN_H;
MRough = Geo.MRough;
% Kres_Rock = 1./Slo_top*24; %%[h] Bedrock aquifer constant
% Kres_Rock(Kres_Rock>600)=600;
%Kres_Rock = 500; %%[h] Bedrock aquifer constant
Kres_Rock = Geo.Kres_Rock; %%[h] Bedrock aquifer constant

% Kres_Rock = Geo.Kres_Rock*0.2; %%[h] Bedrock aquifer constant

rho_g = Geo.rho_g;

Slo_top(Yout,Xout)=0.05;
[m_cell,n_cell]=size(DTM);
%num_cell=numel(DTM);
%x_cell=xllcorner:cellsize:(xllcorner+cellsize*(n_cell-1));
%y_cell=yllcorner:cellsize:(yllcorner+cellsize*(m_cell-1));
MASK=~isnan(DTM);
MASKn=reshape(MASK,num_cell,1);
Kinde = find(MASK==1);
Zsmax =Geo.Zsmax ; 
Slo_top(Slo_top<0.0005)=0.0005;
%ksv = reshape(ksv,num_cell,1);
Aspect(isnan(Aspect))=0;
Asur=(1./cos(atan(Slo_top))); %% Effective Area / Projected Area
Asur=reshape(Asur,num_cell,1);
aTop= 1000*ones(m_cell,n_cell).*(cellsize.^2)./cellsize; %%[mm] Area/Contour length ratio
Ared=ones(num_cell,1);
% Area= (cellsize^2)*sum(sum(MASK)); %% Projected area [m^2]
if isscalar(cellsize) 
    Area= cellsize.^2*sum(sum(MASK)); %% Projected area [m^2]
elseif numel(cellsize) > 1
    Area= sum((cellsize.*MASK).^2,'all')*sum(sum(MASK)); %% Projected area [m^2]
end

% Geo.KBOT = 0.45*ones(size(Geo.DTM,1),size(Geo.DTM,2)).*MASK; 
KBOT = Geo.KBOT;
% Out

Xout = Out.Xout; % points to track
Yout = Out.Yout;
npoint = length(Xout);

% load(['Tracked_points_' basin_name '.mat']);
% P_id=[P_id];
clear choose_point

% Param
cc_max = Param.cc_max; %% Number of vegetation
ms_max = Param.ms_max; %% Number of soil layers
md_max = Param.md_max; %% % Number of debris layers

EOS_maize_main = floor(nanmean(Geo.EOS_maize_main,"all"));
EOS_sc = floor(nanmean(Geo.EOS_sc,"all"));
EOS_tc =  floor(nanmean(Geo.EOS_tc,"all"));
EOS_wc=  floor(nanmean(Geo.EOS_wc,"all"));
SOS_maize_main=  floor(nanmean(Geo.SOS_maize_main,"all"));
SOS_sc=  floor(nanmean(Geo.SOS_sc,"all"));
SOS_tc=  floor(nanmean(Geo.SOS_tc,"all"));
SOS_wc=  floor(nanmean(Geo.SOS_wc,"all"));
TITLE_SAVE = Out.OutDir;
OutDir = Out.OutDir;
% EVcode = Out.EVcode;
% Ccrown_OUT = Out.Ccrown_OUT;
mkdir(OutDir)
%% Run

Oa= 210000;% Intercellular Partial Pressure Oxygen [umolO2/mol] -
dt=3600; %%[s] %%%
dth=1; %%[h]

[YE,MO,DA,HO,~,~] = datevec(Date);
Datam(:,1) = YE; Datam(:,2)= MO; Datam(:,3)= DA; Datam(:,4)= HO;

% L_day=zeros(size(Lon,1),size(Lon,2), length(2:24:length(Datam)));
[Lat, Lon] =  projinv(proj, Geo.x, Geo.y);
L_day=zeros(size(Lon,1),size(Lon,2), 366);
for li = 1:size(Lon,1)
    for lj = 1:size(Lon,2)
        for j = 1:366

            D_doy = datenum([2004 1 1])+(j-1);
            [~,~,~,~,~,L_day(li, lj, j)]= SetSunVariables(datevec(D_doy),DeltaGMT,Lon(li,lj),Lat(li,lj),t_bef,t_aft);
        end
    end
end

Lmax_day = max(L_day, [], 'all');
L_day = reshape(L_day,[size(L_day,1)*size(L_day,2) 366]);

%%%%%%%%%%%%%%%%%%%%

%%%%% Flow potential
T_pot=cell(1,ms_max);
for jk=1:ms_max
    T_pot{jk}= T_flow;
end
%%%%%%%%%%%%%%%%%%%%% -- here

% WC=WC.*MASK;
WC = zeros(size(Geo.SN));
WC(SN==1) = Geo.WC(SN==1);
SN(isnan(SN))=0; %% [Stream Identifier]
SNn=logical(reshape(SN,num_cell,1));
SN = logical(SN);


Pss = [800];
Pwp = [3500]; %%% [kPa]
Kfc = 0.2; %% [mm/h]
Phy = 10000; %% [kPa]
[Osat,L,Pe,Ks,O33]=Soil_parameters_spatial(PSAN,PCLA,PORG);
[Ofc,Oss,Owp,Ohy]=Soil_parametersII_spatial(Osat,L,Pe,Ks,O33,Kfc,Pss,Pwp,Phy);
%%%%%%%%%%%%%%%%%%%%
clear Pss Pwp Kfc Phy Osat L Pe Ks O33 Oss Owp
Ohy_OUT =   Ohy.*MASK;
Ohy_OUT=reshape(Ohy_OUT,num_cell,1);
%%%%%%%%%%%%%%%%
PSAN=reshape(PSAN,num_cell,1);
PCLA=reshape(PCLA,num_cell,1);
PORG=reshape(PORG,num_cell,1);


Zs_OUT = 1000*ones(num_cell,1); %%[mm]

%Ohy_OUT = 0.1154*ones(num_cell,1); %%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% SOLAR PARAMETER     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Computation Horizon Angle
%[HZ,Zasp] = Horizon_angle_polar(DTM,cellsize);
[HZ,Zasp] = Horizon_Angle(DTM, mean(cellsize,'all','omitnan'));
%%% HZ Horizon angle array [angular degree]
%%% Z Azimuth directions  [angualr degree] from N
%%% Sky View Factor and Terrain Configuration Factor
[SvF, Ct] = Sky_View_Factor(DTM,atan(Slo_top)*180/pi,Aspect,HZ,Zasp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      TEMPORAL INITIALIZATION  %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% TOTAL INITIAL CONDITION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% PREVIOUS CONDITION %%%%%%%%%%%%%%%%%%%%%
Ta_t=zeros(num_cell,24);
PAR_t=zeros(num_cell,24);
An_H_t=zeros(num_cell,cc_max,24);
An_L_t=zeros(num_cell,cc_max,24);
Rdark_H_t=zeros(num_cell,cc_max,24);
Rdark_L_t=zeros(num_cell,cc_max,24);
Tdp_H_t=zeros(num_cell,cc_max,24);
Tdp_L_t=zeros(num_cell,cc_max,24);
Psi_x_H_t=zeros(num_cell,cc_max,24);
Psi_l_H_t=zeros(num_cell,cc_max,24);
Psi_x_L_t=zeros(num_cell,cc_max,24);
Psi_l_L_t=zeros(num_cell,cc_max,24);
Qi_out_Rout=zeros(num_cell,cc_max,24);
%%%
Tdp_t=zeros(num_cell,ms_max,24);
O_t=zeros(num_cell,ms_max,24);
V_t=zeros(num_cell,ms_max,24);
Pr_sno_t=zeros(num_cell,24);
%%%

% ij  --> Cells [num_cell]
% t, tday --> Time step [hour] [day]
Fstep=strcat('Final_reached_step_',OutDir);
NutrientResetLog=zeros(0,15);
NutrientResetLogColumns={'iteration','datenum','year','month','day','hour',...
    'cell_linear','cell_row','cell_col','vegetation_layer','pft_index',...
    'reason_code','rNc_at_reset','rPc_at_reset','rKc_at_reset'};
NutrientResetLayerCodes={'1 = high vegetation (H)','2 = low vegetation (L)'};
NutrientResetReasonCodes={'1 = rNc < 0.64','2 = rPc < 0.64','3 = rKc < 0.64',...
    '4 = non-finite nutrient ratio or e_relN','5 = e_relN > 100',...
    '6 = rNc > 40','7 = rPc > 40','8 = rKc > 40',...
    '9 = initial blacklist reset'};
NutrientResetEver_H=false(num_cell,cc_max);
NutrientResetEver_L=false(num_cell,cc_max);
NutrientResetUniqueCellCount=0;
NutrientResetUniquePftCount=0;
NutrientResetLogFile=fullfile(OutDir,'Nutrient_reset_log.mat');
save(NutrientResetLogFile,'NutrientResetLog','NutrientResetLogColumns',...
    'NutrientResetLayerCodes','NutrientResetReasonCodes','NutrientResetEver_H','NutrientResetEver_L',...
    'NutrientResetUniqueCellCount','NutrientResetUniquePftCount');

InitialNutrientResetBlacklist=[];
if isfield(Param,'InitialNutrientResetBlacklist')
    InitialNutrientResetBlacklist=Param.InitialNutrientResetBlacklist;
else
    if isnumeric(basin_name)
        basin_name_text=num2str(basin_name);
    else
        basin_name_text=char(basin_name);
    end
    if strcmp(strtrim(basin_name_text),'01054200')
        InitialNutrientResetBlacklist=[...
            2633 1 4;...
            2743 1 1;...
            2865 1 2;...
            2865 1 4;...
            2865 1 5;...
            2866 1 1;...
            2866 1 4;...
            2866 1 5;...
            2981 1 1;...
            2981 1 4;...
            2984 1 1;...
            2984 1 4;...
            2984 1 5;...
            2985 1 4;...
            2985 1 5;...
            3092 1 1;...
            3209 1 1;...
            3326 1 1;...
            3326 1 2;...
            3326 1 4;...
            3457 1 4;...
            3573 1 2;...
            3574 1 2;...
            3574 1 4;...
            6038 1 1;...
            6038 1 2;...
            6038 1 4;...
            6039 1 2;...
            6039 1 4;...
            6041 1 3;...
            6041 1 4;...
            6042 1 1;...
            6042 1 2;...
            6042 1 4;...
            6495 1 2;...
            6495 1 4;...
            7086 1 1;...
            7086 1 2;...
            7086 1 4;...
            8500 1 5];
    end
end

RoutingTracerDiagnosticCell=[];
if isfield(Param,'RoutingTracerDiagnosticCell')
    RoutingTracerDiagnosticCell=Param.RoutingTracerDiagnosticCell;
else
    if isnumeric(basin_name)
        routing_diagnostic_basin=num2str(basin_name);
    else
        routing_diagnostic_basin=char(basin_name);
    end
    if strcmp(strtrim(routing_diagnostic_basin),'01054200')
        RoutingTracerDiagnosticCell=7086;
    end
end
RoutingTracerDiagnosticFid=-1;
if ~isempty(RoutingTracerDiagnosticCell)
    if ~isscalar(RoutingTracerDiagnosticCell) || RoutingTracerDiagnosticCell < 1 || ...
            RoutingTracerDiagnosticCell > num_cell || RoutingTracerDiagnosticCell ~= floor(RoutingTracerDiagnosticCell)
        error('Param.RoutingTracerDiagnosticCell must be a valid scalar linear cell index.')
    end
    RoutingTracerNames={'nh4','no3','p','k','don','dop','doc'};
    RoutingPoolNames={'nh4','no3','p','k','don','dop','doc_labile','doc_recalcitrant'};
    RoutingTracerDiagnosticColumns={'iteration','datenum','year','month','day','hour','cell_linear',...
        'pr_input','pr_liq','q_runon_in','q_run_after_e','wis','infiltration','rh','rd','eg',...
        'v_sum','vice_sum','qi_in_sum','qi_out_sum','q_runon_out'};
    for routingTracerIndex=1:length(RoutingTracerNames)
        routingTracerName=RoutingTracerNames{routingTracerIndex};
        RoutingTracerDiagnosticColumns=[RoutingTracerDiagnosticColumns,...
            {['tr_sur_in_' routingTracerName],['tr_if_' routingTracerName],['tr_avg_pre_' routingTracerName],...
            ['tr_avg_post_' routingTracerName],['tr_sur_out_' routingTracerName],['mass_infil_' routingTracerName],...
            ['mass_et_' routingTracerName],['mass_rd_loss_' routingTracerName],['mass_leak_' routingTracerName],...
            ['mass_ice_' routingTracerName],['mass_sub_out_' routingTracerName],['mass_sub_in_' routingTracerName]}];
    end
    for routingPoolIndex=1:length(RoutingPoolNames)
        routingPoolName=RoutingPoolNames{routingPoolIndex};
        RoutingTracerDiagnosticColumns=[RoutingTracerDiagnosticColumns,...
            {['p_bg_before_' routingPoolName],['p_bg_after_' routingPoolName],['p_bg_net_' routingPoolName],...
            ['db_flux_sum_' routingPoolName],['p_initial_' routingPoolName],['p_final_' routingPoolName],...
            ['p_routing_net_' routingPoolName]}];
    end
    RoutingTracerDiagnosticColumns=[RoutingTracerDiagnosticColumns,{'ice_volume_change'}];
    RoutingTracerDiagnosticFile=fullfile(OutDir,['Routing_tracer_diagnostic_cell' num2str(RoutingTracerDiagnosticCell) '.csv']);
    RoutingTracerDiagnosticFid=fopen(RoutingTracerDiagnosticFile,'w');
    if RoutingTracerDiagnosticFid < 0
        error('Cannot create routing tracer diagnostic file: %s',RoutingTracerDiagnosticFile)
    end
    fprintf(RoutingTracerDiagnosticFid,'%s\n',strjoin(RoutingTracerDiagnosticColumns,','));
end

tday=0;
%%%%%%%%%%%%%%%% NUMERICAL METHODS OPTIONS
Opt_CR = optimset('TolFun',1);%,'UseParallel','always');
Opt_ST = optimset('TolFun',0.1);%,'UseParallel','always');
Opt_ST2 = NaN;%optimset('TolFun',0.1,'Display','off');
OPT_SM=  odeset('AbsTol',0.05,'MaxStep',dth);
OPT_VD=  odeset('AbsTol',0.05);
OPT_STh = odeset('AbsTol',5e+3);
OPT_PH= odeset('AbsTol',0.01);
OPT_min_SPD = Inf; %% [m] minimum snow pack depth to have a multilayer snow
OPT_VegSnow = 1;
OPT_SoilTemp = 1;
OPT_FR_SOIL =0;
OPT_PlantHydr = 0;
OPT_EnvLimitGrowth = 0;
OPT_VCA=0;
OPT_ALLOME=0;
%%%
OPT_HEAD = 0;
OPT_SoilBiogeochemistry = 1;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% BG added %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if OPT_SoilBiogeochemistry == 1
    P=zeros(num_cell,55);
    Ptm1=zeros(num_cell,55);
    soil_ph_250=load([pwd filesep 'geophysical' filesep 'soil_ph_250.mat'],'soil_ph_250');   % 读入 lai_t_250, lai_time_250
    PHs=nanmean(nanmean(soil_ph_250.soil_ph_250));
    Upl=0.01;
    N_nt=7;
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% SOIL MOISTURE
Param.INI()
DemoNutrientResetInjected=false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic ;
%profile on
%bau = waitbar(0,'Waiting...');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for t=2:N_time_step

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %waitbar(t/N_time_step,bau)
    disp('Iter:'); disp(t);
    %%%%%%%%%%%%%%%% DATE and METEO SPATIAL;
    Datam_S=Datam(t,:);

    % Load current month - Should overwrite older nc files
    % if Datam_S(3)==1 && Datam_S(4)==0 && t>2
    % 
    %     Load_EE_ERA5(Date_dv(t,1), Date_dv(t,2), Geo.shapefile_path, Meteo.Meteo_sc, WorldClim, Meteo.FolderLoc, t_bef, t_aft, DeltaGMT )
    %     % Perform radiation partitioning 
    %     ix_rd = 1; %netcdf Reading index
    % end

    if (Datam_S(4)==1)
        tday = tday+1;
    end
    HRL = Datam_S(4);
    [jDay]= JULIAN_DAY(Datam_S);

    if isfield(Param,'DemoNutrientReset') && Param.DemoNutrientReset.enabled && ...
            ~DemoNutrientResetInjected && HRL==1
        demoReset=Param.DemoNutrientReset;
        if demoReset.cell < 1 || demoReset.cell > num_cell || ...
                demoReset.pft < 1 || demoReset.pft > cc_max || ~MASKn(demoReset.cell)
            error('DemoNutrientReset specifies an invalid cell or PFT.');
        end
        if ~isfield(demoReset,'reserve_value')
            demoReset.reserve_value=1e9;
        end
        if strcmpi(demoReset.layer,'H')
            Nreserve_Htm1(demoReset.cell,demoReset.pft)=demoReset.reserve_value;
        elseif strcmpi(demoReset.layer,'L')
            Nreserve_Ltm1(demoReset.cell,demoReset.pft)=demoReset.reserve_value;
        else
            error('DemoNutrientReset.layer must be H or L.');
        end
        DemoNutrientResetInjected=true;
        fprintf('Demo nutrient reset injected at Iter %d: cell %d, layer %s, PFT %d.\n',...
            t,demoReset.cell,upper(demoReset.layer),demoReset.pft);
    end

    % ASSUMPTION THAT NEEDS TO BE RELAXED FOR BIG SIMULATIONS
    % For big sims, we neeed to change to a local DEM for different time
    % zones
    [h_S,~,zeta_S,~,~,~]= SetSunVariables(Datam_S,DeltaGMT,mean(Lon(MASK),'all'),mean(Lat(MASK),'all'),t_bef,t_aft);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [ShF] = Shadow_Effect(DTM,h_S,zeta_S,HZ,Zasp);

    %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% INPUT  Vector [num_cell]
    %Datam_S Pr_S,Ta_S,Ds_S,Ws_S ,ea_S,N_S,Pre_S, Tdew_S
    % SAB1_S , SAB2_S,SAD1_S,SAD2_S,PARB_S,PARD_S,SvF
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%% SPATIAL DISTRIBUTION OF INPUT
    %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PARD       Rsw_Rif    SAD1       SD         U_Rif      esat_Rif   qta_T
    % Ds_Rif     Pr         SAB1       SAD2       Ta_Rif     Ws_Rif     qta_P
    % PARB       Pre_Rif    SAB2       SB         Tdew_Rif   ea_Rif
    % qta_Pr
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if exist([Meteo.FolderLoc, 'IrD_S.nc'], 'file')
        IrD_S = ncread( [Meteo.FolderLoc, 'IrD_S.nc'], 'Ir', [1, 1, ix_rd], [Inf, Inf, 1] );
    else
        IrD_S =  MASK*0; IrD_S=reshape(IrD_S,num_cell,1); %don't
    end

    if isfield(Geo, 'Salt_S')
        Salt_S = Geo.Salt_S;Salt_S=reshape(Salt_S,num_cell,1);
    else
        Salt_S =  MASK*0; Salt_S=reshape(Salt_S,num_cell,1); %don't
    end

    %% LOAD AND INTERPOLATE DATA TO PROPER GRID

    if Meteo.Interp == 0% no interpolation needed

        Ta_S = ncread( [Meteo.FolderLoc, 'Ta.nc'], 'Ta', [1, 1, ix_rd], [Inf, Inf, 1] );
        U_S = ncread( [Meteo.FolderLoc, 'U.nc'], 'U', [1, 1, ix_rd], [Inf, Inf, 1] );
        Pr_S = ncread( [Meteo.FolderLoc, 'Pr.nc'], 'Pr', [1, 1, ix_rd], [Inf, Inf, 1] );
        Pre_S = ncread( [Meteo.FolderLoc, 'Pre.nc'], 'Pre', [1, 1, ix_rd], [Inf, Inf, 1] );
        N_S = ncread( [Meteo.FolderLoc, 'N.nc'], 'N', [1, 1, ix_rd], [Inf, Inf, 1] );

        esat_S = (610.94*exp( (17.625*Ta_S./(Ta_S+243.04)) ));

        ea_S = esat_S.*U_S;
        %ea_S = ncread( [Meteo.FolderLoc, 'ea_S.nc'], 'ea', [1, 1, ix_rd], [Inf, Inf, 1] );
        Ds_S = esat_S - ea_S;
        Ws_S = ncread( [Meteo.FolderLoc, 'Ws.nc'], 'Ws', [1, 1, ix_rd], [Inf, Inf, 1] );

        Tdew_S = 243.04*( log(U_S) + 17.625.*Ta_S./(243.04+Ta_S) )./(( 17.625 - log(U_S) - 17.625.*Ta_S./(243.04+Ta_S)));
        %Tdew_S = ncread( [Meteo.FolderLoc, 'Tdew_S.nc'], 'Tdew', [1, 1, ix_rd], [Inf, Inf, 1] );

        SAB1_S = ncread( [Meteo.FolderLoc, 'SAB1.nc'], 'SAB1', [1, 1, ix_rd], [Inf, Inf, 1] );
        SAB2_S = ncread( [Meteo.FolderLoc, 'SAB2.nc'], 'SAB2', [1, 1, ix_rd], [Inf, Inf, 1] );
        SAD1_S = ncread( [Meteo.FolderLoc, 'SAD1.nc'], 'SAD1', [1, 1, ix_rd], [Inf, Inf, 1] );
        SAD2_S = ncread( [Meteo.FolderLoc, 'SAD2.nc'], 'SAD2', [1, 1, ix_rd], [Inf, Inf, 1] );

        PARB_S = ncread( [Meteo.FolderLoc, 'PARB.nc'], 'PARB', [1, 1, ix_rd], [Inf, Inf, 1] );
        PARD_S = ncread( [Meteo.FolderLoc, 'PARD.nc'], 'PARD', [1, 1, ix_rd], [Inf, Inf, 1] );
        Ca_S = ncread( [Meteo.FolderLoc, 'Ca.nc'], 'Ca', [1, 1, ix_rd], [Inf, Inf, 1] );

    elseif Meteo.Interp == 1

        % perform a nearst neighbor interpolation from data to grid.
        % Compute only once the indices
        if t == 2
            I_interp = knnsearch([Meteo.x(:), Meteo.y(:)], [Geo.x(:), Geo.y(:)]);
        end
        % =========================================================
        % read a new block only when current ix_rd is outside cache
        % =========================================================
        if ix_rd < blk_start_ix || ix_rd > blk_end_ix

            blk_start_ix = ix_rd;
            blk_end_ix   = min(ix_rd + blockSize - 1, x1 + N_time_step - 1);
            blk_count    = blk_end_ix - blk_start_ix + 1;

            start3 = [1, 1, blk_start_ix];
            count3 = [Inf, Inf, blk_count];

            Ta_blk   = ncread(fnTa,   'Ta',   start3, count3);
            U_blk    = ncread(fnU,    'U',    start3, count3);
            N_blk    = ncread(fnN,    'N',    start3, count3);
            Pr_blk   = ncread(fnPr,   'Pr',   start3, count3);
            Pre_blk  = ncread(fnPre,  'Pre',  start3, count3);
            Ws_blk   = ncread(fnWs,   'Ws',   start3, count3);

            SAB1_blk = ncread(fnSAB1, 'SAB1', start3, count3);
            SAB2_blk = ncread(fnSAB2, 'SAB2', start3, count3);
            SAD1_blk = ncread(fnSAD1, 'SAD1', start3, count3);
            SAD2_blk = ncread(fnSAD2, 'SAD2', start3, count3);

            PARB_blk = ncread(fnPARB, 'PARB', start3, count3);
            PARD_blk = ncread(fnPARD, 'PARD', start3, count3);
            Ca_blk   = ncread(fnCa,   'Ca',   start3, count3);
        end

        % local index within current block
        ib = ix_rd - blk_start_ix + 1;

        % extract one time slice from cached block
        Ta_S_   = Ta_blk(:, :, ib);
        U_S_    = U_blk(:, :, ib);
        N_S_    = N_blk(:, :, ib);
        Pr_S_   = Pr_blk(:, :, ib);
        Pre_S_  = Pre_blk(:, :, ib);
        Ws_S_   = Ws_blk(:, :, ib);

        SAB1_S_ = SAB1_blk(:, :, ib);
        SAB2_S_ = SAB2_blk(:, :, ib);
        SAD1_S_ = SAD1_blk(:, :, ib);
        SAD2_S_ = SAD2_blk(:, :, ib);

        PARB_S_ = PARB_blk(:, :, ib);
        PARD_S_ = PARD_blk(:, :, ib);
        Ca_S_   = Ca_blk(:, :, ib);

        if nanmean(nanmean(Pr_S_)) == 0
            Pr_S_(:) = Pr_sel(t);
        else
            Pr_cr_factor = Pr_sel(t) / nanmean(nanmean(Pr_S_));
            Pr_S_(:) = Pr_S_(:) .* Pr_cr_factor;
        end

        esat_S_ = (610.94 * exp((17.625 * Ta_S_ ./ (Ta_S_ + 243.04))));
        ea_S_   = esat_S_ .* U_S_;
        Tdew_S_ = 243.04 * (log(U_S_) + 17.625 .* Ta_S_ ./ (243.04 + Ta_S_)) ./ ...
                 ((17.625 - log(U_S_) - 17.625 .* Ta_S_ ./ (243.04 + Ta_S_)));

        Pre_S_ = fillmissing(Pre_S_, 'nearest');
        Ws_S_  = fillmissing(Ws_S_, 'nearest');
        Pr_S_  = fillmissing(Pr_S_, 'constant', 0);
        Pr_S_(Pr_S_ < 0.01) = 0;
        U_S_   = fillmissing(U_S_, 'nearest');
        Ta_S_  = fillmissing(Ta_S_, 'nearest');


        Ta_S = zeros(size(DTM));
        %U_S = zeros(size(DTM));
        Pr_S = zeros(size(DTM));
        Pre_S = zeros(size(DTM));
        esat_S = zeros(size(DTM));
        ea_S = zeros(size(DTM));
        Ws_S = zeros(size(DTM));
        Tdew_S = zeros(size(DTM));

        SAB1_S = zeros(size(DTM));
        SAB2_S = zeros(size(DTM));
        SAD1_S = zeros(size(DTM));
        SAD2_S = zeros(size(DTM));

        PARB_S = zeros(size(DTM));
        PARD_S = zeros(size(DTM));
        Ca_S = zeros(size(DTM));
        N_S = zeros(size(DTM));

        nml = numel(DTM);
        Ta_S(1:nml) = Ta_S_(I_interp);
        %Ta_S = Ta_S;
        U_S(1:nml) = U_S_(I_interp);
        Pr_S(1:nml) = Pr_S_(I_interp);
        Pre_S(1:nml) = Pre_S_(I_interp);
        esat_S(1:nml) = esat_S_(I_interp);
        Ws_S(1:nml) = Ws_S_(I_interp);
        Tdew_S(1:nml) = Tdew_S_(I_interp);
        N_S(1:nml) = N_S_(I_interp);
        ea_S(1:nml) = ea_S_(I_interp);
        Ds_S = esat_S - ea_S;

        SAB1_S(1:nml) = SAB1_S_(I_interp);
        SAB2_S(1:nml) = SAB2_S_(I_interp);
        SAD1_S(1:nml) = SAD1_S_(I_interp);
        SAD2_S(1:nml) = SAD2_S_(I_interp);

        PARB_S(1:nml) = PARB_S_(I_interp);
        PARD_S(1:nml) = PARD_S_(I_interp);
        Ca_S(1:nml) = Ca_S_(I_interp);

    elseif Meteo.Interp == 2

        Ta_S_ = ncread( [Meteo.FolderLoc, filesep 'Ta.nc'], 'Ta', [1, 1, ix_rd], [Inf, Inf, 1] );
        U_S_ = ncread( [Meteo.FolderLoc, filesep 'U.nc'], 'U', [1, 1, ix_rd], [Inf, Inf, 1] );
        N_S_ = ncread( [Meteo.FolderLoc, filesep 'N.nc'], 'N', [1, 1, ix_rd], [Inf, Inf, 1] );
        Pr_S_ = ncread( [Meteo.FolderLoc, filesep 'Pr.nc'], 'Pr', [1, 1, ix_rd], [Inf, Inf, 1] );
        Pre_S_ = ncread( [Meteo.FolderLoc, filesep 'Pre.nc'], 'Pre', [1, 1, ix_rd], [Inf, Inf, 1] );
        %esat_S_ = (610.94*exp( (17.625*Ta_S_./(Ta_S_+243.04)) ));
        %ea_S_ = esat_S_.*U_S_;
        %Tdew_S_ = 243.04*( log(U_S_) + 17.625.*Ta_S_./(243.04+Ta_S_) )./(( 17.625 - log(U_S_) - 17.625.*Ta_S_./(243.04+Ta_S_)));
        Ws_S_ = ncread( [Meteo.FolderLoc, filesep 'Ws.nc'], 'Ws', [1, 1, ix_rd], [Inf, Inf, 1] );

        SAB1_S_ = ncread( [Meteo.FolderLoc, filesep 'SAB1.nc'], 'SAB1', [1, 1, ix_rd], [Inf, Inf, 1] );
        SAB2_S_ = ncread( [Meteo.FolderLoc, filesep 'SAB2.nc'], 'SAB2', [1, 1, ix_rd], [Inf, Inf, 1] );
        SAD1_S_ = ncread( [Meteo.FolderLoc, filesep 'SAD1.nc'], 'SAD1', [1, 1, ix_rd], [Inf, Inf, 1] );
        SAD2_S_ = ncread( [Meteo.FolderLoc, filesep 'SAD2.nc'], 'SAD2', [1, 1, ix_rd], [Inf, Inf, 1] );

        PARB_S_ = ncread( [Meteo.FolderLoc, filesep 'PARB.nc'], 'PARB', [1, 1, ix_rd], [Inf, Inf, 1] );
        PARD_S_ = ncread( [Meteo.FolderLoc, filesep 'PARD.nc'], 'PARD', [1, 1, ix_rd], [Inf, Inf, 1] );
        Ca_S_ = ncread( [Meteo.FolderLoc, filesep 'Ca.nc'], 'Ca', [1, 1, ix_rd], [Inf, Inf, 1] );



        Fint = scatteredInterpolant(Meteo.x(:), Meteo.y(:), Ta_S_(:),'natural');
        Ta_S = Fint(Geo.x, Geo.y);
        Ta_S = Ta_S+LS;

        Fint = scatteredInterpolant(Meteo.x(:), Meteo.y(:), U_S_(:),'natural');
        U_S = Fint(Geo.x, Geo.y);

        esat_S = (610.94*exp( (17.625*Ta_S./(Ta_S+243.04)) ));
        ea_S = esat_S.*U_S;
        Tdew_S = 243.04*( log(U_S) + 17.625.*Ta_S./(243.04+Ta_S) )./(( 17.625 - log(U_S) - 17.625.*Ta_S./(243.04+Ta_S)));
        Ds_S = esat_S - ea_S;

        Fint = scatteredInterpolant(Meteo.x(:), Meteo.y(:), Pr_S_(:),'nearest');
        Pr_S = Fint(Geo.x, Geo.y);

        Fint = scatteredInterpolant(Meteo.x(:), Meteo.y(:), N_S_(:),'natural');
        N_S = Fint(Geo.x, Geo.y);

        Fint = scatteredInterpolant(Meteo.x(:), Meteo.y(:), Pre_S_(:),'natural');
        Pre_S = Fint(Geo.x, Geo.y);

        Fint = scatteredInterpolant(Meteo.x(:), Meteo.y(:), Ws_S_(:),'natural');
        Ws_S = Fint(Geo.x, Geo.y);

        Fint = scatteredInterpolant(Meteo.x(:), Meteo.y(:), SAB1_S_(:),'natural');
        SAB1_S = Fint(Geo.x, Geo.y);

        Fint = scatteredInterpolant(Meteo.x(:), Meteo.y(:), SAB2_S_(:),'natural');
        SAB2_S = Fint(Geo.x, Geo.y);

        Fint = scatteredInterpolant(Meteo.x(:), Meteo.y(:), SAD1_S_(:),'natural');
        SAD1_S = Fint(Geo.x, Geo.y);

        Fint = scatteredInterpolant(Meteo.x(:), Meteo.y(:), SAD2_S_(:),'natural');
        SAD2_S = Fint(Geo.x, Geo.y);

        Fint = scatteredInterpolant(Meteo.x(:), Meteo.y(:), PARB_S_(:),'natural');
        PARB_S = Fint(Geo.x, Geo.y);

        Fint = scatteredInterpolant(Meteo.x(:), Meteo.y(:), PARD_S_(:),'natural');
        PARD_S = Fint(Geo.x, Geo.y);

        Fint = scatteredInterpolant(Meteo.x(:), Meteo.y(:), Ca_S_(:),'natural');
        Ca_S = Fint(Geo.x, Geo.y);


    elseif Meteo.Interp == 3

    end

    N_S(N_S > 1) = 1;
    N_S(N_S < 0) = 0;
    %
    % U_S(U_S > 1) = 1;
    % U_S(U_S < 0) = 0;

    Ws_S(Ws_S < 0.1) = 0.1;

    %%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cos_fst = cos(atan(Slo_top))*sin(h_S) + sin(atan(Slo_top)).*cos(h_S).*cos(zeta_S-Aspect*pi/180);
    cos_fst(cos_fst<0)=0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if sin(h_S) <= 0.10 %%[5.73°] Numerical problems
        SAB1_S  =  0*MASKn;
        SAB2_S  =  0*MASKn;
        SAD1_S  =  0*MASKn;
        SAD2_S  =  0*MASKn;
        PARB_S  =  0*MASKn;
        PARD_S  =  0*MASKn;
    else
        SAD1_S = SAD1_S.*SvF + Ct.*rho_g.*((SAB1_S./sin(h_S)).*cos_fst + (1-SvF).*SAD1_S);
        SAD2_S = SAD2_S.*SvF + Ct.*rho_g.*((SAB2_S./sin(h_S)).*cos_fst + (1-SvF).*SAD2_S);
        PARD_S = PARD_S.*SvF + Ct*rho_g.*((PARB_S./sin(h_S)).*cos_fst + (1-SvF).*PARD_S);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        SAB1_S =(SAB1_S./sin(h_S)).*cos_fst.*ShF;
        SAB2_S =(SAB2_S./sin(h_S)).*cos_fst.*ShF;
        PARB_S = (PARB_S./sin(h_S)).*cos_fst.*ShF;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SAB1_S  =  reshape(SAB1_S,num_cell,1);
        % SAB2_S  =  reshape(SAB2_S,num_cell,1);
        % SAD1_S  =  reshape(SAD1_S,num_cell,1);
        % SAD2_S  =  reshape(SAD2_S,num_cell,1);
        % PARB_S  =  reshape(PARB_S,num_cell,1);
        % PARD_S  =  reshape(PARD_S,num_cell,1);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% SPATIAL INITIALIZATION VECTOR PREDEFINING %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if t == 2
        %%%%%%%% VEGETATION

        NPP_H=zeros(num_cell,cc_max);
        ANPP_H=zeros(num_cell,cc_max);
        Rg_H=zeros(num_cell,cc_max);
        RA_H=zeros(num_cell,cc_max);
        Rms_H=zeros(num_cell,cc_max);
        Rmr_H=zeros(num_cell,cc_max);
        Rmc_H=zeros(num_cell,cc_max);
        LAIdead_H=zeros(num_cell,cc_max);
        Sr_H=zeros(num_cell,cc_max);
        Slf_H=zeros(num_cell,cc_max);
        Sfr_H=zeros(num_cell,cc_max);
        Swm_H=zeros(num_cell,cc_max);
        Sll_H=zeros(num_cell,cc_max);
        Rexmy_H= zeros(num_cell,cc_max,3);
        Rrootl_H=zeros(num_cell,cc_max);
        Bfac_dayH=zeros(num_cell,cc_max);
        Bfac_weekH=Bfac_weekHtm1;
        NPPI_H= NPPI_Htm1;
        RB_H = zeros(num_cell,cc_max,7);
        rNc_H =zeros(num_cell,cc_max);
        rPc_H =zeros(num_cell,cc_max);
        rKc_H = zeros(num_cell,cc_max);
        ManIH = zeros(num_cell,cc_max);
        TBio_H= zeros(num_cell,cc_max);
        %%%
        NPP_L=zeros(num_cell,cc_max);
        ANPP_L=zeros(num_cell,cc_max);
        Rg_L=zeros(num_cell,cc_max);
        RA_L=zeros(num_cell,cc_max);
        Rms_L=zeros(num_cell,cc_max);
        Rmr_L=zeros(num_cell,cc_max);
        Rmc_L=zeros(num_cell,cc_max);
        LAIdead_L=zeros(num_cell,cc_max);
        Sr_L=zeros(num_cell,cc_max);
        Slf_L=zeros(num_cell,cc_max);
        Sfr_L=zeros(num_cell,cc_max);
        Swm_L=zeros(num_cell,cc_max);
        Sll_L=zeros(num_cell,cc_max);
        Rexmy_L= zeros(num_cell,cc_max,3);
        Rrootl_L=zeros(num_cell,cc_max);
        Bfac_dayL =zeros(num_cell,cc_max);
        Bfac_weekL=Bfac_weekLtm1;
        NPPI_L=NPPI_Ltm1;
        RB_L  = zeros(num_cell,cc_max,7);
        rNc_L=zeros(num_cell,cc_max);
        rPc_L=zeros(num_cell,cc_max);
        rKc_L=zeros(num_cell,cc_max);
        ManIL=zeros(num_cell,cc_max);
        % %%%%
        % hc_Htm1(1,:) = hc_H;
        % hc_Ltm1(1,:) = hc_L;
        e_rel_H=e_rel_Htm1;
        e_relN_H = e_relN_Htm1;
        LAI_H=LAI_Htm1;
        B_H= B_Htm1;
        PHE_S_H=PHE_S_Htm1;
        dflo_H= dflo_Htm1;
        AgeL_H= AgeL_Htm1;
        AgeDL_H=AgeDL_Htm1;
        SAI_H= SAI_Htm1;
        hc_H= hc_Htm1;
        NBLeaf_H =NBLeaf_Htm1;
        Nreserve_H = Nreserve_Htm1;
        Preserve_H = Preserve_Htm1;
        Kreserve_H = Kreserve_Htm1;
        FNC_H =  FNC_Htm1;
        NupI_H = NupI_Htm1;
        TdpI_H= TdpI_Htm1;
        PARI_H =PARI_Htm1 ;
        NBLI_H  =NBLI_Htm1;
        %%%
        e_rel_L=e_rel_Ltm1;
        e_relN_L = e_relN_Ltm1;
        LAI_L= LAI_Ltm1;
        B_L= B_Ltm1;
        PHE_S_L= PHE_S_Ltm1;
        dflo_L= dflo_Ltm1 ;
        AgeL_L= AgeL_Ltm1 ;
        AgeDL_L =AgeDL_Ltm1;
        SAI_L= SAI_Ltm1;
        hc_L= hc_Ltm1;
        NBLeaf_L =NBLeaf_Ltm1;
        Nreserve_L = Nreserve_Ltm1;
        Preserve_L = Preserve_Ltm1;
        Kreserve_L = Kreserve_Ltm1;
        FNC_L =  FNC_Ltm1;
        NupI_L = NupI_Ltm1;
        TdpI_L= TdpI_Ltm1;
        PARI_L =PARI_Ltm1 ;
        NBLI_L  =NBLI_Ltm1;
        TBio_L= zeros(num_cell,cc_max);
        %%%%%
        Nuptake_H= zeros(num_cell,cc_max);
        Puptake_H= zeros(num_cell,cc_max);
        Kuptake_H= zeros(num_cell,cc_max);
        Nuptake_L= zeros(num_cell,cc_max);
        Puptake_L= zeros(num_cell,cc_max);
        Kuptake_L= zeros(num_cell,cc_max);
        NavlI = zeros(num_cell,3);
        Bam=zeros(num_cell,1);
        Bem=zeros(num_cell,1);
        %%%%
        TexC_H=zeros(num_cell,cc_max);
        TexN_H=zeros(num_cell,cc_max);
        TexP_H=zeros(num_cell,cc_max);
        TexK_H=zeros(num_cell,cc_max);
        TNIT_H=zeros(num_cell,cc_max);
        TPHO_H=zeros(num_cell,cc_max);
        TPOT_H=zeros(num_cell,cc_max);
        SupN_H=zeros(num_cell,cc_max);
        SupP_H=zeros(num_cell,cc_max);
        SupK_H=zeros(num_cell,cc_max);
        ISOIL_H= zeros(num_cell,cc_max,18);
        TexC_L=zeros(num_cell,cc_max);
        TexN_L=zeros(num_cell,cc_max);
        TexP_L=zeros(num_cell,cc_max);
        TexK_L=zeros(num_cell,cc_max);
        TNIT_L=zeros(num_cell,cc_max);
        TPHO_L=zeros(num_cell,cc_max);
        TPOT_L=zeros(num_cell,cc_max);
        SupN_L=zeros(num_cell,cc_max);
        SupP_L=zeros(num_cell,cc_max);
        SupK_L=zeros(num_cell,cc_max);
        ISOIL_L= zeros(num_cell,cc_max,18);
        BA_H =zeros(num_cell,cc_max);
        Tden_H=zeros(num_cell,cc_max);
        AgePl_H=zeros(num_cell,cc_max);
        BA_L =zeros(num_cell,cc_max);
        Tden_L =zeros(num_cell,cc_max);
        AgePl_L =zeros(num_cell,cc_max);
        Ccrown_t =zeros(num_cell,cc_max);
        Ccrown_t_tm1 = Ccrown_t;
        %%%
        NuLit_H= NuLit_Htm1;
        NuLit_L= NuLit_Ltm1;
        BLit=zeros(num_cell,cc_max);
        %%%%%%% HYDROLOGY
        V=zeros(num_cell,ms_max);
        O=zeros(num_cell,ms_max);
        ZWT=zeros(num_cell,1);
        POT=zeros(num_cell,ms_max);
        OF=zeros(num_cell,1);
        OS=zeros(num_cell,1);
        OH=zeros(num_cell,cc_max);
        OL=zeros(num_cell,cc_max);
        Rd=zeros(num_cell,1);
        Qi_out=zeros(num_cell,ms_max);
        Rh=zeros(num_cell,1);
        Lk=zeros(num_cell,1);
        f=zeros(num_cell,1);
        q_run_after_E=zeros(num_cell,1);
        WIS=zeros(num_cell,1);
        Ts=zeros(num_cell,1);
        Csno=zeros(num_cell,1);
        Cice=zeros(num_cell,1);
        Pr_sno=zeros(num_cell,1);
        Pr_liq=zeros(num_cell,1);
        rb_H=zeros(num_cell,cc_max);
        rb_L=zeros(num_cell,cc_max);
        rs_sunH=zeros(num_cell,cc_max);
        rs_sunL=zeros(num_cell,cc_max);
        rs_shdH=zeros(num_cell,cc_max);
        rs_shdL=zeros(num_cell,cc_max);
        rap_H=zeros(num_cell,cc_max);
        rap_L=zeros(num_cell,cc_max);
        r_soil=zeros(num_cell,1);
        b_soil=zeros(num_cell,1);
        alp_soil=zeros(num_cell,1);
        ra=zeros(num_cell,1);
        r_litter=zeros(num_cell,cc_max);
        WR_SP=zeros(num_cell,1);
        U_SWE=zeros(num_cell,1);
        NIn_SWE=zeros(num_cell,1);
        dQ_S=zeros(num_cell,1);
        DQ_S=zeros(num_cell,1);
        DT_S=zeros(num_cell,1);
        Dr_H=zeros(num_cell,cc_max);
        Dr_L=zeros(num_cell,cc_max);
        SE_rock=zeros(num_cell,1);
        SE_urb=zeros(num_cell,1);
        An_L=zeros(num_cell,cc_max);
        An_H=zeros(num_cell,cc_max);
        Rdark_L=zeros(num_cell,cc_max);
        Rdark_H=zeros(num_cell,cc_max);
        Ci_sunH=zeros(num_cell,cc_max);
        Ci_sunL=zeros(num_cell,cc_max);
        Ci_shdH=zeros(num_cell,cc_max);
        Ci_shdL=zeros(num_cell,cc_max);
        Rn=zeros(num_cell,1);
        H=zeros(num_cell,1);
        QE=zeros(num_cell,1);
        Qv=zeros(num_cell,1);
        Lpho=zeros(num_cell,1);
        T_H=zeros(num_cell,cc_max);
        T_L=zeros(num_cell,cc_max);
        EIn_H=zeros(num_cell,cc_max);
        EIn_L=zeros(num_cell,cc_max);
        EG=zeros(num_cell,1);
        ELitter=zeros(num_cell,1);
        ESN=zeros(num_cell,1);
        ESN_In=zeros(num_cell,1);
        EWAT=zeros(num_cell,1);
        EICE=zeros(num_cell,1);
        EIn_urb=zeros(num_cell,1);
        EIn_rock=zeros(num_cell,1);
        dw_SNO=zeros(num_cell,1);
        G=zeros(num_cell,1);
        Gfin=zeros(num_cell,1);
        Tdp=zeros(num_cell,ms_max);
        Tdpsnow=zeros(num_cell,ms_max);
        Tdeb=zeros(num_cell,md_max);
        Vice=zeros(num_cell,ms_max);
        Oice=zeros(num_cell,ms_max);
        Tice=zeros(num_cell,1);
        Imelt=zeros(num_cell,1);
        Smelt=zeros(num_cell,1);
        Tdamp=zeros(num_cell,1);
        Tdp_H=zeros(num_cell,cc_max);
        Tdp_L=zeros(num_cell,cc_max);
        SWE=zeros(num_cell,1);
        SND=zeros(num_cell,1);
        ros=zeros(num_cell,1);
        In_SWE=zeros(num_cell,1);
        SP_wc=zeros(num_cell,1);
        WAT=zeros(num_cell,1);
        ICE=zeros(num_cell,1);
        ICE_D=zeros(num_cell,1);
        IP_wc=zeros(num_cell,1);
        WR_IP=zeros(num_cell,1);
        NIce=zeros(num_cell,1);
        Cicew=zeros(num_cell,1);
        Csnow=zeros(num_cell,1);
        FROCK=zeros(num_cell,1);
        Lk_wat=zeros(num_cell,1);
        Lk_rock=zeros(num_cell,1);
        Qfm=zeros(num_cell,1);
        t_sls=zeros(num_cell,1);
        In_H=In_Htm1;
        In_L=In_Ltm1;
        In_Litter=In_Littertm1;
        In_urb=zeros(num_cell,1);
        In_rock=zeros(num_cell,1);
        er=zeros(num_cell,1);
        snow_albedo=zeros(num_cell,4);
        tau_sno=zeros(num_cell,1);
        e_sno=zeros(num_cell,1);
        EK=zeros(num_cell,1);
        dQVEG=zeros(num_cell,1);
        TsVEG=zeros(num_cell,1);
        Ts_under=zeros(num_cell,1);
        Psi_s_H=zeros(num_cell,cc_max);
        Psi_s_L=zeros(num_cell,cc_max);
        gsr_H=zeros(num_cell,cc_max);
        Psi_x_H=Psi_x_Htm1;
        Psi_l_H=Psi_l_Htm1;
        Jsx_H=zeros(num_cell,cc_max);
        Jxl_H=zeros(num_cell,cc_max);
        Kleaf_H=zeros(num_cell,cc_max);
        Kx_H=zeros(num_cell,cc_max);
        Vx_H=Vx_Htm1;
        Vl_H=Vl_Htm1;
        gsr_L=zeros(num_cell,cc_max);
        Psi_x_L=Psi_x_Ltm1;
        Psi_l_L=Psi_l_Ltm1;
        Jsx_L=zeros(num_cell,cc_max);
        Jxl_L=zeros(num_cell,cc_max);
        Kleaf_L=zeros(num_cell,cc_max);
        Kx_L=zeros(num_cell,cc_max);
        Vx_L=Vx_Ltm1;
        Vl_L=Vl_Ltm1;
        fapar_H=zeros(num_cell,cc_max);
        fapar_L=zeros(num_cell,cc_max);
        SIF_H=zeros(num_cell,cc_max);
        SIF_L=zeros(num_cell,cc_max);
        NDVI=zeros(num_cell,1);
        %CK2=zeros(num_cell,1);
        Ws_under=zeros(num_cell,1);
        HV=zeros(num_cell,1);
        QEV=zeros(num_cell,1);
        CK1=zeros(num_cell,1);
        %%%%%%%%%%%%%%%%%% BG added %%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % if OPT_SoilBiogeochemistry == 1
            % Ptm1=zeros(num_cell,55);
            Lk_t=zeros(num_cell,24);
            T_H_t=zeros(num_cell,cc_max,24);
            T_L_t=zeros(num_cell,cc_max,24);
            NavlItm1 = zeros(num_cell,3);
            dB_flux=zeros(num_cell,55,8);
            Bam=Ptm1(:,20);
            Bem=Ptm1(:,21);
            SupN_Htm1=zeros(num_cell,cc_max);
            SupP_Htm1=zeros(num_cell,cc_max);
            SupK_Htm1=zeros(num_cell,cc_max);
            SupN_Ltm1=zeros(num_cell,cc_max);
            SupP_Ltm1=zeros(num_cell,cc_max);
            SupK_Ltm1=zeros(num_cell,cc_max);
            RexmyItm1=zeros(num_cell,3);
            if ~isempty(InitialNutrientResetBlacklist)
                if size(InitialNutrientResetBlacklist,2) ~= 3
                    error('Param.InitialNutrientResetBlacklist must have columns [cell_linear, vegetation_layer, pft_index].')
                end
                InitialNutrientResetRows=zeros(0,15);
                for resetInitialIndex=1:size(InitialNutrientResetBlacklist,1)
                    resetCell=InitialNutrientResetBlacklist(resetInitialIndex,1);
                    resetLayer=InitialNutrientResetBlacklist(resetInitialIndex,2);
                    resetPft=InitialNutrientResetBlacklist(resetInitialIndex,3);
                    if resetCell < 1 || resetCell > num_cell || resetPft < 1 || resetPft > cc_max
                        warning('Skipping invalid InitialNutrientResetBlacklist row: cell=%d, layer=%d, pft=%d.',resetCell,resetLayer,resetPft)
                        continue
                    end
                    if resetLayer == 1
                        B_Htm1(resetCell,resetPft,:)=0;
                        LAI_Htm1(resetCell,resetPft)=0; PHE_S_Htm1(resetCell,resetPft)=0;
                        dflo_Htm1(resetCell,resetPft)=0; AgeL_Htm1(resetCell,resetPft)=0; AgeDL_Htm1(resetCell,resetPft)=0;
                        NPP_Htm1(resetCell,resetPft)=0; NPPI_Htm1(resetCell,resetPft)=0; TdpI_Htm1(resetCell,resetPft)=0;
                        Bfac_weekHtm1(resetCell,resetPft)=0; Nreserve_Htm1(resetCell,resetPft)=0;
                        Preserve_Htm1(resetCell,resetPft)=0; Kreserve_Htm1(resetCell,resetPft)=0; FNC_Htm1(resetCell,resetPft)=1;
                        SAI_Htm1(resetCell,resetPft)=0; hc_Htm1(resetCell,resetPft)=0; Tden_Htm1(resetCell,resetPft)=0; AgePl_Htm1(resetCell,resetPft)=0;
                        NupI_Htm1(resetCell,resetPft,:)=0; NuLit_Htm1(resetCell,resetPft,:)=0;
                        NBLeaf_Htm1(resetCell,resetPft)=0; PARI_Htm1(resetCell,resetPft,:)=0; NBLI_Htm1(resetCell,resetPft)=0;
                        SupN_Htm1(resetCell,resetPft)=1; SupP_Htm1(resetCell,resetPft)=1; SupK_Htm1(resetCell,resetPft)=1;
                        B_H(resetCell,resetPft,:)=0; LAI_H(resetCell,resetPft)=0; PHE_S_H(resetCell,resetPft)=0;
                        dflo_H(resetCell,resetPft)=0; AgeL_H(resetCell,resetPft)=0; AgeDL_H(resetCell,resetPft)=0;
                        NPP_H(resetCell,resetPft)=0; ANPP_H(resetCell,resetPft)=0; Rg_H(resetCell,resetPft)=0;
                        RA_H(resetCell,resetPft)=0; Rms_H(resetCell,resetPft)=0; Rmr_H(resetCell,resetPft)=0; Rmc_H(resetCell,resetPft)=0;
                        e_rel_H(resetCell,resetPft)=0; e_relN_H(resetCell,resetPft)=0; LAIdead_H(resetCell,resetPft)=0;
                        Sr_H(resetCell,resetPft)=0; Slf_H(resetCell,resetPft)=0; Sfr_H(resetCell,resetPft)=0; Swm_H(resetCell,resetPft)=0; Sll_H(resetCell,resetPft)=0;
                        Rexmy_H(resetCell,resetPft,:)=0; Rrootl_H(resetCell,resetPft)=0; Bfac_dayH(resetCell,resetPft)=0; Bfac_weekH(resetCell,resetPft)=0;
                        NPPI_H(resetCell,resetPft)=0; TdpI_H(resetCell,resetPft)=0; RB_H(resetCell,resetPft,:)=0;
                        Nreserve_H(resetCell,resetPft)=0; Preserve_H(resetCell,resetPft)=0; Kreserve_H(resetCell,resetPft)=0; FNC_H(resetCell,resetPft)=1;
                        NupI_H(resetCell,resetPft,:)=0; NuLit_H(resetCell,resetPft,:)=0; NBLeaf_H(resetCell,resetPft)=0; PARI_H(resetCell,resetPft,:)=0; NBLI_H(resetCell,resetPft)=0;
                        SAI_H(resetCell,resetPft)=0; hc_H(resetCell,resetPft)=0; BA_H(resetCell,resetPft)=0; Tden_H(resetCell,resetPft)=0; AgePl_H(resetCell,resetPft)=0;
                        rNc_H(resetCell,resetPft)=0; rPc_H(resetCell,resetPft)=0; rKc_H(resetCell,resetPft)=0; ManIH(resetCell,resetPft)=0;
                        TexC_H(resetCell,resetPft)=0; TexN_H(resetCell,resetPft)=0; TexP_H(resetCell,resetPft)=0; TexK_H(resetCell,resetPft)=0;
                        TNIT_H(resetCell,resetPft)=0; TPHO_H(resetCell,resetPft)=0; TPOT_H(resetCell,resetPft)=0;
                        SupN_H(resetCell,resetPft)=1; SupP_H(resetCell,resetPft)=1; SupK_H(resetCell,resetPft)=1; ISOIL_H(resetCell,resetPft,:)=0;
                        Ccrown_t(resetCell,resetPft)=0; Ccrown_t_tm1(resetCell,resetPft)=0;
                        NutrientResetEver_H(resetCell,resetPft)=true;
                    elseif resetLayer == 2
                        B_Ltm1(resetCell,resetPft,:)=0;
                        LAI_Ltm1(resetCell,resetPft)=0; PHE_S_Ltm1(resetCell,resetPft)=0;
                        dflo_Ltm1(resetCell,resetPft)=0; AgeL_Ltm1(resetCell,resetPft)=0; AgeDL_Ltm1(resetCell,resetPft)=0;
                        NPP_Ltm1(resetCell,resetPft)=0; NPPI_Ltm1(resetCell,resetPft)=0; TdpI_Ltm1(resetCell,resetPft)=0;
                        Bfac_weekLtm1(resetCell,resetPft)=0; Nreserve_Ltm1(resetCell,resetPft)=0;
                        Preserve_Ltm1(resetCell,resetPft)=0; Kreserve_Ltm1(resetCell,resetPft)=0; FNC_Ltm1(resetCell,resetPft)=1;
                        SAI_Ltm1(resetCell,resetPft)=0; hc_Ltm1(resetCell,resetPft)=0; Tden_Ltm1(resetCell,resetPft)=0; AgePl_Ltm1(resetCell,resetPft)=0;
                        NupI_Ltm1(resetCell,resetPft,:)=0; NuLit_Ltm1(resetCell,resetPft,:)=0;
                        NBLeaf_Ltm1(resetCell,resetPft)=0; PARI_Ltm1(resetCell,resetPft,:)=0; NBLI_Ltm1(resetCell,resetPft)=0;
                        SupN_Ltm1(resetCell,resetPft)=1; SupP_Ltm1(resetCell,resetPft)=1; SupK_Ltm1(resetCell,resetPft)=1;
                        B_L(resetCell,resetPft,:)=0; LAI_L(resetCell,resetPft)=0; PHE_S_L(resetCell,resetPft)=0;
                        dflo_L(resetCell,resetPft)=0; AgeL_L(resetCell,resetPft)=0; AgeDL_L(resetCell,resetPft)=0;
                        NPP_L(resetCell,resetPft)=0; ANPP_L(resetCell,resetPft)=0; Rg_L(resetCell,resetPft)=0;
                        RA_L(resetCell,resetPft)=0; Rms_L(resetCell,resetPft)=0; Rmr_L(resetCell,resetPft)=0; Rmc_L(resetCell,resetPft)=0;
                        e_rel_L(resetCell,resetPft)=0; e_relN_L(resetCell,resetPft)=0; LAIdead_L(resetCell,resetPft)=0;
                        Sr_L(resetCell,resetPft)=0; Slf_L(resetCell,resetPft)=0; Sfr_L(resetCell,resetPft)=0; Swm_L(resetCell,resetPft)=0; Sll_L(resetCell,resetPft)=0;
                        Rexmy_L(resetCell,resetPft,:)=0; Rrootl_L(resetCell,resetPft)=0; Bfac_dayL(resetCell,resetPft)=0; Bfac_weekL(resetCell,resetPft)=0;
                        NPPI_L(resetCell,resetPft)=0; TdpI_L(resetCell,resetPft)=0; RB_L(resetCell,resetPft,:)=0;
                        Nreserve_L(resetCell,resetPft)=0; Preserve_L(resetCell,resetPft)=0; Kreserve_L(resetCell,resetPft)=0; FNC_L(resetCell,resetPft)=1;
                        NupI_L(resetCell,resetPft,:)=0; NuLit_L(resetCell,resetPft,:)=0; NBLeaf_L(resetCell,resetPft)=0; PARI_L(resetCell,resetPft,:)=0; NBLI_L(resetCell,resetPft)=0;
                        SAI_L(resetCell,resetPft)=0; hc_L(resetCell,resetPft)=0; BA_L(resetCell,resetPft)=0; Tden_L(resetCell,resetPft)=0; AgePl_L(resetCell,resetPft)=0;
                        rNc_L(resetCell,resetPft)=0; rPc_L(resetCell,resetPft)=0; rKc_L(resetCell,resetPft)=0; ManIL(resetCell,resetPft)=0;
                        TexC_L(resetCell,resetPft)=0; TexN_L(resetCell,resetPft)=0; TexP_L(resetCell,resetPft)=0; TexK_L(resetCell,resetPft)=0;
                        TNIT_L(resetCell,resetPft)=0; TPHO_L(resetCell,resetPft)=0; TPOT_L(resetCell,resetPft)=0;
                        SupN_L(resetCell,resetPft)=1; SupP_L(resetCell,resetPft)=1; SupK_L(resetCell,resetPft)=1; ISOIL_L(resetCell,resetPft,:)=0;
                        Ccrown_t(resetCell,resetPft)=0; Ccrown_t_tm1(resetCell,resetPft)=0;
                        NutrientResetEver_L(resetCell,resetPft)=true;
                    else
                        warning('Skipping invalid InitialNutrientResetBlacklist layer: cell=%d, layer=%d, pft=%d.',resetCell,resetLayer,resetPft)
                        continue
                    end
                    [resetRow,resetCol]=ind2sub([m_cell,n_cell],resetCell);
                    InitialNutrientResetRows=[InitialNutrientResetRows;...
                        1,Date(1),Datam(1,1:4),resetCell,resetRow,resetCol,resetLayer,resetPft,9,nan,nan,nan];
                end
                if ~isempty(InitialNutrientResetRows)
                    NutrientResetLog=[NutrientResetLog;InitialNutrientResetRows];
                    NutrientResetUniqueCellCount=nnz(any(NutrientResetEver_H | NutrientResetEver_L,2));
                    NutrientResetUniquePftCount=nnz(NutrientResetEver_H)+nnz(NutrientResetEver_L);
                    save(NutrientResetLogFile,'NutrientResetLog','NutrientResetLogColumns',...
                        'NutrientResetLayerCodes','NutrientResetReasonCodes','NutrientResetEver_H','NutrientResetEver_L',...
                        'NutrientResetUniqueCellCount','NutrientResetUniquePftCount');
                    fprintf('Initial nutrient blacklist reset: %d PFT(s) cleared before Iter 2.\n',size(InitialNutrientResetRows,1));
                end
            end
        % end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        aR_ = zeros(num_cell,1);
        Zbio_ = zeros(num_cell,1);
        rsd_ = zeros(num_cell, ms_max);

        Zs_= zeros(num_cell, ms_max+1);       
        EvL_Zs_ = zeros(num_cell, ms_max);
        Inf_Zs_= zeros(num_cell, ms_max);
        Bio_Zs_= zeros(num_cell, ms_max);
        Zinf_= zeros(num_cell,1);
        RfH_Zs_= zeros(num_cell, cc_max, ms_max);
        RfL_Zs_= zeros(num_cell, cc_max, ms_max);
        dz_= zeros(num_cell, ms_max);
        Dz_= zeros(num_cell, ms_max); 
        ms_ = zeros(num_cell,1);
        Kbot_ = zeros(num_cell,1);
        Krock_ = zeros(num_cell,1);
        zatm_ = zeros(num_cell,1);
        Ccrown_= zeros(num_cell, cc_max);
        Cbare_t = zeros(num_cell,1);
        Crock_t = zeros(num_cell,1);
        Curb_t = zeros(num_cell,1);
        Cwat_t = zeros(num_cell,1);
        Color_Class_ = zeros(num_cell,1);
        OM_H_ = zeros(num_cell,1);
        OM_L_ = zeros(num_cell,1);
        d_leaf_H_= zeros(num_cell, cc_max);
        d_leaf_L_= zeros(num_cell, cc_max);
        SPAR_ = zeros(num_cell,1);
        Phy_ = zeros(num_cell,1);
        fpr_ = zeros(num_cell,1);
        aSE_H_ = zeros(num_cell, cc_max); 
        aSE_L_= zeros(num_cell, cc_max);
        fab_H_= zeros(num_cell, cc_max);
        fbe_H_= zeros(num_cell, cc_max);
        fab_L_= zeros(num_cell, cc_max);
        fbe_L_= zeros(num_cell, cc_max);
        ZR95_H_= zeros(num_cell, cc_max);
        ZR95_L_= zeros(num_cell, cc_max);
        In_max_urb_= zeros(num_cell,1);
        In_max_rock_ = zeros(num_cell,1);
        K_usle_ = zeros(num_cell,1);
        % Zs_deb_= zeros(num_cell, cc_max, ms_max);
        
        Sllit_ = zeros(num_cell,1);
        Kct_ = zeros(num_cell,1);
        ExEM_= zeros(num_cell, cc_max);
        a_dis_ = zeros(num_cell,1);
        pow_dis_ = zeros(num_cell,1);
        % a = cell(17,num_cell);
        PFT_opt_H_cell = cell(num_cell,1); 
        PFT_opt_L_cell = cell(num_cell,1); 
        Soil_Param_cell = cell(num_cell,1); 
        Interc_Param_cell = cell(num_cell,1); 
        SnowIce_Param_cell = cell(num_cell,1); 
        VegH_Param_cell = cell(num_cell,1); 
        VegL_Param_cell = cell(num_cell,1); 
        Stoich_H_cell = cell(num_cell,1); 
        Stoich_L_cell = cell(num_cell,1); 
        VegH_Param_Dyn_cell = cell(num_cell,1); 
        VegL_Param_Dyn_cell = cell(num_cell,1); 
        ParEx_H_cell = cell(num_cell,1); 
        Mpar_H_cell = cell(num_cell,1); 
        ParEx_L_cell = cell(num_cell,1); 
        Mpar_L_cell = cell(num_cell,1); 
        Deb_Par_cell = cell(num_cell,1); 
        Zs_deb_cell = cell(num_cell,1); 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    Cwat_all = Geo.Cwat;
    Curb_all = Geo.Curb;
    Crock_all = Geo.Crock;
    Cbare_all = Geo.Cbare;
    PFT_frac_all = Geo.PFT_def;
    % Cwat_all = Geo.Rwat;
    % Curb_all = Geo.Rurb;
    % Crock_all = Geo.Rrock;
    % Cbare_all = Geo.Rbare;
    % PFT_frac_all = Geo.RPFT;

    PFT_frac_all = reshape(PFT_frac_all, [size(PFT_frac_all,1)*size(PFT_frac_all,2) size(PFT_frac_all,3)]);
    PF = Param.MOD_PARAM;

    %tic
     if t ==2
        if isfile(Param.param_sname)
            load(Param.param_sname);
         else
            fprintf('Parameter file not found. Computing parameters...\n');
            pft_info_opts = detectImportOptions("agg_try30_pred6.csv");
            pft_info = readtable("agg_try30_pred6.csv"); 
            pft_info  = pft_info(~isnan(pft_info.sla), :);
    
            site_info_opts = detectImportOptions("Table_Sites_.xlsx");
            site_info_opts.DataRange = "A3:BS253";
            site_info_opts.VariableTypes = ["double","char","char","double","double","double","double","double","char","double","double"...
                "double",repmat("double",1,59)];
            site_info = readtable("Table_Sites_.xlsx",site_info_opts);
            % Fill missing values in numeric columns with NaN
            site_info(:, vartype('numeric')) = fillmissing(site_info(:, vartype('numeric')), 'constant', NaN);
            % Fill missing values in string or character columns with an empty string
            site_info(:, vartype('string')) = fillmissing(site_info(:, vartype('string')), 'constant', "");
            par_names = site_info.Properties.VariableNames(1,10:end);
            %par_con = table2cell(readtable("Table_Sites_AP.xlsx",Range="J1:BS2"));
            par_con = table2cell(readtable("Table_Sites_modified.xlsx",Range="K1:BT2"));
            % Prmean = nansum(squeeze(ncread( [Meteo.FolderLoc, filesep 'Pr_cr_era5_wc.nc'], 'Pr_cr_era5_wc', [1, 1, 1], [1, 1, Inf]) ))/4
            nyear=round(length(squeeze(ncread( [Meteo.FolderLoc, filesep 'Pr_' time_index '.nc'], 'Pr', [1, 1, 1], [1, 1, Inf]) ))/8760);
            Prmean = nansum(squeeze(ncread( [Meteo.FolderLoc, filesep 'Pr_' time_index '.nc'], 'Pr', [1, 1, 1], [1, 1, Inf]) ))/nyear; %% divide by year.
            Tamean = nanmean(squeeze(ncread( [Meteo.FolderLoc, filesep 'Ta_' time_index '.nc'], 'Ta',  [1, 1, 1], [Inf, Inf, Inf])),'all');
            Tamin = nanmin(squeeze(ncread( [Meteo.FolderLoc, filesep 'Ta_' time_index '.nc'], 'Ta',  [1, 1, 1], [Inf, Inf, Inf])),[],'all');
    
            for ij=1:num_cell
                if MASKn(ij)
                    Cwat_ = Cwat_all(ij);
                    Curb_ = Curb_all(ij);
                    Crock_ = Crock_all(ij);
                    Cbare_ = Cbare_all(ij);
                    PFT_frac_ = PFT_frac_all(ij, :);
    
                    [Zbio_(ij),rsd_(ij,:),aR_(ij),Zs_(ij,:),...
                        EvL_Zs_(ij,:),Inf_Zs_(ij,:),Bio_Zs_(ij,:),Zinf_(ij),RfH_Zs_(ij,:,:),RfL_Zs_(ij,:,:),dz_(ij,:),~,Dz_(ij,:),...
                        ms_(ij),Kbot_(ij),Krock_(ij),zatm_(ij),...
                        Ccrown_(ij,:),Cbare_t(ij),Crock_t(ij),Curb_t(ij),Cwat_t(ij),...
                        Color_Class_(ij),OM_H_(ij),OM_L_(ij),PFT_opt_H_,PFT_opt_L_,d_leaf_H_(ij,:),d_leaf_L_(ij,:),...
                        SPAR_(ij),Phy_(ij),Soil_Param_,Interc_Param_,SnowIce_Param_,VegH_Param_,VegL_Param_,fpr_(ij),...
                        VegH_Param_Dyn_,VegL_Param_Dyn_,...
                        Stoich_H_,aSE_H_(ij,:),Stoich_L_,aSE_L_(ij,:),fab_H_(ij,:),fbe_H_(ij,:),fab_L_(ij,:),fbe_L_(ij,:),...
                        ZR95_H_(ij,:),ZR95_L_(ij,:),In_max_urb_(ij),In_max_rock_(ij),K_usle_(ij),...
                        Deb_Par_,Zs_deb_,...
                        Sllit_(ij),Kct_(ij),ExEM_(ij,:),ParEx_H_,Mpar_H_,ParEx_L_,Mpar_L_,a_dis_(ij),pow_dis_(ij)]=...
                        feval(PF,ij,ZATM(ij),PSAN(ij),PCLA(ij),PORG(ij),PFT_frac_, Cwat_,Curb_,Crock_,Cbare_,...
                        KBOT(ij),KROCK(ij),md_max,DTM,PFTdef,eco_idx,leaf_onset_T,leaf_sncs_T,leaf_sncs_day,leaf_onset_day,...
                        lat,IIpar,site_info,par_names,par_con,pft_info,Prmean, Tamean,Tamin,Zsmax, EOS_maize_main,EOS_sc,EOS_tc ,EOS_wc,SOS_maize_main,SOS_sc,SOS_tc,SOS_wc);
                    
                    %%%%%%%%%%calibration added%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%% 01054200 %%%%%%%%%%%%%%%%%%%%%%%%
                    VegH_Param_Dyn_.Tlo(3)=3;
                    VegH_Param_Dyn_.Tlo(1)=5;  
                    VegH_Param_Dyn_.Tlo(4)=3;
                    VegH_Param_Dyn_.Tcold(3)=6;
                    VegH_Param_Dyn_.age_cr(2)=170;
                    VegH_Param_Dyn_.age_cr(3)=170;
                    VegH_Param_.Vmax_H(1)=40;
                    VegH_Param_.Vmax_H(2)=62;
                    VegH_Param_.Vmax_H(4)=72;
                    VegH_Param_Dyn_.Sl(3)=0.026;
                    VegH_Param_Dyn_.Sl(1)=0.01;
                    VegH_Param_Dyn_.Sl(2)=0.025;
                    VegH_Param_.KnitH(3)=0.21;
                    VegH_Param_Dyn_.LDay_cr(1)=10.6;
                    %%%% newly add
                    VegH_Param_Dyn_.PsiG50(5)=-1;                 
                    VegH_Param_Dyn_.PsiG99(5)=-7.5;
                    VegH_Param_.Vmax_H(5)=55;
                    VegH_Param_.Psi_sto_50_H(5)=-7.5;
                    VegH_Param_.PsiL50_H(5)=-8;
                    VegH_Param_.PsiL00_H(5)=-2;
                    VegH_Param_Dyn_.PsiG99(5)=-7.5;
                    VegH_Param_.PsiX50_H(5)=-10;
                    %%%%%%%%%%%%%%%%%%%%% 01466500 %%%%%%%%%%%%%%%%%%%%%%%%
                    % VegH_Param_.Vmax_H(3)=58;
                    % VegH_Param_Dyn_.Sl(1)=0.008;
                    % VegH_Param_Dyn_.Sl(2)=0.034;
                    % VegH_Param_Dyn_.Sl(3)=0.02;
                    % VegH_Param_.KnitH(1)=0.38;
                    % VegH_Param_.KnitH(2)=0.16;
                    % VegH_Param_.KnitH(3)=0.15;
                    %%%%%%%%%%%%%%%%%%%%% 07362587 %%%%%%%%%%%%%%%%%%%%%%%%
                    % VegH_Param_Dyn_.Sl(1)=0.009;
                    % VegH_Param_Dyn_.Sl(2)=0.016;
                    % VegH_Param_Dyn_.Sl(3)=0.022;
                    % VegH_Param_.Vmax_H(1)=45;
                    % VegH_Param_.Vmax_H(2)=65;
                    % VegH_Param_.Vmax_H(3)=55;
                    % VegH_Param_.Vmax_H(4)=63;
                    % VegH_Param_Dyn_.Tcold(3)=8;
                    % VegH_Param_.KnitH(1)= 0.6; 
                    % VegH_Param_.KnitH(2)= 0.12;
                    % VegH_Param_.KnitH(3)=0.55;
                    % VegH_Param_.KnitH(4)=0.2;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    PFT_opt_H_cell{ij} = PFT_opt_H_;
                    PFT_opt_L_cell{ij} = PFT_opt_L_;
                    Soil_Param_cell{ij} = Soil_Param_;
                    Interc_Param_cell{ij} = Interc_Param_;
                    SnowIce_Param_cell{ij} = SnowIce_Param_;
                    VegH_Param_cell{ij} = VegH_Param_;
                    VegL_Param_cell{ij} = VegL_Param_;
                    Stoich_H_cell{ij} = Stoich_H_;
                    Stoich_L_cell{ij} = Stoich_L_;
                    VegH_Param_Dyn_cell{ij} = VegH_Param_Dyn_;
                    VegL_Param_Dyn_cell{ij} = VegL_Param_Dyn_;
                    ParEx_H_cell{ij} = ParEx_H_;
                    Mpar_H_cell{ij} = Mpar_H_;
                    ParEx_L_cell{ij} = ParEx_L_;
                    Mpar_L_cell{ij} = Mpar_L_;
                    Deb_Par_cell{ij} = Deb_Par_;
                    Zs_deb_cell{ij} = Zs_deb_;
    
                end 
            end
            save(Param.param_sname,'Zbio_','rsd_','aR_','Zs_',...
                        'EvL_Zs_','Inf_Zs_','Bio_Zs_','Zinf_','RfH_Zs_','RfL_Zs_','dz_','Dz_',...
                        'ms_','Kbot_','Krock_','zatm_',...
                        'Ccrown_','Cbare_t','Crock_t','Curb_t','Cwat_t',...
                        'Color_Class_','OM_H_','OM_L_','PFT_opt_H_cell','PFT_opt_L_cell','d_leaf_H_','d_leaf_L_',...
                        'SPAR_','Phy_','Soil_Param_cell','Interc_Param_cell','SnowIce_Param_cell','VegH_Param_cell','VegL_Param_cell','fpr_',...
                        'VegH_Param_Dyn_cell','VegL_Param_Dyn_cell',...
                        'Stoich_H_cell','aSE_H_','Stoich_L_cell','aSE_L_','fab_H_','fbe_H_','fab_L_','fbe_L_',...
                        'ZR95_H_','ZR95_L_','In_max_urb_','In_max_rock_','K_usle_',...
                        'Deb_Par_cell','Zs_deb_cell',...
                        'Sllit_','Kct_','ExEM_','ParEx_H_cell','Mpar_H_cell','ParEx_L_cell','Mpar_L_cell','a_dis_','pow_dis_')
         end
    end 
    % load('matlab.mat'); 
    %toc 

    NutrientReset_H=false(num_cell,cc_max);
    NutrientReset_L=false(num_cell,cc_max);
    NutrientResetInfo_H=nan(num_cell,cc_max,4);
    NutrientResetInfo_L=nan(num_cell,cc_max,4);
    parfor ij=1:num_cell
     % for ij=1:num_cell
        if MASKn(ij)
            %[il, jl] = ind2sub([m_cell,n_cell],ij);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% BOUNDARY CONDITION  %%% INTRODUCED SOIL AND VEG. for ij

            % Cwat = PFTdef.Cwat;
            % Curb = PFTdef.Curb;
            % Crock = PFTdef.Crock;
            % Cbare = PFTdef.Cbare;
            % Ccrown = PFTdef.PFT_frac;

            %PFTdef_c =struct;

            % 
            % Cwat_ = Cwat_all(ij);
            % Curb_ = Curb_all(ij);
            % Crock_ = Crock_all(ij);
            % Cbare_ = Cbare_all(ij);
            % PFT_frac_ = PFT_frac_all(ij, :);
            % 
            % if t ==2
            %     [aR_(ij),Zs_(ij,:),...
            %         EvL_Zs_(ij,:),Inf_Zs_(ij,:),Bio_Zs_(ij,:),Zinf_(ij),RfH_Zs_(ij,:,:),RfL_Zs_(ij,:,:),dz_(ij,:),~,Dz_(ij,:),...
            %         ms_(ij),Kbot_(ij),Krock_(ij),zatm_(ij),...
            %         Ccrown_(ij,:),Cbare_t(ij),Crock_t(ij),Curb_t(ij),Cwat_t(ij),...
            %         Color_Class_(ij),OM_H_(ij),OM_L_(ij),PFT_opt_H_,PFT_opt_L_,d_leaf_H_(ij,:),d_leaf_L_(ij,:),...
            %         SPAR_(ij),Phy_(ij),Soil_Param_,Interc_Param_,SnowIce_Param_,VegH_Param_,VegL_Param_,fpr_(ij),...
            %         VegH_Param_Dyn_,VegL_Param_Dyn_,...
            %         Stoich_H_,aSE_H_(ij,:),Stoich_L_,aSE_L_(ij,:),fab_H_(ij,:),fbe_H_(ij,:),fab_L_(ij,:),fbe_L_(ij,:),...
            %         ZR95_H_(ij,:),ZR95_L_(ij,:),In_max_urb_(ij),In_max_rock_(ij),K_usle_(ij),...
            %         Deb_Par_,Zs_deb_(ij,:,:),...
            %         Sllit_(ij),Kct_(ij),ExEM_(ij,:),ParEx_H_,Mpar_H_,ParEx_L_,Mpar_L_,a_dis_(ij),pow_dis_(ij)]=...
            %         feval(PF,ij,ZATM(ij),PSAN(ij),PCLA(ij),PORG(ij),PFT_frac_, Cwat_,Curb_,Crock_,Cbare_,...
            %         KBOT(ij),KROCK(ij),md_max,DTM,PFTdef,eco_idx,leaf_onset_T,leaf_sncs_T,leaf_sncs_day,leaf_onset_day,lat);
            % 
            % 
            %         a{1,ij} = PFT_opt_H_;
            %         a{2,ij} = PFT_opt_L_;
            %         a{3,ij} = Soil_Param_;
            %         a{4,ij} = Interc_Param_;
            %         a{5,ij} = SnowIce_Param_;
            %         a{6,ij} = VegH_Param_;
            %         a{7,ij} = VegL_Param_;
            %         a{8,ij} = Stoich_H_; 
            %         a{9,ij} = Stoich_L_;
            %         a{10,ij} = VegH_Param_Dyn_;
            %         a{11,ij} = VegL_Param_Dyn_;
            %         a{12,ij} = ParEx_H_;
            %         a{13,ij} = Mpar_H_;
            %         a{14,ij} = ParEx_L_;
            %         a{15,ij} = Mpar_L_;
            %         a{16,ij} = Deb_Par_;
            % end 

            Zbio = Zbio_(ij);
            rsd = rsd_(ij,:);
            aR = aR_(ij);
            Zs =  Zs_(ij,:);
            EvL_Zs = EvL_Zs_(ij,:);
            Inf_Zs=Inf_Zs_(ij,:);
            Bio_Zs=Bio_Zs_(ij,:);
            Zinf=Zinf_(ij);
            RfH_Zs = squeeze(RfH_Zs_(ij,:,:));
            RfL_Zs = squeeze(RfL_Zs_(ij,:,:));
            dz = dz_(ij,:);
            Dz = Dz_(ij,:);
            ms = ms_(ij);
            Kbot = Kbot_(ij);
            Krock = Krock_(ij);
            zatm = zatm_(ij);
            Ccrown = Ccrown_(ij,:);
            Cbare = Cbare_t(ij);
            Crock = Crock_t(ij);
            Curb = Curb_t(ij);
            Cwat = Cwat_t(ij);
            Color_Class = Color_Class_(ij);
            OM_H = OM_H_(ij);
            OM_L = OM_L_(ij);
            PFT_opt_H = PFT_opt_H_cell{ij};
            PFT_opt_L = PFT_opt_L_cell{ij};
            d_leaf_H = d_leaf_H_(ij,:);
            d_leaf_L = d_leaf_L_(ij,:);
            SPAR = SPAR_(ij);
            Phy = Phy_(ij);
            Soil_Param = Soil_Param_cell{ij};
            Interc_Param =  Interc_Param_cell{ij};
            SnowIce_Param = SnowIce_Param_cell{ij};
            VegH_Param = VegH_Param_cell{ij};
            VegL_Param = VegL_Param_cell{ij};
            fpr = fpr_(ij);
            VegH_Param_Dyn =  VegH_Param_Dyn_cell{ij};
            VegL_Param_Dyn =  VegL_Param_Dyn_cell{ij};
            Stoich_H = Stoich_H_cell{ij};
            aSE_H = aSE_H_(ij,:);
            Stoich_L = Stoich_L_cell{ij};
            aSE_L = aSE_L_(ij,:);
            fab_H = fab_H_(ij,:);
            fbe_H = fbe_H_(ij,:);
            fab_L = fab_L_(ij,:);
            fbe_L = fbe_L_(ij,:);
            ZR95_H = ZR95_H_(ij,:);
            ZR95_L = ZR95_L_(ij,:);
            In_max_urb = In_max_urb_(ij);
            In_max_rock = In_max_rock_(ij);
            K_usle = K_usle_(ij);
            Deb_Par = Deb_Par_cell{ij};
            Sllit = Sllit_(ij);
            Kct = Kct_(ij);
            ExEM = ExEM_(ij,:);
            if length(unique(ExEM))>1
                print('Warning:ExEM not always 0 or 1')
            else
                ExEM=mode(ExEM);
            end
            ParEx_H =ParEx_H_cell{ij};
            Mpar_H = Mpar_H_cell{ij};
            ParEx_L =ParEx_L_cell{ij};
            Mpar_L =Mpar_L_cell{ij};
            a_dis = a_dis_(ij);
            pow_dis = pow_dis_(ij);
            Zs_deb = Zs_deb_cell{ij};
            
            Urb_Par = struct;
            Urb_Par.alb= UALB(ij);
            Urb_Par.e_sur= Ue_sur(ij);
            Urb_Par.zom= BuildH(ij)*0.123;
            if Curb==1 && BuildH(ij)==0
                Urb_Par.zom = 0.2;
            end

            if (HRL==1)

                % [hc_H(ij,:),SAI_H(ij,:),B_Htm1(ij,:,:),Nreserve_Htm1(ij,:),Preserve_Htm1(ij,:),Kreserve_Htm1(ij,:),...
                %     hc_L(ij,:),SAI_L(ij,:),B_Ltm1(ij,:,:),Nreserve_Ltm1(ij,:),Preserve_Ltm1(ij,:),Kreserve_Ltm1(ij,:)] = ...
                %     ...
                %     CropHeightPAR...
                %     ...
                %     (Ccrown, hc_Htm1(ij,:),SAI_Htm1(ij,:),aSE_H,LAI_H(ij,:),LAIdead_H(ij,:),B_Htm1(ij,:,:),Nreserve_Htm1(ij,:),Preserve_Htm1(ij,:),Kreserve_Htm1(ij,:),ManIH(ij,:),Mpar_H,VegH_Param_Dyn,...
                %     hc_Ltm1(ij,:),SAI_Ltm1(ij,:),aSE_L,LAI_L(ij,:),LAIdead_L(ij,:),B_Ltm1(ij,:,:),Nreserve_Ltm1(ij,:),Preserve_Ltm1(ij,:),Kreserve_Ltm1(ij,:),ManIL(ij,:),Mpar_L,VegL_Param_Dyn,OPT_SoilBiogeochemistry);
                % 

                %%%%%%%%%%%%%
                %%%%% SOIL BIOGEOCHEMISTRY MODULE
                [Se_bio,Se_fc,Psi_bio,Tdp_bio,VSUM,VTSUM]=Biogeo_environment([squeeze(Tdp_t(ij,:,:))]',[squeeze(O_t(ij,:,:))]',[squeeze(V_t(ij,:,:))]',...
                    Soil_Param,Phy,SPAR,Bio_Zs);%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if OPT_SoilBiogeochemistry == 1
                    [B_IO]=Biogeochemistry_IO(Zs,mean(Lat(MASK),'all'),mean(Lon(MASK),'all'),Upl);
                    [dB_flux(ij,:,:), P(ij,:),BLit(ij,:),Bam(ij),Bem(ij),...
                        Nuptake_H(ij,:),Puptake_H(ij,:),Kuptake_H(ij,:),Nuptake_L(ij,:),Puptake_L(ij,:),Kuptake_L(ij,:),RexmyI(ij,:),...
                        R_litter(ij),R_microbe(ij),R_litter_sur(ij),R_ew(ij),VOL(ij),N2flx(ij),Min_N(ij),Min_P(ij),...
                        R_bacteria(ij),RmycAM(ij),RmycEM(ij),Prod_B(ij),Prod_F(ij),BfixN(ij),NavlI(ij,:),LitFirEmi(ij,:)]=Biogeo_MODULE_PAR(ISOIL_L(ij,:,:),ISOIL_H(ij,:,:),Rexmy_L(ij,:,:),Rexmy_H(ij,:,:),ManIH(ij,:),ManIL(ij,:),Ptm1(ij,:),Zbio,sum(Bio_Zs.*rsd),PHs,Tdp_bio,mean(Ta_t(ij,:)),Psi_bio,Se_bio,Se_fc,VSUM,VTSUM,...
                        Ccrown,Bio_Zs,RfH_Zs,RfL_Zs,sum(Lk_t(ij,:)),sum(T_H_t(ij,:,:),3),sum(T_L_t(ij,:,:),3),B_Htm1(ij,:,3),B_Ltm1(ij,:,3),LAI_H(ij,:),LAI_L(ij,:),... % should be LAI_H not LAI_Htm1
                        SupN_Htm1(ij,:),SupP_Htm1(ij,:),SupK_Htm1(ij,:),SupN_Ltm1(ij,:),SupP_Ltm1(ij,:),SupK_Ltm1(ij,:),...
                        RexmyItm1(ij,:),ExEM,NavlItm1(ij,:),PCLA(ij),PSAN(ij),...
                        B_IO,jDay,0,cc_max);
                else
                    % Keep prior local pools so only the later routing step can update them.
                    P(ij,:)=Ptm1(ij,:);
                    dB_flux(ij,:,:)=0;
                    BLit(ij,:)=0;
                    Nuptake_H(ij,:)= 0.0;
                    Puptake_H(ij,:)= 0.0;
                    Kuptake_H(ij,:)= 0.0; %% [gK/m^2 day]
                    %%%%
                    Nuptake_L(ij,:)= 0.0;  %% [gN/m^2 day]
                    Puptake_L(ij,:)= 0.0;
                    Kuptake_L(ij,:)= 0.0;
                    %%%%
                    NavlI(ij,:)=[1 1 1];
                    Bam(ij)=0; Bem(ij)=0;
                    RexmyI(ij,:)=0;
                    R_litter(ij)=0; R_microbe(ij)=0; R_litter_sur(ij)=0; R_ew(ij)=0;
                    VOL(ij)=0; N2flx(ij)=0; Min_N(ij)=0; Min_P(ij)=0;
                    R_bacteria(ij)=0; RmycAM(ij)=0; RmycEM(ij)=0;
                    Prod_B(ij)=0; Prod_F(ij)=0; BfixN(ij)=0; LitFirEmi(ij,:)=0;
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% Biogeochemistry Unit
                %[I,J] = ind2sub(size(DTM), ij);

                % L_day = L_day_f.L_day(ij, j);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                 [LAI_H(ij,:),B_H(ij,:,:),NPP_H(ij,:),ANPP_H(ij,:),Rg_H(ij,:),RA_H(ij,:),...
                    Rms_H(ij,:),Rmr_H(ij,:),Rmc_H(ij,:),PHE_S_H(ij,:),...
                    dflo_H(ij,:),AgeL_H(ij,:),e_rel_H(ij,:),e_relN_H(ij,:),...
                    LAI_L(ij,:),B_L(ij,:,:),NPP_L(ij,:),ANPP_L(ij,:),Rg_L(ij,:),RA_L(ij,:),...
                    Rms_L(ij,:),Rmr_L(ij,:),Rmc_L(ij,:),PHE_S_L(ij,:),...
                    dflo_L(ij,:),AgeL_L(ij,:),e_rel_L(ij,:),e_relN_L(ij,:),...
                    SAI_H(ij,:),hc_H(ij,:),SAI_L(ij,:),hc_L(ij,:),...
                    LAIdead_H(ij,:),NBLeaf_H(ij,:),Sr_H(ij,:),Slf_H(ij,:),Sfr_H(ij,:),Sll_H(ij,:),Swm_H(ij,:),Rexmy_H(ij,:,:),NupI_H(ij,:,:),NuLit_H(ij,:,:),...
                    LAIdead_L(ij,:),NBLeaf_L(ij,:),Sr_L(ij,:),Slf_L(ij,:),Sfr_L(ij,:),Sll_L(ij,:),Swm_L(ij,:),Rexmy_L(ij,:,:),NupI_L(ij,:,:),NuLit_L(ij,:,:),...
                    Rrootl_H(ij,:),AgeDL_H(ij,:),Bfac_dayH(ij,:),Bfac_weekH(ij,:),NPPI_H(ij,:),TdpI_H(ij,:),PARI_H(ij,:,:),NBLI_H(ij,:),RB_H(ij,:,:),FNC_H(ij,:),...
                    Nreserve_H(ij,:),Preserve_H(ij,:),Kreserve_H(ij,:),rNc_H(ij,:),rPc_H(ij,:),rKc_H(ij,:),ManIH(ij,:),...
                    Rrootl_L(ij,:),AgeDL_L(ij,:),Bfac_dayL(ij,:),Bfac_weekL(ij,:),NPPI_L(ij,:),TdpI_L(ij,:),PARI_L(ij,:,:),NBLI_L(ij,:),RB_L(ij,:,:),FNC_L(ij,:),...
                    Nreserve_L(ij,:),Preserve_L(ij,:),Kreserve_L(ij,:),rNc_L(ij,:),rPc_L(ij,:),rKc_L(ij,:),ManIL(ij,:),...
                    TexC_H(ij,:),TexN_H(ij,:),TexP_H(ij,:),TexK_H(ij,:),TNIT_H(ij,:),TPHO_H(ij,:),TPOT_H(ij,:),...
                    SupN_H(ij,:),SupP_H(ij,:),SupK_H(ij,:),ISOIL_H(ij,:,:),...
                    TexC_L(ij,:),TexN_L(ij,:),TexP_L(ij,:),TexK_L(ij,:),TNIT_L(ij,:),TPHO_L(ij,:),TPOT_L(ij,:),...
                    SupN_L(ij,:),SupP_L(ij,:),SupK_L(ij,:),ISOIL_L(ij,:,:),...
                    BA_H(ij,:),Tden_H(ij,:),AgePl_H(ij,:),BA_L(ij,:),Tden_L(ij,:),AgePl_L(ij,:),Ccrown_t(ij,:),...
                    NutrientReset_H(ij,:),NutrientReset_L(ij,:),NutrientResetInfo_H_cell,NutrientResetInfo_L_cell]=...
                    ...
                    VEGETATION_MODULE_PAR_Nutrient...
                    ...
                    (cc_max,Ccrown,ZR95_H,ZR95_L,B_Htm1(ij,:,:),...
                    PHE_S_Htm1(ij,:),dflo_Htm1(ij,:),AgeL_Htm1(ij,:),AgeDL_Htm1(ij,:),...
                    Ta_t(ij,:),PAR_t(ij,:),Tdp_H_t(ij,:,:),Psi_x_H_t(ij,:,:),Psi_l_H_t(ij,:,:),An_H_t(ij,:,:),Rdark_H_t(ij,:,:),NPP_Htm1(ij,:),jDay,Datam_S,...
                    NPPI_Htm1(ij,:),TdpI_Htm1(ij,:),Bfac_weekHtm1(ij,:),...
                    Stoich_H,aSE_H,VegH_Param_Dyn,...
                    Nreserve_Htm1(ij,:),Preserve_Htm1(ij,:),Kreserve_Htm1(ij,:),Nuptake_H(ij,:),Puptake_H(ij,:),Kuptake_H(ij,:),FNC_Htm1(ij,:),Tden_Htm1(ij,:),AgePl_Htm1(ij,:),...
                    fab_H,fbe_H,ParEx_H,Mpar_H,TBio_H(ij,:),SAI_Htm1(ij,:),hc_Htm1(ij,:),...
                    B_Ltm1(ij,:,:),PHE_S_Ltm1(ij,:),dflo_Ltm1(ij,:),AgeL_Ltm1(ij,:),AgeDL_Ltm1(ij,:),...
                    Tdp_L_t(ij,:,:),Psi_x_L_t(ij,:,:),Psi_l_L_t(ij,:,:),An_L_t(ij,:,:),Rdark_L_t(ij,:,:),NPP_Ltm1(ij,:),...
                    NPPI_Ltm1(ij,:),TdpI_Ltm1(ij,:),Bfac_weekLtm1(ij,:),...
                    NupI_Htm1(ij,:,:),NupI_Ltm1(ij,:,:),NuLit_Htm1(ij,:,:),NuLit_Ltm1(ij,:,:),NBLeaf_Htm1(ij,:),NBLeaf_Ltm1(ij,:),...
                    PARI_Htm1(ij,:,:),NBLI_Htm1(ij,:),PARI_Ltm1(ij,:,:),NBLI_Ltm1(ij,:),...
                    Stoich_L,aSE_L,VegL_Param_Dyn,...
                    NavlI(ij,:),Bam(ij),Bem(ij),Ccrown_t_tm1(ij,:),...
                    Nreserve_Ltm1(ij,:),Preserve_Ltm1(ij,:),Kreserve_Ltm1(ij,:),Nuptake_L(ij,:),Puptake_L(ij,:),Kuptake_L(ij,:),FNC_Ltm1(ij,:),Tden_Ltm1(ij,:),AgePl_Ltm1(ij,:),...
                    fab_L,fbe_L,ParEx_L,Mpar_L,TBio_L(ij,:),SAI_Ltm1(ij,:),hc_Ltm1(ij,:),...
                    ExEM,Lmax_day,L_day(ij,jDay),Se_bio,Tdp_bio,OPT_EnvLimitGrowth,OPT_VD,OPT_VCA,OPT_ALLOME,OPT_SoilBiogeochemistry);

                NutrientResetInfo_H(ij,:,:)=reshape(NutrientResetInfo_H_cell,[1,cc_max,4]);
                NutrientResetInfo_L(ij,:,:)=reshape(NutrientResetInfo_L_cell,[1,cc_max,4]);

                if (sum(Ccrown_t(ij,:)) + Cwat + Curb + Crock + Cbare) < 1 % adjust for harvest
                    Cbare = 1 - (sum(Ccrown_t(ij,:)) + Cwat + Curb + Crock);
                end

                %%%%%%
                BLit(ij,:)= 0.0 ; % %% %%[kg DM / m2]
            end

            %%%%%%%%%%%%%%%%%%% HYDROLOGY MODULE


            [q_run_after_E(ij),V(ij,:),O(ij,:),Vice(ij,:),Oice(ij,:),ZWT(ij),OF(ij),OS(ij),OH(ij,:),OL(ij,:),Psi_s_H(ij,:),Psi_s_L(ij,:),Rd(ij),Qi_out(ij,:),...
                Rh(ij),Lk(ij),f(ij),WIS(ij),Ts(ij),Csno(ij),Cice(ij),NDVI(ij),...
                Pr_sno(ij),Pr_liq(ij),rb_H(ij,:),rb_L(ij,:),rs_sunH(ij,:),rs_sunL(ij,:),rs_shdH(ij,:),rs_shdL(ij,:),...
                rap_H(ij,:),rap_L(ij,:),r_soil(ij),b_soil(ij),alp_soil(ij),ra(ij),r_litter(ij,:),...
                WR_SP(ij),U_SWE(ij),NIn_SWE(ij),dQ_S(ij),DQ_S(ij),DT_S(ij),...
                WAT(ij),ICE(ij),ICE_D(ij),IP_wc(ij),WR_IP(ij),NIce(ij),Cicew(ij),Csnow(ij),FROCK(ij),...
                Dr_H(ij,:),Dr_L(ij,:),SE_rock(ij),SE_urb(ij),Lk_wat(ij),Lk_rock(ij),...
                An_L(ij,:),An_H(ij,:),Rdark_L(ij,:),Rdark_H(ij,:),Ci_sunH(ij,:),Ci_sunL(ij,:),Ci_shdH(ij,:),Ci_shdL(ij,:),Rn(ij),...
                H(ij),QE(ij),Qv(ij),Lpho(ij),T_H(ij,:),T_L(ij,:),EIn_H(ij,:),EIn_L(ij,:),EG(ij),ELitter(ij),ESN(ij),ESN_In(ij),...
                EWAT(ij),EICE(ij),EIn_urb(ij),EIn_rock(ij),dw_SNO(ij),Imelt(ij),Smelt(ij),...
                G(ij),Gfin(ij),Tdp(ij,:),Tdpsnow(ij,:),Tdeb(ij,:),Tice(ij,:),Tdamp(ij),Tdp_H(ij,:),Tdp_L(ij,:),SWE(ij),SND(ij),ros(ij),In_SWE(ij),SP_wc(ij),Qfm(ij),t_sls(ij),...
                In_H(ij,:),In_L(ij,:),In_Litter(ij),In_urb(ij),In_rock(ij),...
                gsr_H(ij,:),Psi_x_H(ij,:),Psi_l_H(ij,:),Jsx_H(ij,:),Jxl_H(ij,:),Kleaf_H(ij,:),Kx_H(ij,:),Vx_H(ij,:),Vl_H(ij,:),...
                gsr_L(ij,:),Psi_x_L(ij,:),Psi_l_L(ij,:),Jsx_L(ij,:),Jxl_L(ij,:),Kleaf_L(ij,:),Kx_L(ij,:),Vx_L(ij,:),Vl_L(ij,:),...
                fapar_H(ij,:),fapar_L(ij,:),SIF_H(ij,:),SIF_L(ij,:),...
                Ws_under(ij),er(ij),snow_albedo(ij,:),tau_sno(ij),e_sno(ij),...
                HV(ij),QEV(ij),dQVEG(ij),TsVEG(ij),Ts_under(ij),EK(ij),POT(ij,:),CK1(ij)]=...
                ...
                HYDROLOGY_MODULE_PAR...
                ...
                (Vtm1(ij,:),Oicetm1(ij,:),...
                aR,Zs,...
                EvL_Zs,Inf_Zs,Zinf,RfH_Zs,RfL_Zs,dz,Dz,ms,Kbot,Pr_S(ij),Ta_S(ij),Ds_S(ij),Ws_S(ij),zatm,Tstm1(ij),dt,dth,ea_S(ij),N_S(ij),Pre_S(ij),Tstm0(ij),...
                LAI_H(ij,:),SAI_H(ij,:),LAI_L(ij,:),SAI_L(ij,:),LAIdead_H(ij,:),LAIdead_L(ij,:),...
                Rrootl_H(ij,:),Rrootl_L(ij,:),BLit(ij,:),Sllit,Kct,...
                Datam_S,DeltaGMT,Lon(ij),Lat(ij),t_bef,t_aft,...
                Ccrown,Cbare,Crock,Curb,Cwat,...
                SAB1_S(ij),SAB2_S(ij),SAD1_S(ij),SAD2_S(ij),PARB_S(ij),PARD_S(ij),SvF(ij),SNDtm1(ij),...
                snow_albedotm1(ij,:),...
                Color_Class,OM_H,OM_L,PFT_opt_H,PFT_opt_L,hc_H(ij,:),hc_L(ij,:),d_leaf_H,d_leaf_L,...
                Soil_Param,Interc_Param,SnowIce_Param,VegH_Param,VegL_Param,...
                Ca_S(ij),Oa,Citm1_sunH(ij,:),Citm1_shdH(ij,:),Citm1_sunL(ij,:),Citm1_shdL(ij,:),...
                e_rel_H(ij,:),e_relN_H(ij,:),e_rel_L(ij,:),e_relN_L(ij,:),...
                e_snotm1(ij),In_Htm1(ij,:),In_Ltm1(ij,:),In_Littertm1(ij),In_urbtm1(ij),In_rocktm1(ij),SWEtm1(ij),In_SWEtm1(ij),....
                Tdebtm1(ij,:),Ticetm1(ij),Tdptm1(ij,:),Tdpsnowtm1(ij,:),Tdamptm1(ij),Ts_undertm1(ij),...
                WATtm1(ij),ICEtm1(ij),IP_wctm1(ij),ICE_Dtm1(ij),Cicewtm1(ij),...
                Vx_Htm1(ij,:),Vl_Htm1(ij,:),Vx_Ltm1(ij,:),Vl_Ltm1(ij,:),Psi_x_Htm1(ij,:),Psi_l_Htm1(ij,:),Psi_x_Ltm1(ij,:),Psi_l_Ltm1(ij,:),...
                ZR95_H,ZR95_L,...
                FROCKtm1(ij),Krock,...
                Urb_Par,Deb_Par,Zs_deb,...
                Tdew_S(ij),t_slstm1(ij),rostm1(ij),SP_wctm1(ij),fpr,IrD_S(ij),...
                In_max_urb,In_max_rock,K_usle,tau_snotm1(ij),Ta_t(ij,:),...
                Slo_top(ij),Slo_head(ij,:),Asur(ij),Ared(ij),aTop(ij),EKtm1(ij),q_runon(ij),Qi_in(ij,:),...
                Ws_undertm1(ij),Pr_sno_t(ij,:),...
                pow_dis,a_dis,Salt_S(ij),...
                SPAR,SNn(ij),OPT_min_SPD,OPT_VegSnow,OPT_SoilTemp,OPT_PlantHydr,Opt_CR,Opt_ST,Opt_ST2,OPT_SM,OPT_STh,OPT_FR_SOIL,OPT_PH);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            if OPT_SoilBiogeochemistry==1
                NEE(ij) = -(NPP_H(ij,:)+NPP_L(ij,:))*Ccrown' + R_litter(ij) + R_microbe(ij) + R_ew(ij);%Net ecocystem exchange
            end
        end
    end

    [resetCellH,resetPftH]=find(NutrientReset_H);
    [resetCellL,resetPftL]=find(NutrientReset_L);
    if ~isempty(resetCellH) || ~isempty(resetCellL)
        resetCell=[resetCellH;resetCellL];
        resetPft=[resetPftH;resetPftL];
        resetLayer=[ones(length(resetCellH),1);2*ones(length(resetCellL),1)];
        resetIndexH=sub2ind([num_cell,cc_max],resetCellH,resetPftH);
        resetIndexL=sub2ind([num_cell,cc_max],resetCellL,resetPftL);
        resetInfoH=zeros(length(resetCellH),4);
        resetInfoL=zeros(length(resetCellL),4);
        for resetInfoColumn=1:4
            resetInfoPlane=NutrientResetInfo_H(:,:,resetInfoColumn);
            resetInfoH(:,resetInfoColumn)=resetInfoPlane(resetIndexH);
            resetInfoPlane=NutrientResetInfo_L(:,:,resetInfoColumn);
            resetInfoL(:,resetInfoColumn)=resetInfoPlane(resetIndexL);
        end
        resetInfo=[resetInfoH;resetInfoL];
        [resetRow,resetCol]=ind2sub([m_cell,n_cell],resetCell);
        nReset=length(resetCell);
        resetRows=[repmat(t,nReset,1),repmat(Date(t),nReset,1),...
            repmat(Datam(t,1:4),nReset,1),resetCell,resetRow,resetCol,resetLayer,resetPft,resetInfo];
        NutrientResetLog=[NutrientResetLog;resetRows];
        NutrientResetEver_H=NutrientResetEver_H | NutrientReset_H;
        NutrientResetEver_L=NutrientResetEver_L | NutrientReset_L;
        NutrientResetUniqueCellCount=nnz(any(NutrientResetEver_H | NutrientResetEver_L,2));
        NutrientResetUniquePftCount=nnz(NutrientResetEver_H)+nnz(NutrientResetEver_L);
        save(NutrientResetLogFile,'NutrientResetLog','NutrientResetLogColumns',...
            'NutrientResetLayerCodes','NutrientResetReasonCodes','NutrientResetEver_H','NutrientResetEver_L',...
            'NutrientResetUniqueCellCount','NutrientResetUniquePftCount');
        fprintf('Nutrient reset at Iter %d: %d PFT(s) cleared; %d unique cell(s) so far.\n',...
            t,nReset,NutrientResetUniqueCellCount);
    end

    ix_rd = ix_rd + 1; %netcdf Reading index - Update for next step
    P((P<0))=0;%%% sometimes if the previpus step a P value is 0, it becomes -e-15 in the next step. avoid this issue.
    for ii=1:ms_max
        Qi_out_Rout(:,:,ii) =  reshape(Qi_out(:,ii),m_cell,n_cell); %%[mm/h]
    end


    Rd = reshape(Rd,m_cell,n_cell); %%[mm]
    Rh = reshape(Rh,m_cell,n_cell); %%[mm]
    
    %%%%%%%%%%% ROUTING MODULE
    %%%%%%%%%%%%%%%%%%%%%%%% BG added %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%
    if OPT_SoilBiogeochemistry==1
        P_rout_old(:,:)=[P(:,31),P(:,32),P(:,43),P(:,52),P(:,33),P(:,47),P(:,12),P(:,13)];% NH4,NO3,P,K,DON,DOP,DOC
        if RoutingTracerDiagnosticFid > 0
            RoutingTracerDiagnosticPoolIndices=[31,32,43,52,33,47,12,13];
            RoutingTracerDiagnosticPBgBefore=Ptm1(RoutingTracerDiagnosticCell,RoutingTracerDiagnosticPoolIndices);
            RoutingTracerDiagnosticPBgAfter=P_rout_old(RoutingTracerDiagnosticCell,:);
            RoutingTracerDiagnosticDbFluxSum=reshape(sum(dB_flux(RoutingTracerDiagnosticCell,RoutingTracerDiagnosticPoolIndices,:),3),1,[]);
            RoutingTracerDiagnosticRunonIn=q_runon(RoutingTracerDiagnosticCell);
            RoutingTracerDiagnosticQiIn=sum(Qi_in(RoutingTracerDiagnosticCell,:));
            RoutingTracerDiagnosticQiOut=sum(Qi_out(RoutingTracerDiagnosticCell,:));
        end
        [TR_ch,TR_ice,TR_sur,P_rout,LEAK_nt,N_check,q_runon,Q_channel,Qi_in_Rout,Slo_pot,...
            Q_exit,Qsub_exit,...
            TR_exit,...
            TR_exitS,...
            TRQexit,...
            T_pot,QpointH,QpointC,UpointH,UpointC,RoutingTracerDiagnostic]= ROUTING_MODULE_BG_v4(TR_ch,TR_ice,Oicetm1,Vice,dz,N_nt,TR_sur,q_run_after_E,Asur,Ared,EG,f,Pr_S,P_rout_old,Vtm1,Lk,dt,dth,Rd,Rh,Qi_out_Rout,Q_channel,...
            cellsize,Area,DTM,NMAN_H,NMAN_C,MRough,WC,SN,T_flow,T_pot,Slo_top,ms_max,POT,ZWT,OPT_HEAD,Xout,Yout,RoutingTracerDiagnosticCell);
         %%%% get the P_rout diff
        P_rout_dif=P_rout-P_rout_old;
        % reallocate back to the P pool
        P(:,31)=P_rout(:,1);
        P(:,32)=P_rout(:,2);
        P(:,43)=P_rout(:,3);
        P(:,52)=P_rout(:,4);
        P(:,33)=P_rout(:,5);
        P(:,47)=P_rout(:,6);
        P(:,12)=P_rout(:,7);
        P(:,13)=P_rout(:,8);
    else
         [q_runon,Q_channel,Qi_in_Rout,Slo_pot,Q_exit,Qsub_exit,T_pot,...
            QpointH,QpointC,UpointH,UpointC]= ROUTING_MODULE(dt,dth,Rd,Rh,Qi_out_Rout,Q_channel,...
            cellsize,Area,DTM,NMAN_H,NMAN_C,MRough,WC,SN,T_flow,T_pot,Slo_top,ms_max,POT,ZWT,OPT_HEAD,Xout,Yout);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%
    % [q_runon,Q_channel,Qi_in_Rout,Slo_pot,Q_exit,Qsub_exit,T_pot,...
    %     QpointH,QpointC,UpointH,UpointC]= ROUTING_MODULE(dt,dth,Rd,Rh,Qi_out_Rout,Q_channel,...
    %     cellsize,Area,DTM,NMAN_H,NMAN_C,MRough,WC,SN,T_flow,T_pot,Slo_top,ms_max,POT,ZWT,OPT_HEAD,Xout,Yout);
    disp(Q_exit)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Qi_in = NaN(num_cell, ms_max);
    Slo_head = NaN(num_cell, ms_max);
    for ii=1:ms_max
        Qi_in(:,ii) =  reshape(Qi_in_Rout(:,:,ii),num_cell,1); %%[mm]
        Slo_head(:,ii) = reshape(Slo_pot(:,:,ii),num_cell,1);
    end
    Rd = reshape(Rd,num_cell,1); %%[mm]
    Rh = reshape(Rh,num_cell,1); %%[mm]
    q_runon = reshape(q_runon,num_cell,1); %%[mm]
    %%%%%%%%%%%%%
    Qi_in=Qi_in/dth;%%% [mm/h]
    q_runon = q_runon/dth; %%% [mm/h]
    if RoutingTracerDiagnosticFid > 0
        RoutingTracerDiagnosticBase=[t,Date(t),Datam(t,1:4),RoutingTracerDiagnosticCell,...
            Pr_S(RoutingTracerDiagnosticCell),Pr_liq(RoutingTracerDiagnosticCell),RoutingTracerDiagnosticRunonIn,...
            q_run_after_E(RoutingTracerDiagnosticCell),WIS(RoutingTracerDiagnosticCell),f(RoutingTracerDiagnosticCell),...
            Rh(RoutingTracerDiagnosticCell),Rd(RoutingTracerDiagnosticCell),EG(RoutingTracerDiagnosticCell),...
            sum(V(RoutingTracerDiagnosticCell,:)),sum(Vice(RoutingTracerDiagnosticCell,:)),...
            RoutingTracerDiagnosticQiIn,RoutingTracerDiagnosticQiOut,q_runon(RoutingTracerDiagnosticCell)];
        RoutingTracerDiagnosticTracerValues=reshape([RoutingTracerDiagnostic.TR_sur_in;RoutingTracerDiagnostic.TR_If;...
            RoutingTracerDiagnostic.TR_avg_pre;RoutingTracerDiagnostic.TR_avg_post;RoutingTracerDiagnostic.TR_sur_out;...
            RoutingTracerDiagnostic.infiltration_mass;RoutingTracerDiagnostic.et_mass;RoutingTracerDiagnostic.rd_loss_mass;...
            RoutingTracerDiagnostic.leakage_mass;RoutingTracerDiagnostic.ice_mass;RoutingTracerDiagnostic.subsurface_out_mass;...
            RoutingTracerDiagnostic.subsurface_in_mass],1,[]);
        RoutingTracerDiagnosticPoolValues=reshape([RoutingTracerDiagnosticPBgBefore;RoutingTracerDiagnosticPBgAfter;...
            RoutingTracerDiagnosticPBgAfter-RoutingTracerDiagnosticPBgBefore;RoutingTracerDiagnosticDbFluxSum;...
            RoutingTracerDiagnostic.P_initial;RoutingTracerDiagnostic.P_final;...
            RoutingTracerDiagnostic.P_final-RoutingTracerDiagnostic.P_initial],1,[]);
        RoutingTracerDiagnosticRow=[RoutingTracerDiagnosticBase,RoutingTracerDiagnosticTracerValues,...
            RoutingTracerDiagnosticPoolValues,RoutingTracerDiagnostic.ice_volume_change];
        fprintf(RoutingTracerDiagnosticFid,[repmat('%.12g,',1,length(RoutingTracerDiagnosticRow)-1),'%.12g\n'],RoutingTracerDiagnosticRow);
        if mod(t,24)==0
            fclose(RoutingTracerDiagnosticFid);
            RoutingTracerDiagnosticFid=fopen(RoutingTracerDiagnosticFile,'a');
            if RoutingTracerDiagnosticFid < 0
                error('Cannot append routing tracer diagnostic file: %s',RoutingTracerDiagnosticFile)
            end
        end
    end
    
    %%% Q_exit Qsub_exit [mm] over the entire domain

    if not(isreal(sum(sum(q_runon))))
        disp('The Program fails because of runon numerical instability')
        %break
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% RESERVOIR COMPONENT %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if exist("RES_ID_List")
        if not(isempty(RES_ID_List))
            [WAT,Q_channel,q_runon,Q_out_Res,H_Res,VOL_Res]= RESERVOIRS(DTM,cellsize,t,dth,dt,SN,WAT,Q_channel,q_runon,...
                RES_ID_List,RES_ID,RES_Outlet,Res_prop,Res_TS);
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% AVALANCHES COMPONENT %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SND=reshape(SND,m_cell,n_cell);
    SWE=reshape(SWE,m_cell,n_cell);
    ros=reshape(ros,m_cell,n_cell);
    [SND,SWE,ros,Swe_exit]= AVALANCHES(DTM,cellsize,Area,reshape(Asur,m_cell,n_cell),Slo_top,SND,SWE,ros);
    SND= reshape(SND,num_cell,1);
    SWE = reshape(SWE,num_cell,1);
    ros = reshape(ros,num_cell,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% FRACTURED ROCK COMPONENT %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %[Q_channel,FROCK,Qflow_rock]= FRACTURED_ROCK(Q_channel,FROCK,SPRINGn,dth,m_cell,n_cell,num_cell,Kres_Rock);
    % [Q_channel,FROCK,Qflow_rock]= FRACTURED_ROCK(Q_channel,FROCK,SNn,dth,m_cell,n_cell,num_cell,Kres_Rock);
    if OPT_SoilBiogeochemistry==1
        [TR_rock,TR_ch,Q_channel,FROCK,Qflow_rock]= FRACTURED_ROCK_BG(N_nt,TR_rock,TR_ch,LEAK_nt,Q_channel,FROCK,FROCKtm1,SNn,dth,m_cell,n_cell,num_cell,Kres_Rock);
    else
        [Q_channel,FROCK,Qflow_rock]= FRACTURED_ROCK(Q_channel,FROCK,SNn,dth,m_cell,n_cell,num_cell,Kres_Rock);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if t==2
        V_tgtm1 =sum(Ared.*Asur.*sum(Vtm1,2)).*(cellsize.^2)/Area;
        Vice_tgtm1 =sum(Ared.*Asur.*sum(Vicetm1,2)).*(cellsize.^2)/Area;
        SWE_tgtm1 = sum(SWEtm1).*(cellsize.^2)/Area; %%
        In_tgtm1 = (sum(sum(In_Htm1)) + sum(sum(In_Ltm1)) + sum(sum(In_Littertm1)) + sum(SP_wctm1) + ...
            sum(In_SWEtm1) + sum(In_urbtm1) + sum(In_rocktm1) +sum(IP_wctm1) ).*(cellsize.^2)/Area ;
        ICE_tgtm1 =  sum(ICEtm1).*(cellsize.^2)/Area; %%
        WAT_tgtm1 =  sum(WATtm1).*(cellsize.^2)/Area; %%
        FROCK_tgtm1 =  sum(FROCKtm1).*(cellsize.^2)/Area; %%
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% INITIAL CONDITION FOR THE NEXT STEP        %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Tstm0 =2*Ts-Tstm1;
    Vtm1=V;
    Tstm1 =Ts;
    SNDtm1 = SND ;
    snow_albedotm1 = snow_albedo ;
    Citm1_sunH = Ci_sunH;
    Citm1_sunL = Ci_sunL;
    Citm1_shdH = Ci_shdH;
    Citm1_shdL = Ci_shdL;
    e_snotm1 = e_sno;
    In_Htm1 = In_H ;
    In_Ltm1 = In_L;
    In_Littertm1 = In_Litter;
    In_urbtm1 = In_urb;
    In_rocktm1 = In_rock;
    SWEtm1 = SWE ;
    In_SWEtm1 = In_SWE;
    Tdamptm1 = Tdamp;
    Tdptm1 = Tdp;
    Vx_Htm1 =  Vx_H;
    Vl_Htm1 = Vl_H;
    Vx_Ltm1 =  Vx_L;
    Vl_Ltm1 =  Vl_L;
    Psi_x_Htm1 = Psi_x_H;
    Psi_l_Htm1 = Psi_l_H;
    Psi_x_Ltm1 = Psi_x_L;
    Psi_l_Ltm1 = Psi_l_L;
    WATtm1=WAT;
    ICEtm1=ICE;
    IP_wctm1=IP_wc;
    ICE_Dtm1=ICE_D;
    Cicewtm1=Cicew;
    FROCKtm1=FROCK;
    t_slstm1= t_sls;
    rostm1=ros;
    SP_wctm1=SP_wc ;
    tau_snotm1 = tau_sno;
    EKtm1 = EK;
    Ws_undertm1 = Ws_under;
    Tdebtm1= Tdeb;
    Ticetm1= Tice;
    Ts_undertm1=Ts_under;
    Tdpsnowtm1=Tdpsnow;
    Oicetm1 = Oice;
    %%%%%%%%%%%%%%%%%% BG added %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if OPT_SoilBiogeochemistry==1
        Ptm1=P; %% now Ptm1 should be updated every hour
    end    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (Datam_S(4)==1)
        if OPT_SoilBiogeochemistry==1
            %%% BG added
            SupN_Htm1=SupN_H;
            SupP_Htm1=SupP_H;
            SupK_Htm1=SupK_H;
            SupN_Ltm1=SupN_L;
            SupP_Ltm1=SupP_L;
            SupK_Ltm1=SupK_L;
            RexmyItm1=RexmyI;
            NavlItm1=NavlI;
        end
        B_Htm1= B_H;
        PHE_S_Htm1=PHE_S_H;
        dflo_Htm1= dflo_H;
        AgeL_Htm1= AgeL_H;
        AgeDL_Htm1 = AgeDL_H;
        NPP_Htm1 =  NPP_H;
        NPPI_Htm1 = NPPI_H;
        TdpI_Htm1 = TdpI_H;
        Bfac_weekHtm1 = Bfac_weekH ;
        Nreserve_Htm1 = Nreserve_H;
        Preserve_Htm1 = Preserve_H;
        Kreserve_Htm1 = Kreserve_H;
        FNC_Htm1 = FNC_H;
        SAI_Htm1= SAI_H;
        hc_Htm1= hc_H;
        Tden_Htm1 =Tden_H;
        AgePl_Htm1 = AgePl_H;
        NBLeaf_Htm1=NBLeaf_H;
        PARI_Htm1 = PARI_H;
        NBLI_Htm1 = NBLI_H;
        NuLit_Htm1=NuLit_H;
        %%%
        B_Ltm1= B_L;
        PHE_S_Ltm1= PHE_S_L;
        dflo_Ltm1= dflo_L ;
        AgeL_Ltm1= AgeL_L ;
        AgeDL_Ltm1 = AgeDL_L;
        NPP_Ltm1 =  NPP_L;
        NPPI_Ltm1 = NPPI_L;
        TdpI_Ltm1 = TdpI_L;
        Bfac_weekLtm1 = Bfac_weekL ;
        Nreserve_Ltm1 = Nreserve_L;
        Preserve_Ltm1 = Preserve_L;
        Kreserve_Ltm1 = Kreserve_L;
        FNC_Ltm1 = FNC_L;
        SAI_Ltm1= SAI_L;
        hc_Ltm1= hc_L;
        Tden_Ltm1=  Tden_L;
        AgePl_Ltm1=  AgePl_L;
        NBLeaf_Ltm1=NBLeaf_L;
        PARI_Ltm1 = PARI_L;
        NBLI_Ltm1 = NBLI_L;
        NuLit_Ltm1=NuLit_L;
        %%%
        Ccrown_t_tm1 = Ccrown_t;
        NupI_Htm1 =  NupI_H;
        NupI_Ltm1 =  NupI_L;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% MEMORY CONDITION FOR VEGETATION  MODEL    %%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% 1 day
    if t > 24
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%% BG added %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if OPT_SoilBiogeochemistry==1
            Lk_t(:,1:23)=Lk_t(:,2:24);
            Lk_t(:,24)=Lk;
            T_H_t(:,:,1:23)=T_H_t(:,:,2:24);
            T_H_t(:,:,24)=T_H;
            T_L_t(:,:,1:23)=T_L_t(:,:,2:24);
            T_L_t(:,:,24)=T_L;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Ta_t(:,1:23)=Ta_t(:,2:24);
        Ta_t(:,24)=Ta_S(:);
        PAR_t(:,1:23)= PAR_t(:,2:24);
        PAR_t(:,24)=PARB_S(:) + PARD_S(:);
        Tdp_t(:,:,1:23)=Tdp_t(:,:,2:24);
        Tdp_t(:,:,24)=Tdp;
        O_t(:,:,1:23)=O_t(:,:,2:24);
        O_t(:,:,24)=O;
        V_t(:,:,1:23)=V_t(:,:,2:24);
        V_t(:,:,24)=V;
        Pr_sno_t(:,1:23)=Pr_sno_t(:,2:24);
        Pr_sno_t(:,24)=Pr_sno;
        %%%
        An_H_t(:,:,1:23)=An_H_t(:,:,2:24);
        An_H_t(:,:,24)=An_H;
        An_L_t(:,:,1:23)=An_L_t(:,:,2:24);
        An_L_t(:,:,24)=An_L;
        Rdark_H_t(:,:,1:23)=Rdark_H_t(:,:,2:24);
        Rdark_H_t(:,:,24)=Rdark_H;
        Rdark_L_t(:,:,1:23)=Rdark_L_t(:,:,2:24);
        Rdark_L_t(:,:,24)=Rdark_L;
        Psi_x_H_t(:,:,1:23)=Psi_x_H_t(:,:,2:24);
        Psi_x_H_t(:,:,24)=Psi_x_H;
        Psi_l_H_t(:,:,1:23)=Psi_l_H_t(:,:,2:24);
        Psi_l_H_t(:,:,24)=Psi_l_H;
        Tdp_H_t(:,:,1:23)=Tdp_H_t(:,:,2:24);
        Tdp_H_t(:,:,24)=Tdp_H;
        Psi_x_L_t(:,:,1:23)=Psi_x_L_t(:,:,2:24);
        Psi_x_L_t(:,:,24)=Psi_x_L;
        Psi_l_L_t(:,:,1:23)=Psi_l_L_t(:,:,2:24);
        Psi_l_L_t(:,:,24)=Psi_l_L;
        Tdp_L_t(:,:,1:23)=Tdp_L_t(:,:,2:24);
        Tdp_L_t(:,:,24)=Tdp_L;
    else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%% BG added %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if OPT_SoilBiogeochemistry==1
            Lk_t(:,t)=Lk;
            T_H_t(:,:,t)=T_H;
            T_L_t(:,:,t)=T_L;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Ta_t(:,t)=Ta_S(:);
        PAR_t(:,24)=PARB_S(:) + PARD_S(:);
        Tdp_t(:,:,t)=Tdp;
        O_t(:,:,t)=O;
        V_t(:,:,t)=V;
        An_H_t(:,:,t)=An_H;
        An_L_t(:,:,t)=An_L;
        Rdark_H_t(:,:,t)=Rdark_H;
        Rdark_L_t(:,:,t)=Rdark_L;
        Psi_x_H_t(:,:,t)=Psi_x_H;
        Psi_l_H_t(:,:,t)=Psi_l_H;
        Tdp_H_t(:,:,t)=Tdp_H;
        Psi_x_L_t(:,:,t)=Psi_x_L;
        Psi_l_L_t(:,:,t)=Psi_l_L;
        Tdp_L_t(:,:,t)=Tdp_L;
        Pr_sno_t(:,t)=Pr_sno;
    end

    %%%%%%%%%%%%%%%%% OUTPUT WRITING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if OPT_SoilBiogeochemistry==1
        % run('OUTPUT_MANAGER_PAR_BG_pix')
       % run('OUTPUT_MANAGER_bin_BG');
        % run('OUTPUT_MANAGER_PAR_BG');
        run('OUTPUT_MANAGER_PAR_BG_flexible');
    else   
        run('OUTPUT_MANAGER_bin');% to improve with nc files
    end
    %%%%%%%%%
  
    if  ismember(t, Out.SaveStep)  ||  t==N_time_step

        Fstep2= strcat(Fstep,'_',num2str(t));

        if strcmp(Out.SaveVariables, 'all')
            save([Out.OutDir, filesep, Fstep2]);%to add variable list to be saved
        else
            eval([' save("', Out.OutDir, filesep, Fstep2,'.mat",' Out.SaveVariables, ');']);%to add variable list to be saved
        end
    end

    if  jDay == 1 && Datam(t, 4) == 0 ||  t==N_time_step
        % Save workspace for restarts - Every start of year
        %save([Out.OutDir, '\State_', Fstep2]);%to add variable list to be saved
        save([Out.OutDir, '/State_', Fstep2]);%to add variable list to be saved

    end

end
if RoutingTracerDiagnosticFid > 0
    fclose(RoutingTracerDiagnosticFid);
end
Computational_Time =toc;
disp('COMPUTATIONAL TIME [h] ')
disp(Computational_Time/3600)
disp(' COMPUTATIONAL TIME [s/cycle] ')
disp(Computational_Time/N_time_step)
end

