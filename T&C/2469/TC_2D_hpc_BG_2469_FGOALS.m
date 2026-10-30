clear 
clc
 
%% Set the python environment with the required packages - Figure out how to get this done on the HPC - Should be done with -venv
% myCluster = parcluster('Processes')
% delete(myCluster.Jobs)
% %% Open threads for parallel computing
%parpool('local', 128) %- If enabled you need to change in the RUN2D file the main loop from "for" to "parfor"

%% Code directory
addpath([pwd filesep 'Model_run_camels_BG_2469' filesep 'T&C_code_EE_BG']);
addpath(genpath(fullfile(pwd,'geophysical','Topotoolbox_2D','TopoToolbox')))


%% Set up the simulation
future_climat='FGOALS-g3';
basin_name = '2469'; 
TITLE_SAVE = '2469_AFGOALS-g3'; % in which folder to save the results
Date = datenum([2015 1 1 0 0 0]):(1/24):datenum([2020 12 31 23 0 0]);

Meteo.fn_future=['Factor_' basin_name '_' future_climat];
%% Load weather data - this example assumes uniform weather. It can be changed to 2D fields

Meteo.Interp = 1;%1 near neighbour, 2: linear interpolation; 0 same grid; 3 time series
Meteo.FolderLoc = [pwd filesep 'meteo' filesep '2469']; 
% pr = ncread( [fir_dir,'\meteo\shp_12285001_110\Pr_cr_era5_wc.nc'],'Pr_cr_era5_wc');
% %%
% for i =1:35064
%     imagesc(pr(:,:,i))
% end 
Meteo.Meteo_sc = 250;

Meteo.x1 = 1;% which time steps to run x1: start, x2:end
Meteo.x2 = 365*24*6;% 
% Meteo.DeltaGMT = -5;% Only if time stamps are in local times and not GMT
% Meteo.t_bef = 0;
% Meteo.t_aft = 1;
Meteo.proj = false;%if false, Meteo is in lon, lat. If not, in x, y. Assumption same proj with Geo
Meteo.Tls = 5.5; %lapse rate [K/km] - Could be extended to get data from WorldClim (currently not there)
%% Load spatial domain - See preprocessing step for how to set up the variables inside the Geo structure
Geo = load([pwd filesep 'geophysical' filesep 'distributed_scale_',basin_name,'_reso250_T.mat']);
%Geo = load("TC_Earth_Engine/Vil_4000.mat");
Geo.lai_m =double(Geo.lai_m);
%%%%%%%%%%%%%%%%% if no need to change any landcover
% Geo.PFT_def=permute(Geo.PFT_def,[2 3 1]);
% Geo.IIpar = logical(Geo.II); 
% Geo.PFT_def = Geo.PFT_def(:,:,Geo.IIpar); 
 
%%%%%%% mannuly remove crop, replace with bare land
Geo.PFT_def=permute(Geo.PFT_def,[2 3 1]);
Geo.IIpar = logical(Geo.II); 
Geo.IIpar_crop=Geo.IIpar;
Geo.IIpar_crop(1:7)=0;
Geo.crop_def=Geo.PFT_def(:,:,Geo.IIpar_crop);
Geo.crop_def=sum(Geo.crop_def,3);
Geo.Cbare=Geo.Cbare+Geo.crop_def;
Geo.IIpar(8:end)=0;
Geo.II(8:end)=0;
Geo.PFT_def = Geo.PFT_def(:,:,Geo.IIpar);  %%% this is a m*n*II map
 

%%% load Zsmax information
% Geo.Zsmax = round(T3_sub.soil_depth_statsgo(T3_sub.gauge_id_str==basin_name),1)*1000;
Geo.Zsmax = 1500; %% further information can be used

Geo.DeltaGMT = 1;  %%% mannual?
Geo.KBOT = Geo.KBOT*0.3;  %%% mannual?
Geo.Kres_Rock = Geo.Kres_Rock;
%% Parameters load
%addpath('ModelParams') % folder where the model parameters are defined
addpath(['Model_run_camels_BG_2469' filesep 'MOD_PARAM_BG_2469' filesep])

%addpath('TC_Earth_Engine/') % folder where the model parameters are defined

%addpath('TC_Earth_Engine/MOD_PARAM generalisation') % folder where the model parameters are defined

Param.cc_max = sum(Geo.II==1); % max number of PFTs in each grid
Param.ms_max = 7; % number of soil layers
Param.md_max = 1; % debris layers
Param.MOD_PARAM = @PARAMETERS_2D_PFTs; % Function including the parameters for all PFTs in the area
Param.INI = @INITIAL_CONDITIONS_2Dnc_BG; % Script for initializing the simulation
% rmpath('ModelParams')
% Prmean = nansum(squeeze(ncread( [Meteo.FolderLoc, '\Pr_hr_20202023yr.nc'], 'Pr_hr', [1, 1, 1], [1, 1, Inf]) ))/4;
% Tamean = nanmean(squeeze(ncread( [Meteo.FolderLoc, '\Ta_20202023yr.nc'], 'Ta',  [1, 1, 1], [1, 1, Inf])));
Param.param_sname = [pwd filesep 'Parameters_BG' filesep 'param_2469_zs1500_kbot004_VegCal.mat']; %%%kbot 0.01
Param.basin_name=basin_name;

%%% to enable RF spin up
Param.spin=1;
%% Outout setup
%%% why not Geo.Xout, sometimes Geo.Xout is different from real measurement
%%% point
Out.Yout =[136]; % i=indices of points to track [i MATLAB coordinate]
Out.Xout = [75]; % j=indices of points to track [j MATLAB coordinate]
% Out.Yout = [114 69]; % i=indices of points to track [i MATLAB coordinate]
% Out.Xout = [33 39]; % j=indices of points to track [j MATLAB coordinate]

Out.OutDir = TITLE_SAVE;
Out.SaveStep = [12:24:24*365*6]; %how frequently to save fields
%Out.SaveVariables = 'all'; %select variables to save
Out.SaveVariables = '"RA_H","OH","OL","LAI_H","LAI_L","NPP_H","NPP_L","QE","H","Rn","SND","Rd","Rh","Csno","f","G","GPP_H_spatial","GPP_L_spatial","O_space","Q_channel","T_H","T_L","EIn_H", "EIn_H", "EIn_L","EG","ESN", "ESN_In","EWAT","EICE","Dr_H","Dr_L", "EIn_urb","EIn_rock","Qi_in","Qi_out","q_runon"'; %select variables to save ... expand when necessary


%% Run 
% parpool('local',16)
%%%% final
RUN_TC_2D_hpc_BG_future_template_CH(Meteo, Geo, Param, Out,Date)  


