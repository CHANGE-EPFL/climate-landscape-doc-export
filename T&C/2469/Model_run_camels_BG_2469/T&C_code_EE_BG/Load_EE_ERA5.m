function Load_EE_ERA5(year, month, shapefile_path, Meteo_sc, WorldClim, FolderLoc, t_bef, t_aft, DeltaGMT)
% Load 1 month worth of earth engine meteorological data
DeltaGMT = 0;

% mkdir(Meteo.FolderLoc)
%% Step 0: Delete previous month's Data: Must be a better way. That's not clean
files = fullfile(FolderLoc, '*.nc'); % Change to whatever pattern you need.
fn = dir(files);
for k = 1 : length(fn)
    baseFileName = fn(k).name;
    fullFileName = fullfile(FolderLoc, baseFileName);
    delete(fullFileName);
end

%% Step 1: Get the Data from the EE

py.ModuleERAmonthly.ERA5Land_mon_retrieve_save(year, month, shapefile_path, Meteo_sc, WorldClim, FolderLoc )
disp(year); disp(month)

%% Step 2: split radiation into bands + add Ca based on observations
%FolderLoc = 'WRF_preproccessing/TC_Vil_test/'; 
% Zbas = ncread([FolderLoc, '\DEM.nc'], 'DEM');
Zbas = ncread([FolderLoc, '/DEM.nc'], 'DEM');
t_y = (datenum(year, month, 1, 0, 0, 0):(1/24):datenum(year, month, eomday(year,month), 23, 0, 0))';

% SW = ncread([FolderLoc, '\S.nc'], 'S');%shortwave down [W/m2]
% Ta = ncread([FolderLoc, '\Ta.nc'], 'Ta');% air temperature
% U = ncread([FolderLoc, '\U.nc'], 'U');%relative humidity
% Pr = ncread([FolderLoc, '\Pr.nc'], 'Pr');%precip
SW = ncread([FolderLoc, '/S.nc'], 'S');%shortwave down [W/m2]
Ta = ncread([FolderLoc, '/Ta.nc'], 'Ta');% air temperature
U = ncread([FolderLoc, '/U.nc'], 'U');%relative humidity
Pr = ncread([FolderLoc, '/Pr.nc'], 'Pr');%precip
Pr =fillmissing(Pr,'constant',0);  
Pr(Pr<0.01) = 0; 

U  = fillmissing(U,'nearest');
Ta =fillmissing(Ta,'nearest');
SW = fillmissing(SW,'nearest');

% Calculate dew point
Tdew = 243.04*( log(U) + 17.625.*Ta./(243.04+Ta) )./(( 17.625 - log(U) - 17.625.*Ta./(243.04+Ta)));

% info = ncinfo([FolderLoc, '\S.nc']);
info = ncinfo([FolderLoc, '/S.nc']);
projM = projcrs(int32(str2num(info.Attributes.Value)));

% x = ncread([FolderLoc, '\S.nc'], 'x');
% y = ncread([FolderLoc, '\S.nc'], 'y');
x = ncread([FolderLoc, '/S.nc'], 'x');
y = ncread([FolderLoc, '/S.nc'], 'y');
[x, y] = meshgrid(x, y);
% Reproject to Geo System. Should be ideantical, but as a fail safe
[lat, lon] = projinv(projM, x, y);

SAD1 = NaN(size(SW));
SAD2 = NaN(size(SW));
SAB1 = NaN(size(SW));
SAB2 = NaN(size(SW));
PARB = NaN(size(SW));
PARD = NaN(size(SW));
N = NaN(size(SW));

% vars = ["Pr","Pre","S","U","Ws","Ta"]

for ii = 1:size(lon, 1)
    for jj = 1:size(lon, 2)

        % fprintf(['ii-%d jj-%d \n'],ii,jj)

        Tdewloop = Tdew(ii, jj, :);
        Prloop = Pr(ii, jj, :);

        % fix the time here
        SWt = SW(ii, jj, :);

        % tic
        if ~isnan(Zbas(ii, jj))
        % [~,~,SAD1t,SAD2t,SAB1t,SAB2t,PARBt,PARDt,Nt,~,~,~]...
        %     =...
        %     B_Automatic_Radiation_Partition(t_y ,lat(ii, jj), lon(ii, jj) ,Zbas(ii, jj),...
        %     DeltaGMT, squeeze(Prloop(1:length(t_y))), squeeze(Tdewloop(1:length(t_y))), squeeze(SWt(1:length(t_y))), 0, t_bef, t_aft);
        [~,~,SAD1t,SAD2t,SAB1t,SAB2t,PARBt,PARDt,Nt,~,~,~]...
            =...
            C_Automatic_Radiation_Partition(t_y ,lat(ii, jj), lon(ii, jj) ,Zbas(ii, jj),...
            DeltaGMT, squeeze(Prloop(1:length(t_y))), squeeze(Tdewloop(1:length(t_y))), squeeze(SWt(1:length(t_y))), 0, t_bef, t_aft);

        % toc
        SAD1(ii, jj, 1:length(t_y)) = SAD1t;
        SAD2(ii, jj, 1:length(t_y)) = SAD2t;
        SAB1(ii, jj, 1:length(t_y)) = SAB1t;
        SAB2(ii, jj, 1:length(t_y)) = SAB2t;
        PARB(ii, jj, 1:length(t_y)) = PARBt;
        PARD(ii, jj, 1:length(t_y)) = PARDt;

        N(ii, jj, 1:length(t_y)) = Nt;
        end
    end
end
% a = squeeze(SW(70,20,:));
% plot(a)
% figure
% b = squeeze(PARB(70,20,:));
% plot(b)


Ca = zeros(size(N));
Ca(:) = CO2ann(year);


nccreate([FolderLoc, '/N.nc'], "N", "Dimensions", {"x", size(x, 1), "y", size(y, 2), "time", length(t_y)}, "DeflateLevel",9 );
nccreate([FolderLoc, '/SAB1.nc'], "SAB1", "Dimensions", {"x", size(x, 1), "y", size(y, 2), "time", length(t_y)}, "DeflateLevel",9 );
nccreate([FolderLoc, '/SAB2.nc'], "SAB2", "Dimensions", {"x", size(x, 1), "y", size(y, 2), "time", length(t_y)}, "DeflateLevel",9 );
nccreate([FolderLoc, '/SAD1.nc'], "SAD1", "Dimensions", {"x", size(x, 1), "y", size(y, 2), "time", length(t_y)}, "DeflateLevel",9 );
nccreate([FolderLoc, '/SAD2.nc'], "SAD2", "Dimensions", {"x", size(x, 1), "y", size(y, 2), "time", length(t_y)}, "DeflateLevel",9 );
nccreate([FolderLoc, '/PARB.nc'], "PARB", "Dimensions", {"x", size(x, 1), "y", size(y, 2), "time", length(t_y)}, "DeflateLevel",9 );
nccreate([FolderLoc, '/PARD.nc'], "PARD", "Dimensions", {"x", size(x, 1), "y", size(y, 2), "time", length(t_y)}, "DeflateLevel",9 );
nccreate([FolderLoc, '/Ca.nc'], "Ca", "Dimensions", {"x", size(x, 1), "y", size(y, 2), "time", length(t_y)}, "DeflateLevel",9 );

V = {'Ca.nc','N.nc',...
    'SAB1.nc','SAB2.nc','SAD1.nc','SAD2.nc','PARB.nc',...
    'PARD.nc'};

for vv = 1:length(V)
    sn = [FolderLoc, '/'];
    nccreate([sn, V{vv}], "Date", "Dimensions", {"time", length(t_y)} );
    nccreate([sn, V{vv}], "x", "Dimensions", {"x", size(lon, 1), "y", size(lon, 2)} )
    nccreate([sn, V{vv}], "y", "Dimensions", {"x", size(lon, 1), "y", size(lon, 2)} )

    ncwrite([sn, V{vv}], "x", x)
    ncwrite([sn, V{vv}], "y", y)
    ncwrite([sn, V{vv}], "Date", t_y)
end

ncwrite([sn, 'SAB1.nc'], "SAB1", SAB1(:, :, 1:length(t_y)), [1 1 1]);
ncwrite([sn, 'SAB2.nc'], "SAB2", SAB2(:, :, 1:length(t_y)), [1 1 1]);
ncwrite([sn, 'SAD1.nc'], "SAD1", SAD1(:, :, 1:length(t_y)), [1 1 1]);
ncwrite([sn, 'SAD2.nc'], "SAD2", SAD2(:, :, 1:length(t_y)), [1 1 1]);
ncwrite([sn, 'PARB.nc'], "PARB", PARB(:, :, 1:length(t_y)), [1 1 1]);
ncwrite([sn, 'PARD.nc'], "PARD", PARD(:, :, 1:length(t_y)), [1 1 1]);
ncwrite([sn, 'N.nc'], "N", N(:, :, 1:length(t_y)), [1 1 1]);
ncwrite([sn, 'Ca.nc'], "Ca", Ca(:, :, 1:length(t_y)), [1 1 1]);


end

function Ca_c = CO2ann(year)

Data = ...
    [1979   336.85     0.11
    1980   338.91     0.07
    1981   340.11     0.09
    1982   340.86     0.03
    1983   342.53     0.06
    1984   344.07     0.08
    1985   345.54     0.07
    1986   346.97     0.07
    1987   348.68     0.10
    1988   351.16     0.07
    1989   352.78     0.07
    1990   354.05     0.07
    1991   355.39     0.07
    1992   356.09     0.06
    1993   356.83     0.07
    1994   358.33     0.08
    1995   360.17     0.05
    1996   361.93     0.04
    1997   363.05     0.05
    1998   365.70     0.04
    1999   367.80     0.05
    2000   368.96     0.06
    2001   370.57     0.05
    2002   372.59     0.04
    2003   375.15     0.04
    2004   376.95     0.06
    2005   378.98     0.05
    2006   381.15     0.05
    2007   382.90     0.04
    2008   385.02     0.05
    2009   386.50     0.04
    2010   388.76     0.06
    2011   390.63     0.05
    2012   392.65     0.06
    2013   395.40     0.06
    2014   397.34     0.05
    2015   399.65     0.05
    2016   403.06     0.06
    2017   405.22     0.07
    2018   407.61     0.07
    2019   410.07     0.07
    2020   412.44     0.06
    2021   414.70     0.07
    2022   417.07     0.07
    2023   419.31     0.15];% To be exapnded with climate projections

if year > max(Data(:,1))
    Ca_c = 420;
elseif year < min(Data(:,1))
    Ca_c = 330;
else
    ix = Data(:,1)==year;
    Ca_c = Data(ix,2);
end
end