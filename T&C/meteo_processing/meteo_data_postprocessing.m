function meteo_data_postprocessing(site_ids, overwrite)
%METEO_DATA_POSTPROCESSING Create derived T&C meteorological inputs.
%
% Run from any directory. By default, all study sites are processed:
%   meteo_data_postprocessing
%
% Process one or more selected sites:
%   meteo_data_postprocessing("01434025")
%   meteo_data_postprocessing(["01434025", "2203"])
%
% Replace existing derived files:
%   meteo_data_postprocessing("01434025", true)

if nargin < 1 || isempty(site_ids)
    site_ids = ["01054200", "01434025", "01545600", "07362587", ...
        "09352900", "10343500", "12056500", "12447390", ...
        "2203", "2420", "2469", "2481"];
end
if nargin < 2
    overwrite = false;
end

site_ids = string(site_ids(:));
script_dir = fileparts(mfilename("fullpath"));
repository_root = fileparts(script_dir);
addpath(script_dir);

for site_index = 1:numel(site_ids)
    site_id = site_ids(site_index);
    [start_year, end_year] = site_period(site_id);
    process_site(repository_root, site_id, start_year, end_year, overwrite);
end
end


function process_site(repository_root, site_id, start_year, end_year, overwrite)
period = sprintf("%d%dyr", start_year, end_year);
meteo_dir = fullfile(repository_root, site_id, "meteo", site_id);

input_files = ["DEM.nc", ...
    "Pr_" + period + ".nc", "Pre_" + period + ".nc", ...
    "S_" + period + ".nc", "U_" + period + ".nc", ...
    "Ws_" + period + ".nc", "Ta_" + period + ".nc"];
for file_index = 1:numel(input_files)
    input_path = fullfile(meteo_dir, input_files(file_index));
    if ~isfile(input_path)
        error("Missing input for site %s: %s", site_id, input_path)
    end
end

output_variables = ["N", "SAB1", "SAB2", "SAD1", "SAD2", ...
    "PARB", "PARD", "Ca"];
output_files = fullfile(meteo_dir, output_variables + ".nc");
output_exists = arrayfun(@isfile, output_files);

if all(output_exists) && ~overwrite
    fprintf("Site %s already processed; skipping.\n", site_id);
    return
elseif any(output_exists) && ~overwrite
    error(["Site %s has an incomplete set of derived files. " ...
        "Run with overwrite=true to replace them."], site_id)
end

fprintf("\n=== Postprocessing site %s: %d-%d ===\n", ...
    site_id, start_year, end_year);
fprintf("Meteo directory: %s\n", meteo_dir);

dem_path = fullfile(meteo_dir, "DEM.nc");
shortwave_path = fullfile(meteo_dir, "S_" + period + ".nc");
temperature_path = fullfile(meteo_dir, "Ta_" + period + ".nc");
humidity_path = fullfile(meteo_dir, "U_" + period + ".nc");
precipitation_path = fullfile(meteo_dir, "Pr_" + period + ".nc");

elevation = ncread(dem_path, "DEM");
shortwave = ncread(shortwave_path, "S");
temperature = ncread(temperature_path, "Ta");
relative_humidity = ncread(humidity_path, "U");
precipitation = ncread(precipitation_path, "Pr");

precipitation = fillmissing(precipitation, "constant", 0);
precipitation(precipitation < 0.01) = 0;
relative_humidity = fillmissing(relative_humidity, "nearest");
temperature = fillmissing(temperature, "nearest");
shortwave = fillmissing(shortwave, "nearest");

time_datetime = (datetime(start_year, 1, 1, 0, 0, 0):hours(1): ...
    datetime(end_year, 12, 31, 23, 0, 0)).';
time_datenum = datenum(time_datetime);
number_of_hours = numel(time_datenum);

if size(shortwave, 3) ~= number_of_hours
    error(["Unexpected time dimension for site %s: S contains %d hours, " ...
        "but %d are expected for %d-%d."], site_id, size(shortwave, 3), ...
        number_of_hours, start_year, end_year)
end
if ~isequal(size(shortwave), size(temperature), size(relative_humidity), ...
        size(precipitation))
    error("The S, Ta, U and Pr arrays have different dimensions for site %s.", ...
        site_id)
end

dew_point = 243.04 .* (log(relative_humidity) + ...
    17.625 .* temperature ./ (243.04 + temperature)) ./ ...
    (17.625 - log(relative_humidity) - ...
    17.625 .* temperature ./ (243.04 + temperature));

crs_code = str2double(string(ncreadatt(shortwave_path, "/", "crs")));
projected_crs = projcrs(crs_code);
x_coordinate = ncread(shortwave_path, "x");
y_coordinate = ncread(shortwave_path, "y");
[x_grid, y_grid] = meshgrid(x_coordinate, y_coordinate);
[latitude, longitude] = projinv(projected_crs, x_grid, y_grid);

spatial_size = size(shortwave, [1, 2]);
if ~isequal(size(elevation), spatial_size) || ...
        ~isequal(size(latitude), spatial_size)
    error("The DEM and meteorological grids do not align for site %s.", site_id)
end

SAD1 = nan(size(shortwave));
SAD2 = nan(size(shortwave));
SAB1 = nan(size(shortwave));
SAB2 = nan(size(shortwave));
PARB = nan(size(shortwave));
PARD = nan(size(shortwave));
N = nan(size(shortwave));

t_bef = 0;
t_aft = 1;
delta_gmt = 0;

for row = 1:size(longitude, 1)
    fprintf("  Radiation partition: row %d of %d\n", row, size(longitude, 1));
    for column = 1:size(longitude, 2)
        if isnan(elevation(row, column))
            continue
        end

        [~, ~, SAD1_cell, SAD2_cell, SAB1_cell, SAB2_cell, ...
            PARB_cell, PARD_cell, N_cell, ~, ~, ~] = ...
            C_Automatic_Radiation_Partition(time_datenum, ...
            latitude(row, column), longitude(row, column), ...
            elevation(row, column), delta_gmt, ...
            squeeze(precipitation(row, column, :)), ...
            squeeze(dew_point(row, column, :)), ...
            squeeze(shortwave(row, column, :)), 0, t_bef, t_aft);

        SAD1(row, column, :) = SAD1_cell;
        SAD2(row, column, :) = SAD2_cell;
        SAB1(row, column, :) = SAB1_cell;
        SAB2(row, column, :) = SAB2_cell;
        PARB(row, column, :) = PARB_cell;
        PARD(row, column, :) = PARD_cell;
        N(row, column, :) = N_cell;
    end
end

Ca = zeros(size(N));
hour_year = year(time_datetime);
for current_year = start_year:end_year
    Ca(:, :, hour_year == current_year) = historical_co2(current_year);
end

if overwrite
    for file_index = 1:numel(output_files)
        if isfile(output_files(file_index))
            delete(output_files(file_index));
        end
    end
end

write_forcing(meteo_dir, "N", N, x_grid, y_grid, time_datenum);
write_forcing(meteo_dir, "SAB1", SAB1, x_grid, y_grid, time_datenum);
write_forcing(meteo_dir, "SAB2", SAB2, x_grid, y_grid, time_datenum);
write_forcing(meteo_dir, "SAD1", SAD1, x_grid, y_grid, time_datenum);
write_forcing(meteo_dir, "SAD2", SAD2, x_grid, y_grid, time_datenum);
write_forcing(meteo_dir, "PARB", PARB, x_grid, y_grid, time_datenum);
write_forcing(meteo_dir, "PARD", PARD, x_grid, y_grid, time_datenum);
write_forcing(meteo_dir, "Ca", Ca, x_grid, y_grid, time_datenum);

fprintf("Completed site %s.\n", site_id);
end


function write_forcing(meteo_dir, variable_name, values, x_grid, y_grid, time)
output_path = fullfile(meteo_dir, variable_name + ".nc");
number_of_rows = size(values, 1);
number_of_columns = size(values, 2);
number_of_hours = size(values, 3);

nccreate(output_path, variable_name, "Dimensions", ...
    {"x", number_of_rows, "y", number_of_columns, "time", number_of_hours}, ...
    "DeflateLevel", 9);
nccreate(output_path, "Date", "Dimensions", {"time", number_of_hours});
nccreate(output_path, "x", "Dimensions", ...
    {"x", number_of_rows, "y", number_of_columns});
nccreate(output_path, "y", "Dimensions", ...
    {"x", number_of_rows, "y", number_of_columns});

ncwrite(output_path, variable_name, values);
ncwrite(output_path, "Date", time);
ncwrite(output_path, "x", x_grid);
ncwrite(output_path, "y", y_grid);
end


function [start_year, end_year] = site_period(site_id)
switch site_id
    case {"01434025", "01545600", "10343500", "12447390"}
        start_year = 2009;
        end_year = 2014;
    case {"01054200", "07362587", "09352900", "12056500"}
        start_year = 2013;
        end_year = 2018;
    case {"2203", "2420", "2469", "2481"}
        start_year = 2015;
        end_year = 2020;
    otherwise
        error("Unknown site ID: %s", site_id)
end
end


function concentration = historical_co2(year_value)
years = (1979:2023).';
concentrations = [336.85; 338.91; 340.11; 340.86; 342.53; 344.07; ...
    345.54; 346.97; 348.68; 351.16; 352.78; 354.05; 355.39; ...
    356.09; 356.83; 358.33; 360.17; 361.93; 363.05; 365.70; ...
    367.80; 368.96; 370.57; 372.59; 375.15; 376.95; 378.98; ...
    381.15; 382.90; 385.02; 386.50; 388.76; 390.63; 392.65; ...
    395.40; 397.34; 399.65; 403.06; 405.22; 407.61; 410.07; ...
    412.44; 414.70; 417.07; 419.31];

index = find(years == year_value, 1);
if isempty(index)
    error("Historical CO2 is unavailable for year %d.", year_value)
end
concentration = concentrations(index);
end
