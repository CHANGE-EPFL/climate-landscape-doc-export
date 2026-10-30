%% Create future CO2 Ca.nc file by remapping current annual Ca values
% The current Ca.nc already has the correct grid, Date, and hourly structure.
% This script copies that file and replaces the six unique annual Ca values
% with future SSP585 values for 2088-2093.

clear; clc;

%% User settings
site_ids = {'01054200', '01434025', '01545600', '07362587', ...
    '09352900', '10343500', '12056500', '12447390', ...
    '2203', '2420', '2469', '2481'};

% Resolve all paths relative to this script:
% <repository>/<site>/meteo/<site>/Ca.nc
script_dir = fileparts(mfilename('fullpath'));
repository_root = fileparts(script_dir);

future_years = [2088 2089 2090 2091 2092 2093];
future_ca = [971.10 984.68 998.32 1012.02 1025.74 1039.45];

out_name = 'Ca_585_88_93.nc';

%% Process sites
summary = table('Size', [0 4], ...
    'VariableTypes', {'string', 'string', 'double', 'string'}, ...
    'VariableNames', {'site_id', 'target_file', 'n_unique_current_ca', 'status'});

for ss = 1:numel(site_ids)
    site_id = site_ids{ss};
    target_dir = fullfile(repository_root, site_id, 'meteo', site_id);
    source_file = fullfile(target_dir, 'Ca.nc');
    target_file = fullfile(target_dir, out_name);

    fprintf('\n=== %s ===\n', site_id);

    try
        if ~isfolder(target_dir)
            mkdir(target_dir);
        end

        if ~isfile(source_file)
            error('Source Ca.nc not found: %s', source_file);
        end

        [current_unique, counts] = create_one_future_ca_file( ...
            source_file, target_file, future_years, future_ca);

        fprintf('Saved future Ca file:\n%s\n', target_file);
        fprintf('Source current values -> future SSP585 values:\n');
        for ii = 1:numel(current_unique)
            fprintf('  %.2f -> %d %.2f ppm, count %d\n', ...
                current_unique(ii), future_years(ii), future_ca(ii), counts(ii));
        end

        summary(end+1, :) = {string(site_id), string(target_file), numel(current_unique), "ok"}; %#ok<SAGROW>
    catch ME
        warning('Failed to create future Ca file for %s: %s', site_id, ME.message);
        summary(end+1, :) = {string(site_id), string(target_file), NaN, string(ME.message)}; %#ok<SAGROW>
    end
end

fprintf('\nBatch summary:\n');
disp(summary);

function [current_unique, counts] = create_one_future_ca_file(source_file, target_file, future_years, future_ca)
copyfile(source_file, target_file, 'f');

Ca = ncread(source_file, 'Ca');
current_unique = unique(Ca(:));
current_unique = current_unique(isfinite(current_unique));
current_unique = sort(current_unique(:));

if numel(current_unique) ~= numel(future_ca)
    error('Expected %d unique Ca values but found %d in %s.', ...
        numel(future_ca), numel(current_unique), source_file);
end

Ca_future = Ca;
counts = zeros(numel(current_unique), 1);
for ii = 1:numel(current_unique)
    mask = Ca == current_unique(ii);
    counts(ii) = nnz(mask);
    Ca_future(mask) = future_ca(ii);
end

ncwrite(target_file, 'Ca', Ca_future);

% Add lightweight provenance attributes. Ignore gracefully if the file format
% does not allow writing attributes in a given environment.
try
    ncwriteatt(target_file, 'Ca', 'future_scenario', 'ssp585');
    ncwriteatt(target_file, 'Ca', 'future_years', sprintf('%d-%d', future_years(1), future_years(end)));
    ncwriteatt(target_file, 'Ca', 'source_current_file', 'Ca.nc in the same directory');
    ncwriteatt(target_file, 'Ca', 'note', ...
        'Current six annual Ca values were remapped in sorted order to SSP585 2088-2093 annual CO2 values.');
catch
end

u_future = unique(Ca_future(:));
u_future = u_future(isfinite(u_future));
if numel(u_future) ~= numel(future_ca) || any(abs(sort(u_future(:)) - future_ca(:)) > 1e-6)
    error('Output future Ca values do not match requested SSP585 values in %s.', target_file);
end
end
