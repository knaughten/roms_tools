% Create ROMS file for surface fluxes of heat, freshwater, and momentum
% for the circumpolar quarter-degree configuration of ROMS. 
% Use monthly averages of the daily Tamura data for heat and freshwater 
% fluxes, and the monthly ERA-Interim renalysis for winds.

% Code by Kaitlin Alexander, adapted from the code by David Gwyther:
% in /ds/projects/iomp/totten/ana/dgwyther/netcdf/tisom008/sbc, the files
% make_sbc.m, read_tamura_daily.m, ERA_interim_misom_grid_stress_annual.m,
% and do_ISOM_sbc_nc_daily_1file_v1.m.

function caisom_sbc

% User-defined parameters

% Path to grid file
grid_file = 'caisom001_OneQuartergrd.nc';
% Desired path to output file
output_file = 'caisom001_surfacefluxes.nc';
% Years of data to select (exists from 1992-2012)
min_year = 1992;
max_year = 2012;
% If year_avg=true, a monthly climatology will be created by averaging
% each month's fluxes over the selected years.
year_avg = true;

% Parameters for conversions
rhoAir = 1.3; % density of air
Cd = 1.4e-3; % drag coefficient
ref_salt = 34.4; % reference sea surface salinity

% End user-defined parameters


% Load some ROMS scripts
addpath(genpath('/ds/projects/iomp/matlab_scripts'));

% Path to wind file
wind_file = '/ds/projects/iomp/obs/ERA_Interim/ERA_Interim_1992_2012_monthly.nc';

% Paths to Tamura data
tamura_base = '/ds/projects/iomp/obs/Tamura_air_sea_fluxes/daily_latest/TSDM2hb_';
tamura_grid = '/ds/projects/iomp/obs/Tamura_air_sea_fluxes/daily/latlon.data';
month_names = ['jan'; 'feb'; 'mar'; 'apr'; 'may'; 'jun'; 'jul'; 'aug'; 'sep'; ...
    'oct'; 'nov'; 'dec'];
days_per_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
leap_years = [1992:4:max_year];

% Read wind data
ncload(wind_file, 'u10', 'v10', 'longitude', 'latitude');
% Re-order the axes for u and v
u10 = permute(u10, [3 2 1]);
v10 = permute(v10, [3 2 1]);
% Mask out missing values
u10(u10 == -32767) = NaN;
v10(v10 == -32767) = NaN;
% Deal with the scale and offset which exist for some reason
[u10_scale, tmp] = attnc(wind_file, 'u10', 'scale_factor', 0);
[u10_offset, tmp] = attnc(wind_file, 'u10', 'add_offset', 0);
[v10_scale, tmp] = attnc(wind_file, 'v10', 'scale_factor', 0);
[v10_offset, tmp] = attnc(wind_file, 'v10', 'add_offset', 0);
u10 = u10*u10_scale + u10_offset;
v10 = v10*v10_scale + v10_offset;
% Throw away data north of 50S
u10 = u10(:,187:end,:);
v10 = v10(:,187:end,:);
[lon_wind, lat_wind] = meshgrid(longitude, latitude(187:end));
% Convert to wind stress
taux_raw = rhoAir*Cd*u10.^2.*sign(u10);
tauy_raw = rhoAir*Cd*v10.^2.*sign(v10);

% Find lengths of latitude and longitude axes in ROMS
ncload(grid_file, 'lon_rho', 'lat_rho', 'mask_rho', 'mask_u', 'mask_v');
num_lon_rho = size(lon_rho, 2);
num_lon_u = num_lon_rho - 1;
num_lon_v = num_lon_rho;
num_lat_rho = size(lon_rho, 1);
num_lat_u = num_lat_rho;
num_lat_v = num_lat_rho - 1;
% Convert ROMS longitude values from range (0, 360) to range (-180, 180)
lon_rho_fix = lon_rho(:,:);
index = lon_rho_fix > 180.0;
lon_rho_fix(index) = lon_rho_fix(index) - 360.0;

% Read Tamura grid
grid_id = fopen(tamura_grid, 'r');
grid_obs = reshape(fread(grid_id, 721*721*2, 'float32=>double'), 721, 721, 2);
lon_tam = squeeze(grid_obs(:,:,2));
lat_tam = squeeze(grid_obs(:,:,1));

% Set up arrays to store final surface fluxes
if year_avg
    taux = zeros(num_lon_u, num_lat_u, 12);
    tauy = zeros(num_lon_v, num_lat_v, 12);
    shflux = zeros(num_lon_rho, num_lat_rho, 12);
    ssflux = zeros(num_lon_rho, num_lat_rho, 12);    
else
    taux = zeros(num_lon_u, num_lat_u, 12*(max_year-min_year+1));
    tauy = zeros(num_lon_v, num_lat_v, 12*(max_year-min_year+1));
    shflux = zeros(num_lon_rho, num_lat_rho, 12*(max_year-min_year+1));
    ssflux = zeros(num_lon_rho, num_lat_rho, 12*(max_year-min_year+1));    
end

% Loop through each month to interpolate fluxes
for year = min_year:max_year
    for month = 1:12
        
        disp(['Processing month ', num2str(month), ' of ', num2str(year)]);
        
        % Select wind stress data for current month
        index = (year-min_year)*12 + month;  
        taux_raw_curr = squeeze(taux_raw(:,:,index))';
        tauy_raw_curr = squeeze(tauy_raw(:,:,index))';
        % Interpolate to ROMS rho grid
        taux_interp = inpaint_nans(interp2(lon_wind, lat_wind, taux_raw_curr, ...
            lon_rho, lat_rho, 'linear'), 2);
        tauy_interp = inpaint_nans(interp2(lon_wind, lat_wind, tauy_raw_curr, ...
            lon_rho, lat_rho, 'linear'), 2);
        % Interpolate from ROMS rho grid to u and v grids
        taux_interp = (taux_interp(:,1:end-1,:) + taux_interp(:,2:end,:))/2.0;        
        tauy_interp = (tauy_interp(1:end-1,:,:) + tauy_interp(2:end,:,:))/2.0;
        % Mask out land points
        taux_interp(mask_u == 0) = 0.0;
        tauy_interp(mask_v == 0) = 0.0;
        
        % Read heat/salt flux data for this month
        if month == 2 && any(year==leap_years)
            days = 29;
        else
            days = days_per_month(month);
        end        
        data_file = [tamura_base, num2str(year), '_', month_names(month,:), '.data'];
        id = fopen(data_file, 'r');
        % Number of variables in file changed in 2002
        if year < 2002
            data = reshape(fread(id, 721*721*days*9, 'float32=>double'), ...
                721, 721, 9, days);
        else
            data = reshape(fread(id, 721*721*days*6, 'float32=>double'), ...
                721, 721, 6, days);
        end
        % Mask out missing values
        data(data > 1e9) = NaN;  
        % Select the heat and salt flux data, depending on number of
        % variables for that year
        if year < 2002
            shflux_raw = squeeze(data(:,:,4,:));
            ssflux_raw = squeeze(data(:,:,6,:));
        else
            shflux_raw = squeeze(data(:,:,1,:));
            ssflux_raw = squeeze(data(:,:,3,:));
        end
        
        % Calculate monthly average
        shflux_raw = mean(shflux_raw, 3);
        ssflux_raw = mean(ssflux_raw, 3);
        % Interpolate to ROMS rho grid
        shflux_interp = inpaint_nans(griddata(lon_tam, lat_tam, shflux_raw, ...
            lon_rho_fix, lat_rho, 'cubic'), 2);
        ssflux_interp = inpaint_nans(griddata(lon_tam, lat_tam, ssflux_raw, ...
            lon_rho_fix, lat_rho, 'cubic'), 2);
        % Mask out land points
        shflux_interp(mask_rho == 0) = 0.0;
        ssflux_interp(mask_rho == 0) = 0.0;
        
        % Save results for this month
        if year_avg
            taux(:,:,month) = taux(:,:,month) + ...
                taux_interp'/(max_year-min_year+1);
            tauy(:,:,month) = tauy(:,:,month) + ...
                tauy_interp'/(max_year-min_year+1);
            shflux(:,:,month) = shflux(:,:,month) + ...
                shflux_interp'/(max_year-min_year+1);
            ssflux(:,:,month) = ssflux(:,:,month) + ...
                ssflux_interp'/(max_year-min_year+1);
        else
            taux(:,:,index) = taux_interp';
            tauy(:,:,index) = tauy_interp';
            shflux(:,:,index) = shflux_interp';
            ssflux(:,:,index) = ssflux_interp';
        end      
        
    end
    
end

% Convert any NaN values to default NetCDF missing value
taux(isnan(taux)) = -1e34;
tauy(isnan(tauy)) = -1e34;
shflux(isnan(shflux)) = -1e34;
ssflux(isnan(ssflux)) = -1e34;

% Convert from salt flux to freshwater flux
swflux = -ssflux*100.0/ref_salt;

% Calculate time axis in days, centered in the middle of each month
if year_avg
    time = [0.5:11.5]*365.25/12;
else
    time = [0.5:12*(max_year-min_year+1)-0.5]*365.25/12;
end

% Save results in NetCDF file
% This will overwrite any existing file with the same name
id = netcdf.create(output_file, 'clobber');

% Set up dimensions
xi_u_dim = netcdf.defDim(id, 'xi_u', num_lon_u);
eta_u_dim = netcdf.defDim(id, 'eta_u', num_lat_u);
xi_v_dim = netcdf.defDim(id, 'xi_v', num_lon_v);
eta_v_dim = netcdf.defDim(id, 'eta_v', num_lat_v);
xi_rho_dim = netcdf.defDim(id, 'xi_rho', num_lon_rho);
eta_rho_dim = netcdf.defDim(id, 'eta_rho', num_lat_rho);
sms_time_dim = netcdf.defDim(id, 'sms_time', length(time));
shf_time_dim = netcdf.defDim(id, 'shf_time', length(time));
swf_time_dim = netcdf.defDim(id, 'swf_time', length(time));

% Set up variables
sms_time_id = netcdf.defVar(id, 'sms_time', 'double', sms_time_dim);
netcdf.putAtt(id, sms_time_id, 'long_name', 'surface momentum stress time');
netcdf.putAtt(id, sms_time_id, 'units', 'days');
if year_avg
    netcdf.putAtt(id, sms_time_id, 'cycle_length', 365.25);
end
shf_time_id = netcdf.defVar(id, 'shf_time', 'double', shf_time_dim);
netcdf.putAtt(id, shf_time_id, 'long_name', 'surface heat flux time');
netcdf.putAtt(id, shf_time_id, 'units', 'days');
if year_avg
    netcdf.putAtt(id, shf_time_id, 'cycle_length', 365.25);
end
swf_time_id = netcdf.defVar(id, 'swf_time', 'double', swf_time_dim);
netcdf.putAtt(id, swf_time_id, 'long_name', 'surface freshwater flux time');
netcdf.putAtt(id, swf_time_id, 'units', 'days');
if year_avg
    netcdf.putAtt(id, swf_time_id, 'cycle_length', 365.25);
end
sustr_id = netcdf.defVar(id, 'sustr', 'double', [xi_u_dim, eta_u_dim, sms_time_dim]);
netcdf.putAtt(id, sustr_id, 'long_name', 'surface u-momentum stress');
netcdf.putAtt(id, sustr_id, 'units', 'Newton meter-2');
netcdf.putAtt(id, sustr_id, 'missing_value', -1e34);
netcdf.putAtt(id, sustr_id, '_FillValue', -1e34);
svstr_id = netcdf.defVar(id, 'svstr', 'double', [xi_v_dim, eta_v_dim, sms_time_dim]);
netcdf.putAtt(id, svstr_id, 'long_name', 'surface v-momentum stress');
netcdf.putAtt(id, svstr_id, 'units', 'Newton meter-2');
netcdf.putAtt(id, svstr_id, 'missing_value', -1e34);
netcdf.putAtt(id, svstr_id, '_FillValue', -1e34);
shflux_id = netcdf.defVar(id, 'shflux', 'double', [xi_rho_dim, eta_rho_dim, shf_time_dim]);
netcdf.putAtt(id, shflux_id, 'long_name', 'surface net heat flux');
netcdf.putAtt(id, shflux_id, 'units', 'Watts meter-2');
netcdf.putAtt(id, shflux_id, 'missing_value', -1e34);
netcdf.putAtt(id, shflux_id, '_FillValue', -1e34);
swflux_id = netcdf.defVar(id, 'swflux', 'double', [xi_rho_dim, eta_rho_dim, swf_time_dim]);
netcdf.putAtt(id, swflux_id, 'long_name', 'surface freshwater flux (E-P)');
netcdf.putAtt(id, swflux_id, 'units', 'centimeter day-1');
netcdf.putAtt(id, swflux_id, 'missing_value', -1e34);
netcdf.putAtt(id, swflux_id, '_FillValue', -1e34);
netcdf.putAtt(id, netcdf.getConstant('NC_GLOBAL'), 'title', 'CAISOM surface heat, salt fluxes and winds');
netcdf.putAtt(id, netcdf.getConstant('NC_GLOBAL'), 'date', date);
netcdf.putAtt(id, netcdf.getConstant('NC_GLOBAL'), 'grd_file', grid_file);
netcdf.putAtt(id, netcdf.getConstant('NC_GLOBAL'), 'type', 'ROMS forcing file');
netcdf.endDef(id);

% Write variables
netcdf.putVar(id, sms_time_id, time);
netcdf.putVar(id, shf_time_id, time);
netcdf.putVar(id, swf_time_id, time);
netcdf.putVar(id, sustr_id, taux);
netcdf.putVar(id, svstr_id, tauy);
netcdf.putVar(id, shflux_id, shflux);
netcdf.putVar(id, swflux_id, swflux);

netcdf.close(id);

end
