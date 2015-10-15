% Create ROMS file for surface atmospheric conditions (pressure, temperature, 
% relative humidity, cloud fraction, winds, rainfall, snowfall) to be used
% in the computation of bulk fluxes for the circumpolar quarter-degree 
% configuration of ROMS. All data is from ERA-Interim monthly averages.
% Code by Kaitlin Alexander

function caisom_bulkfluxes

% User-defined parameters

% Path to grid file
grid_file = 'caisom001_OneQuartergrd.nc';
% Desired path to output file
output_file = 'caisom001_surface.nc';
% Years of data to select (exists from 1992-2012)
min_year = 1992;
max_year = 2012;
% If year_avg=true, a monthly climatology will be created by averaging
% each month's fluxes over the selected years.
year_avg = true;

% End user-defined parameters


% Path to file containing monthly averages of all data except precipitation and snow
data_file = '/ds/projects/iomp/obs/ERA_Interim/ERA_Interim_atmosphere_1992_2012_monthly.nc';
% Path to file containing monthly averages of total precipitation
precip_file = '/ds/projects/iomp/obs/ERA_Interim/ERA_Interim_precip_1992_2012_monthly.nc';
% Path to file containing monthly averages of snowfall
snow_file = '/ds/projects/iomp/obs/ERA_Interim/ERA_Interim_snow_1992_2012_monthly.nc';

Lv = 2.5e6; % latent heat of vapourisation
Rv = 461.5; % gas constant for moist air

% Load some ROMS scripts
addpath(genpath('/ds/projects/iomp/matlab_scripts'));

% Read first batch of data
ncload(data_file, 'sp', 't2m', 'd2m', 'tcc', 'u10', 'v10', 'longitude', 'latitude');

% Surface pressure
% Re-order axes
sp = permute(sp, [3 2 1]);
% Mask out missing values
sp(sp == -32767) = NaN;
% Apply scale and offset
[sp_scale, tmp] = attnc(data_file, 'sp', 'scale_factor', 0);
[sp_offset, tmp] = attnc(data_file, 'sp', 'add_offset', 0);
sp = sp*sp_scale + sp_offset;
% Throw away data north of 50S
sp = sp(:,187:end,:);

% Temperature at 2 metres
t2m = permute(t2m, [3 2 1]);
t2m(t2m == -32767) = NaN;
[t2m_scale, tmp] = attnc(data_file, 't2m', 'scale_factor', 0);
[t2m_offset, tmp] = attnc(data_file, 't2m', 'add_offset', 0);
t2m = t2m*t2m_scale + t2m_offset;
t2m = t2m(:,187:end,:);

% Dew point at 2 metres
d2m = permute(d2m, [3 2 1]);
d2m(d2m == -32767) = NaN;
[d2m_scale, tmp] = attnc(data_file, 'd2m', 'scale_factor', 0);
[d2m_offset, tmp] = attnc(data_file, 'd2m', 'add_offset', 0);
d2m = d2m*d2m_scale + d2m_offset;
d2m = d2m(:,187:end,:);

% Total cloud cover
tcc = permute(tcc, [3 2 1]);
tcc(tcc == -32767) = NaN;
[tcc_scale, tmp] = attnc(data_file, 'tcc', 'scale_factor', 0);
[tcc_offset, tmp] = attnc(data_file, 'tcc', 'add_offset', 0);
tcc = tcc*tcc_scale + tcc_offset;
tcc = tcc(:,187:end,:);

% U-winds at 10 metres
u10 = permute(u10, [3 2 1]);
u10(u10 == -32767) = NaN;
[u10_scale, tmp] = attnc(data_file, 'u10', 'scale_factor', 0);
[u10_offset, tmp] = attnc(data_file, 'u10', 'add_offset', 0);
u10 = u10*u10_scale + u10_offset;
u10 = u10(:,187:end,:);

% V-winds at 10 metres
v10 = permute(v10, [3 2 1]);
v10(v10 == -32767) = NaN;
[v10_scale, tmp] = attnc(data_file, 'v10', 'scale_factor', 0);
[v10_offset, tmp] = attnc(data_file, 'v10', 'add_offset', 0);
v10 = v10*v10_scale + v10_offset;
v10 = v10(:,187:end,:);

% Calculate relative humidity from temperature and dew point
rh = exp(Lv/Rv*((t2m).^(-1) - (d2m).^(-1)));

% Set up ERA-Interim grid
[lon_era, lat_era] = meshgrid(longitude, latitude(187:end));

% Read total precipitation
ncload(precip_file, 'tp');
tp = permute(tp, [3 2 1]);
tp(tp == -32767) = NaN;
[tp_scale, tmp] = attnc(precip_file, 'tp', 'scale_factor', 0);
[tp_offset, tmp] = attnc(precip_file, 'tp', 'add_offset', 0);
tp = tp*tp_scale + tp_offset;
tp = tp(:,187:end,:);
% Convert from m per 12 h to kg/m^2/s
tp = tp*1e3/(12*60*60);

% Read snowfall
ncload(snow_file, 'sf');
sf = permute(sf, [3 2 1]);
sf(sf == -32767) = NaN;
[sf_scale, tmp] = attnc(snow_file, 'sf', 'scale_factor', 0);
[sf_offset, tmp] = attnc(snow_file, 'sf', 'add_offset', 0);
sf = sf*sf_scale + sf_offset;
sf = sf(:,187:end,:);
% Convert from m (water equivalent) per 12 h to kg/m^2/s
sf = sf*1e3/(12*60*60);

% Read ROMS grid
ncload(grid_file, 'lon_rho', 'lat_rho', 'mask_rho', 'lon_u', 'lat_u', ...
    'mask_u', 'lon_v', 'lat_v', 'mask_v');
% Find lengths of latitude and longitude axes
num_lon_rho = size(lon_rho, 2);
num_lon_u = size(lon_u, 2);
num_lon_v = size(lon_v, 2);
num_lat_rho = size(lon_rho, 1);
num_lat_u = size(lon_u, 1);
num_lat_v = size(lon_v, 1);

% Set up arrays to store final surface fluxes
if year_avg
    pair = zeros(num_lon_rho, num_lat_rho, 12);
    tair = zeros(num_lon_rho, num_lat_rho, 12);
    qair = zeros(num_lon_rho, num_lat_rho, 12);
    cloud = zeros(num_lon_rho, num_lat_rho, 12);
    uwind = zeros(num_lon_u, num_lat_u, 12);
    vwind = zeros(num_lon_v, num_lat_v, 12);
    rain = zeros(num_lon_rho, num_lat_rho, 12);
    snow = zeros(num_lon_rho, num_lat_rho, 12);
else
    pair = zeros(num_lon_rho, num_lat_rho, 12*(max_year-min_year+1));
    tair = zeros(num_lon_rho, num_lat_rho, 12*(max_year-min_year+1));
    qair = zeros(num_lon_rho, num_lat_rho, 12*(max_year-min_year+1));
    cloud = zeros(num_lon_rho, num_lat_rho, 12*(max_year-min_year+1));
    uwind = zeros(num_lon_u, num_lat_u, 12*(max_year-min_year+1));
    vwind = zeros(num_lon_v, num_lat_v, 12*(max_year-min_year+1));
    rain = zeros(num_lon_rho, num_lat_rho, 12*(max_year-min_year+1));   
    snow = zeros(num_lon-rho, num_lat_rho, 12*(max_year-min_year+1));
end

% Loop through each month to interpolate data
for year = min_year:max_year
    for month = 1:12
        
        disp(['Processing month ', num2str(month), ' of ', num2str(year)]);
        
        % Select data for this month
        index = (year-min_year)*12 + month;
        pair_raw_curr = squeeze(sp(:,:,index))';
        % Convert temperature from K to C
        tair_raw_curr = squeeze(t2m(:,:,index))' - 273.15;
        qair_raw_curr = squeeze(rh(:,:,index))';
        cloud_raw_curr = squeeze(tcc(:,:,index))';
        uwind_raw_curr = squeeze(u10(:,:,index))';
        vwind_raw_curr = squeeze(v10(:,:,index))';
        precip_raw_curr = squeeze(tp(:,:,index))';
        snow_raw_curr = squeeze(sf(:,:,index))';        
        
        % Interpolate to ROMS grid
        pair_interp = inpaint_nans(interp2(lon_era, lat_era, pair_raw_curr, ...
            lon_rho, lat_rho, 'linear'), 2);
        tair_interp = inpaint_nans(interp2(lon_era, lat_era, tair_raw_curr, ...
            lon_rho, lat_rho, 'linear'), 2);
        qair_interp = inpaint_nans(interp2(lon_era, lat_era, qair_raw_curr, ...
            lon_rho, lat_rho, 'linear'), 2);
        cloud_interp = inpaint_nans(interp2(lon_era, lat_era, cloud_raw_curr, ...
            lon_rho, lat_rho, 'linear'), 2);
        uwind_interp = inpaint_nans(interp2(lon_era, lat_era, uwind_raw_curr, ...
            lon_u, lat_u, 'linear'), 2);
        vwind_interp = inpaint_nans(interp2(lon_era, lat_era, vwind_raw_curr, ...
            lon_v, lat_v, 'linear'), 2);
        rain_interp = inpaint_nans(interp2(lon_era, lat_era, precip_raw_curr, ...
            lon_rho, lat_rho, 'linear'), 2);
        snow_interp = inpaint_nans(interp2(lon_era, lat_era, snow_raw_curr, ...
            lon_rho, lat_rho, 'linear'), 2);
        
        % mask out land points
        pair_interp(mask_rho == 0) = NaN;
        tair_interp(mask_rho == 0) = NaN;
        cloud_interp(mask_rho == 0) = NaN;
        uwind_interp(mask_u == 0) = NaN;
        vwind_interp(mask_v == 0) = NaN;
        qair_interp(mask_rho == 0) = NaN;
        rain_interp(mask_rho == 0) = NaN;
        snow_interp(mask_rho == 0) = NaN;
        
        % Save results for this month
        if year_avg
            pair(:,:,month) = pair(:,:,month) + ...
                pair_interp'/(max_year-min_year+1);
            tair(:,:,month) = tair(:,:,month) + ...
                tair_interp'/(max_year-min_year+1);
            qair(:,:,month) = qair(:,:,month) + ...
                qair_interp'/(max_year-min_year+1);
            cloud(:,:,month) = cloud(:,:,month) + ...
                cloud_interp'/(max_year-min_year+1);
            uwind(:,:,month) = uwind(:,:,month) + ...
                uwind_interp'/(max_year-min_year+1);
            vwind(:,:,month) = vwind(:,:,month) + ...
                vwind_interp'/(max_year-min_year+1);           
            rain(:,:,month) = rain(:,:,month) + ...
                rain_interp'/(max_year-min_year+1);
            snow(:,:,month) = snow(:,:,month) + ...
                snow_interp'/(max_year-min_year+1);
        else
            pair(:,:,index) = pair_interp';
            tair(:,:,index) = tair_interp';
            qair(:,:,index) = qair_interp';
            cloud(:,:,index) = cloud_interp';
            uwind(:,:,index) = uwind_interp';
            vwind(:,:,index) = vwind_interp';            
            rain(:,:,index) = rain_interp';
            snow(:,:,index) = snow_interp';
        end
        
    end
    
end

% Convert any NaN values to default NetCDF missing value
pair(isnan(pair)) = -1e34;
tair(isnan(tair)) = -1e34;
qair(isnan(qair)) = -1e34;
cloud(isnan(cloud)) = -1e34;
uwind(isnan(uwind)) = -1e34;
vwind(isnan(vwind)) = -1e34;
rain(isnan(rain)) = -1e34;
snow(isnan(snow)) = -1e34;

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
pair_time_dim = netcdf.defDim(id, 'pair_time', length(time));
tair_time_dim = netcdf.defDim(id, 'tair_time', length(time));
qair_time_dim = netcdf.defDim(id, 'qair_time', length(time));
srf_time_dim = netcdf.defDim(id, 'srf_time', length(time));
wind_time_dim = netcdf.defDim(id, 'wind_time', length(time));
rain_time_dim = netcdf.defDim(id, 'rain_time', length(time));
snow_time_dim = netcdf.defDim(id, 'snow_time', length(time));

% Set up variables
pair_time_id = netcdf.defVar(id, 'pair_time', 'double', pair_time_dim);
netcdf.putAtt(id, pair_time_id, 'long_name', 'surface air pressure time');
netcdf.putAtt(id, pair_time_id, 'units', 'days');
if year_avg
    netcdf.putAtt(id, pair_time_id, 'cycle_length', 365.25);
end
tair_time_id = netcdf.defVar(id, 'tair_time', 'double', tair_time_dim);
netcdf.putAtt(id, tair_time_id, 'long_name', 'surface air temperature time');
netcdf.putAtt(id, tair_time_id, 'units', 'days');
if year_avg
    netcdf.putAtt(id, tair_time_id, 'cycle_length', 365.25);
end
qair_time_id = netcdf.defVar(id, 'qair_time', 'double', qair_time_dim);
netcdf.putAtt(id, qair_time_id, 'long_name', 'surface relative humidity time');
netcdf.putAtt(id, qair_time_id, 'units', 'days');
if year_avg
    netcdf.putAtt(id, qair_time_id, 'cycle_length', 365.25);
end
srf_time_id = netcdf.defVar(id, 'srf_time', 'double', srf_time_dim);
netcdf.putAtt(id, srf_time_id, 'long_name', 'surface time');
netcdf.putAtt(id, srf_time_id, 'units', 'days');
if year_avg
    netcdf.putAtt(id, srf_time_id, 'cycle_length', 365.25);
end
wind_time_id = netcdf.defVar(id, 'wind_time', 'double', wind_time_dim);
netcdf.putAtt(id, wind_time_id, 'long_name', 'surface wind time');
netcdf.putAtt(id, wind_time_id, 'units', 'days');
if year_avg
    netcdf.putAtt(id, wind_time_id, 'cycle_length', 365.25);
end
rain_time_id = netcdf.defVar(id, 'rain_time', 'double', rain_time_dim);
netcdf.putAtt(id, rain_time_id, 'long_name', 'rain time');
netcdf.putAtt(id, rain_time_id, 'units', 'days');
if year_avg
    netcdf.putAtt(id, rain_time_id, 'cycle_length', 365.25);
end
snow_time_id = netcdf.defVar(id, 'snow_time', 'double', snow_time_dim);
netcdf.putAtt(id, snow_time_id, 'long_name', 'snow time');
netcdf.putAtt(id, snow_time_id, 'units', 'days');
if year_avg
    netcdf.putAtt(id, snow_time_id, 'cycle_length', 365.25);
end
pair_id = netcdf.defVar(id, 'Pair', 'double', [xi_rho_dim, eta_rho_dim, pair_time_dim]);
netcdf.putAtt(id, pair_id, 'long_name', 'surface air pressure');
netcdf.putAtt(id, pair_id, 'units', 'Pascal');
netcdf.putAtt(id, pair_id, 'missing_value', -1e34);
netcdf.putAtt(id, pair_id, '_FillValue', -1e34);
tair_id = netcdf.defVar(id, 'Tair', 'double', [xi_rho_dim, eta_rho_dim, tair_time_dim]);
netcdf.putAtt(id, tair_id, 'long_name', 'surface air temperature');
netcdf.putAtt(id, tair_id, 'units', 'Celsius');
netcdf.putAtt(id, tair_id, 'missing_value', -1e34);
netcdf.putAtt(id, tair_id, '_FillValue', -1e34);
qair_id = netcdf.defVar(id, 'Qair', 'double', [xi_rho_dim, eta_rho_dim, qair_time_dim]);
netcdf.putAtt(id, qair_id, 'long_name', 'surface air relative humidity');
netcdf.putAtt(id, qair_id, 'units', 'kg/kg');
netcdf.putAtt(id, qair_id, 'missing_value', -1e34);
netcdf.putAtt(id, qair_id, '_FillValue', -1e34);
cloud_id = netcdf.defVar(id, 'cloud', 'double', [xi_rho_dim, eta_rho_dim, srf_time_dim]);
netcdf.putAtt(id, cloud_id, 'long_name', 'cloud fraction');
netcdf.putAtt(id, cloud_id, 'units', 'nondimensional');
netcdf.putAtt(id, cloud_id, 'missing_value', -1e34);
netcdf.putAtt(id, cloud_id, '_FillValue', -1e34);
uwind_id = netcdf.defVar(id, 'Uwind', 'double', [xi_u_dim, eta_u_dim, wind_time_dim]);
netcdf.putAtt(id, uwind_id, 'long_name', 'surface u-wind component');
netcdf.putAtt(id, uwind_id, 'units', 'meter second-1');
netcdf.putAtt(id, uwind_id, 'missing_value', -1e34);
netcdf.putAtt(id, uwind_id, '_FillValue', -1e34);
vwind_id = netcdf.defVar(id, 'Vwind', 'double', [xi_v_dim, eta_v_dim, wind_time_dim]);
netcdf.putAtt(id, vwind_id, 'long_name', 'surface v-wind component');
netcdf.putAtt(id, vwind_id, 'units', 'meter second-1');
netcdf.putAtt(id, vwind_id, 'missing_value', -1e34);
netcdf.putAtt(id, vwind_id, '_FillValue', -1e34);
rain_id = netcdf.defVar(id, 'rain', 'double', [xi_rho_dim, eta_rho_dim, rain_time_dim]);
netcdf.putAtt(id, rain_id, 'long_name', 'rain fall rate');
netcdf.putAtt(id, rain_id, 'units', 'kilogram meter-2 second-1');
netcdf.putAtt(id, rain_id, 'missing_value', -1e34);
netcdf.putAtt(id, rain_id, '_FillValue', -1e34);
snow_id = netcdf.defVar(id, 'snow', 'double', [xi_rho_dim, eta_rho_dim, snow_time_dim]);
netcdf.putAtt(id, snow_id, 'long_name', 'snow fall rate');
netcdf.putAtt(id, snow_id, 'units', 'kilogram meter-2 second-1');
netcdf.putAtt(id, snow_id, 'missing_value', -1e34);
netcdf.putAtt(id, snow_id, '_FillValue', -1e34);
netcdf.putAtt(id, netcdf.getConstant('NC_GLOBAL'), 'title', 'CAISOM surface atmospheric conditions for bulk fluxes');
netcdf.putAtt(id, netcdf.getConstant('NC_GLOBAL'), 'date', date);
netcdf.putAtt(id, netcdf.getConstant('NC_GLOBAL'), 'grd_file', grid_file);
netcdf.putAtt(id, netcdf.getConstant('NC_GLOBAL'), 'type', 'ROMS forcing file');
netcdf.endDef(id);

% Write variables
netcdf.putVar(id, pair_time_id, time);
netcdf.putVar(id, tair_time_id, time);
netcdf.putVar(id, qair_time_id, time);
netcdf.putVar(id, srf_time_id, time);
netcdf.putVar(id, wind_time_id, time);
netcdf.putVar(id, rain_time_id, time);
netcdf.putVar(id, snow_time_id, time);
netcdf.putVar(id, pair_id, pair);
netcdf.putVar(id, tair_id, tair);
netcdf.putVar(id, qair_id, qair);
netcdf.putVar(id, cloud_id, cloud);
netcdf.putVar(id, uwind_id, uwind);
netcdf.putVar(id, vwind_id, vwind);
netcdf.putVar(id, rain_id, rain);
netcdf.putVar(id, snow_id, snow);

netcdf.close(id);

end
