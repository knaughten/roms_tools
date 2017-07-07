% Read ROMS grid
id = netcdf.open('circ30S_quarterdegree.nc','NOWRITE');
lat_id = netcdf.inqVarID(id, 'lat_rho');
lat = netcdf.getVar(id, lat_id);
lon_id = netcdf.inqVarID(id, 'lon_rho');
lon = netcdf.getVar(id, lon_id);
h_id = netcdf.inqVarID(id, 'h');
h = netcdf.getVar(id, h_id);
zice_id = netcdf.inqVarID(id, 'zice');
zice = netcdf.getVar(id, zice_id);
netcdf.close(id);

% Make sure longitude is in the range (-180, 180)
index = lon > 180;
lon(index) = lon(index)-360;
num_lon = size(lon,1);
num_lat = size(lon,2);

% Get water column thickness
wct = h + zice;

% Call TMD for the amplitudes of tidal transports U and V in m^2/s
% for each tidal constituent
[amp_U,~,~,~] = tmd_extract_HC('DATA/Model_CATS2008a_opt', lat, lon, 'U');
[amp_V,~,~,~] = tmd_extract_HC('DATA/Model_CATS2008a_opt', lat, lon, 'V');
% Get directionless amplitude with sqrt(amp_U^2 + amp_V^2)
% Divide by sqrt(2) for RMS of sinusoid
rms_components = sqrt(amp_U.^2 + amp_V.^2)/sqrt(2);
% Sum the tidal consitutents and divide by water column thickness for
% result in m/s
rms_tide_vel_tmp = squeeze(sum(rms_components, 1))./wct;
% Add a time axis of size 1
rms_tide_vel = zeros([1, size(rms_tide_vel_tmp,1), size(rms_tide_vel_tmp,2)]);
rms_tide_vel(1,:,:) = rms_tide_vel_tmp;

% Remove NaNs
index = isnan(rms_tide_vel);
rms_tide_vel(index) = 0.0;
% Set to zero outside ice shelf cavities
index = zice==0;
rms_tide_vel(index) = 0;

% Write to file
id = netcdf.create('rms_tides.nc', 'CLOBBER');
time_dim = netcdf.defDim(id, 'rms_time', netcdf.getConstant('NC_UNLIMITED'));
eta_dim = netcdf.defDim(id, 'eta_rho', num_lat);
xi_dim = netcdf.defDim(id, 'xi_rho', num_lon);
time_id = netcdf.defVar(id, 'rms_time', 'double', time_dim);
netcdf.putAtt(id, time_id, 'calendar', 'none');
tide_id = netcdf.defVar(id, 'tide_RMSvel', 'NC_DOUBLE', [xi_dim, eta_dim, time_dim]);
netcdf.putAtt(id, tide_id, 'units', 'm/s');
netcdf.endDef(id);
netcdf.putVar(id, time_id, 0, 1, [0.]);
start = [0, 0, 0];
count = [num_lon, num_lat, 1];
netcdf.putVar(id, tide_id, start, count, rms_tide_vel);
netcdf.close(id);
