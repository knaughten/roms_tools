% Create ROMS file for LBCs at the northern boundary (50S) of the circumpolar,
% quarter-degree configuration of ROMS.
% Use ECCO2 output for temperature, salinity, velocity, and barotropic
% (vertically averaged) velocity. Set zeta (sea surface height) to zero.

% Code by Kaitlin Alexander, adapted from the code by David Gwyther:
% in /ds/projects/iomp/totten/ana/dgwyther/netcdf/tisom008/lbc, the files
% make_lbc.m, do_interp_ecco2_theta_cube92.m, do_ISOM_lbc_cube92.m, and
% do_ISOM_lbc_nc_cube92.m.

function caisom_lbc

% User-defined parameters

% Path to grid file
grid_file = 'caisom001_OneQuartergrd.nc';
% Path to desired LBC file
output_file = 'caisom001_bry.nc';
% Years of ECCO2 to select (exists from 1992-2012):
min_year = 1992;
max_year = 2012;
% If avg=true, a monthly climatology will be created by averaging each
% month's LBCs over the selected years.
avg = true;

% Parameters for s-coordinates (check grid file for values if unsure)
Vtransform = 2;
Vstretching = 2;
Tcline = 20;
theta_s = 5;
theta_b = 0.4;
hc = 20;
% Number of vertical levels
N = 31;
% Radius of earth in m (for spherical coordinates transformations)
r = 6.371e6;

% End user-defined parameters


% Load some ROMS scripts
addpath(genpath('/ds/projects/iomp/matlab_scripts'));

% Paths to ECCO2 output
ecco2_dir = '/ds/projects/iomp/obs/ECCO2/cube92_real/';
theta_base = [ecco2_dir, 'THETA_monthly.nc/THETA.1440x720x50.'];
salt_base = [ecco2_dir, 'SALT_monthly.nc/SALT.1440x720x50.'];
u_base = [ecco2_dir, 'UVEL_monthly.nc/UVEL.1440x720x50.'];
v_base = [ecco2_dir, 'VVEL_monthly.nc/VVEL.1440x720x50.'];
sample_tail = '199201.nc';
% Read the ECCO2 grid
ncload([theta_base, sample_tail], 'LATITUDE_T', 'LONGITUDE_T', 'DEPTH_T');
% Select indices of ECCO2 horizontal grid for circumpolar domain to 50S
min_i = 1;
max_i = 1440;
min_j = 1;
max_j = 161;

% Calculate horizontal cell thickness along the northern boundary
% of the ROMS v-grid
ncload(grid_file, 'lon_rho', 'lat_rho', 'lon_u', 'lat_u', 'lon_v', 'lat_v');
dlat = lat_v(end,:) - lat_v(end-1,:);
% dy = r*dlat where dlat is converted to radians
% Copy this into a 2D matrix of longitude vs depth
dy = repmat(pi*r/180.0*dlat, N, 1);
% dlon is constant at this latitude
dlon = 0.25;
lat = lat_v(end,:);
% dx = r*cos(lat)*dlon where lat and dlon are converted to radians
% Copy this into a 2D matrix of longitude vs depth
dx = repmat(pi*r*cos(pi*lat/180.0)/180.0*dlon, N, 1);

% Find the length of the longitude axis in ROMS
num_lon_rho = size(lon_rho, 2); % for rho-grid, same as v-grid
num_lon_u = num_lon_rho - 1; % for u-grid

% Set up arrays to store final LBCs
% Dimensions are longitude x depth x time
if avg
    temp_north = zeros(num_lon_rho, N, 12);
    salt_north = zeros(num_lon_rho, N, 12);
    u_north = zeros(num_lon_u, N, 12);
    v_north = zeros(num_lon_rho, N, 12);
    ubar_north = zeros(num_lon_u, 12);
    vbar_north = zeros(num_lon_rho, 12);
else
    temp_north = zeros(num_lon_rho, N, 12*(max_year-min_year+1));
    salt_north = zeros(num_lon_rho, N, 12*(max_year-min_year+1));
    u_north = zeros(num_lon_u, N, 12*(max_year-min_year+1));
    v_north = zeros(num_lon_rho, N, 12*(max_year-min_year+1));
    ubar_north = zeros(num_lon_u, 12*(max_year-min_year+1));
    vbar_north = zeros(num_lon_rho, 12*(max_year-min_year+1));
end

% Loop through each month; read and interpolate ECCO2 output
for year=min_year:max_year
    for month=1:12   
        
        disp(['Processing month ', num2str(month), ' of ', num2str(year)]);
        % Complete path to ECCO2 files
        if month < 10
            tail = [num2str(year), '0', num2str(month), '.nc'];
        else
            tail = [num2str(year), num2str(month), '.nc'];
        end
        
        % Set up arrays to store ECCO2 output over selected domain
        theta_raw = nan(length(DEPTH_T), max_j-min_j+1, max_i-min_i+1);
        salt_raw = nan(length(DEPTH_T), max_j-min_j+1, max_i-min_i+1);
        u_raw = nan(length(DEPTH_T), max_j-min_j+1, max_i-min_i+1);
        v_raw = nan(length(DEPTH_T), max_j-min_j+1, max_i-min_i+1);
        
        % Read ECCO2 output for T, S, u, v
        ncload([theta_base, tail], 'THETA');
        theta_raw(:,:,:) = THETA(:, min_j:max_j, min_i:max_i);
        clear THETA;
        ncload([salt_base, tail], 'SALT');
        salt_raw(:,:,:) = SALT(:, min_j:max_j, min_i:max_i);
        clear SALT;
        ncload([u_base, tail], 'UVEL');
        u_raw(:,:,:) = UVEL(:, min_j:max_j, min_i:max_i);
        clear UVEL;
        ncload([v_base, tail], 'VVEL');
        v_raw(:,:,:) = VVEL(:, min_j:max_j, min_i:max_i);
        clear VVEL;       
                            
        % Replace missing values with NaN
        mask = theta_raw < -1e20;
        theta_raw(mask) = NaN;
        salt_raw(mask) = NaN;
        u_raw(mask) = NaN;
        v_raw(mask) = NaN;
        
        % Calculate s-coordinates in ROMS grid
        ncload(grid_file, 'h', 'zice', 'mask_zice', 'mask_rho');
        % Mask out bathymetry (h) and ice shelf draft (zice)
        h = h.*mask_rho;
        zice = zice.*mask_zice;
        % Call script to calculate s-levels
        [tmp, s_r, Cs_r] = scoord_zice(h', zice', lon_rho', lat_rho', Vtransform, ...
            Vstretching, theta_s, theta_b, hc, N, 0, 0, 1, 0);
        % Calculate depth of each s-point in Cartesian coordinates
        wct = h(end,:); % water column thickness (no ice shelves at this boundary)
        % Linearly interpolate wct to u and v grids
        wct_u = (h(end,1:end-1) + h(end,2:end))/2;
        wct_v = (h(end,:) + h(end-1,:))/2;
        % Flip s values so starting from surface, not bottom
        s_flip = fliplr(s_r);
        % Calculate z = s*wct at each point on the boundary
        z_u = zeros(length(s_flip), length(wct_u));
        z_v = zeros(length(s_flip), length(wct_v));
        for k=1:length(s_flip)
            z_u(k,:) = -wct_u*s_flip(k);
            z_v(k,:) = -wct_v*s_flip(k);
        end
        % These are the depths of the centres of each cell; now interpolate
        % to find the depths of the edges
        z_edges_u = zeros(size(z_u, 1)+1, size(z_u, 2));
        z_edges_u(2:end-1,:) = (z_u(1:end-1,:) + z_u(2:end,:))/2;
        z_edges_u(end,:) = z_u(end,:);
        % Calculate thickness dz of each cell
        dz_u = z_edges_u(2:end,:) - z_edges_u(1:end-1,:);
        % Repeat for v-grid
        z_edges_v = zeros(size(z_v, 1)+1, size(z_v, 2));
        z_edges_v(2:end-1,:) = (z_v(1:end-1,:) + z_v(2:end,:))/2;
        z_edges_v(end,:) = z_v(end,:);
        dz_v = z_edges_v(2:end,:) - z_edges_v(1:end-1,:);
        
        % Create 2D arrays showing depth, latitude, and longitude values
        % along northern boundary of ECCO2 grid (longitude vs depth)        
        depth_bry_rho = flipud(repmat(wct, length(Cs_r), 1).*repmat(Cs_r', ...
            1, size(wct, 2)));        
        depth_bry_u = flipud(repmat(wct_u, length(Cs_r), 1).*repmat(Cs_r', ...
            1, size(wct_u, 2)));
        depth_bry_v = flipud(repmat(wct_v, length(Cs_r), 1).*repmat(Cs_r', ...
            1, size(wct_v, 2)));
        lat_bry_rho = repmat(lat_rho(end,:), length(Cs_r), 1);
        lon_bry_rho = repmat(lon_rho(end,:), length(Cs_r), 1);
        lat_bry_u = repmat(lat_u(end,:), length(Cs_r), 1);
        lon_bry_u = repmat(lon_u(end,:), length(Cs_r), 1);
        lat_bry_v = repmat(lat_v(end,:), length(Cs_r), 1);
        lon_bry_v = repmat(lon_v(end,:), length(Cs_r), 1);
        
        % Interpolate T, S, u, v from ECCO2 grid to correct ROMS grid
        [y, z, x] = meshgrid(LATITUDE_T(min_j:max_j), -DEPTH_T, LONGITUDE_T(min_i:max_i));        
        theta_interp = inpaint_nans(interp3(y, z, x, theta_raw(:,:,:), ...
            lat_bry_rho, depth_bry_rho, lon_bry_rho, 'linear'), 2);
        salt_interp = inpaint_nans(interp3(y, z, x, salt_raw(:,:,:), ...
            lat_bry_rho, depth_bry_rho, lon_bry_rho, 'linear'), 2);
        u_interp = inpaint_nans(interp3(y, z, x, u_raw(:,:,:), ...
            lat_bry_u, depth_bry_u, lon_bry_u, 'linear'), 2);
        v_interp = inpaint_nans(interp3(y, z, x, v_raw(:,:,:), ...
            lat_bry_v, depth_bry_v, lon_bry_v, 'linear'), 2);
        
        % Calculate average northward velocity at the northern boundary
        v_int = squeeze(nansum(squeeze(nansum(v_interp.*dx.*dy.*dz_v, 1)), 2));
        volume = squeeze(nansum(squeeze(nansum(dx.*dy.*dz_v, 1)), 2));
        v_avg = v_int/volume;
        % Subtract this from v so LBCs conserve volume
        v_interp = v_interp - v_avg;
        
        % Calculate barotropic u and v
        ubar_interp = squeeze(nansum(u_interp.*dz_u, 1))./wct_u;
        vbar_interp = squeeze(nansum(v_interp.*dz_v, 1))./wct_v;
        ubar_interp(wct_u == 0) = 0.0;
        vbar_interp(wct_v == 0) = 0.0;
        
        % Flip along the depth axes and transpose so dimensions are in the
        % right order for NetCDF
        theta_interp = flipdim(theta_interp, 1)';
        salt_interp = flipdim(salt_interp, 1)';
        u_interp = flipdim(u_interp, 1)';
        v_interp = flipdim(v_interp, 1)';
        ubar_interp = ubar_interp';
        vbar_interp = vbar_interp';

        % Save results for this month
        if avg
            % Divide by number of years and add to existing values
            % so that by the end it will be an average
            temp_north(:,:,month) = temp_north(:,:,month) + ...
                theta_interp/(max_year-min_year+1);
            salt_north(:,:,month) = salt_north(:,:,month) + ...
                salt_interp/(max_year-min_year+1);
            u_north(:,:,month) = u_north(:,:,month) + ...
                u_interp/(max_year-min_year+1);
            v_north(:,:,month) = v_north(:,:,month) + ...
                v_interp/(max_year-min_year+1);
            ubar_north(:,month) = ubar_north(:,month) + ...
                ubar_interp/(max_year-min_year+1);
            vbar_north(:,month) = vbar_north(:,month) + ...
                vbar_interp/(max_year-min_year+1);
        else
            index = (year-min_year)*12 + month;
            temp_north(:,:,index) = theta_interp;
            salt_north(:,:,index) = salt_interp;
            u_north(:,:,index) = u_interp;
            v_north(:,:,index) = v_interp;
            ubar_north(:,index) = ubar_interp;
            vbar_north(:,index) = vbar_interp;
        end
        
    end
end

% Replace NaNs with the default missing value for NetCDF
temp_north(isnan(temp_north)) = -1e34;
salt_north(isnan(salt_north)) = -1e34;
u_north(isnan(u_north)) = -1e34;
v_north(isnan(v_north)) = -1e34;
ubar_north(isnan(ubar_north)) = -1e34;
vbar_north(isnan(vbar_north)) = -1e34;

% Calculate time axis in days, centered in the middle of each month
if avg
    time = [0.5:11.5]*365.25/12;
else
    time = [0.5:12*(max_year-min_year)-0.5]*365.25/12;
end

% Save results in NetCDF file
id = netcdf.create(output_file, 'clobber');

% Set up dimensions
xi_u_dim = netcdf.defDim(id, 'xi_u', num_lon_u);
xi_rho_dim = netcdf.defDim(id, 'xi_rho', num_lon_rho);
s_rho_dim = netcdf.defDim(id, 's_rho', N);
time_dim = netcdf.defDim(id, 'ocean_time', length(time));
one_dim = netcdf.defDim(id, 'one', 1);

% Set up variables
theta_s_id = netcdf.defVar(id, 'theta_s', 'double', one_dim);
netcdf.putAtt(id, theta_s_id, 'long_name', 'S-coordinate surface control parameter');
netcdf.putAtt(id, theta_s_id, 'units', 'nondimensional');
theta_b_id = netcdf.defVar(id, 'theta_b', 'double', one_dim);
netcdf.putAtt(id, theta_b_id, 'long_name', 'S-coordinate bottom control parameter');
netcdf.putAtt(id, theta_b_id, 'units', 'nondimensional');
tcline_id = netcdf.defVar(id, 'Tcline', 'double', one_dim);
netcdf.putAtt(id, tcline_id, 'long_name', 'S-coordinate surface/bottom layer width');
netcdf.putAtt(id, tcline_id, 'units', 'meter');
hc_id = netcdf.defVar(id, 'hc', 'double', one_dim);
netcdf.putAtt(id, hc_id, 'long_name', 'S-coordinate parameter, critical depth');
netcdf.putAtt(id, hc_id, 'units', 'meter');
sc_r_id = netcdf.defVar(id, 'sc_r', 'double', s_rho_dim);
netcdf.putAtt(id, sc_r_id, 'long_name', 'S-coordinate at RHO-points');
netcdf.putAtt(id, sc_r_id, 'units', 'nondimensional');
cs_r_id = netcdf.defVar(id, 'Cs_r', 'double', s_rho_dim);
netcdf.putAtt(id, cs_r_id, 'long_name', 'S-coordinate stretching curves at RHO-points');
netcdf.putAtt(id, cs_r_id, 'units', 'nondimensional');
netcdf.putAtt(id, cs_r_id, 'valid_min', -1);
netcdf.putAtt(id, cs_r_id, 'valid_max', 0);
time_id = netcdf.defVar(id, 'ocean_time', 'double', time_dim);
netcdf.putAtt(id, time_id, 'long_name', 'time for climatology');
netcdf.putAtt(id, time_id, 'units', 'day');
if avg
    netcdf.putAtt(id, time_id, 'cycle_length', 365.25);
end
temp_north_id = netcdf.defVar(id, 'temp_north', 'double', [xi_rho_dim, s_rho_dim, time_dim]);
netcdf.putAtt(id, temp_north_id, 'long_name', 'northern boundary potential temperature');
netcdf.putAtt(id, temp_north_id, 'units', 'Celsius');
netcdf.putAtt(id, temp_north_id, 'missing_value', -1e34);
netcdf.putAtt(id, temp_north_id, '_FillValue', -1e34);
salt_north_id = netcdf.defVar(id, 'salt_north', 'double', [xi_rho_dim, s_rho_dim, time_dim]);
netcdf.putAtt(id, salt_north_id, 'long_name', 'northern boundary salinity');
netcdf.putAtt(id, salt_north_id, 'units', 'PSU');
netcdf.putAtt(id, salt_north_id, 'missing_value', -1e34);
netcdf.putAtt(id, salt_north_id, '_FillValue', -1e34);
u_north_id = netcdf.defVar(id, 'u_north', 'double', [xi_u_dim, s_rho_dim, time_dim]);
netcdf.putAtt(id, u_north_id, 'long_name', 'northern boundary u-monentum component');
netcdf.putAtt(id, u_north_id, 'units', 'meter second-1');
netcdf.putAtt(id, u_north_id, 'missing_value', -1e34);
netcdf.putAtt(id, u_north_id, '_FillValue', -1e34);
v_north_id = netcdf.defVar(id, 'v_north', 'double', [xi_rho_dim, s_rho_dim, time_dim]);
netcdf.putAtt(id, v_north_id, 'long_name', 'northern boundary v-momentum component');
netcdf.putAtt(id, v_north_id, 'units', 'meter second-1');
netcdf.putAtt(id, v_north_id, 'missing_value', -1e34);
netcdf.putAtt(id, v_north_id, '_FillValue', -1e34);
ubar_north_id = netcdf.defVar(id, 'ubar_north', 'double', [xi_u_dim, time_dim]);
netcdf.putAtt(id, ubar_north_id, 'long_name', 'northern boundary vertically integrated u-momentum component');
netcdf.putAtt(id, ubar_north_id, 'units', 'meter second-1');
netcdf.putAtt(id, ubar_north_id, 'missing_value', -1e34);
netcdf.putAtt(id, ubar_north_id, '_FillValue', -1e34);
vbar_north_id = netcdf.defVar(id, 'vbar_north', 'double', [xi_rho_dim, time_dim]);
netcdf.putAtt(id, vbar_north_id, 'long_name', 'northern boundary vertically integrated v-momentum component');
netcdf.putAtt(id, vbar_north_id, 'units', 'meter second-1');
netcdf.putAtt(id, vbar_north_id, 'missing_value', -1e34);
netcdf.putAtt(id, vbar_north_id, '_FillValue', -1e34);
zeta_north_id = netcdf.defVar(id, 'zeta_north', 'double', [xi_rho_dim, time_dim]);
netcdf.putAtt(id, zeta_north_id, 'long_name', 'northern boundary sea surface height');
netcdf.putAtt(id, zeta_north_id, 'units', 'meter');
netcdf.putAtt(id, netcdf.getConstant('NC_GLOBAL'), 'title', 'Lateral Boundaries');
netcdf.putAtt(id, netcdf.getConstant('NC_GLOBAL'), 'date', date);
netcdf.putAtt(id, netcdf.getConstant('NC_GLOBAL'), 'clim_file', output_file);
netcdf.putAtt(id, netcdf.getConstant('NC_GLOBAL'), 'grd_file', grid_file);
netcdf.putAtt(id, netcdf.getConstant('NC_GLOBAL'), 'type', 'BOUNDARY file');
netcdf.putAtt(id, netcdf.getConstant('NC_GLOBAL'), 'title', 'ROMS');
netcdf.endDef(id);

% Write variables
netcdf.putVar(id, theta_s_id, theta_s);
netcdf.putVar(id, theta_b_id, theta_b);
netcdf.putVar(id, tcline_id, Tcline);
netcdf.putVar(id, hc_id, hc);
netcdf.putvar(id, sc_r_id, s_r);
netcdf.putVar(id, cs_r_id, Cs_r);
netcdf.putVar(id, time_id, time);
netcdf.putVar(id, temp_north_id, temp_north);
netcdf.putVar(id, salt_north_id, salt_north);
netcdf.putVar(id, u_north_id, u_north);
netcdf.putVar(id, v_north_id, v_north);
netcdf.putVar(id, ubar_north_id, ubar_north);
netcdf.putVar(id, vbar_north_id, vbar_north);
netcdf.putVar(id, zeta_north_id, zeros(num_lon_rho, length(time)));
        
netcdf.close(id);
        
        
        
        
        











        
        


