
function OI(argo_float_data, wod_ctd_data, FM_fname, AdditionalFields, output_fname,  lon_grid, lat_grid, option);

%DERNIERE VERSION 9 FEVRIER

% climatology from float data - maybe do this basin wide at a time ... ...
% especially if this data set get's combined with CTD data.

% mex make_dist_matrix_mex_exp_scal.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
% gcc 4-2 or newer
% check mex opts file!!!


% -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
% -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
% VARIABLES YOU MAY TUNE
% -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
% -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

[x,y]=meshgrid(lon_grid, lat_grid);	 % LON, LAT

if length(option)~=30
	disp('Option input has wrong dimensions')
	return
end
eps100			=option(1); 
eps999 			=option(2); 
create_new 		=option(3); 
compute_interior  	=option(4); 
drawit 			=option(5); 
useargo 		=option(6); 
upperocean 		=option(7); 
upperocean_limit 	=option(8); 
wfkt_multiplier  	=option(9); 
max_profs_inmem 	=option(10); 
min_CTD_profs		=option(11); 
max_num_profs 		=option(12); 
min_num_profs 		=option(13); 
num_of_rand_profiles 	=option(14); 
min_allowed_weight 	=option(15); 
rnd_within 		=option(16); 
iqr_mult 		=option(17); 
continue_at 		=option(18); 
continue_upto		=option(19); 
continue_lon_at		=option(20); 
continue_lon_upto	=option(21);  
OMP_NUM_THREADS_DAY 	=option(22); 
OMP_NUM_THREADS_NIGHT 	=option(23); 
signal_to_noise 	=option(24); 
horiz_scale_set		=option(25); 
time_scale 		=option(26); 
tim 			=option(27); 
centeryear		=option(28); 
dec_scale		=option(29); 
dec_noise_offset 	=option(30); 



%%
% -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
% -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
% VARIABLE YOU MAY   = NOT =   TUNE
% -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
% -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%% preliminary constants and initialization
%expweight = 1.50; % FIXED, DO NOT CHANGE % use something that can be expressed by x * sqrt(x)
% or similar, since x.^y is VERY time consuming for rem(y,0.5) ~= 0
% (DOUBLE the time for other exponents!!! (25% more time for the whole
% project!!! expweight = 0.5,1,1.5,2 ... the ones to choose from but to implement manually.)

expweight = 2; %!!!!

numsortedprofs = max_num_profs-num_of_rand_profiles;

% Matrix inverse options:
opts.SYM = true;  % COVARIANCE MATRIX ALWAYS symetrical
opts.POSDEF = true; % covariance matrix always positive defined

[m,n]=size(x)  % grid to compute

horiz_scale = (horiz_scale_set / 110)^2;


% find max dist for nearest neighbor
maxdist_finegrid = 0.15 * sqrt(2);
maxdist_coarsegrid = 0.51 * sqrt(2);


%% INITIALIZE ACCESS TO ISOPYCNAL CTD AND FLOAT DATA:    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Argo_file.nc
if useargo
    ncid_AFD = netcdf.open(argo_float_data,'NOWRITE');
    disp('loading Argo float profile locations')
    AFD_id_LON = netcdf.inqVarID(ncid_AFD,'lon');
    AFD_id_LAT = netcdf.inqVarID(ncid_AFD,'lat');
    AFD_id_DYR = netcdf.inqVarID(ncid_AFD,'dyr');
    AFD_id_SIGMA = netcdf.inqVarID(ncid_AFD,'sigma');
    AFD_id_TEMP = netcdf.inqVarID(ncid_AFD,'temp');
    AFD_id_SAL = netcdf.inqVarID(ncid_AFD,'sal');
    AFD_id_PRES = netcdf.inqVarID(ncid_AFD,'pres');
    AFD_id_ML_PRES = netcdf.inqVarID(ncid_AFD,'mld_pres');
    AFD_id_ML_DENS = netcdf.inqVarID(ncid_AFD,'mld_dens');
    AFD_id_ML_SALT = netcdf.inqVarID(ncid_AFD,'mld_salt');
    AFD_id_ML_TEMP = netcdf.inqVarID(ncid_AFD,'mld_temp');
    for ifield=1:length(AdditionalFields);
	field=AdditionalFields{ifield};
	eval(['AFD_id_ML_' field ' = netcdf.inqVarID(ncid_AFD,''' field ''');']);%=============================MODIF
    end
    AFD_LAT = netcdf.getVar(ncid_AFD,AFD_id_LAT);
    AFD_LON = netcdf.getVar(ncid_AFD,AFD_id_LON);
    AFD_DYR = netcdf.getVar(ncid_AFD,AFD_id_DYR);
    sg0_pleth = double(netcdf.getVar(ncid_AFD,AFD_id_SIGMA));
    
    % create index vector to find the right profiles in netcdf file.
    AFD_IDVEC = (0:numel(AFD_LAT)-1)';% negative values for CTD file, positive for Argo
    disp(' finished reading Argo float profile locations')
    
else
    AFD_LAT = 190;  %#ok<UNRCH> % arbitary 'off globe' location for a single non-existant profile
    AFD_LON = 0;
    AFD_DYR = -100;
    sg0_pleth = [];
    % create index vector to find the right profiles in netcdf file.
    AFD_IDVEC = 100000;
end



ncid_CTD = netcdf.open(wod_ctd_data,'NOWRITE');	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CTD_file.nc
CTD_id_LON = netcdf.inqVarID(ncid_CTD,'lon');
CTD_id_LAT = netcdf.inqVarID(ncid_CTD,'lat');
CTD_id_DYR = netcdf.inqVarID(ncid_CTD,'dyr');
CTD_id_SIGMA = netcdf.inqVarID(ncid_CTD,'sigma');
CTD_id_TEMP = netcdf.inqVarID(ncid_CTD,'temp');
CTD_id_SAL = netcdf.inqVarID(ncid_CTD,'sal');
CTD_id_PRES = netcdf.inqVarID(ncid_CTD,'pres');
CTD_id_ML_PRES = netcdf.inqVarID(ncid_CTD,'mld_pres');
CTD_id_ML_DENS = netcdf.inqVarID(ncid_CTD,'mld_dens');
CTD_id_ML_SALT = netcdf.inqVarID(ncid_CTD,'mld_salt');
CTD_id_ML_TEMP = netcdf.inqVarID(ncid_CTD,'mld_temp');
for ifield=1:length(AdditionalFields);
	field=AdditionalFields{ifield};
	eval(['CTD_id_ML_' field ' = netcdf.inqVarID(ncid_CTD,''' field ''');']);%=============================MODIF
end

CTD_LAT = netcdf.getVar(ncid_CTD,CTD_id_LAT);
CTD_LON = netcdf.getVar(ncid_CTD,CTD_id_LON);
CTD_DYR = netcdf.getVar(ncid_CTD,CTD_id_DYR);

if isempty(sg0_pleth)
    sg0_pleth = netcdf.getVar(ncid_CTD,CTD_id_SIGMA);
end


% create index vector to find the right profiles in netcdf file.
CTD_IDVEC = (-numel(CTD_LAT):-1)'; % negative values for CTD file, positive for Argo
CTD_offset =numel(CTD_LAT);
disp(' finished reading WOC CTD profile locations')

%%

disp(' pre allocating matrices')


% pick out levels on which to map things
sg_lev=1:1:length(sg0_pleth); % << full length as set in interpolated profiles.
% or subset: ONLY CONTINIOUS SUBSETS POSSIBLE SO FAR!!!
%sg_lev = find(sg0_pleth>24.0,1,'first'):find(sg0_pleth<28.9,1,'last');
%sg_lev = find(sg0_pleth>24.3,1,'first'):find(sg0_pleth<26.65,1,'last');

%sg_l_orig=numel(sg_lev);

% just use those levels - remove others
sg0_pleth=sg0_pleth(sg_lev);


sg_l=numel(sg0_pleth); % number of sigma levels


% initialize time - sigma matrices
E_Z = NaN(max_num_profs);

dummy_T_sig = NaN([tim sg_l],'double');

ones_sig = ones(1,sg_l);

if compute_interior
    disp('initialize time')
    sg_sa_gr = dummy_T_sig;
    sg_pr_gr = dummy_T_sig;
    sg_th_gr = dummy_T_sig;
    sg_yr_gr = dummy_T_sig;
    
    sg_sa_gr_wmean = dummy_T_sig;
    sg_pr_gr_wmean = dummy_T_sig;
    sg_th_gr_wmean = dummy_T_sig;
    % sg_yr_gr_wmean = dummy_T_sig;
    
    sg_weight = dummy_T_sig;
    sg_r_dist = complex(dummy_T_sig);
    sg_num_used = dummy_T_sig;
    sg_area = dummy_T_sig;
    sg_error= dummy_T_sig;
    sg_dZ_scale = dummy_T_sig;
    disp('end initialize time')
end

dummy_T = NaN([tim 1],'double');
ML_dens_gr= dummy_T;	
ML_pres_gr= dummy_T;
ML_salt_gr= dummy_T;
ML_temp_gr= dummy_T;
ML_year_gr= dummy_T;
ML_salt_gr_wmean= dummy_T;
ML_temp_gr_wmean= dummy_T;
ML_dens_gr_wmean= dummy_T;
ML_pres_gr_wmean= dummy_T;
for ifield=1:length(AdditionalFields);
	field=AdditionalFields{ifield};
	eval(['ML_' field '_gr = dummy_T;']);%=============================MODIF
	eval(['ML_' field '_gr_wmean = dummy_T;']);
end
ML_area= dummy_T;
ML_error= dummy_T;
ML_weight=dummy_T;
ML_num_used=dummy_T;
ML_r_dist = dummy_T;
ML_Drho = dummy_T;

sg_pr = NaN([max_profs_inmem,sg_l],'single');
sg_sa = sg_pr;
sg_th = sg_pr;

MLpres = NaN([max_profs_inmem,1],'single');
MLdens = MLpres;
MLsalt = MLpres;
MLtemp = MLpres;
for ifield=1:length(AdditionalFields);
	field=AdditionalFields{ifield};
	eval(['ML' field ' = MLpres;']);%=============================MODIF
	eval(['ML' field ' = MLpres;']);
end

profs_inmem = (-100000:-100000+max_profs_inmem-1)-9e12; % set to 10,000 unique values not within IDs from floats or CTD ... ...

%% INITIALIZE ACCESS TO FAST MARCHING DATA:
disp(' loading Fast Marching gridpoints')

ncid_FM = netcdf.open(FM_fname,'NOWRITE');
FM_id_LON = netcdf.inqVarID(ncid_FM,'lon');
FM_id_LAT = netcdf.inqVarID(ncid_FM,'lat');
FM_id_DPT = netcdf.inqVarID(ncid_FM,'depth');
%FM_id_NUM = netcdf.inqVarID(ncid_FM,'num_subpoints');
FM_id_SUBDIS = netcdf.inqVarID(ncid_FM,'sub_distance');
FM_id_SUBANG = netcdf.inqVarID(ncid_FM,'sub_angle');
FM_id_SUBLAT = netcdf.inqVarID(ncid_FM,'sub_lat');
FM_id_SUBLON = netcdf.inqVarID(ncid_FM,'sub_lon');

FM_id_SUBDIS_C = netcdf.inqVarID(ncid_FM,'sub_distance_coarse');
FM_id_SUBANG_C = netcdf.inqVarID(ncid_FM,'sub_angle_coarse');
FM_id_SUBLAT_C = netcdf.inqVarID(ncid_FM,'sub_lat_coarse');
FM_id_SUBLON_C = netcdf.inqVarID(ncid_FM,'sub_lon_coarse');

dimid = netcdf.inqDimID(ncid_FM,'subgrid');
[~,FM_ALL] = netcdf.inqDim(ncid_FM,dimid);
dimid = netcdf.inqDimID(ncid_FM,'subgrid_coarse');
[~,FM_ALL_C] = netcdf.inqDim(ncid_FM,dimid);

FM_LAT = double(netcdf.getVar(ncid_FM,FM_id_LAT));
FM_LON = double(netcdf.getVar(ncid_FM,FM_id_LON));
FM_DPT = double(netcdf.getVar(ncid_FM,FM_id_DPT));
%FM_NUM = double(netcdf.getVar(ncid_FM,FM_id_NUM));
% cut down data set in memory to only data that might get needed:

disp(' wrapping data around globe')

% remove data that won't be used from LAT and LON vectors
% doing this PRIOR merging CTD and Argo data since both vectors already can be VERY large.
found_it = find(FM_LON < 0);

FM_LON(found_it) = FM_LON(found_it)+360;
%%
found_it = find(AFD_LON < 0);

AFD_LON(found_it) = AFD_LON(found_it)+360;

found_it = find(CTD_LON < 0);

CTD_LON(found_it) = CTD_LON(found_it)+360;


% Grid to compute adapted to FM grid: 
FM_LAT=round(FM_LAT.*100)/100;
FM_LON=round(FM_LON.*100)/100;

yy=unique(FM_LAT);
xx=unique(FM_LON);
xx_ok=find(xx>=min(lon_grid) & xx<=max(lon_grid));
yy_ok=find(yy>=min(lat_grid) & yy<=max(lat_grid));
lon_grid=xx(xx_ok); 
lat_grid=yy(yy_ok);
[x,y]=meshgrid(lon_grid, lat_grid);
[m,n]=size(x)  % grid to compute

%% join all data sets

lon = double([AFD_LON; CTD_LON]);
lat = double([AFD_LAT; CTD_LAT]);
dyr = double([AFD_DYR; CTD_DYR]);

dyr_map = dyr;

dec_dyr = double(dyr - centeryear);
%dec_dyr(dec_dyr > 0) = 0;
dec_dyr = (dec_noise_offset + signal_to_noise) - dec_noise_offset*exp(-(abs(dec_dyr./dec_scale)).^expweight) ;

dyr = double(rem(dyr,1));
idv = [AFD_IDVEC; CTD_IDVEC];

dist_dummy = NaN(size(lat));

%free some memory
clear AFD_LAT AFD_LON AFD_DYR CTD_LAT CTD_LON CTD_DYR CTD_IDVEC AFD_IDVEC found_it

%% INITIALIZE FINAL OUTPUT NETCDF FILE

% write to one LARGE netcdf file:::::
if ~exist(output_fname,'file') || create_new
    disp('si il nexiste pas: create new file Clim.nc')
    ncid_CLIM = netcdf.create(output_fname,'NC_64BIT_OFFSET');
    
    % Define the dimensions of the variables.
    dimid_gridLON    = netcdf.defDim(ncid_CLIM,'LON_VEC',n);
    dimid_gridLAT    = netcdf.defDim(ncid_CLIM,'LAT_VEC',m);
    dimid_months     = netcdf.defDim(ncid_CLIM,'MONTH_VEC',tim);
    dimid_isopycnals = netcdf.defDim(ncid_CLIM,'PDENS_VEC',sg_l);
    dimid_scalar     = netcdf.defDim(ncid_CLIM,'SCALAR',1);
    
    
    % Define new variables in file.
    latID  = netcdf.defVar(ncid_CLIM,'LATITUDE','float',dimid_gridLAT);
    lonID  = netcdf.defVar(ncid_CLIM,'LONGITUDE','float',dimid_gridLON);
    sigID  = netcdf.defVar(ncid_CLIM,'NEUTRAL_DENSITY','float',dimid_isopycnals);
    timID  = netcdf.defVar(ncid_CLIM,'MONTH','float',dimid_months);
    
    iqrID  = netcdf.defVar(ncid_CLIM,'IQR_RANGE','float',dimid_scalar);
    max_numID    = netcdf.defVar(ncid_CLIM,'MAX_NUM_PROFS','short',dimid_scalar);
    max_num2ID    = netcdf.defVar(ncid_CLIM,'MAX_NUM_PROFS_IN_MEMORY_GRIDPOINT','short',dimid_scalar);
    min_numID    = netcdf.defVar(ncid_CLIM,'MIN_NUM_PROFS','short',dimid_scalar);
    rnd_numID    = netcdf.defVar(ncid_CLIM,'RND_NUM_PROFS','short',dimid_scalar);
    ctd_numID    = netcdf.defVar(ncid_CLIM,'REQ_CTD_PROFS','short',dimid_scalar);
    horiz_scaleID    = netcdf.defVar(ncid_CLIM,'HORIZONTAL_SCALE','float',dimid_scalar);
    temporal_scaleID    = netcdf.defVar(ncid_CLIM,'TEMPORAL_SCALE','float',dimid_scalar);
    noise_signalID    = netcdf.defVar(ncid_CLIM,'NOISE_SIGNAL_RATIO','float',dimid_scalar);
    decadal_scaleID    = netcdf.defVar(ncid_CLIM,'DECADAL_SCALE','float',dimid_scalar);
    decadal_noiseID    = netcdf.defVar(ncid_CLIM,'DECADAL_NOISE','float',dimid_scalar);
    disp('end: create new file Clim.nc')
    
    if compute_interior
        
        data_pres_ID = netcdf.defVar(ncid_CLIM,'PRESSURE','float',[dimid_isopycnals dimid_months dimid_gridLAT dimid_gridLON]);
        netcdf.putAtt(ncid_CLIM,data_pres_ID,'_FillValue',NaN('single'))
        data_salt_ID = netcdf.defVar(ncid_CLIM,'SALINITY','float',[dimid_isopycnals dimid_months dimid_gridLAT dimid_gridLON]);
        netcdf.putAtt(ncid_CLIM,data_salt_ID,'_FillValue',NaN('single'))
        data_temp_ID = netcdf.defVar(ncid_CLIM,'POTENTIAL_TEMPERATURE','float',[dimid_isopycnals dimid_months dimid_gridLAT dimid_gridLON]);
        netcdf.putAtt(ncid_CLIM,data_temp_ID,'_FillValue',NaN('single'))
        data_pres_wm_ID = netcdf.defVar(ncid_CLIM,'PRESSURE_WEIGHTED_MEAN','float',[dimid_isopycnals dimid_months dimid_gridLAT dimid_gridLON]);
        netcdf.putAtt(ncid_CLIM,data_pres_wm_ID,'_FillValue',NaN('single'))
        data_salt_wm_ID = netcdf.defVar(ncid_CLIM,'SALINITY_WEIGHTED_MEAN','float',[dimid_isopycnals dimid_months dimid_gridLAT dimid_gridLON]);
        netcdf.putAtt(ncid_CLIM,data_salt_wm_ID,'_FillValue',NaN('single'))
        data_temp_wm_ID = netcdf.defVar(ncid_CLIM,'POTENTIAL_TEMPERATURE_WEIGHTED_MEAN','float',[dimid_isopycnals dimid_months dimid_gridLAT dimid_gridLON]);
        netcdf.putAtt(ncid_CLIM,data_temp_wm_ID,'_FillValue',NaN('single'))
        data_error_ID= netcdf.defVar(ncid_CLIM,'ERROR_OBJECTIVE_ANALYSIS','float',[dimid_isopycnals dimid_months dimid_gridLAT dimid_gridLON]);
        netcdf.putAtt(ncid_CLIM,data_error_ID,'_FillValue',NaN('single'))
        data_area_ID = netcdf.defVar(ncid_CLIM,'AREA_OBJECTIVE_ANALYSIS','float',[dimid_isopycnals dimid_months dimid_gridLAT dimid_gridLON]);
        netcdf.putAtt(ncid_CLIM,data_area_ID,'_FillValue',NaN('single'))
        data_nwgt_ID = netcdf.defVar(ncid_CLIM,'WEIGHT_SUMMED','float',[dimid_isopycnals dimid_months dimid_gridLAT dimid_gridLON]);
        netcdf.putAtt(ncid_CLIM,data_nwgt_ID,'_FillValue',NaN('single'))
        data_dist_Z_ID = netcdf.defVar(ncid_CLIM,'DISTANCE_WEIGHTED_ZONAL','short',[dimid_isopycnals dimid_months dimid_gridLAT dimid_gridLON]);
        netcdf.putAtt(ncid_CLIM,data_dist_Z_ID,'_FillValue',int16(-32000))
        data_dist_M_ID = netcdf.defVar(ncid_CLIM,'DISTANCE_WEIGHTED_MERIDIONAL','short',[dimid_isopycnals dimid_months dimid_gridLAT dimid_gridLON]);
        netcdf.putAtt(ncid_CLIM,data_dist_M_ID,'_FillValue',int16(-32000))
        data_num_data_ID = netcdf.defVar(ncid_CLIM,'NUMBER_DATA_POINTS','short',[dimid_isopycnals dimid_months dimid_gridLAT dimid_gridLON]);
        netcdf.putAtt(ncid_CLIM,data_num_data_ID,'_FillValue',int16(-32000))
        data_z_weight_ID = netcdf.defVar(ncid_CLIM,'DELTA_Z_SCALE','short',[dimid_isopycnals dimid_months dimid_gridLAT dimid_gridLON]);
        netcdf.putAtt(ncid_CLIM,data_z_weight_ID,'_FillValue',int16(-32000))
        
        data_year_ID= netcdf.defVar(ncid_CLIM,'YEAR_OF_DATA','short',[dimid_isopycnals dimid_months dimid_gridLAT dimid_gridLON]);
        netcdf.putAtt(ncid_CLIM,data_year_ID,'_FillValue',int16(-32000))
        %  data_year_wm_ID = netcdf.defVar(ncid_CLIM,'YEAR_OF_DATA_WEIGHTED_MEAN','short',[dimid_isopycnals dimid_months dimid_gridLAT dimid_gridLON]);
        %  netcdf.putAtt(ncid_CLIM,data_year_wm_ID,'_FillValue',int16(-32000))
        
    end
    
    disp('Define new variables for ML')
    ML_pres_ID  = netcdf.defVar(ncid_CLIM,'ML_MAX_PRESSURE','float',[dimid_months dimid_gridLAT dimid_gridLON]);
    netcdf.putAtt(ncid_CLIM,ML_pres_ID,'_FillValue',NaN('single'))
    ML_dens_ID  = netcdf.defVar(ncid_CLIM,'ML_DENSITY','float',[dimid_months dimid_gridLAT dimid_gridLON]);
    netcdf.putAtt(ncid_CLIM,ML_dens_ID,'_FillValue',NaN('single'))
    ML_salt_ID  = netcdf.defVar(ncid_CLIM,'ML_SALINITY','float',[dimid_months dimid_gridLAT dimid_gridLON]);
    netcdf.putAtt(ncid_CLIM,ML_salt_ID,'_FillValue',NaN('single'))
    ML_temp_ID  = netcdf.defVar(ncid_CLIM,'ML_TEMPERATURE','float',[dimid_months dimid_gridLAT dimid_gridLON]);
    netcdf.putAtt(ncid_CLIM,ML_temp_ID,'_FillValue',NaN('single'))
    ML_pres_wm_ID  = netcdf.defVar(ncid_CLIM,'ML_PRESSURE_WEIGHTED_MEAN','float',[dimid_months dimid_gridLAT dimid_gridLON]);
    netcdf.putAtt(ncid_CLIM,ML_pres_wm_ID,'_FillValue',NaN('single'))
    ML_dens_wm_ID  = netcdf.defVar(ncid_CLIM,'ML_DENSITY_WEIGHTED_MEAN','float',[dimid_months dimid_gridLAT dimid_gridLON]);
    netcdf.putAtt(ncid_CLIM,ML_dens_wm_ID,'_FillValue',NaN('single'))
    ML_salt_wm_ID  = netcdf.defVar(ncid_CLIM,'ML_SALINITY_WEIGHTED_MEAN','float',[dimid_months dimid_gridLAT dimid_gridLON]);
    netcdf.putAtt(ncid_CLIM,ML_salt_wm_ID,'_FillValue',NaN('single'))
    ML_temp_wm_ID  = netcdf.defVar(ncid_CLIM,'ML_TEMPERATURE_WEIGHTED_MEAN','float',[dimid_months dimid_gridLAT dimid_gridLON]);
    netcdf.putAtt(ncid_CLIM,ML_temp_wm_ID,'_FillValue',NaN('single'))
    ML_error_ID= netcdf.defVar(ncid_CLIM,'ML_ERROR','float',[dimid_months dimid_gridLAT dimid_gridLON]);
    netcdf.putAtt(ncid_CLIM,ML_error_ID,'_FillValue',NaN('single'))
    ML_area_ID = netcdf.defVar(ncid_CLIM,'ML_AREA','float',[dimid_months dimid_gridLAT dimid_gridLON]);
    netcdf.putAtt(ncid_CLIM,ML_area_ID,'_FillValue',NaN('single'))
    ML_nwgt_ID = netcdf.defVar(ncid_CLIM,'ML_WEIGHT_SUMMED','float',[dimid_months dimid_gridLAT dimid_gridLON]);
    netcdf.putAtt(ncid_CLIM,ML_nwgt_ID,'_FillValue',NaN('single'))
    ML_dist_Z_ID = netcdf.defVar(ncid_CLIM,'ML_DISTANCE_AVG_ZONAL','short',[dimid_months dimid_gridLAT dimid_gridLON]);
    netcdf.putAtt(ncid_CLIM,ML_dist_Z_ID,'_FillValue',int16(-32000))
    ML_dist_M_ID = netcdf.defVar(ncid_CLIM,'ML_DISTANCE_AVG_MERIDIONAL','short',[dimid_months dimid_gridLAT dimid_gridLON]);
    netcdf.putAtt(ncid_CLIM,ML_dist_M_ID,'_FillValue',int16(-32000))
    ML_num_data_ID = netcdf.defVar(ncid_CLIM,'ML_NUM_PROFS','short',[dimid_months dimid_gridLAT dimid_gridLON]);
    netcdf.putAtt(ncid_CLIM,ML_num_data_ID,'_FillValue',int16(0))
    ML_drho_ID = netcdf.defVar(ncid_CLIM,'ML_DELTA_RHO_SCALE','float',[dimid_months dimid_gridLAT dimid_gridLON]);
    netcdf.putAtt(ncid_CLIM,ML_drho_ID,'_FillValue',NaN('single'))
    ML_year_ID= netcdf.defVar(ncid_CLIM,'ML_YEAR_OF_DATA','short',[dimid_months dimid_gridLAT dimid_gridLON]);
    netcdf.putAtt(ncid_CLIM,ML_year_ID,'_FillValue',int16(-32000))
    for ifield=1:length(AdditionalFields);
	field=AdditionalFields{ifield};
	eval(['ML_' field '_ID  = netcdf.defVar(ncid_CLIM,''ML_' field ''',''float'',[dimid_months dimid_gridLAT dimid_gridLON]);']); %=======================	MODIF
	eval(['netcdf.putAtt(ncid_CLIM,ML_' field '_ID,''_FillValue'',NaN(''single''))']);
	eval(['ML_' field '_wm_ID  = netcdf.defVar(ncid_CLIM,''ML_' field '_WEIGHTED_MEAN'',''float'',[dimid_months dimid_gridLAT dimid_gridLON]);']); %=======================	MODIF
	eval(['netcdf.putAtt(ncid_CLIM,ML_' field '_wm_ID,''_FillValue'',NaN(''single''))']);
    end
    % ML_year_wm_ID = netcdf.defVar(ncid_CLIM,'ML_YEAR_OF_DATA_WEIGHTED_MEAN','short',[dimid_months dimid_gridLAT dimid_gridLON]);
    % netcdf.putAtt(ncid_CLIM,ML_year_wm_ID,'_FillValue',int16(-32000))
    disp('end: Define new variables for ML')
    
    % Leave define mode and enter data mode to write data.
    netcdf.endDef(ncid_CLIM)
    
    % Write constants and grid:
    disp('Write constants and grid')
    netcdf.putVar(ncid_CLIM,latID,y(:,1));
    netcdf.putVar(ncid_CLIM,lonID,x(1,:));
    netcdf.putVar(ncid_CLIM,sigID,sg0_pleth);
    netcdf.putVar(ncid_CLIM,timID,(1:tim)/tim*12);
    netcdf.putVar(ncid_CLIM,iqrID,iqr_mult);
    netcdf.putVar(ncid_CLIM,max_numID,max_num_profs);
    
    netcdf.putVar(ncid_CLIM,max_num2ID,max_profs_inmem);
    netcdf.putVar(ncid_CLIM,min_numID,min_num_profs);
    netcdf.putVar(ncid_CLIM,rnd_numID,num_of_rand_profiles);
    netcdf.putVar(ncid_CLIM,ctd_numID,min_CTD_profs);
    netcdf.putVar(ncid_CLIM,horiz_scaleID,horiz_scale_set);
    netcdf.putVar(ncid_CLIM,temporal_scaleID,time_scale);
    
    
    netcdf.putVar(ncid_CLIM,noise_signalID,signal_to_noise-1);
    netcdf.putVar(ncid_CLIM,decadal_scaleID,dec_scale);
    netcdf.putVar(ncid_CLIM,decadal_noiseID,dec_noise_offset);
    disp('end: Write constants and grid')
    
else
    disp('si Clim.nc existe: ouverture')
    ncid_CLIM = netcdf.open(output_fname,'WRITE');
    % get variable IDs in file.
    latID  = netcdf.inqVarID(ncid_CLIM,'LATITUDE');
    lonID  = netcdf.inqVarID(ncid_CLIM,'LONGITUDE');
    sigID  = netcdf.inqVarID(ncid_CLIM,'NEUTRAL_DENSITY');
    timID  = netcdf.inqVarID(ncid_CLIM,'MONTH');
    iqrID  = netcdf.inqVarID(ncid_CLIM,'IQR_RANGE');
    max_numID    = netcdf.inqVarID(ncid_CLIM,'MAX_NUM_PROFS');
    if compute_interior
        data_pres_ID = netcdf.inqVarID(ncid_CLIM,'PRESSURE');
        data_salt_ID = netcdf.inqVarID(ncid_CLIM,'SALINITY');
        data_temp_ID = netcdf.inqVarID(ncid_CLIM,'TEMPERATURE');
        data_pres_wm_ID = netcdf.inqVarID(ncid_CLIM,'PRESSURE_WEIGTHED_MEAN');
        data_salt_wm_ID = netcdf.inqVarID(ncid_CLIM,'SALINITY_WEIGTHED_MEAN');
        data_temp_wm_ID = netcdf.inqVarID(ncid_CLIM,'TEMPERATURE_WEIGTHED_MEAN');
        data_error_ID= netcdf.inqVarID(ncid_CLIM,'ERROR');
        data_area_ID = netcdf.inqVarID(ncid_CLIM,'AREA');
        data_nwgt_ID = netcdf.inqVarID(ncid_CLIM,'WEIGHT_SUMMED');
        data_dist_Z_ID = netcdf.inqVarID(ncid_CLIM,'DISTANCE_AVG_ZONAL');
        data_dist_M_ID = netcdf.inqVarID(ncid_CLIM,'DISTANCE_AVG_MERIDIONAL');
        data_num_data_ID = netcdf.inqVarID(ncid_CLIM,'NUM_PROFS');
    end
    ML_pres_ID  = netcdf.inqVarID(ncid_CLIM,'ML_MAX_PRESSURE');
    ML_dens_ID  = netcdf.inqVarID(ncid_CLIM,'ML_DENSITY');
    ML_salt_ID  = netcdf.inqVarID(ncid_CLIM,'ML_SALINITY');
    ML_temp_ID  = netcdf.inqVarID(ncid_CLIM,'ML_TEMPERATURE');
    ML_pres_wm_ID  = netcdf.inqVarID(ncid_CLIM,'ML_PRESSURE_WEIGHTED_MEAN');
    ML_dens_wm_ID  = netcdf.inqVarID(ncid_CLIM,'ML_DENSITY_WEIGHTED_MEAN');
    ML_error_ID= netcdf.inqVarID(ncid_CLIM,'ML_ERROR');
    ML_area_ID = netcdf.inqVarID(ncid_CLIM,'ML_AREA');
    ML_nwgt_ID = netcdf.inqVarID(ncid_CLIM,'ML_WEIGHT_SUMMED');
    ML_dist_Z_ID = netcdf.inqVarID(ncid_CLIM,'ML_DISTANCE_AVG_ZONAL');
    ML_dist_M_ID = netcdf.inqVarID(ncid_CLIM,'ML_DISTANCE_AVG_MERIDIONAL');
    ML_num_data_ID = netcdf.inqVarID(ncid_CLIM,'ML_NUM_PROFS');
    ML_drho_ID  = netcdf.inqVarID(ncid_CLIM,'ML_DELTA_RHO_SCALE');
    ML_salt_wm_ID  = netcdf.inqVarID(ncid_CLIM,'ML_SALINITY_WEIGHTED_MEAN');
    ML_temp_wm_ID  = netcdf.inqVarID(ncid_CLIM,'ML_TEMPERATURE_WEIGHTED_MEAN');
    ML_year_ID  = netcdf.inqVarID(ncid_CLIM,'ML_YEAR_OF_DATA');

    for ifield=1:length(AdditionalFields);
	field=AdditionalFields{ifield};
	eval(['ML_' field '_ID     = netcdf.inqVarID(ncid_CLIM, ''ML_' field ''');']); 
	eval(['ML_' field '_wm_ID  = netcdf.inqVarID(ncid_CLIM, ''ML_' field '_WEIGHTED_MEAN'');']); 

    end
    disp('end de l ouverture de Clim.nc')
end


[~, seed_rndprofs] = sort(rand(1,max_profs_inmem - (numsortedprofs)));
seed_rndprofs = seed_rndprofs + numsortedprofs;
time_scale = time_scale^2;

%% ALL INITIALIZED FOR WARP ... hmmm COMPUTATION  -   engage

% take time of start

tstart = now;
%disp([' Initializing complete at ' datestr(tstart,'HH:MM:SS') ', starting computation.'])%sunke version
warning('off','MATLAB:lscov:RankDefDesignMat')

for istep = 1:m 	%Latitude loop 140 times / 280 times if 0.5? ...
    tic
    if y(istep,1) >= continue_at && y(istep,1) <= continue_upto
        time_now = rem(now,1);
        if time_now > 0.25 && time_now <0.75
            setenv('OMP_NUM_THREADS', num2str(OMP_NUM_THREADS_DAY));
        else
            setenv('OMP_NUM_THREADS', num2str(OMP_NUM_THREADS_NIGHT));
        end
        %  fprintf(' calculating longitude:                     ')%sunke version
        for jstep=1:n %longitude loop 360 times / 720 times if 0.5? ...
            if x(istep,jstep) >= continue_lon_at && x(istep,jstep) <= continue_lon_upto
                % keeping track of current executions:
                %fprintf(['\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b', num2str(x(istep,jstep),'%05.1f') ' after ' datestr(now-tstart,'HH:MM:SS')]) %sunke version
                
                % Find Fast marching grid point
                FM_pos = find(FM_LAT == y(istep,jstep) & ( FM_LON == x(istep,jstep) | FM_LON == x(istep,jstep)-360));
                %disp(['FM: ' num2str(y(istep,jstep)) ';' num2str(x(istep,jstep))])
                if ~isempty(FM_pos)
                    FM_pos = FM_pos(1);
                    %disp('CA PASSE')
                    if upperocean
                        local_depth = min([upperocean_limit, abs(FM_DPT(FM_pos))]);
                    else
                        local_depth = abs(FM_DPT(FM_pos)); %#ok<UNRCH>
                    end
                    
                    % load fast marching subgrid
                    fm_dis = double(netcdf.getVar(ncid_FM,FM_id_SUBDIS,[0 FM_pos-1],[FM_ALL 1]));
                    
                    % keyboard
                    
                    gefFM = find(fm_dis<2000);
                    
                    if numel(gefFM) > 10
                        
                        % get lon from FM grid
                        fm_lat = double(netcdf.getVar(ncid_FM,FM_id_SUBLAT,[0 FM_pos-1],[FM_ALL 1]));
                        fm_lon = double(netcdf.getVar(ncid_FM,FM_id_SUBLON,[0 FM_pos-1],[FM_ALL 1]));
                        fm_ang = double(netcdf.getVar(ncid_FM,FM_id_SUBANG,[0 FM_pos-1],[FM_ALL 1]));
                        
                        % load coarse grid:
                        fm_lat_c = double(netcdf.getVar(ncid_FM,FM_id_SUBLAT_C,[0 FM_pos-1],[FM_ALL_C 1]));
                        fm_lon_c = double(netcdf.getVar(ncid_FM,FM_id_SUBLON_C,[0 FM_pos-1],[FM_ALL_C 1]));
                        fm_dis_c = double(netcdf.getVar(ncid_FM,FM_id_SUBDIS_C,[0 FM_pos-1],[FM_ALL_C 1]));
                        fm_ang_c = double(netcdf.getVar(ncid_FM,FM_id_SUBANG_C,[0 FM_pos-1],[FM_ALL_C 1]));
                        
                        % use only valid data:
                        fm_dis = fm_dis(gefFM);
                        fm_ang = fm_ang(gefFM);
                        fm_lat = fm_lat(gefFM);
                        fm_lon = fm_lon(gefFM);
                        fm_lon(fm_lon<0) = fm_lon(fm_lon<0) +360;
                        max_small_FM = max(fm_dis);
                        gefFMc = find(fm_dis_c < 2800 & fm_dis_c >= min([max_small_FM 300]));
                        fm_dis_c = fm_dis_c(gefFMc);
                        fm_ang_c = fm_ang_c(gefFMc);
                        fm_lat_c = fm_lat_c(gefFMc);
                        fm_lon_c = fm_lon_c(gefFMc);
                        fm_lon_c(fm_lon_c<0) = fm_lon_c(fm_lon_c<0) +360;
                        %keyboard
                        
                        % take subset of all float/CTD data within range of FM data points
                        maxfmlat = max([max(fm_lat_c) max(fm_lat)]);
                        minfmlat = min([min(fm_lat_c) min(fm_lat)]);
                        maxfmlon = max([max(fm_lon_c) max(fm_lon)]);
                        minfmlon = min([min(fm_lon_c) min(fm_lon)]);
                        
                        fm_fit = find(abs(lat-(maxfmlat+minfmlat)/2) < (maxfmlat-minfmlat)/2 & ...
                            abs(lon-(maxfmlon+minfmlon)/2) < (maxfmlon-minfmlon)/2);
                        
                        % check if enough data is found
                        if ~isempty(fm_fit) && numel(fm_fit) > min_num_profs
                               %keyboard
                            
                            tmp_fm_lon = fm_lon-x(istep,jstep);
                            tmp_fm_lon(tmp_fm_lon>180) = tmp_fm_lon(tmp_fm_lon>180) - 360;
                            tmp_fm_lon(tmp_fm_lon<-180) = tmp_fm_lon(tmp_fm_lon<-180) + 360;
                            tmp_fm_lon = tmp_fm_lon.*cosd(fm_lat);
                            
                            
                            tmp_fm_lon_c = fm_lon_c-x(istep,jstep);
                            tmp_fm_lon_c(tmp_fm_lon_c>180) = tmp_fm_lon_c(tmp_fm_lon_c>180) - 360;
                            tmp_fm_lon_c(tmp_fm_lon_c<-180) = tmp_fm_lon_c(tmp_fm_lon_c<-180) + 360;
                            tmp_fm_lon_c = tmp_fm_lon_c.*cosd(fm_lat_c);
                            
                            % initialize irregular grid interpolant
                            DT = DelaunayTri(double(tmp_fm_lon),double(fm_lat));
                            
                            [PI,D] = nearestNeighbor(DT,tmp_fm_lon_c,fm_lat_c);
                            to_close = find(D<maxdist_finegrid);
                            tmp_fm_lon_c(to_close) = [];
                            fm_lat_c(to_close) = [];
                            fm_dis_c(to_close) = [];
                            fm_ang_c(to_close) = [];
                            
                            if numel(fm_lat_c) > 0
                                % extend Delaunay Triangulation by coarse gridpoints
                                DT.X(end+1 : end+numel(fm_lat_c) ,:) = [tmp_fm_lon_c fm_lat_c];
                            end
                            
                            %
                            if ~isempty(DT.Triangulation)
                                tmp_lon = lon(fm_fit)-x(istep,jstep);
                                tmp_lon(tmp_lon>180) = tmp_lon(tmp_lon>180) - 360;
                                tmp_lon(tmp_lon<-180) = tmp_lon(tmp_lon<-180) + 360;
                                tmp_lon = tmp_lon.*cosd(lat(fm_fit));
                                
                                % find nearest neighbour indices
                                [PI,D] = nearestNeighbor(DT,tmp_lon,lat(fm_fit));
                                
                                fm_dis_vec_combined = [fm_dis; fm_dis_c];
                                
                                % accept divergences up to half FM grid point
                                D_sub = find((D<=maxdist_finegrid+1 & fm_dis_vec_combined(PI) <= max_small_FM+50) | (D<=maxdist_coarsegrid+1+maxdist_finegrid & fm_dis_vec_combined(PI) >= max_small_FM-50));
                                D_subset = fm_fit(D_sub);
                                DD_sub = PI(D_sub);
                                % GRB CLOSE DATA HERE >>>>
                                if numel(D_subset) > min_num_profs &&  numel(D_sub) >= 10
                                    
                                    
                                    % do a more accurate distance computation with subset: (cannot be done earlier due to interpolation errors at bays, inlets, ridges ...
                                    TSI = TriScatteredInterp(DT,fm_dis_vec_combined);
                                    
                                    % in future get rid of one of these 2 dis vectors ... only keep compex one
                                    % put all data into complex dist_vec
                                    dist_vec = dist_dummy;
                                    dist_vec_complex = complex(dist_vec,0);
                                    
                                    % fill distance vector with distance
                                    tmp =  TSI(double(tmp_lon(D_sub)),lat(fm_fit(D_sub)))./110;
                                    gefu = find(isnan(tmp));
                                    if ~isempty(tmp)
                                        dist_vec(D_subset)  = tmp;
                                        gef = isnan(dist_vec(D_subset));
                                        D_subset(gef) = [];
                                        DD_sub(gef) = [];
                                        
                                        if ~isempty(D_subset)
                                            
                                            complex_angle_sin = sind([fm_ang; fm_ang_c]);
                                            complex_angle_cos = cosd([fm_ang; fm_ang_c]);
                                            TSI_sin = TriScatteredInterp(DT,complex_angle_sin);
                                            TSI_cos = TriScatteredInterp(DT,complex_angle_cos);
                                            tmp_sin = TSI_sin(double(tmp_lon(D_sub)),lat(fm_fit(D_sub)));
                                            tmp_cos = TSI_cos(double(tmp_lon(D_sub)),lat(fm_fit(D_sub)));
                                            fm_ang_vec_combined = atan2(tmp_sin,tmp_cos) *180/pi;
                                            fm_ang_vec_combined(gefu) = [];
                                            if numel(fm_ang_vec_combined) ~= numel(D_subset)
                                                disp('this should not happen!!!')
                                                keyboard
                                            end
                                            
                                            % add some noise to direction to prevent colinear data in bays, inlets ...etc
                                            %fm_ang_vec_combined(DD_sub) = fm_ang_vec_combined(DD_sub) + randi(3,numel(DD_sub),1) -2 ;
                                            
                                            %fill complex distance vector with data
                                            dist_vec_complex(D_subset)  = - sind(fm_ang_vec_combined) .* dist_vec(D_subset) + 1i .* cosd(fm_ang_vec_combined) .* dist_vec(D_subset);
                                            
                                            % find profiles within range
                                            ii_dist = D_subset;
                                            %dist_vec(ii_dist) = dist_vec(ii_dist)./maxrange;
                                            %dist_vec_complex(ii_dist) = dist_vec_complex(ii_dist);%./maxrange;
                                            
                                            
                                            % if more than MAX_PROFS_INMEM profiles are in close range, only use closest XXXX don't load more;
                                            if numel(ii_dist) > max_profs_inmem
                                                
                                                % find the closest profs .... BUT maxe sure in open ocean at least 500 CTD full depth profiles are included ...
                                                [~, kk]=sort(dist_vec(ii_dist));
                                                ii_dist_tmp = ii_dist(kk(1:max_profs_inmem));
                                                
                                                if sum(idv(kk(1:max_profs_inmem))<0) < min_CTD_profs
                                                    gef_CTD = find(idv(ii_dist) < 0);
                                                    numCTD = numel(gef_CTD);
                                                    
                                                    [~, kk2]=sort(dist_vec(ii_dist(gef_CTD)));
                                                    ii_dist2 = ii_dist(gef_CTD(kk2(1:min([min_CTD_profs numCTD]))));
                                                    
                                                    [~,IA,~] =setxor(ii_dist_tmp, ii_dist2);
                                                    ii_dist_tmp = ii_dist_tmp(IA);
                                                    
                                                    ii_dist_tmp = [ii_dist2; ii_dist_tmp(1:max_profs_inmem-numel(ii_dist2))];
                                                end
                                                
                                                ii_dist = ii_dist_tmp;
                                            end
                                            
                                            if drawit
                                                figure(1) %#ok<UNRCH>
                                                clf
                                                m_proj('Azimuthal Equal-area','lat',y(istep,jstep),'long',x(istep,jstep),'rad',20,'rect','on');
                                                m_plot(DT.X(1:4:end,1)./cosd(DT.X(1:4:end,2))+x(istep,jstep),DT.X(1:4:end,2),'.k')
                                                hold on
                                                m_plot(lon(ii_dist(1:10:end)),lat(ii_dist(1:10:end)),'+r')
                                                
                                                m_plot(x(istep,jstep),y(istep,jstep),'og')
                                                title('Every 4th FM point and every 10th data point (out of the up to 2000 closest)')
                                                m_coast
                                                m_grid
                                                drawnow
                                            end
                                            % gather all information and to load profiles not yet in memory...
                                            profs_needed = idv(ii_dist);
                                            
                                            % what to read and what not to read
                                            [~, AI, BI] = setxor(profs_needed,profs_inmem);	%
                                            
                                            % these we do need to read, since not yet in memory
                                            profs_toload = profs_needed(AI);
                                            
                                            % these are in memory but no longer needed
                                            profs_redundant = BI;
                                            
                                            % gather all indices that are right next to each other to load at once and not in multiple steps
                                            multi_load_vec = single(diff(profs_toload) == 1) ;
                                            
                                            % load all profiles from file
                                            while ~isempty(profs_toload)
                                                
                                                % do final summarize of data that is adjecent to each other to load
                                                nc_load_width = sum(cumprod(multi_load_vec))+1;
                                                
                                                if numel(profs_redundant)>=nc_load_width
                                                    % assign id in memory block to overwrite
                                                    ID_overwrite = profs_redundant(1:nc_load_width);
                                                    profs_redundant(1:nc_load_width) = [];
                                                else
                                                    % assign  IDs at end of 'known' memory block
                                                    ID_overwrite = (numel(profs_inmem)+1) : (numel(profs_inmem)+nc_load_width);
                                                end
                                                
                                                
                                                switch sign(profs_toload(1)) + sign(profs_toload(nc_load_width)) % determine if CTD or Argo float data is needed:
                                                    case -2 % load data from CTD data base
                                                        if compute_interior
                                                            sg_th(ID_overwrite,:) = netcdf.getVar(ncid_CTD,CTD_id_TEMP,[sg_lev(1)-1 CTD_offset+profs_toload(1)],[sg_l nc_load_width])';
                                                            tmp = ones(numel(ID_overwrite),sg_l) ; %sg_th(ID_overwrite,:).*0+1;
                                                            tmp(sg_th(ID_overwrite,:)>100) = NaN;
                                                            sg_th(ID_overwrite,:) = sg_th(ID_overwrite,:) .* tmp;
                                                            sg_pr(ID_overwrite,:) = netcdf.getVar(ncid_CTD,CTD_id_PRES,[sg_lev(1)-1 CTD_offset+profs_toload(1)],[sg_l nc_load_width])'.*tmp;
                                                            sg_sa(ID_overwrite,:) = netcdf.getVar(ncid_CTD,CTD_id_SAL,[sg_lev(1)-1 CTD_offset+profs_toload(1)],[sg_l nc_load_width])'.*tmp;
                                                        end
                                                        ttmp = netcdf.getVar(ncid_CTD,CTD_id_ML_DENS,CTD_offset+profs_toload(1),nc_load_width);
                                                        if ttmp > 100
                                                            MLdens(ID_overwrite) = ttmp -1000; %================== MODIF
                                                        else
                                                            MLdens(ID_overwrite) = ttmp;
                                                        end
                                                        MLpres(ID_overwrite) = netcdf.getVar(ncid_CTD,CTD_id_ML_PRES,CTD_offset+profs_toload(1),nc_load_width);
                                                        MLsalt(ID_overwrite) = netcdf.getVar(ncid_CTD,CTD_id_ML_SALT,CTD_offset+profs_toload(1),nc_load_width);
                                                        MLtemp(ID_overwrite) = netcdf.getVar(ncid_CTD,CTD_id_ML_TEMP,CTD_offset+profs_toload(1),nc_load_width);
    							for ifield=1:length(AdditionalFields);
								field=AdditionalFields{ifield};
								eval(['ML' field '(ID_overwrite)  = netcdf.getVar(ncid_CTD,CTD_id_ML_' field ',CTD_offset+profs_toload(1),nc_load_width);']); %=======================	MODIF
    							end

                                                    case 2 % load data from Argo float data base
                                                        if compute_interior
                                                            sg_th(ID_overwrite,:) = netcdf.getVar(ncid_AFD,AFD_id_TEMP,[sg_lev(1)-1 profs_toload(1)],[sg_l nc_load_width])';
                                                            tmp = ones(numel(ID_overwrite),sg_l) ; %tmp = sg_th(ID_overwrite,:).*0+1;
                                                            tmp(sg_th(ID_overwrite,:)>100) = NaN;
                                                            sg_th(ID_overwrite,:) = sg_th(ID_overwrite,:) .* tmp;
                                                            sg_pr(ID_overwrite,:) = netcdf.getVar(ncid_AFD,AFD_id_PRES,[sg_lev(1)-1 profs_toload(1)],[sg_l nc_load_width])'.*tmp;
                                                            sg_sa(ID_overwrite,:) = netcdf.getVar(ncid_AFD,AFD_id_SAL,[sg_lev(1)-1 profs_toload(1)],[sg_l nc_load_width])'.*tmp;
                                                        end
                                                        ttmp = netcdf.getVar(ncid_AFD,AFD_id_ML_DENS,profs_toload(1),nc_load_width);
                                                        if ttmp > 100
                                                            MLdens(ID_overwrite) = ttmp-1000;
                                                        else
                                                            MLdens(ID_overwrite) = ttmp;
                                                        end
                                                        MLpres(ID_overwrite) = netcdf.getVar(ncid_AFD,AFD_id_ML_PRES,profs_toload(1),nc_load_width);
                                                        MLsalt(ID_overwrite) = netcdf.getVar(ncid_AFD,AFD_id_ML_SALT,profs_toload(1),nc_load_width);
                                                        MLtemp(ID_overwrite) = netcdf.getVar(ncid_AFD,AFD_id_ML_TEMP,profs_toload(1),nc_load_width);
    							for ifield=1:length(AdditionalFields);
								field=AdditionalFields{ifield};
								eval(['ML' field '(ID_overwrite)  = netcdf.getVar(ncid_AFD,AFD_id_ML_' field ',profs_toload(1),nc_load_width);']); %=======================	MODIF
    							end
                                                    otherwise % load data from both data sets
                                                        for loopstep = 1:nc_load_width % read each profile individually - ignore multi-load
                                                            if profs_toload(loopstep)<0
                                                                if compute_interior
                                                                    sg_th(ID_overwrite(loopstep),:) = netcdf.getVar(ncid_CTD,CTD_id_TEMP,[sg_lev(1)-1 CTD_offset+profs_toload(loopstep)],[sg_l 1])';
                                                                    tmp = ones_sig ;
                                                                    tmp(sg_th(ID_overwrite(loopstep),:)>100) = NaN;
                                                                    sg_th(ID_overwrite(loopstep),:) = sg_th(ID_overwrite(loopstep),:) .* tmp;
                                                                    sg_pr(ID_overwrite(loopstep),:) = netcdf.getVar(ncid_CTD,CTD_id_PRES,[sg_lev(1)-1 CTD_offset+profs_toload(loopstep)],[sg_l 1])'.*tmp;
                                                                    sg_sa(ID_overwrite(loopstep),:) = netcdf.getVar(ncid_CTD,CTD_id_SAL,[sg_lev(1)-1 CTD_offset+profs_toload(loopstep)],[sg_l 1])'.*tmp;
                                                                end
                                                                ttmp =  netcdf.getVar(ncid_CTD,CTD_id_ML_DENS,CTD_offset+profs_toload(loopstep),1);
                                                                if ttmp > 100
                                                                    MLdens(ID_overwrite(loopstep)) = ttmp-1000;
                                                                else
                                                                    MLdens(ID_overwrite(loopstep)) = ttmp;
                                                                end
                                                                MLpres(ID_overwrite(loopstep)) = netcdf.getVar(ncid_CTD,CTD_id_ML_PRES,CTD_offset+profs_toload(loopstep),1);
                                                                MLsalt(ID_overwrite(loopstep)) = netcdf.getVar(ncid_CTD,CTD_id_ML_SALT,CTD_offset+profs_toload(loopstep),1);
                                                                MLtemp(ID_overwrite(loopstep)) = netcdf.getVar(ncid_CTD,CTD_id_ML_TEMP,CTD_offset+profs_toload(loopstep),1);
    								for ifield=1:length(AdditionalFields);
									field=AdditionalFields{ifield};
									eval(['ML' field '(ID_overwrite(loopstep))  = netcdf.getVar(ncid_CTD,CTD_id_ML_' field ',CTD_offset+profs_toload(loopstep),1);']); %=======================	MODIF
    								end
                                                             else
                                                                if compute_interior
                                                                    sg_th(ID_overwrite(loopstep),:) = netcdf.getVar(ncid_AFD,AFD_id_TEMP,[sg_lev(1)-1 profs_toload(loopstep)],[sg_l 1])';
                                                                    tmp = ones_sig ;%tmp = sg_th(ID_overwrite(loopstep),:).*0+1;
                                                                    tmp(sg_th(ID_overwrite(loopstep),:)>100) = NaN;
                                                                    sg_th(ID_overwrite(loopstep),:) = sg_th(ID_overwrite(loopstep),:) .* tmp;
                                                                    sg_pr(ID_overwrite(loopstep),:) = netcdf.getVar(ncid_AFD,AFD_id_PRES,[sg_lev(1)-1 profs_toload(loopstep)],[sg_l 1])'.*tmp;
                                                                    sg_sa(ID_overwrite(loopstep),:) = netcdf.getVar(ncid_AFD,AFD_id_SAL,[sg_lev(1)-1 profs_toload(loopstep)],[sg_l 1])'.*tmp;
                                                                end
                                                                ttmp = netcdf.getVar(ncid_AFD,AFD_id_ML_DENS,profs_toload(loopstep),1);
                                                                if ttmp > 100
                                                                    MLdens(ID_overwrite(loopstep)) = ttmp-1000; %==================MODIF
                                                                else
                                                                    MLdens(ID_overwrite(loopstep)) = ttmp; %==================MODIF
                                                                end
                                                                MLpres(ID_overwrite(loopstep)) = netcdf.getVar(ncid_AFD,AFD_id_ML_PRES,profs_toload(loopstep),1);
                                                                MLsalt(ID_overwrite(loopstep)) = netcdf.getVar(ncid_AFD,AFD_id_ML_SALT,profs_toload(loopstep),1);
                                                                MLtemp(ID_overwrite(loopstep)) = netcdf.getVar(ncid_AFD,AFD_id_ML_TEMP,profs_toload(loopstep),1);
    								for ifield=1:length(AdditionalFields);
									field=AdditionalFields{ifield};
									eval(['ML' field '(ID_overwrite(loopstep))  = netcdf.getVar(ncid_AFD,AFD_id_ML_' field ',profs_toload(loopstep),1);']); %=======================	MODIF
    								end
                                                            end
                                                        end
                                                end
                                                
                                                gef = find(MLdens(ID_overwrite) < 0 | MLdens(ID_overwrite) > 10000);
                                                MLdens(ID_overwrite(gef)) = NaN;
                                                MLpres(ID_overwrite(gef)) = NaN;
                                                MLsalt(ID_overwrite(gef)) = NaN;
                                                MLtemp(ID_overwrite(gef)) = NaN;
   						for ifield=1:length(AdditionalFields);
							field=AdditionalFields{ifield};
							eval(['ML' field '(ID_overwrite(gef)) = NaN;']); %=======================	MODIF
    						end
                                                
                                                % require ML for data to be used
                                                sg_pr(ID_overwrite(gef),:) = NaN;
                                                
                                                % which are now newly in mem?
                                                profs_inmem(ID_overwrite) = profs_toload(1:nc_load_width);
                                                
                                                % they no longer need to be loaded
                                                profs_toload(1:nc_load_width) = [];
                                                
                                                % remove multiload by same ammount and take care of overflow at end of vector
                                                if numel(multi_load_vec) > nc_load_width
                                                    multi_load_vec(1:nc_load_width) = [];
                                                else
                                                    multi_load_vec = 0;
                                                end
                                                
                                            end
                                            
                                            %  disp('done loading') --->>>  now sort! to have same indices for all
                                            i2_dist = find(ismember(profs_inmem,profs_needed));
                                            [~, kk] = sort(profs_needed);
                                            ii_dist = ii_dist(kk);
                                            [~, kk] = sort(profs_inmem(i2_dist));
                                            i2_dist = i2_dist(kk);
                                            
                                            sg_pr = double(sg_pr);
                                            sg_sa = double(sg_sa);
                                            sg_th = double(sg_th);
                                            MLdens = double(MLdens);
                                            MLpres = double(MLpres);
                                            MLsalt = double(MLsalt);
                                            MLtemp = double(MLtemp);
   					    for ifield=1:length(AdditionalFields);
						field=AdditionalFields{ifield};
						eval(['ML' field '= double(ML' field ');']); %=======================	MODIF
    					    end
                                            
                                            dyrTMP = dyr(ii_dist);
                                            dyrMAP = dyr_map(ii_dist);
                                            
                                            dist_vec_complexTMP = dist_vec_complex(ii_dist);
                                            dec_dyrTMP = dec_dyr(ii_dist);
                                            %compute the huge covariancematrix
                                            %   E_all2 = ((double(((repmat(real(dist_vec_complex(ii_dist)),1,numel(ii_dist)) - repmat(real(dist_vec_complex(ii_dist))',numel(ii_dist),1)).^2) ...
                                            %                   + ((repmat(imag(dist_vec_complex(ii_dist)),1,numel(ii_dist)) - repmat(imag(dist_vec_complex(ii_dist))',numel(ii_dist),1)).^2))));
                                            %               E_all2 =( sqrt(E_all2).^(expweight))./horiz_scale.^expweight;
                                            %     E_all2 = E_all2+ double(abs(repmat(dec_dyr(ii_dist),1,numel(ii_dist)) - repmat(dec_dyr(ii_dist)',numel(ii_dist),1)).^expweight ) ./ dec_scale.^expweight ;
                                            
                                            
                                            % C_all =   (sqrt(abs(real(dist_vec_complexTMP).*real(dist_vec_complexTMP)) ...
                                            %     + abs(imag(dist_vec_complexTMP).*imag(dist_vec_complexTMP))));
                                            % identical to:
                                            C_all =abs(dist_vec_complexTMP);
                                            
                                            % EXP 1.5
                                            %C_all = C_all .* sqrt(C_all) ./ (horiz_scale.*sqrt(horiz_scale));
                                            % EXP 2
                                            C_all = C_all .* C_all ./ horiz_scale;
                                            
                                            %EXP 1.5
                                            %E_all_T = make_dist_matrix_mex_exp_scal(...
                                            %    real(dist_vec_complexTMP), ...
                                            %    imag(dist_vec_complexTMP), ...
                                            %    horiz_scale.*sqrt(horiz_scale)) ...
                                            %    + make_time_matrix_mex_exp_scal(dyrTMP,time_scale.*sqrt(time_scale));
                                            %
                                            % EXP 2
                                            E_all_T = make_dist_matrix_mex_exp_scal_gauss(...
                                                real(dist_vec_complexTMP), ...
                                                imag(dist_vec_complexTMP), ...
                                                horiz_scale) ...
                                                + make_time_matrix_mex_exp_scal_gauss(dyrTMP,time_scale);
                                            
                                            
                                            % -------------------------------------
                                            % loop over months/seasons
                                            % -------------------------------------
                                            
                                            for tstep =  1:tim % about 12 times ... or 4 times ...
                                                
                                                switch tstep
                                                    case 1
                                                        i2_dist_tmp2 = i2_dist;
                                                        ii_dist_tmp2 = ii_dist;
                                                    otherwise
                                                        i2_dist = i2_dist_tmp2;
                                                        ii_dist = ii_dist_tmp2;
                                                end
                                                
                                                % current month is:
                                                month_step= (tstep/tim)*12; % if tim == 12 then this is equal to month_step = tstep
                                                
                                                % compute fraction of year from month and center all data
                                                % in time around this month set to 0
                                                time_dist=( dyrTMP- (month_step-0.5)/12 ...
                                                    -((dyrTMP-(month_step-0.5)/12) >  0.5)...
                                                    +((dyrTMP-(month_step-0.5)/12) < -0.5) )';
                                                
                                                C_all_T2 = abs(time_dist);
                                                
                                                %EXP 1.5
                                                %C_all_T2 = C_all_T2 .* sqrt(C_all_T2) ./ (time_scale.*sqrt(time_scale));
                                                % EXP 2
                                                C_all_T2 = C_all_T2 .* C_all_T2 ./ (time_scale);
                                                
                                                C_all_T2 = C_all' + C_all_T2;
                                                
                                                C_all_T = exp(-C_all_T2);
                                                C_all_T = C_all_T .* eps999 + eps100;
                                                
                                                % only continue if atleast 'min_num_profs' profiles still with us
                                                if numel(ii_dist) > min_num_profs
                                                    
                                                    
                                                    % -------------------------------------------------------------------------
                                                    % -------------------------------------------------------------------------
                                                    % -------------------------------------------------------------------------
                                                    % -------------------------------------------------------------------------
                                                    % COMPUTE MIXED LAYED DEPTH
                                                    % -------------------------------------------------------------------------
                                                    % -------------------------------------------------------------------------
                                                    % -------------------------------------------------------------------------
                                                    % -------------------------------------------------------------------------
                                                    
                                                    [~, kk ] = sort(C_all_T,'descend');
                                                    
                                                    
                                                    i2_dist = i2_dist(kk);
                                                    %ii_dist = ii_dist(kk);
                                                    
                                                    
                                                    % all ready for computation, do mixed layer and underlying isopycnals:
						    sum_str='MLdens(i2_dist) + MLpres(i2_dist) +  MLsalt(i2_dist) + MLtemp(i2_dist)';
						    for ifield=1:length(AdditionalFields);
							field=AdditionalFields{ifield};
							eval(['sum_str= [sum_str  ''+ ML' field '(i2_dist)''];']);%=======================	MODIF
						    end
						    eval(['index_isfin=find(isfinite(' sum_str '));']);                                            
						    %index_isfin=find(isfinite(MLdens(i2_dist) + MLatan030(i2_dist) + MLatan0200(i2_dist) + MLpres(i2_dist) +  MLsalt(i2_dist) + MLtemp(i2_dist)));% & thirdofyear');%=============MODIF
                                                    
                                                    ll=numel(index_isfin);
                                                    if ll > min_num_profs
                                                        
                                                        % get number of data points
                                                        % grab closest (weightspace) data profile + random, discard others
                                                        if ll>max_num_profs
                                                            if ll > rnd_within
                                                                ll = rnd_within;
                                                            end
                                                            tmp = seed_rndprofs(seed_rndprofs<=ll);
                                                            index_isfin = index_isfin([1:numsortedprofs tmp(1:num_of_rand_profiles)]);
                                                        end
                                                        
                                                        % compute IQR on MLD density
                                                        index_goodIQR  = findgoodIQR(MLdens(i2_dist(index_isfin)),iqr_mult);
                                                        index_goodIQR2 = findgoodIQR(MLpres(i2_dist(index_isfin(index_goodIQR))),iqr_mult);
                                                        index_goodIQR  = index_goodIQR(index_goodIQR2);
                                                        
                                                        %  this data is good and can be used
                                                        i2_final = i2_dist(index_isfin(index_goodIQR));
                                                        % ii_final = ii_dist(index_isfin(index_goodIQR));
                                                        i3_final = kk(index_isfin(index_goodIQR));
                                                        % find subset of data for computation:
                                                        
                                                        ll = numel(i2_final);
                                                        
                                                        % set temporary variables
                                                        CaTTMP = C_all_T(i3_final);
                                                        denTMP = MLdens(i2_final);
                                                        preTMP = MLpres(i2_final);
                                                        salTMP = MLsalt(i2_final);
                                                        temTMP = MLtemp(i2_final);
                                                        yeaTMP = dyrMAP(i3_final);
   					   		for ifield=1:length(AdditionalFields);
								field=AdditionalFields{ifield};
								eval([field 'TMP = ML' field '(i2_final);']); %=======================	MODIF
    					    		end

                                                        
                                                        mi = (CaTTMP * denTMP) / sum(CaTTMP);
                                                        
                                                        
                                                        % take 100 closest or less
                                                        minleng = min([150; ll]);
                                                        % to compute weighted STD of density
                                                        wfkt = sqrt(var( denTMP(1:minleng),CaTTMP(1:minleng))) .* wfkt_multiplier ;
                                                        ML_Drho(tstep) = wfkt;
                                                        wfkt_sq = wfkt*wfkt;
                                                        % if STD lower than 0.05 use 0.05
                                                        %if wfkt < 0.05
                                                        %    wfkt = 0.05;
                                                        %elseif wfkt > 1.05
                                                        %    wfkt = 1.05;
                                                        %end
                                                        
                                                        denTMP2 = abs(denTMP - mi);
                                                        % add density weighting to weight vector and finally take the negative exponential
                                                        % EXP 1.5
                                                        %CCC = exp(-(C_all_T2(i3_final) + ((abs(denTMP).*sqrt(abs(denTMP))./(wfkt.*sqrt(wfkt)))')));
                                                        % EXP 2
                                                        CCC = exp(-(C_all_T2(i3_final) + (((denTMP2.*denTMP2)./wfkt_sq)')));
                                                        CCC = CCC .* eps999 + eps100;
                                                        




                                                        % compute weighted mean
                                                        tmpsw = sum(CCC);
                                                        mi = (CCC * denTMP) /tmpsw;
                                                        mp = (CCC * preTMP) /tmpsw;
                                                        ms = (CCC * salTMP) /tmpsw;
                                                        mt = (CCC * temTMP) /tmpsw;
                                                        my = (CCC * yeaTMP) /tmpsw;
  					   		for ifield=1:length(AdditionalFields);
								field=AdditionalFields{ifield};
								eval(['mi' field ' = (CCC * ' field 'TMP) /tmpsw; ']); %=======================	MODIF
    					    		end
                                                        
                                                        % remove  mean
                                                        denTMP = (denTMP - mi);
                                                        preTMP = (preTMP - mp);
                                                        salTMP = (salTMP - ms);
                                                        temTMP = (temTMP - mt);
                                                        yeaTMP = (yeaTMP - my);
  					   		for ifield=1:length(AdditionalFields);
								field=AdditionalFields{ifield};
								eval([field 'TMP = (' field 'TMP - mi' field ');']); %=======================	MODIF
    					    		end

                                                        % write weighted mean of dens and pres to file
                                                        ML_dens_gr_wmean(tstep) = mi;
                                                        ML_pres_gr_wmean(tstep) = mp;
                                                        ML_salt_gr_wmean(tstep) = ms;
                                                        ML_temp_gr_wmean(tstep) = mt;
  					   		for ifield=1:length(AdditionalFields);
								field=AdditionalFields{ifield};
								eval(['ML_' field '_gr_wmean(tstep) = mi' field ';']); %=======================	MODIF
    					    		end
                                                        
                                                        CCC = exp(-(C_all_T2(i3_final) + (((denTMP.*denTMP)./wfkt_sq)')));
                                                        CCC = CCC .* eps999 + eps100;
                                                        
                                                        tmpsw = sum(CCC);
                                                        
                                                        % save summed weight to file
                                                        ML_weight(tstep) = tmpsw;
                                                        
                                                        %   E_Z(1:ll,1:ll) = exp(-(E_all_T(i3_final,i3_final) + (abs(double(repmat(double(denTMP)',ll,1)-repmat(double(denTMP),1,ll))).^expweight)./wfkt.^expweight));
                                                        %E_Z(1:max_num_profs+1:end) = E_Z(1:max_num_profs+1:end) + signal_to_noise;
                                                        
                                                        % add density covariance to covariance matrix and finally take the negative exponential (IN C and parallel)
                                                        % further more add diagonal signal/noise
                                                        
                                                        % EXP 1.5
                                                        %E_Z = plus_uptri_and_exp_mex_par2(E_all_T(i3_final,i3_final),denTMP,wfkt.*sqrt(wfkt),dec_dyrTMP(i3_final));
                                                        % EXP 2
                                                        %disp('plus_uptri_and_exp_mex_par2_gauss');
                                                        E_Z = plus_uptri_and_exp_mex_par2_gauss(E_all_T(i3_final,i3_final),denTMP,wfkt_sq,dec_dyrTMP(i3_final));
                                                        E_Z = E_Z .* eps999 + eps100;
                                                        
                                                        % do OA with anomalies and add mean
							an_tosolve='denTMP preTMP CCC'' ones(numel(i2_final),1) salTMP temTMP yeaTMP';
							wm_tosolve='mi mp 0 0  ms mt my';
  					   		for ifield=1:length(AdditionalFields);
								field=AdditionalFields{ifield};
								eval(['an_tosolve= [an_tosolve '' ' field 'TMP''];']); %=======================	MODIF
								eval(['wm_tosolve= [wm_tosolve  '' mi' field '''];']); %=======================	MODIF
    					    		end

							eval(['b_ML = CCC*linsolve(E_Z,[' an_tosolve '],opts) + [' wm_tosolve '];']); %=======================	MODIF
                                                        %b_ML = CCC*linsolve(E_Z,[denTMP preTMP CCC' ones(numel(i2_final),1) salTMP temTMP yeaTMP atan30TMP atan200TMP],opts) + [mi mp 0 0  ms mt my mi30 mi200];
							
							
                                                        % keyboard
                                                        ML_salt_gr(tstep) = b_ML(5);
                                                        ML_temp_gr(tstep) = b_ML(6);
                                                        ML_dens_gr(tstep) = b_ML(1);
                                                        ML_pres_gr(tstep) = b_ML(2);
                                                        ML_year_gr(tstep) = b_ML(7);
  					   		for ifield=1:length(AdditionalFields);
								field=AdditionalFields{ifield};
								eval(['ML_' field '_gr(tstep) = b_ML( ' int2str(7+ifield) ');']); %=======================	MODIF
    					    		end                                                        
                                                        %save all errors into grid point TIME matrix
                                                        ML_area(tstep) = b_ML(4);
                                                        ML_error(tstep) = b_ML(3);
                                                        
                                                        ML_num_used(tstep) = ll;
                                                        ML_r_dist(tstep) = (CaTTMP * dist_vec_complexTMP(i3_final)) *110;
                                                        
                                                        old_depth = ML_pres_gr(tstep);
                                                        %keyboard
                                                    else
                                                        old_depth = 0;
                                                    end
                                                    
                                                    
                                                    if isfinite(ML_dens_gr(tstep) + ML_pres_gr(tstep)) && ML_dens_gr(tstep) > 2
                                                        sg_ind_ML = find(sg0_pleth < ML_dens_gr(tstep),1,'last'); % start at this isopycnal IN ML
                                                        sg_ind = sg_ind_ML;
                                                        sgindmax = 1;
                                                    else
                                                        sg_ind_ML = 1;
                                                        sg_ind=1;
                                                        sgindmax = 1;
                                                    end
                                                    
                                                    
                                                    % -------------------------------------------------------------------------
                                                    % -------------------------------------------------------------------------
                                                    % -------------------------------------------------------------------------
                                                    % -------------------------------------------------------------------------
                                                    % COMPUTE ALL ISOPYCNALS
                                                    % -------------------------------------------------------------------------
                                                    % -------------------------------------------------------------------------
                                                    % -------------------------------------------------------------------------
                                                    % -------------------------------------------------------------------------
                                                    
                                                    if compute_interior
                                                        %take subsample (makes code more readable, despite being slightly slower (0.5%)
                                                        sg_saTMP_r = sg_sa(i2_dist,:);
                                                        sg_prTMP_r = sg_pr(i2_dist,:);
                                                        
                                                        % loop over all desired sigma levels ... up to hundreds of times
                                                        while sg_ind<=numel(sg0_pleth);
                                                            
                                                            % again take subsample of data (more readbable but slower)
                                                            sg_saTMP = sg_saTMP_r(:,sg_ind);
                                                            sg_prTMP = sg_prTMP_r(:,sg_ind);
                                                            
                                                            % last check that we get rid of all NaN, Inf, in raw data ...
                                                            index_isfin=find(isfinite(sg_saTMP + sg_prTMP));
                                                            
                                                            ll = numel(index_isfin);
                                                            % only keep going if still enough data is available
                                                            if ll > min_num_profs
                                                                
                                                                % get number of data points
                                                                % grab closest (weightspace) data profile + random, discard others
                                                                if ll>max_num_profs
                                                                    if ll > rnd_within
                                                                        ll = rnd_within;
                                                                    end
                                                                    tmp = seed_rndprofs(seed_rndprofs<=ll);
                                                                    index_isfin = index_isfin([1:numsortedprofs tmp(1:num_of_rand_profiles)]);
                                                                end
                                                                
                                                                
                                                                % compute IQR on salinity and pressure surface
                                                                disp('findgood IQR')
                                                                [index_goodIQR2]= findgoodIQR(sg_saTMP(index_isfin),iqr_mult);
                                                                [index_goodIQR ] = findgoodIQR(sg_prTMP(index_isfin(index_goodIQR2)),iqr_mult);
                                                                
                                                                index_goodIQR  = index_goodIQR2(index_goodIQR);
                                                                
                                                                %  this data is good and can be used
                                                                i2_final = i2_dist(index_isfin(index_goodIQR));
                                                                i3_final = kk(index_isfin(index_goodIQR));
                                                                
                                                                ll = numel(i2_final);
                                                                
                                                                % save number of profiles available and used for this data point:
                                                                sg_num_used(tstep,sg_ind)=numel(index_goodIQR);
                                                                
                                                                % take subsamples of data to minimize indicing and further make code more readable
                                                                saTMP = sg_sa(i2_final,sg_ind) ;
                                                                thTMP = sg_th(i2_final,sg_ind) ;
                                                                yrTMP = dyrMAP(i3_final);
                                                                prTMP = sg_pr(i2_final,sg_ind) ;
                                                                
                                                                CaTTMP = C_all_T(i3_final);
                                                                
                                                                % compute weighted mean
                                                                tmpsw = sum(CaTTMP);
                                                                mp = (CaTTMP * prTMP) /tmpsw;
                                                                prTMP2 = abs(prTMP - mp);
                                                                
                                                                %compute delta-Z weights
                                                                %minleng = min([150; ll]);
                                                                %wfkt = sqrt(var(prTMP(1:minleng),CaTTMP(1:minleng)));
                                                                wfkt = sqrt( CaTTMP*(prTMP2.*prTMP2)/tmpsw ) .* wfkt_multiplier;
                                                                
                                                                sg_dZ_scale(tstep,sg_ind)=wfkt;
                                                                if wfkt < 10
                                                                    wfkt = 10;
                                                                    wfkt_sq = 100;
                                                                    %elseif wfkt > 500
                                                                    %    wfkt = 500;
                                                                    %    wfkt_sq = 250000;
                                                                else
                                                                    wfkt_sq = wfkt*wfkt;
                                                                end
                                                                
                                                                %keyboard
                                                                
                                                                %if weighted mean of depth is shallower than last depth- correct:
                                                                if mp < old_depth && sg_ind > sg_ind_ML
                                                                    mp = old_depth;
                                                                    prTMP2 = abs(prTMP - mp);
                                                                elseif  sg_ind == sg_ind_ML && mp > old_depth
                                                                    mp = max([old_depth-5;5]) ;
                                                                    prTMP2 = abs(prTMP - mp);
                                                                end
                                                                
                                                                
                                                                % add dZ covariance to weightvector
                                                                % EXP 1.5
                                                                %CCC = exp(-(C_all_T2(i3_final) + ((abs(prTMP).*sqrt(abs(prTMP))./(wfkt.*sqrt(wfkt)))')));
                                                                % EXP 2
                                                                CCC = exp(-(C_all_T2(i3_final) + (((prTMP2.*prTMP2)./wfkt_sq)')));
                                                                CCC = CCC .* eps999 + eps100;
                                                                
                                                                tmpsw = sum(CCC);
                                                                
                                                                if tmpsw > min_allowed_weight
                                                                    
                                                                    mp = (CCC * prTMP) /tmpsw;
                                                                    %if weighted mean of depth is shallower than last depth- correct:
                                                                    if mp < old_depth && sg_ind > sg_ind_ML
                                                                        mp = old_depth;
                                                                    elseif  sg_ind == sg_ind_ML && mp > old_depth
                                                                        mp = max([old_depth-5;5]) ;
                                                                    end
                                                                    
                                                                    ms = (CCC * saTMP) /tmpsw;
                                                                    mt = (CCC * thTMP) /tmpsw;
                                                                    my = (CCC * yrTMP) /tmpsw;
                                                                    
                                                                    
                                                                    % remove mean from raw data
                                                                    prTMP = prTMP - mp;
                                                                    saTMP = saTMP - ms;
                                                                    thTMP = thTMP - mt;
                                                                    yrTMP = yrTMP - my;
                                                                    
                                                                    CCC = exp(-(C_all_T2(i3_final) + (((prTMP.*prTMP)./wfkt_sq)')));
                                                                    CCC = CCC .* eps999 + eps100;
                                                                    
                                                                    tmpsw = sum(CCC);
                                                                    
                                                                    if tmpsw > min_allowed_weight
                                                                        
                                                                        % save weighted means
                                                                        sg_pr_gr_wmean(tstep,sg_ind) = mp;
                                                                        sg_sa_gr_wmean(tstep,sg_ind) = ms;
                                                                        sg_th_gr_wmean(tstep,sg_ind) = mt;
                                                                        %sg_yr_gr_wmean(tstep,sg_ind) = my;
                                                                        
                                                                        % sum weights and save them
                                                                        sg_weight(tstep,sg_ind)=tmpsw;
                                                                        
                                                                        % save r_dist% this is distance in cumulative DISTANCE km
                                                                        sg_r_dist(tstep,sg_ind)=  (CCC * dist_vec_complexTMP(i3_final)) * 110;
                                                                        
                                                                        % create dZ covariance matrix, add to existing covariance matrix, take negative exponential, add noise to diagonal
                                                                        % EXP 1.5
                                                                        %E_Z = plus_uptri_and_exp_mex_par2(E_all_T(i3_final,i3_final),prTMP,wfkt.*sqrt(wfkt),dec_dyrTMP(i3_final));
                                                                        % EXP 2
                                                                        %disp('plus_uptri_and_exp_mex_par2_gauss')
                                                                        E_Z = plus_uptri_and_exp_mex_par2_gauss(E_all_T(i3_final,i3_final),prTMP,wfkt_sq,dec_dyrTMP(i3_final));
                                                                        
                                                                        E_Z = E_Z .* eps999 + eps100;
                                                                        
                                                                        % do the Optimal Interpolation and add back the weighted means.
                                                                        b=    CCC ...
                                                                            * linsolve(E_Z,  ...
                                                                            [prTMP saTMP thTMP yrTMP CCC' ones(ll,1)],...
                                                                            opts) ...
                                                                            + [mp ms mt my 0 0];
                                                                        %keyboard
                                                                        
                                                                        % save OI outputs
                                                                        sg_sa_gr(tstep,sg_ind)=b(2);
                                                                        sg_pr_gr(tstep,sg_ind)=b(1);
                                                                        sg_th_gr(tstep,sg_ind)=b(3);
                                                                        sg_yr_gr(tstep,sg_ind)=b(4);
                                                                        
                                                                        % check if we hit the bottom yet and if we are deeper than before
                                                                        if b(1) > old_depth && sg_ind > sg_ind_ML
                                                                            old_depth = b(1);
                                                                            % DO SO IN CLEANUP !!!!
                                                                            %if  old_depth >= local_depth
                                                                            %    % set to bottom pressure if we ended up to deep
                                                                            %    sg_pr_gr(tstep,sg_ind) = local_depth;
                                                                            %end
                                                                        end
                                                                        
                                                                        %save all errors into grid point TIME-DENS matrix
                                                                        sg_area(tstep,sg_ind)= b(6); %se_b(1,2);
                                                                        sg_error(tstep,sg_ind)= b(5); %se_b(1,1);
                                                                    end
                                                                end
                                                                
                                                            end
                                                            
                                                            if sg_ind > sgindmax
                                                                sgindmax = sg_ind;
                                                            end
                                                            % if we hit the bottom or were not able to compute this isopycnal go on to next month/grid point
                                                            if isnan(sg_pr_gr(tstep,sg_ind)) || old_depth >= local_depth
                                                                
                                                                break
                                                            end
                                                            
                                                            % go on to next denser isopycnal
                                                            sg_ind = sg_ind + 1;
                                                            
                                                        end
                                                    end
                                                end
                                            end
                                            
                                            % save all data (all months all isopycnals) from this gridpoint to netcdf file
                                            if compute_interior
                                                sg_w = numel(sg_ind_ML:sgindmax);
                                                if sg_w > 0
                                                    netcdf.putVar(ncid_CLIM,data_temp_ID,[sg_ind_ML-1 0 istep-1 jstep-1],[sg_w 12 1 1],single(sg_th_gr(:,sg_ind_ML:sgindmax)'));
                                                    netcdf.putVar(ncid_CLIM,data_salt_ID,[sg_ind_ML-1 0 istep-1 jstep-1],[sg_w 12 1 1],single(sg_sa_gr(:,sg_ind_ML:sgindmax)'));
                                                    netcdf.putVar(ncid_CLIM,data_pres_ID,[sg_ind_ML-1 0 istep-1 jstep-1],[sg_w 12 1 1],single(sg_pr_gr(:,sg_ind_ML:sgindmax)'));
                                                    
                                                    netcdf.putVar(ncid_CLIM,data_temp_wm_ID,[sg_ind_ML-1 0 istep-1 jstep-1],[sg_w 12 1 1],single(sg_th_gr_wmean(:,sg_ind_ML:sgindmax)'));
                                                    netcdf.putVar(ncid_CLIM,data_salt_wm_ID,[sg_ind_ML-1 0 istep-1 jstep-1],[sg_w 12 1 1],single(sg_sa_gr_wmean(:,sg_ind_ML:sgindmax)'));
                                                    netcdf.putVar(ncid_CLIM,data_pres_wm_ID,[sg_ind_ML-1 0 istep-1 jstep-1],[sg_w 12 1 1],single(sg_pr_gr_wmean(:,sg_ind_ML:sgindmax)'));
                                                    
                                                    
                                                    netcdf.putVar(ncid_CLIM,data_year_ID,[sg_ind_ML-1 0 istep-1 jstep-1],[sg_w 12 1 1],int16(round((sg_yr_gr(:,sg_ind_ML:sgindmax)'-2000).*1000)));
                                                    %  netcdf.putVar(ncid_CLIM,data_year_wm_ID,[sg_ind_ML-1 0 istep-1 jstep-1],[sg_w 12 1 1],int16(round((sg_yr_gr_wmean(:,sg_ind_ML:sgindmax)'-2000).*1000)));
                                                    
                                                    netcdf.putVar(ncid_CLIM,data_nwgt_ID,[sg_ind_ML-1 0 istep-1 jstep-1],[sg_w 12 1 1],single(sg_weight(:,sg_ind_ML:sgindmax)'));
                                                    netcdf.putVar(ncid_CLIM,data_dist_M_ID,[sg_ind_ML-1 0 istep-1 jstep-1],[sg_w 12 1 1],int16(imag(sg_r_dist(:,sg_ind_ML:sgindmax))'));
                                                    netcdf.putVar(ncid_CLIM,data_dist_Z_ID,[sg_ind_ML-1 0 istep-1 jstep-1],[sg_w 12 1 1],int16(real(sg_r_dist(:,sg_ind_ML:sgindmax))'));
                                                    netcdf.putVar(ncid_CLIM,data_error_ID,[sg_ind_ML-1 0 istep-1 jstep-1],[sg_w 12 1 1],single(sg_error(:,sg_ind_ML:sgindmax)'));
                                                    netcdf.putVar(ncid_CLIM,data_area_ID,[sg_ind_ML-1 0 istep-1 jstep-1],[sg_w 12 1 1],single(sg_area(:,sg_ind_ML:sgindmax)'));
                                                    netcdf.putVar(ncid_CLIM,data_num_data_ID,[sg_ind_ML-1 0 istep-1 jstep-1],[sg_w 12 1 1],int16(round(sg_num_used(:,sg_ind_ML:sgindmax))'));
                                                    netcdf.putVar(ncid_CLIM,data_z_weight_ID,[sg_ind_ML-1 0 istep-1 jstep-1],[sg_w 12 1 1],int16(round(sg_dZ_scale(:,sg_ind_ML:sgindmax))'));
                                                    
                                                    % clean temporary data and set to NaN
                                                    sg_sa_gr = dummy_T_sig;
                                                    sg_pr_gr = dummy_T_sig;
                                                    sg_th_gr = dummy_T_sig;
                                                    sg_sa_gr_wmean = dummy_T_sig;
                                                    sg_pr_gr_wmean = dummy_T_sig;
                                                    sg_th_gr_wmean = dummy_T_sig;
                                                    sg_yr_gr = dummy_T_sig;
                                                    %sg_yr_gr_wmean = dummy_T_sig;
                                                    sg_weight = dummy_T_sig;
                                                    sg_r_dist = complex(dummy_T_sig);
                                                    sg_num_used = dummy_T_sig;
                                                    sg_area = dummy_T_sig;
                                                    sg_error= dummy_T_sig;
                                                    sg_dZ_scale = dummy_T_sig;
                                                end
                                            end
                                            netcdf.putVar(ncid_CLIM,ML_dens_ID,[0 istep-1 jstep-1],[12 1 1],single(ML_dens_gr));
                                            netcdf.putVar(ncid_CLIM,ML_pres_ID,[0 istep-1 jstep-1],[12 1 1],single(ML_pres_gr));
                                            netcdf.putVar(ncid_CLIM,ML_dens_wm_ID,[0 istep-1 jstep-1],[12 1 1],single(ML_dens_gr_wmean));
                                            netcdf.putVar(ncid_CLIM,ML_pres_wm_ID,[0 istep-1 jstep-1],[12 1 1],single(ML_pres_gr_wmean));
                                            netcdf.putVar(ncid_CLIM,ML_salt_wm_ID,[0 istep-1 jstep-1],[12 1 1],single(ML_salt_gr_wmean));
                                            netcdf.putVar(ncid_CLIM,ML_temp_wm_ID,[0 istep-1 jstep-1],[12 1 1],single(ML_temp_gr_wmean));
                                            netcdf.putVar(ncid_CLIM,ML_salt_ID,[0 istep-1 jstep-1],[12 1 1],single(ML_salt_gr));
                                            netcdf.putVar(ncid_CLIM,ML_temp_ID,[0 istep-1 jstep-1],[12 1 1],single(ML_temp_gr));
                                            netcdf.putVar(ncid_CLIM,ML_nwgt_ID,[0 istep-1 jstep-1],[12 1 1],single(ML_weight));
 					    for ifield=1:length(AdditionalFields);
					    	field=AdditionalFields{ifield};
					    	eval(['netcdf.putVar(ncid_CLIM,ML_' field '_ID,[0 istep-1 jstep-1],[12 1 1],single(ML_' field '_gr));']);
						eval(['netcdf.putVar(ncid_CLIM,ML_' field '_wm_ID,[0 istep-1 jstep-1],[12 1 1],single(ML_' field '_gr_wmean));']);	
    					    end                                                        


                                            
                                            netcdf.putVar(ncid_CLIM,ML_year_ID,[0 istep-1 jstep-1],[12 1 1],int16((ML_year_gr-2000).*1000));
                                            %netcdf.putVar(ncid_CLIM,ML_year_wm_ID,[0 istep-1 jstep-1],[12 1 1],int16((ML_year_gr_wmean-2000).*1000));
                                            
                                            netcdf.putVar(ncid_CLIM,ML_dist_M_ID,[0 istep-1 jstep-1],[12 1 1],int16(imag(ML_r_dist)));
                                            netcdf.putVar(ncid_CLIM,ML_dist_Z_ID,[0 istep-1 jstep-1],[12 1 1],int16(real(ML_r_dist)));
                                            netcdf.putVar(ncid_CLIM,ML_error_ID,[0 istep-1 jstep-1],[12 1 1],single(ML_error));
                                            netcdf.putVar(ncid_CLIM,ML_area_ID,[0 istep-1 jstep-1],[12 1 1],single(ML_area));
                                            netcdf.putVar(ncid_CLIM,ML_num_data_ID,[0 istep-1 jstep-1],[12 1 1],int16(ML_num_used));
                                            netcdf.putVar(ncid_CLIM,ML_drho_ID,[0 istep-1 jstep-1],[12 1 1],single(ML_Drho));
                                            
                                            ML_dens_gr = dummy_T;
                                            ML_pres_gr = dummy_T;
                                            ML_year_gr = dummy_T;
                                            ML_temp_gr = dummy_T;
                                            ML_salt_gr = dummy_T;
                                            ML_dens_gr_wmean = dummy_T;
                                            ML_salt_gr_wmean = dummy_T;
                                            ML_temp_gr_wmean = dummy_T;
                                            ML_pres_gr_wmean = dummy_T;
                                            ML_weight = dummy_T;
                                            ML_r_dist = complex(dummy_T);
                                            ML_num_used = dummy_T;
                                            ML_area = dummy_T;
                                            ML_error= dummy_T;
                                            ML_Drho = dummy_T;
 					    for ifield=1:length(AdditionalFields);
					    	field=AdditionalFields{ifield};
					    	eval(['ML_' field '_gr_wmean = dummy_T;']);%=======================	MODIF
					    	eval(['ML_' field '_gr_wmean = dummy_T;']);
    					    end                                                        
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        % print out latitude processed and time used since start ...
        fprintf(' \n ')
        disp([' calculated latitude: ', ...
            num2str(y(istep,jstep)), '?N  - running since ',...
            datestr(now-tstart,'HH:MM:SS'), ' HH:MM:SS.'])
        
    end
    toc
end

%%
disp(' Closing all open file handles ... ')
% close all handles
netcdf.close(ncid_FM);
netcdf.close(ncid_CLIM);
if useargo
    netcdf.close(ncid_AFD);
end
netcdf.close(ncid_CTD);

disp(' All done! Have a nice day ... ')




% =-=-=-=-=-=-==-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-==-=-
% =-=-=-=-=-=-==-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-==-=-
%
%           SUB FUNCTIONS
%
% =-=-=-=-=-=-==-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-==-=-
% =-=-=-=-=-=-==-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-==-=-



