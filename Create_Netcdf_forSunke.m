output_fname='Database/Update2017/MLD_Argo_Seals.nc';
OceanBasins={'SO2','SO1'}
Source={'Argo','Seals_noInt'}
disp('Load and Merge...');

rmf={'length_prof','pts_above','pts_below','fit','thrs','grad','T','S'};
Merge_1d={'NT15','NS15','SA','CT','lon','lat','date','P','gap'};
Rename_1d={'NT15','NS15','mlSA','mlCT','lon','lat','date','mlP','gap'};
Merge_2d={'SA','CT','P'};
clear Merge
iloop=0;
for ib=length(OceanBasins):-1:1
 for is=length(Source):-1:1
	disp(['Database/Update2017/MLD_' OceanBasins{ib} '_' Source{is} '.mat'])
	iloop=iloop+1;    	
	load(['Database/Update2017/MLD_' OceanBasins{ib} '_' Source{is} '.mat']);
	for irm=1:length(rmf)
		eval(['MLD=rmfield(MLD,''' rmf{irm} ''');']);
	end	
	iok=find(MLD.perc2a2<0.25);

	if iloop==1;
		for i=1:length(Merge_1d);
			field=Merge_1d{i};
			fieldrn=Rename_1d{i};
			eval(['Merge.' fieldrn '=MLD.' field '(iok);']);				
		end
	else
		for i=1:length(Merge_1d);
			field=Merge_1d{i};
			fieldrn=Rename_1d{i};
			eval(['Merge.' fieldrn '=[Merge.' fieldrn '; MLD.' field '(iok)];']);				
		end
	end
	clear MLD;

	load(['Database/Update2017/MergeSig_' OceanBasins{ib} '_' Source{is} '.mat']);
	if iloop==1;
		Merge.Sigma0=Merge_sig.Sigma0;
		for i=1:length(Merge_2d);
			field=Merge_2d{i};
			eval(['Merge.' field '=Merge_sig.' field '(iok,:);']);				
		end
	else
		for i=1:length(Merge_2d);
			field=Merge_2d{i};
			eval(['Merge.' field '=[Merge.' field '; Merge_sig.' field '(iok,:)];']);				
		end
	end
	clear Merge_sig
 end
end

[nbdata nbdep]=size(Merge.SA);
year=str2num(datestr(Merge.date, 10));
dyr=year+(Merge.date-datenum(year,1,1))./(datenum(year+1,1,1)-datenum(year,1,1));
sigma0=squeeze(nanmean(Merge.Sigma0));
Merge.mlSig0=gsw_sigma0(Merge.mlSA,Merge.mlCT);

AdditionalFields={'NT15','NS15'};

%	CRéATION DU FICHIER Netcdf ENTREE
   ncid_CLIM = netcdf.create(output_fname,'NC_64BIT_OFFSET');

%	Define the dimensions of the variables.
disp('create dim')
    dimid_scalar     = netcdf.defDim(ncid_CLIM,'SCALAR',1);
    dimid_nbdata     = netcdf.defDim(ncid_CLIM,	'NB_DATA',nbdata);
    dimid_level      = netcdf.defDim(ncid_CLIM,'NB_PROF',nbdep);
        
%	Define new variables in file.
disp('create var')
    latID  = netcdf.defVar(ncid_CLIM,'lat','float',[dimid_nbdata]);
    lonID  = netcdf.defVar(ncid_CLIM,'lon','float',[dimid_nbdata]);
    dyrID  = netcdf.defVar(ncid_CLIM,'dyr','float',[dimid_nbdata]);
    mld_presID = netcdf.defVar(ncid_CLIM,'mld_pres','float',[dimid_nbdata]);
    mld_gapID = netcdf.defVar(ncid_CLIM,'mld_gap','float',[dimid_nbdata]);
    mld_saltID = netcdf.defVar(ncid_CLIM,'mld_salt','float',[dimid_nbdata]);
    mld_tempID = netcdf.defVar(ncid_CLIM,'mld_temp','float',[dimid_nbdata]);
    mld_densID = netcdf.defVar(ncid_CLIM,'mld_dens','float',[dimid_nbdata]);
    for ifield=1:length(AdditionalFields);
	field=AdditionalFields{ifield};
	eval(['mld' field 'ID = netcdf.defVar(ncid_CLIM,''' field ''',''float'',[dimid_nbdata]);']);%=============================MODIF
    end

    sigmaID = netcdf.defVar(ncid_CLIM,'sigma','float',[dimid_level]');
    tempID = netcdf.defVar(ncid_CLIM,'temp','float',[dimid_level dimid_nbdata ]); 	
    salID = netcdf.defVar(ncid_CLIM,'sal','float',[dimid_level dimid_nbdata ]);		
    presID = netcdf.defVar(ncid_CLIM,'pres','float',[dimid_level dimid_nbdata ]);	

%	Write netCDF attributes
    disp('create attrib')
    netcdf.putAtt(ncid_CLIM,latID,'_FillValue',NaN('single'))
    netcdf.putAtt(ncid_CLIM,lonID,'_FillValue',NaN('single'))
    netcdf.putAtt(ncid_CLIM,dyrID,'_FillValue',NaN('single'))
    netcdf.putAtt(ncid_CLIM,mld_presID ,'_FillValue',NaN('single'))
    netcdf.putAtt(ncid_CLIM,mld_densID ,'_FillValue',NaN('single'))
    netcdf.putAtt(ncid_CLIM,mld_gapID,'_FillValue',NaN('single'))
    for ifield=1:length(AdditionalFields);
	field=AdditionalFields{ifield};
	eval(['netcdf.putAtt(ncid_CLIM,mld' field 'ID,''_FillValue'',NaN(''single''))']);
    end
    netcdf.putAtt(ncid_CLIM,mld_saltID,'_FillValue',NaN('single'))
    netcdf.putAtt(ncid_CLIM,mld_tempID,'_FillValue',NaN('single'))
    netcdf.putAtt(ncid_CLIM,sigmaID,'_FillValue',NaN('single'))
    netcdf.putAtt(ncid_CLIM,tempID,'_FillValue',NaN('single'))
    netcdf.putAtt(ncid_CLIM,salID,'_FillValue',NaN('single'))
    netcdf.putAtt(ncid_CLIM,presID,'_FillValue',NaN('single'))
   
    
%	Leave define mode and enter data mode to write data.
    netcdf.endDef(ncid_CLIM)
    
%	Write constants and grid: 
    disp('write var')
    netcdf.putVar(ncid_CLIM,latID,Merge.lat);
    netcdf.putVar(ncid_CLIM,lonID,Merge.lon);
    netcdf.putVar(ncid_CLIM,dyrID,dyr);% attention ici on met le vecteur sur lequel on a fait transformation en decimal years soit sur MLD.date soit sur Merge.date

    netcdf.putVar(ncid_CLIM,tempID,Merge.CT); 
    netcdf.putVar(ncid_CLIM,salID,Merge.SA);
    netcdf.putVar(ncid_CLIM,sigmaID,Merge.Sigma0);
    netcdf.putVar(ncid_CLIM,presID,Merge.P);

    netcdf.putVar(ncid_CLIM,mld_saltID,Merge.mlSA); 
    netcdf.putVar(ncid_CLIM,mld_tempID,Merge.mlCT);
    netcdf.putVar(ncid_CLIM,mld_presID,Merge.mlP);
    netcdf.putVar(ncid_CLIM,mld_densID,Merge.mlSig0);

    netcdf.putVar(ncid_CLIM,mld_gapID,Merge.gap);
    for ifield=1:length(AdditionalFields);
	field=AdditionalFields{ifield};
	eval(['netcdf.putVar(ncid_CLIM,mld' field 'ID, Merge.' field ');']);
    end

    disp('write mld_densID')
    disp('tutto è finito! ')                                      
    netcdf.close(ncid_CLIM);

clear all

disp('NOW WORK ON SHIP')

output_fname='Database/Update2017/MLD_Ship.nc';
OceanBasins={'SO2','SO1'}
Source={'Ship'}
disp('Load and Merge...');

rmf={'length_prof','pts_above','pts_below','fit','thrs','grad','T','S'};
Merge_1d={'NT15','NS15','SA','CT','lon','lat','date','P','gap'};
Rename_1d={'NT15','NS15','mlSA','mlCT','lon','lat','date','mlP','gap'};
Merge_2d={'SA','CT','P'};
clear Merge
iloop=0;
for ib=length(OceanBasins):-1:1
 for is=length(Source):-1:1
	disp(['Database/Update2017/MLD_' OceanBasins{ib} '_' Source{is} '.mat'])
	iloop=iloop+1;    	
	load(['Database/Update2017/MLD_' OceanBasins{ib} '_' Source{is} '.mat']);
	for irm=1:length(rmf)
		eval(['MLD=rmfield(MLD,''' rmf{irm} ''');']);
	end	
	iok=find(MLD.perc2a2<0.25);

	if iloop==1;
		for i=1:length(Merge_1d);
			field=Merge_1d{i};
			fieldrn=Rename_1d{i};
			eval(['Merge.' fieldrn '=MLD.' field '(iok);']);				
		end
	else
		for i=1:length(Merge_1d);
			field=Merge_1d{i};
			fieldrn=Rename_1d{i};
			eval(['Merge.' fieldrn '=[Merge.' fieldrn '; MLD.' field '(iok)];']);				
		end
	end
	clear MLD;

	load(['Database/Update2017/MergeSig_' OceanBasins{ib} '_' Source{is} '.mat']);
	if iloop==1;
		Merge.Sigma0=Merge_sig.Sigma0;
		for i=1:length(Merge_2d);
			field=Merge_2d{i};
			eval(['Merge.' field '=Merge_sig.' field '(iok,:);']);				
		end
	else
		for i=1:length(Merge_2d);
			field=Merge_2d{i};
			eval(['Merge.' field '=[Merge.' field '; Merge_sig.' field '(iok,:)];']);				
		end
	end
	clear Merge_sig
 end
end


[nbdata nbdep]=size(Merge.SA);
year=str2num(datestr(Merge.date, 10));
dyr=year+(Merge.date-datenum(year,1,1))./(datenum(year+1,1,1)-datenum(year,1,1));
sigma0=squeeze(nanmean(Merge.Sigma0));
Merge.mlSig0=gsw_sigma0(Merge.mlSA,Merge.mlCT);

AdditionalFields={'NT15','NS15'};


%	CRéATION DU FICHIER Netcdf ENTREE
   ncid_CLIM = netcdf.create(output_fname,'NC_64BIT_OFFSET');

%	Define the dimensions of the variables.
disp('create dim')
    dimid_scalar     = netcdf.defDim(ncid_CLIM,'SCALAR',1);
    dimid_nbdata     = netcdf.defDim(ncid_CLIM,	'NB_DATA',nbdata);
    dimid_level      = netcdf.defDim(ncid_CLIM,'NB_PROF',nbdep);
        
%	Define new variables in file.
disp('create var')
    latID  = netcdf.defVar(ncid_CLIM,'lat','float',[dimid_nbdata]);
    lonID  = netcdf.defVar(ncid_CLIM,'lon','float',[dimid_nbdata]);
    dyrID  = netcdf.defVar(ncid_CLIM,'dyr','float',[dimid_nbdata]);
    mld_presID = netcdf.defVar(ncid_CLIM,'mld_pres','float',[dimid_nbdata]);
    mld_gapID = netcdf.defVar(ncid_CLIM,'mld_gap','float',[dimid_nbdata]);
    mld_saltID = netcdf.defVar(ncid_CLIM,'mld_salt','float',[dimid_nbdata]);
    mld_tempID = netcdf.defVar(ncid_CLIM,'mld_temp','float',[dimid_nbdata]);
    mld_densID = netcdf.defVar(ncid_CLIM,'mld_dens','float',[dimid_nbdata]);
    for ifield=1:length(AdditionalFields);
	field=AdditionalFields{ifield};
	eval(['mld' field 'ID = netcdf.defVar(ncid_CLIM,''' field ''',''float'',[dimid_nbdata]);']);%=============================MODIF
    end

    sigmaID = netcdf.defVar(ncid_CLIM,'sigma','float',[dimid_level]');
    tempID = netcdf.defVar(ncid_CLIM,'temp','float',[dimid_level dimid_nbdata ]); 	
    salID = netcdf.defVar(ncid_CLIM,'sal','float',[dimid_level dimid_nbdata ]);		
    presID = netcdf.defVar(ncid_CLIM,'pres','float',[dimid_level dimid_nbdata ]);	

%	Write netCDF attributes
    disp('create attrib')
    netcdf.putAtt(ncid_CLIM,latID,'_FillValue',NaN('single'))
    netcdf.putAtt(ncid_CLIM,lonID,'_FillValue',NaN('single'))
    netcdf.putAtt(ncid_CLIM,dyrID,'_FillValue',NaN('single'))
    netcdf.putAtt(ncid_CLIM,mld_presID ,'_FillValue',NaN('single'))
    netcdf.putAtt(ncid_CLIM,mld_densID ,'_FillValue',NaN('single'))
    netcdf.putAtt(ncid_CLIM,mld_gapID,'_FillValue',NaN('single'))
    for ifield=1:length(AdditionalFields);
	field=AdditionalFields{ifield};
	eval(['netcdf.putAtt(ncid_CLIM,mld' field 'ID,''_FillValue'',NaN(''single''))']);
    end
    netcdf.putAtt(ncid_CLIM,mld_saltID,'_FillValue',NaN('single'))
    netcdf.putAtt(ncid_CLIM,mld_tempID,'_FillValue',NaN('single'))
    netcdf.putAtt(ncid_CLIM,sigmaID,'_FillValue',NaN('single'))
    netcdf.putAtt(ncid_CLIM,tempID,'_FillValue',NaN('single'))
    netcdf.putAtt(ncid_CLIM,salID,'_FillValue',NaN('single'))
    netcdf.putAtt(ncid_CLIM,presID,'_FillValue',NaN('single'))
   
    
%	Leave define mode and enter data mode to write data.
    netcdf.endDef(ncid_CLIM)
    
%	Write constants and grid: 
    disp('write var')
    netcdf.putVar(ncid_CLIM,latID,Merge.lat);
    netcdf.putVar(ncid_CLIM,lonID,Merge.lon);
    netcdf.putVar(ncid_CLIM,dyrID,dyr);% attention ici on met le vecteur sur lequel on a fait transformation en decimal years soit sur MLD.date soit sur Merge.date

    netcdf.putVar(ncid_CLIM,tempID,Merge.CT); 
    netcdf.putVar(ncid_CLIM,salID,Merge.SA);
    netcdf.putVar(ncid_CLIM,sigmaID,Merge.Sigma0);
    netcdf.putVar(ncid_CLIM,presID,Merge.P);

    netcdf.putVar(ncid_CLIM,mld_saltID,Merge.mlSA); 
    netcdf.putVar(ncid_CLIM,mld_tempID,Merge.mlCT);
    netcdf.putVar(ncid_CLIM,mld_presID,Merge.mlP);
    netcdf.putVar(ncid_CLIM,mld_densID,Merge.mlSig0);

    netcdf.putVar(ncid_CLIM,mld_gapID,Merge.gap);
    for ifield=1:length(AdditionalFields);
	field=AdditionalFields{ifield};
	eval(['netcdf.putVar(ncid_CLIM,mld' field 'ID, Merge.' field ');']);
    end

    disp('write mld_densID')
    disp('tutto è finito! ')                                      
    netcdf.close(ncid_CLIM);

