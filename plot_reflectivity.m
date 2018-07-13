


settings.location='Jacobshavn';
save_en=1;

load(['/cresis/snfs1/scratch/manjish/new_jacobshavn/results/reflectivity_med_att3.mat'])


  %   depth_G_r_corrected = depth_G_r_corrected / 1000;
  %
  %   reference_power = nanmean(ice_bed_power_G_r_corrected(isfinite(ice_bed_power_G_r_corrected)));
  %
  %   %reference_power=25;
  %   %relative_ice_bed_power_G_r_corrected = lp(ice_bed_power_G_r_corrected/reference_power);
  %   relative_ice_bed_power_G_r_corrected = (ice_bed_power_G_r_corrected)-(reference_power);
  %
  %   relative_reflectivty_using_constant_attn = relative_ice_bed_power_G_r_corrected + constant_attenuation ;
  %   relative_reflectivty_using_variable_attn = relative_ice_bed_power_G_r_corrected + variable_attenuation ;
  
 % geotiff_fn='/cresis/snfs1/dataproducts/GIS_data/greenland/Jakobshavn/10300100069AA400_ortho_11b.tif';
  %geotiff_fn='/cresis/snfs1/dataproducts/GIS_data/greenland/plummer_jakobshavn/jak_grid.tif';
  %geotiff_fn='/users/manjish/maps/peterman2.tif';
 % geotiff_fn='/users/manjish/maps/peterman_bed.tif';
 %geotiff_fn='/users/manjish/maps/STILL_jakobshavn2000.tif';
  geotiff_fn = '/cresis/snfs1/dataproducts/GIS_data/greenland/Landsat-7/Greenland_natural.tif';
 % geotiff_fn='/users/manjish/maps/peterman_flow_map2.tif';
% geotiff_fn='/users/manjish/maps/jacob_bed2.tif';

  proj = geotiffinfo(geotiff_fn);
  %proj = geotiffinfo('X:\GIS_data\antarctica\Landsat-7\Antarctica_LIMA_480m.tif');
  
  [A CMAP R]= geotiffread(geotiff_fn);
  
  %% Reflectivity using constant na
  
  figure(1)
  mapshow(rgb2gray(A),CMAP/1e3);
  %plot_geotiff(geotiff_fn);
  xlabel('X (km)');
  ylabel('Y (km)');
  if  strcmp(settings.location,'Peterman')
    xlim([-350 -50]);
    ylim([-1250 -900]);
  else
    xlim([-250 150]);
    ylim([-2450 -2100]);
  end
  
  hold on
  clear gps.x gps.y
  [gps.x,gps.y] = projfwd(proj,out.Latitude,out.Longitude);
  
  gps.x = gps.x / 1000;
  gps.y = gps.y / 1000;
  figure(1)
  hold on;
  
  scatter(gps.x,gps.y,20,out.Refl_const,'fill')
  caxis([-15 15])
  colorbar;
  title('Reflectivity using constant na with att method 3')
  
  
  if save_en
    if strcmp(settings.location,'Peterman')
      save_path=['/cresis/snfs1/scratch/manjish/new_peterman/img_results/att3_bed/const_refl'];
    else
      save_path=['/cresis/snfs1/scratch/manjish/new_jacobshavn/img_results/att3/const_reflectivity'];
    end
    [save_dir] =fileparts(save_path);
    if ~exist(save_dir,'dir')
      
      mkdir(save_dir);
    end
    saveas(figure(1),save_path,'jpg')
    % geotiffwrite(save_path,A,R);
  end
  
  %Histogram
  figure(2), hist(out.Refl_const,30);
  title('Reflectivity using constant na ')
  
  if save_en
    if strcmp(settings.location,'Peterman')
      save_path=['/cresis/snfs1/scratch/manjish/new_peterman/img_results/att3_bed/const_refl_hist'];
    else
      save_path=['/cresis/snfs1/scratch/manjish/new_jacobshavn/img_results/att3/const_reflectivity_hist'];
    end
    [save_dir] =fileparts(save_path);
    if ~exist(save_dir,'dir')
      
      mkdir(save_dir);
    end
    saveas(figure(2),save_path,'jpg')
  end
  
  %% Reflectivity using variable attenuation
  
  clear idx;
  [~,idx]=find(isnan(out.Refl_var));
   lat= out.Latitude;
  lon=out.Longitude;
  Refl_var=out.Refl_var;
  lat(idx)=[];
  lon(idx)=[];
  Refl_var(idx)=[];
  
  figure(3)
  mapshow(rgb2gray(A),CMAP/1e3);
  xlabel('X (km)');
  ylabel('Y (km)');
  if  strcmp(settings.location,'Peterman')
    xlim([-350 -50]);
    ylim([-1250 -900]);
  else
    xlim([-250 150]);
    ylim([-2450 -2100]);
  end
  hold on
  clear gps.x gps.y
  [gps.x,gps.y] = projfwd(proj,lat,lon);
  
  gps.x = gps.x / 1000;
  gps.y = gps.y / 1000;
  hold on;
  
  scatter(gps.x,gps.y,20,Refl_var,'fill')
  caxis([-15 15])
  colorbar;
  title('Reflectivity using variable na ')
  
  if save_en
    if strcmp(settings.location,'Peterman')
      save_path=['/cresis/snfs1/scratch/manjish/new_peterman//img_results/att3_bed/var_refl'];
    else
      save_path=['/cresis/snfs1/scratch/manjish/new_jacobshavn/img_results/att3/var_reflectivity_median_att3'];
    end
    [save_dir] =fileparts(save_path);
    if ~exist(save_dir,'dir')
      
      mkdir(save_dir);
    end
    saveas(figure(3),save_path,'jpg')
  end
  
  %Histogram
  figure(4), hist(out.Refl_var,30);
  title('Reflectivity using variable na ')
  
  
  %  keyboard;
  %% Plot Total Constant Attn
  figure(5)
  mapshow(rgb2gray(A),CMAP/1e3);
  xlabel('X (km)');
  ylabel('Y (km)');
  if  strcmp(settings.location,'Peterman')
    xlim([-350 -50]);
    ylim([-1250 -900]);
  else
    xlim([-250 150]);
    ylim([-2450 -2100]);
  end
  hold on
  clear gps.x gps.y
  [gps.x,gps.y] = projfwd(proj,out.Latitude,out.Longitude);
  gps.x = gps.x / 1000;
  gps.y = gps.y / 1000;
  hold on;
  
  scatter(gps.x,gps.y,20,out.const_attenuation,'fill')
  %caxis([-15 15])
  colorbar;
  title('Total Constant Attenuation')
  
  
  
  if save_en
    if strcmp(settings.location,'Peterman')
      save_path=['/cresis/snfs1/scratch/manjish/new_peterman/img_results/att3_bed/const_att'];
    else
      save_path=['/cresis/snfs1/scratch/manjish/new_jacobshavn/img_results/att3/const_att_median'];
    end
    [save_dir] =fileparts(save_path);
    if ~exist(save_dir,'dir')
      
      mkdir(save_dir);
    end
    saveas(figure(5),save_path,'jpg')
  end
  
  
  %% Plot Total Variable Attn
  figure(6)
  mapshow(rgb2gray(A),CMAP/1e3);
  xlabel('X (km)');
  ylabel('Y (km)');
  if  strcmp(settings.location,'Peterman')
    xlim([-350 -50]);
    ylim([-1250 -900]);
  else
    xlim([-250 150]);
    ylim([-2450 -2100]);  end
  hold on
  clear gps.x gps.y
  [gps.x,gps.y] = projfwd(proj,out.Latitude,out.Longitude);
  gps.x = gps.x / 1000;
  gps.y = gps.y / 1000;
  hold on;
  
  scatter(gps.x,gps.y,20,out.var_attenuation,'fill')
  %caxis([-15 15])
  colorbar;
  title('Total Variable Attenuation')
  
  if save_en
    if strcmp(settings.location,'Peterman')
      save_path=['/cresis/snfs1/scratch/manjish/new_peterman//img_results/att3_bed/var_att'];
    else
      save_path=['/cresis/snfs1/scratch/manjish/new_jacobshavn/img_results/att3/var_att_median_att3'];
    end
    [save_dir] =fileparts(save_path);
    if ~exist(save_dir,'dir')
      
      mkdir(save_dir);
    end
    saveas(figure(6),save_path,'jpg')
  end
  
  
  %% Plot Value of mod Na
 
   clear idx;
  [~,idx]=find(isnan(out.modified_Na));
   lat= out.Latitude;
  lon=out.Longitude;
 modified_Na1= out.modified_Na;
  lat(idx)=[];
  lon(idx)=[];
  modified_Na1(idx)=[];
 
  figure(7)
  mapshow(rgb2gray(A),CMAP/1e3);
  xlabel('X (km)');
  ylabel('Y (km)');
  if  strcmp(settings.location,'Peterman')
    xlim([-350 -50]);
    ylim([-1250 -900]);
  else
    xlim([-250 150]);
    ylim([-2450 -2100]);
  end
  
  hold on
  clear gps.x gps.y
  [gps.x,gps.y] = projfwd(proj,lat,lon);
  gps.x = gps.x / 1000;
  gps.y = gps.y / 1000;
  hold on;
  
  scatter(gps.x,gps.y,20,modified_Na1,'fill')
  %caxis([-15 15])
  colorbar;
  title('Modified Na ')
  
  if save_en
    if strcmp(settings.location,'Peterman')
      save_path=['/cresis/snfs1/scratch/manjish/new_peterman//img_results/att3_bed/mod_NA'];
    else
      save_path=['/cresis/snfs1/scratch/manjish/new_jacobshavn/img_results/att3/mod_Na_median_att3'];
    end
    [save_dir] =fileparts(save_path);
    if ~exist(save_dir,'dir')
      
      mkdir(save_dir);
    end
    saveas(figure(7),save_path,'jpg')
  end

  %% D eptjh
   %% Plot Total Constant Attn
   
    clear idx;
  [~,idx]=find(isnan(out.Depth));
   lat= out.Latitude;
  lon=out.Longitude;
    dpth= out.Depth;
  lat(idx)=[];
  lon(idx)=[];
  dpth(idx)=[];
   
   
  figure(8)
  mapshow(rgb2gray(A),CMAP/1e3);
  xlabel('X (km)');
  ylabel('Y (km)');
  if  strcmp(settings.location,'Peterman')
    xlim([-350 -50]);
    ylim([-1250 -900]);
  else
    xlim([-250 150]);
    ylim([-2450 -2100]);
  end
  hold on
  clear gps.x gps.y
  [gps.x,gps.y] = projfwd(proj,lat,lon);
  gps.x = gps.x / 1000;
  gps.y = gps.y / 1000;
  hold on;
  
  scatter(gps.x,gps.y,20,dpth,'fill')
  %caxis([-15 15])
  colorbar;
  title('Depth')
  
  
  
  if save_en
    if strcmp(settings.location,'Peterman')
      save_path=['/cresis/snfs1/scratch/manjish/new_peterman//img_results/att3_bed/depth'];
    else
      save_path=['/cresis/snfs1/scratch/manjish/new_jacobshavn/img_results/att3/const_att_median_att3'];
    end
    [save_dir] =fileparts(save_path);
    if ~exist(save_dir,'dir')
      
      mkdir(save_dir);
    end
    saveas(figure(5),save_path,'jpg')
  end
  