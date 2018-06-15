


settings.location='Jacobshavn';

load(['/cresis/snfs1/scratch/manjish/new_jacobshavn/reflectivity_median_att3.mat'])


if 1
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
  
  geotiff_fn = '/cresis/snfs1/dataproducts/GIS_data/greenland/Landsat-7/Greenland_natural.tif';
  proj = geotiffinfo(geotiff_fn);
  %proj = geotiffinfo('X:\GIS_data\antarctica\Landsat-7\Antarctica_LIMA_480m.tif');
  
  [A CMAP R]= geotiffread(geotiff_fn);
  
  %% Reflectivity using constant na
  
  figure(1)
  mapshow(rgb2gray(A),CMAP/1e3);
  xlabel('X (km)');
  ylabel('Y (km)');
  if  strcmp(settings.location,'Peterman')
    xlim([-350 -50]);
    ylim([-1250 -900]);
  else
 xlim([-250 -50]);
    ylim([-2400 -2160]); 
  end
  
  hold on
  clear gps.x gps.y
  [gps.x,gps.y] = projfwd(proj,out.Latitude,out.Longitude);
  
  gps.x = gps.x / 1000;
  gps.y = gps.y / 1000;
  hold on;
  
  scatter(gps.x,gps.y,20,out.Refl_const,'fill')
  caxis([-15 15])
  colorbar;
  title('Reflectivity using constant na with att method 3')
  

    
  if strcmp(settings.location,'Peterman')
    save_path=['/cresis/snfs1/scratch/manjish/new_peterman/const_reflectivity_median_att3'];
  else
    save_path=['/cresis/snfs1/scratch/manjish/new_jacobshavn/const_reflectivity_median_att3'];
  end
  [save_dir] =fileparts(save_path);
  if ~exist(save_dir,'dir')
    
    mkdir(save_dir);
  end
  saveas(figure(1),save_path,'jpg')
  
  
  
  %Histogram
  figure(2), hist(out.Refl_const,30);
  title('Reflectivity using constant na ')
  
  if strcmp(settings.location,'Peterman')
    save_path=['/cresis/snfs1/scratch/manjish/new_peterman/const_reflectivity_hist'];
  else
    save_path=['/cresis/snfs1/scratch/manjish/new_jacobshavn/const_reflectivity_hist'];
  end
  [save_dir] =fileparts(save_path);
  if ~exist(save_dir,'dir')
    
    mkdir(save_dir);
  end
  saveas(figure(2),save_path,'jpg')
  
  %% Reflectivity using variable attenuation
  figure(3)
  mapshow(rgb2gray(A),CMAP/1e3);
  xlabel('X (km)');
  ylabel('Y (km)');
  if  strcmp(settings.location,'Peterman')
    xlim([-350 -50]);
    ylim([-1250 -900]);
  else
    xlim([-250 -50]);
    ylim([-2400 -2160]);
  end
  hold on
  clear gps.x gps.y
  [gps.x,gps.y] = projfwd(proj,out.Latitude,out.Longitude);
  
  gps.x = gps.x / 1000;
  gps.y = gps.y / 1000;
  hold on;
  
  scatter(gps.x,gps.y,20,out.Refl_var,'fill')
  caxis([-15 15])
  colorbar;
  title('Reflectivity using variable na ')
  
  if strcmp(settings.location,'Peterman')
    save_path=['/cresis/snfs1/scratch/manjish/new_peterman/var_reflectivity_median_att3'];
  else
    save_path=['/cresis/snfs1/scratch/manjish/new_jacobshavn/var_reflectivity_median_att3'];
  end
  [save_dir] =fileparts(save_path);
  if ~exist(save_dir,'dir')
    
    mkdir(save_dir);
  end
  saveas(figure(3),save_path,'jpg')
  
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
   xlim([-250 -50]);
    ylim([-2400 -2160]);
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
  
  if strcmp(settings.location,'Peterman')
    save_path=['/cresis/snfs1/scratch/manjish/new_peterman/const_att_median_att3'];
  else
    save_path=['/cresis/snfs1/scratch/manjish/new_jacobshavn/const_att_median_att3'];
  end
  [save_dir] =fileparts(save_path);
  if ~exist(save_dir,'dir')
    
    mkdir(save_dir);
  end
  saveas(figure(5),save_path,'jpg')
  
  %% Plot Total Variable Attn
  figure(6)
  mapshow(rgb2gray(A),CMAP/1e3);
  xlabel('X (km)');
  ylabel('Y (km)');
  if  strcmp(settings.location,'Peterman')
    xlim([-350 -50]);
    ylim([-1250 -900]);
  else
  xlim([-250 -50]);
    ylim([-2400 -2160]);
  end
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
  
  if strcmp(settings.location,'Peterman')
    save_path=['/cresis/snfs1/scratch/manjish/new_peterman/var_att_median_att3'];
  else
    save_path=['/cresis/snfs1/scratch/manjish/new_jacobshavn/var_att_median_att3'];
  end
  [save_dir] =fileparts(save_path);
  if ~exist(save_dir,'dir')
    
    mkdir(save_dir);
  end
  saveas(figure(6),save_path,'jpg')
  
  %% Plot Value of DN
  %   figure(7)
  %   mapshow(rgb2gray(A),CMAP/1e3);
  %   xlabel('X (km)');
  %   ylabel('Y (km)');
  %   %xlim([350 650]);
  %   %ylim([-1000 -600]);
  %
  %   hold on
  %   clear gps.x gps.y
  %   [gps.x,gps.y] = projfwd(proj,lat_G_r_corrected,lon_G_r_corrected);
  %   gps.x = gps.x / 1000;
  %   gps.y = gps.y / 1000;
  %   hold on;
  %
  %   scatter(gps.x,gps.y,20,estimated_DN,'fill')
  %   %caxis([-15 15])
  %   colorbar;
  %   title('Value of DN ')
end