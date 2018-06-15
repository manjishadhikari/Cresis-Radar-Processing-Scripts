%Script: estimation_of_relative_reflectivity_values
%  Fit Na and DN for every 1 km
%{
from englacial_att_method1

median_power =
    8.8886
mean_power =
  -15.3333
max_power =
   66.3854
avg_depth =
   1505
%}

%% setup
clear
close all
clc
dbstop error

plots =1;
coh_int=0;
ice_bed_power_G_r_corrected = [];
lat_G_r_corrected = [];
lon_G_r_corrected = [];
depth_G_r_corrected = [];
constant_attenuation = [];
estimated_Na = [];
estimated_DN=[];
variable_attenuation=[];
all_GR=[];

%% loading the data


%%
disp('Englacial Attn Method 2')
for M =1:35
  
  clearvars -except all_GR M range_power coh_int plots ice_bed_power_G_r_corrected lat_G_r_corrected lon_G_r_corrected depth_G_r_corrected cross_lines constant_attenuation estimated_Na estimated_DN variable_attenuation
  clc
  
  param.radar.fs = 195000000;
  if M<21
    cross_lines = 1;
    load(['/cresis/snfs1/scratch/manjish/peterman/radar_w_index/crossline',num2str(M)]);
  else
    M1=M-20;
    cross_lines = 0;
    load(['/cresis/snfs1/scratch/manjish/peterman/radar_w_index/verticalline',num2str(M1)]);
  end
  physical_constants
  %%
  
  Greenland.depth = (Greenland.ice_bed_time - Greenland.surface_time)*c/2/sqrt(er_ice);
  Greenland.surface_height = (Greenland.surface_time)*c/2;
  
  
  geometric_loss = (2*(Greenland.surface_height+Greenland.depth/sqrt(er_ice))).^2;
  geometric_loss_surface = (2*(Greenland.surface_height)).^2;
  
  if plots
    figure(1);plot(Greenland.depth, lp(Greenland.ice_bed_power));
    grid on; title('Depth vs Power')
    figure(2); plot( lp(Greenland.ice_bed_power));
    grid on; title('Along Track vs Power')
    
    if 0 
      geotiff_fn = '/cresis/snfs1/dataproducts/GIS_data/greenland/Landsat-7/Greenland_natural.tif';
      proj = geotiffinfo(geotiff_fn);
      %proj = geotiffinfo('X:\GIS_data\antarctica\Landsat-7\Antarctica_LIMA_480m.tif');
      
      [A CMAP R]= geotiffread(geotiff_fn);
      
      figure(100)
      mapshow(rgb2gray(A),CMAP/1e3);
      xlabel('X (km)');
      ylabel('Y (km)');
      %xlim([350 650]);
      %ylim([-1000 -600]);
      
      hold on
      clear gps.x gps.y
      [gps.x,gps.y] = projfwd(proj,Greenland.Latitude,Greenland.Longitude);
      
      gps.x = gps.x / 1000;
      gps.y = gps.y / 1000;
      hold on;
      
      scatter(gps.x,gps.y,20,lp(Greenland.ice_bed_power),'fill')
      scatter(gps.x(1),gps.y(1),50,'X')
      %caxis([-15 15])
      colorbar;
      title('Radar line')
    close
   end
   
    
  end
  
  
  %Coherent Integration to increase SNR 
  
   %% COHERENT INTEGRATIONS
  if coh_int~=0 
  Greenland=coh_integration(Greenland,coh_int);
  end
  
  Greenland.ice_bed_power=abs(Greenland.ice_bed_power).^2;
  Greenland.ice_surface_power=abs(Greenland.ice_surface_power).^2;


  Greenland.ice_bed_power_cgl =lp(Greenland.ice_bed_power)+lp(geometric_loss);
  if plots
      figure;plot(Greenland.depth,lp(Greenland.ice_bed_power));
  
    figure(4);subplot(2,1,1);plot(Greenland.ice_bed_power_cgl);title('Bed power after Geom corr')
    subplot(2,1,2);plot(lp(geometric_loss)); title('Geom correction')
  end
  
  id = ~(isfinite( Greenland.ice_bed_power_cgl ));
  Greenland.ice_bed_power_cgl(id) = [];
  Greenland.Latitude(id) = [];
  Greenland.Longitude(id) = [];
  Greenland.depth(id) = [];
  
  
  %% attenuation_fitting
  %reference_power = 25 ;
  reference_power = median((Greenland.ice_bed_power_cgl));
  %reference_power=-15;
  relative_ice_bed_power_G_r_corrected = (Greenland.ice_bed_power_cgl)-reference_power;
  
  %relative_ice_bed_power_G_r_corrected= relative_ice_bed_power_G_r_corrected-mean(relative_ice_bed_power_G_r_corrected);
   Greenland.depth = Greenland.depth/1000;
  %relative_depth = mean( Greenland.depth_avg);
%   relative_depth = nanmean(isfinite(Greenland.depth_avg));
  relative_depth =1.505;
  % relative_depth=1.6033;
  if plots
    figure;plot(relative_ice_bed_power_G_r_corrected);
    title('Relative Power')
  end

  if range(relative_ice_bed_power_G_r_corrected)>102
    keyboard
  end
  
  range_power(M)= range(relative_ice_bed_power_G_r_corrected);
  all_GR=cat(2,all_GR,relative_ice_bed_power_G_r_corrected);
%   window1=200;
%   window2=10;
%   power_filtered_long=sgolayfilt(relative_ice_bed_power_G_r_corrected,2,window1+1,gausswin(window1+1));
%   power_filtered_short=sgolayfilt(relative_ice_bed_power_G_r_corrected,2,window2+1,gausswin(window2+1));
%   
%   if plots
%     figure;plot(power_filtered_long); hold on; plot(relative_ice_bed_power_G_r_corrected)
%     legend('filtered','original'); title('Long filter')
%     figure;plot(power_filtered_short); hold on; plot(relative_ice_bed_power_G_r_corrected)
%     legend('filtered','original'); title('short filter')
%   end
%   
%  
%   
%   constant_attenuation =  cat(2,constant_attenuation, const_attenuation);
%   estimated_Na = cat(2,estimated_Na, estimated_na);
%   estimated_DN=cat(2,estimated_DN,estimated_dn);
%   variable_attenuation=cat(2,variable_attenuation,var_attenuation);
%   
%   close all
%   %%
  %ref=Greenland.ice_bed_power_cgl-mean(Greenland.ice_bed_power_cgl)+const_attenuation;
  %ref2=Greenland.ice_bed_power_cgl-mean(Greenland.ice_bed_power_cgl)+var_attenuation;
end


%% Reflectivity using constant na

% figure(1)
% mapshow(rgb2gray(A),CMAP/1e3);
% xlabel('X (km)');
% ylabel('Y (km)');
% %xlim([350 650]);
% %ylim([-1000 -600]);
% 
% hold on
% clear gps.x gps.y
% [gps.x,gps.y] = projfwd(proj,lat_G_r_corrected,lon_G_r_corrected);
% 
% gps.x = gps.x / 1000;
% gps.y = gps.y / 1000;
% hold on;
% 
% scatter(gps.x,gps.y,20,relative_reflectivty_using_constant_attn,'fill')
% caxis([-15 15])
% colorbar;
% title('Reflectivity using constant na ')
% 
% %Histogram
% figure(2), hist(relative_reflectivty_using_constant_attn,30);
% title('Reflectivity using constant na ')
% 
