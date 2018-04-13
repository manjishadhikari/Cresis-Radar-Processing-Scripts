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
dbstop if error

plots =1;
numofCohInt=0;
ice_bed_power_G_r_corrected = [];
lat_G_r_corrected = [];
lon_G_r_corrected = [];
depth_G_r_corrected = [];
constant_attenuation = [];
estimated_Na = [];
estimated_DN=[];
variable_attenuation=[];
c_ref=[];
v_ref=[];
%% loading the data


disp('Englacial Attn Method 2')
for iter=1
  for M =1:20
    
    % clearvars -except M coh_int plots ice_bed_power_G_r_corrected lat_G_r_corrected lon_G_r_corrected depth_G_r_corrected cross_lines constant_attenuation estimated_Na estimated_DN variable_attenuation
    clc
    param.radar.fc = 195000000;  %Center Frequency
    if M<21
      cross_lines = 1;
      M1=0;
      load(['/cresis/snfs1/scratch/manjish/new_peterman/radar_w_idx_new/crossline',num2str(M)]);
    else
      M1=M-20;
      cross_lines = 0;
      load(['/cresis/snfs1/scratch/manjish/new_peterman/radar_w_idx_new/verticalline',num2str(M1)]);
    end
    physical_constants
    %%
    
    Greenland.depth = (Greenland.ice_bed_time - Greenland.surface_time)*c/2/sqrt(er_ice);
    Greenland.surface_height = (Greenland.surface_time)*c/2;
    Greenland.geometric_loss = (2*(Greenland.surface_height+Greenland.depth/sqrt(er_ice))).^2;
    Greenland.geometric_loss_surface = (2*(Greenland.surface_height)).^2;
    
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
        xlim([-350 -50]);
        ylim([-1250 -900]);
        %  xlim([-350 -50]);
        %ylim([-2450 -2150]);
        
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
        keyboard
        close
      end
      
      
    end
    
    
    %Coherent Integration to increase SNR
    
    %% COHERENT INTEGRATIONS
    if numofCohInt~=0
      [Greenland]=coh_integration(Greenland,numofCohInt);
    end
    
    
    Greenland.ice_bed_power=abs(Greenland.ice_bed_power).^2;   %Ice bed power
    Greenland.ice_surface_power=abs(Greenland.ice_surface_power).^2;  %Ice surface power
    Greenland.roll=Greenland.Roll*180/pi;   %Roll in degrees
    
    if ~isempty(find(isnan(Greenland.ice_bed_power),1))
      disp(sprintf('%d Nan values found for bed power \n',length(find(isnan(Greenland.ice_bed_power),1))))
    end
    
    %% compensating reflected bed power for surface roughness
    settings.num_int=1000;
    settings.repeat_after=10;
    settings.type='surface';
    settings.M=M;
    settings.M1=M1;
    % [Greenland,sf_rms]=surf_roughness(Greenland,num_int,repeat_after);
    [Greenland,sf_rms,sf_corr_power,orig_avg_power]=surf_roughness(Greenland,settings);
    
    
    if plots
      figure(3);subplot(3,1,1); plot(lp(orig_avg_power));
      
      hold on; plot(lp( Greenland.ice_bed_power_avg));
      grid on
      title('Ice Bed Power surface roughness corrected')
      subplot(3,1,2); plot(lp(sf_corr_power));title('Sf corrected power')
    end
    
    
    
    
    %% compensating for bed roughness
   if 0
    settings.type='bed';
    settings.iter=iter;
    if settings.iter==2
      figure(201); subplot(2,1,1); plot(lp(Greenland.ice_bed_power_avg));
      Greenland.ice_bed_power_avg=Greenland.ice_bed_power_avg.*(10.^(Attenuation.const_attenuation/10));
      figure(201); subplot(2,1,1);hold on; plot(lp(Greenland.ice_bed_power_avg));
      legend('Before const attn','After const attn')
      figure(201); hold on;subplot(2,1,2); hold on; plot(Attenuation.const_attenuation);
    end
    [Greenland,bed_rms,bed_corr_power]=bed_roughness(Greenland,settings);
    
    if plots
      figure(3);subplot(3,1,1)
      hold on;
      plot(lp( Greenland.ice_bed_power_avg));
      grid on
      legend('Original','Sf corrected ',' Bed roughness corrected')
      figure;plot(lp(orig_avg_power)-lp(Greenland.ice_bed_power_avg));
      title('Power Difference after roughness correction')
      if exist('bed_corr_power','var')
        figure(3); subplot(3,1,3); plot(lp(bed_corr_power));title('Bed Corrected power')
      end
    end
   end
%     if iter==2
%       continue;
%     end
%     
    %%  relative geometrically corrected bed - echo power
    %     for i = 1:length(Greenland.ice_bed_power)
    %         if i < 3
    %             Greenland.ice_bed_power_rgc(i) = Greenland.ice_bed_power(i)-nanmean(Greenland.ice_bed_power(1:6-i));
    %         elseif i+2 > length(Greenland.ice_bed_power)
    %             Greenland.ice_bed_power_rgc(i) = Greenland.ice_bed_power(i)- nanmean(Greenland.ice_bed_power(i-5:end));
    %         else
    %             Greenland.ice_bed_power_rgc(i) = Greenland.ice_bed_power(i)- nanmean(Greenland.ice_bed_power(i-2:i+2));
    %         end
    %     end
    
    if length(Greenland.ice_bed_power_avg) ~= length(Greenland.geometric_loss_avg)
      keyboard
    end
    
    Greenland.ice_bed_power_cgl =lp(Greenland.ice_bed_power_avg)+lp(Greenland.geometric_loss_avg);
    if plots
      figure(4);subplot(3,1,1);plot(Greenland.ice_bed_power_cgl);title('Bed power after Geom corr')
      subplot(3,1,2);plot(lp(Greenland.geometric_loss_avg)); title('Geom correction')
      subplot(3,1,3); plot(Greenland.depth_avg); title('Depth')
    end
    
    %     id = ~(isfinite( Greenland.ice_bed_power_cgl ));
    %     Greenland.ice_bed_power_cgl(id) = [];
    %     Greenland.Latitude_avg(id) = [];
    %     Greenland.Longitude_avg(id) = [];
    %     Greenland.depth_avg(id) = [];
    
    
    ice_bed_power_G_r_corrected = cat(2,ice_bed_power_G_r_corrected,Greenland.ice_bed_power_cgl);
    lat_G_r_corrected =  cat(2,lat_G_r_corrected,Greenland.Latitude_avg);
    lon_G_r_corrected =  cat(2,lon_G_r_corrected,Greenland.Longitude_avg);
    depth_G_r_corrected =  cat(2,depth_G_r_corrected,Greenland.depth_avg);
    
    
    %% attenuation_fitting
    %reference_power = 25 ;
    Greenland.reference_power = nanmean((Greenland.ice_bed_power_cgl));
    %reference_power=-15;
    Greenland.relative_ice_bed_power_G_r_corrected = (Greenland.ice_bed_power_cgl)-Greenland.reference_power;
    
    %relative_ice_bed_power_G_r_corrected= relative_ice_bed_power_G_r_corrected-mean(relative_ice_bed_power_G_r_corrected);
    Greenland.depth_avg = Greenland.depth_avg/1000;
    Greenland.relative_depth = nanmean( Greenland.depth_avg);
    %   relative_depth = nanmean(isfinite(Greenland.depth_avg));
    % relative_depth =1.505;
    % relative_depth=1.6033;
    if plots
      figure;plot(Greenland.relative_ice_bed_power_G_r_corrected);
      title('Relative Power')
    end
    dist=geodetic_to_along_track(lat_G_r_corrected,lon_G_r_corrected);
    
    window1=floor(30000/(dist(100)-dist(99)));  %Long filter of 30 km
    % window1=500;
    window2=floor(1500/(dist(100)-dist(99)));       %Short filter of 1.5 km
    if window1> size(Greenland.relative_ice_bed_power_G_r_corrected)
      window1=floor(size(Greenland.relative_ice_bed_power_G_r_corrected,2)/10);
      disp('Window is small')
    end
    
    if mod(window1,2)==0
      window1=window1-1;  %making odd for sgolayfilt
    end
    if mod(window2,2)==0
      window2=window2-1;  %making odd for sgolayfilt
    end
    %     window1=240;
    %     window2=40;
    
    nanidx=(find(isnan(Greenland.relative_ice_bed_power_G_r_corrected)));
    %  nanidx2=find(isnan(power_filtered_short));
    % nanidx=unique([nanidx1 nanidx2]);
    test=zeros(1,length(Greenland.relative_ice_bed_power_G_r_corrected));
    test(nanidx)=nan;
    notnanidx=find(~isnan(test));
    
    %notnanidx=find(power_filtered_long(~nanidx));
    Greenland.relative_ice_bed_power_G_r_corrected(nanidx)=[];
    % power_filtered_short(nanidx)=[];
    Greenland.depth_avg(nanidx)=[];
    Greenland.Latitude_avg(nanidx)=[];
    Greenland.Longitude_avg(nanidx)=[];
    Greenland.along_track_avg(nanidx)=[];
    Greenland.geometric_loss_avg(nanidx)=[];
    
    
    
    if iter<5
      power_filtered_long=sgolayfilt(Greenland.relative_ice_bed_power_G_r_corrected,2,window1,gausswin(window1));
      power_filtered_short=sgolayfilt(Greenland.relative_ice_bed_power_G_r_corrected,2,window2,gausswin(window2));
    else
      power_filtered_long=Greenland.relative_ice_bed_power_G_r_corrected;
      power_filtered_short=Greenland.relative_ice_bed_power_G_r_corrected;
      
    end
    if plots
      figure;plot(power_filtered_long); hold on; plot(Greenland.relative_ice_bed_power_G_r_corrected)
      legend('filtered','original'); title('Long filter')
      figure;plot(power_filtered_short); hold on; plot(Greenland.relative_ice_bed_power_G_r_corrected)
      legend('filtered','original'); title('short filter')
    end
    
    
    
    
    
    %    notnanidx=find(power_filtered_long(~nanidx));
    %      nanidx=find(isnan(power_filtered_short));
    %     power_filtered_long(nanidx)=[];
    %     power_filtered_short(nanidx)=[];
    %       Greenland.depth_avg(nanidx)=[];
    %     Greenland.Latitude_avg(nanidx)=[];
    %     Greenland.Longitude_avg(nanidx)=[];
    %     Greenland.along_track_avg(nanidx)=[];
    %     Greenland.geometric_loss_avg(nanidx)=[];
    %      Greenland.ice_bed_power_cgl(nanidx)=[];
    
    if plots
      figure;plot(power_filtered_long);title('Power Filtered Long')
      figure;plot(power_filtered_short);title('Power Filtered Short')
    end
    
       %% Attenuation calculation
     %Method 1 fit Na and DN for evry 1 km
   [Attenuation]=attenuation_calculation_method1(Greenland,power_filtered_long,power_filtered_short);

    %%Method 2 use Na from 1 and fit DN
    %  [Attenuation]=attenuation_calculation_method1(Greenland,power_filtered_short);
    
 
    if plots
        figure;plot(Attenuation.const_attenuation);title('Const Attenuation')
        figure;plot(Attenuation.var_attenuation); title('Var attenuation')
        figure;plot(Greenland.depth); title('Depth')
     end 
    
    if ~isempty(nanidx)
      constant_att=zeros(1,size(Greenland.ice_bed_power_cgl,2));
      var_att=zeros(1,size(Greenland.ice_bed_power_cgl,2));
      const_att(nanidx)=nan;
      var_att(nanidx)=nan;
      const_att(notnanidx)=const_attenuation;
      var_att(notnanidx)=var_attenuation;
      const_attenuation=const_att;
      var_attenuation=var_att;
    end
    
    constant_attenuation =  cat(2,constant_attenuation, Attenuation.const_attenuation);
    estimated_Na = cat(2,estimated_Na, Attenuation.estimated_na);
    estimated_DN=cat(2,estimated_DN,Attenuation.estimated_dn);
    variable_attenuation=cat(2,variable_attenuation,Attenuation.var_attenuation);
    
      ref=Greenland.ice_bed_power_cgl-nanmean(Greenland.ice_bed_power_cgl)+Attenuation.const_attenuation;
      ref2=Greenland.ice_bed_power_cgl-nanmean(Greenland.ice_bed_power_cgl)+Attenuation.var_attenuation;
    c_ref=cat(2,c_ref,ref);
    v_ref=cat(2,v_ref,ref2);
    
    if plots
%       ref=Greenland.ice_bed_power_cgl-nanmean(Greenland.ice_bed_power_cgl)+Attenuation.const_attenuation;
%       ref2=Greenland.ice_bed_power_cgl-nanmean(Greenland.ice_bed_power_cgl)+Attenuation.var_attenuation;
      figure;plot(ref);title('Reflectivity using const Na')
      figure;histogram(ref);title('Reflectivity using const Na')
      figure;plot(ref2);title('Reflectivity using var Na')
      figure;histogram(ref2);title('Reflectivity using var Na')
    end
    
    
%     const_att_iter{iter}=Attenuation.const_attenuation;
%     var_att_iter{iter}=Attenuation.var_attenuation;
%     c_ref{iter}=ref;
%     v_ref{iter}=ref2;
    close all
    
  end
end


% for i=1:iter
%   
%   figure(10000);hold on; plot(const_att_iter{i});title('Const Attn')
%   figure(10001);hold on; plot(var_att_iter{i});title('Var Attn')
%   figure(10002);hold on; plot(c_ref{i});title('Reflectivity using const Na')
%   figure(10003);hold on; plot(v_ref{i});title('Reflectivity using var Na')
%   
% end
% figure(10004); histogram(c_ref{i});title('Iter 3 consta att')
% figure(10005); histogram(v_ref{i});title('Iter 3 const att')
% figure;(10006); histogram(c_ref{1}); title('Iter 1 const att')
keyboard
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
  %xlim([350 650]);
  %ylim([-1000 -600]);
  
  hold on
  clear gps.x gps.y
  [gps.x,gps.y] = projfwd(proj,lat_G_r_corrected,lon_G_r_corrected);
  
  gps.x = gps.x / 1000;
  gps.y = gps.y / 1000;
  hold on;
  
  scatter(gps.x,gps.y,20,c_ref,'fill')
  caxis([-15 15])
  colorbar;
  title('Reflectivity using constant na ')
  
  %Histogram
  figure(2), hist(c_ref,30);
  title('Reflectivity using constant na ')
  
  %% Reflectivity using variable attenuation
  figure(3)
  mapshow(rgb2gray(A),CMAP/1e3);
  xlabel('X (km)');
  ylabel('Y (km)');
  %xlim([350 650]);
  %ylim([-1000 -600]);
  
  hold on
  clear gps.x gps.y
  [gps.x,gps.y] = projfwd(proj,lat_G_r_corrected,lon_G_r_corrected);
  
  gps.x = gps.x / 1000;
  gps.y = gps.y / 1000;
  hold on;
  
  scatter(gps.x,gps.y,20,v_ref,'fill')
  caxis([-15 15])
  colorbar;
  title('Reflectivity using variable na ')
  
  %Histogram
  figure(4), hist(v_ref,30);
  title('Reflectivity using variable na ')
  
  keyboard;
  %% Plot Total Constant Attn
  figure(5)
  mapshow(rgb2gray(A),CMAP/1e3);
  xlabel('X (km)');
  ylabel('Y (km)');
  %xlim([350 650]);
  %ylim([-1000 -600]);
  
  hold on
  clear gps.x gps.y
  [gps.x,gps.y] = projfwd(proj,lat_G_r_corrected,lon_G_r_corrected);
  gps.x = gps.x / 1000;
  gps.y = gps.y / 1000;
  hold on;
  
  scatter(gps.x,gps.y,20,constant_attenuation,'fill')
  %caxis([-15 15])
  colorbar;
  title('Total Constant Attenuation')
  
  %% Plot Total Variable Attn
  figure(6)
  mapshow(rgb2gray(A),CMAP/1e3);
  xlabel('X (km)');
  ylabel('Y (km)');
  %xlim([350 650]);
  %ylim([-1000 -600]);
  
  hold on
  clear gps.x gps.y
  [gps.x,gps.y] = projfwd(proj,lat_G_r_corrected,lon_G_r_corrected);
  gps.x = gps.x / 1000;
  gps.y = gps.y / 1000;
  hold on;
  
  scatter(gps.x,gps.y,20,variable_attenuation,'fill')
  %caxis([-15 15])
  colorbar;
  title('Total Variable Attenuation')
  
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