%Script: estimation_of_relative_reflectivity_values
%  Fit Na and DN for every 1 km
%{
from englacial_att_method1
peterman
median_power =
  -59.6864
mean_power =
  -27.6142
max_power =
    2.7618
avg_depth =
   1.5050e+03
Na_reg =
    9.7746
Na_bar =
    9.7500

%%Coherent int 600 peterman
median_power =
  -76.2135
mean_power =
  -41.1230
max_power =
   -7.1976
avg_depth =
   1.5000e+03
Na_reg =
    9.7221
Na_bar =
    9.7000


%%Jacobshavn
median_power =
  -61.5947
mean_power =
  -43.2931
max_power =
   -5.4816
avg_depth =
   1.1496e+03
Na_reg =
   13.1028
Na_bar =
   13.1000
%}

%% setup
clear
close all
clc
dbstop if error

settings.location={'Jacobshavn'};


plots =1;
numofCohInt=0;
ice_bed_power_G_r_corrected = [];
lat_G_r_corrected = [];
lon_G_r_corrected = [];
depth_G_r_corrected = [];
constant_attenuation = [];
estimated_Na = [];
estimated_DN=[];
modified_Na=[];
variable_attenuation=[];
c_ref=[];
v_ref=[];
line_no=[];
  Ice_bed_elevation=[];
%% loading the data

if strcmp(settings.location,'Jacobshavn')
  
  median_power= -61.3157;
  mean_power =  -46.2426;
  max_power =-13.6709;
  avg_depth =1.1874e+03;
  Na_reg =12.9632;
  Na_bar =12.95;
  num_of_lines=73;
  cross_line_no=33;
  
  %Crossline
  median_power=-61.84;
  Na_bar=13.4;
    avg_depth= 1.2426e+03;
  %  Na_bar=13;
   % median_power=-62.99;
  
elseif strcmp(settings.location,'Peterman')
%   median_power =-59.6864;
%   mean_power =-27.6142;
%   max_power =2.7618;
%   avg_depth =1.5050e+03;
%   Na_reg = 9.7746;
%   Na_bar =9.7500;
 mean_power =-46.2754;
   median_power =-60.1428;
   max_power =   -8.37;
   avg_depth =1.5456e+03;
   Na_reg = 8.129;
   Na_bar =8.1;
  num_of_lines=35;
  cross_line_no=20;
end

disp('Englacial Attn Method 2')
for iter=1
  for M =30
    
    % clearvars -except M   Ice_bed_elevation coh_int plots ice_bed_power_G_r_corrected lat_G_r_corrected lon_G_r_corrected depth_G_r_corrected cross_lines constant_attenuation estimated_Na estimated_DN variable_attenuation
    clc
    param.radar.fc = 195000000;  %Center Frequency
    if strcmp(settings.location,'Jacobshavn')
      if M<=cross_line_no
        cross_lines = 1;
        M1=0;
        load(['/cresis/snfs1/scratch/manjish/new_jacobshavn/radar_w_idx_new/crossline',num2str(M)]);
      else
        M1=M-cross_line_no;
        cross_lines = 0;
        load(['/cresis/snfs1/scratch/manjish/new_jacobshavn/radar_w_idx_new/verticalline',num2str(M1)]);
      end
      
    elseif strcmp(settings.location,'Peterman')
      if M<21
        cross_lines = 1;
        M1=0;
        load(['/cresis/snfs1/scratch/manjish/new_peterman/radar_w_idx_new/crossline',num2str(M)]);
      else
        M1=M-20;
        cross_lines = 0;
        load(['/cresis/snfs1/scratch/manjish/new_peterman/radar_w_idx_new/verticalline',num2str(M1)]);
      end
      
      
    end
    
      figure; plot( lp(Greenland.ice_bed_power));
      grid on; title('Along Track vs Power before roll correction')
    
    Greenland.roll=Greenland.Roll*180/pi;
    idx=find(abs(Greenland.roll)>5);
    Greenland.ice_bed_power(idx)=nan;
    
    
    
    physical_constants
    figure; plot( lp(Greenland.ice_bed_power));
      grid on; title('Along Track vs Power before coh integration and roll correction')
      
      %% COHERENT INTEGRATIONS
    if numofCohInt~=0
      [Greenland]=coh_integration(Greenland,numofCohInt);    
    end
    
    %%
    Greenland.depth = (Greenland.ice_bed_time - Greenland.surface_time)*c/2/sqrt(er_ice);
    Greenland.surface_height = (Greenland.surface_time)*c/2;
    Greenland.geometric_loss = (2*(Greenland.surface_height+Greenland.depth/sqrt(er_ice))).^2;
    Greenland.geometric_loss_surface = (2*(Greenland.surface_height)).^2;
    Greenland.ice_bed_elevation=Greenland.Elevation-Greenland.depth-Greenland.depth;
    if plots
      figure(1);plot(Greenland.depth, lp(Greenland.ice_bed_power));
      grid on; title('Depth vs Power')
      figure(2); plot( lp(Greenland.ice_bed_power));
      grid on; title('Along Track vs Power')
      
      %Location
      if 0
        geotiff_fn = '/cresis/snfs1/dataproducts/GIS_data/greenland/Landsat-7/Greenland_natural.tif';
        proj = geotiffinfo(geotiff_fn);
        %proj = geotiffinfo('X:\GIS_data\antarctica\Landsat-7\Antarctica_LIMA_480m.tif');
        
        [A CMAP R]= geotiffread(geotiff_fn);
        
        figure(100)
        mapshow(rgb2gray(A),CMAP/1e3);
        xlabel('X (km)');
        ylabel('Y (km)');
        
        if  strcmp(settings.location,'Peterman')
          xlim([-350 -50]);
          ylim([-1250 -900]);
        else
          xlim([-350 -50]);
          ylim([-2450 -2150]);
        end
        
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
    
    Greenland.ice_bed_power=abs(Greenland.ice_bed_power).^2;   %Ice bed power
    Greenland.ice_surface_power=abs(Greenland.ice_surface_power).^2;  %Ice surface power
    Greenland.roll=Greenland.Roll*180/pi;   %Roll in degrees
    
    if ~isempty(find(isnan(Greenland.ice_bed_power),1))
      disp(sprintf('%d Nan values found for bed power \n',length(find(isnan(Greenland.ice_bed_power)))))
    end
    
    %% compensating reflected bed power for surface roughness
    settings.num_int=1000;
    settings.repeat_after=10;
    settings.type='surface';
    settings.cross_lines=cross_lines;
    settings.M=M;
    settings.M1=M1;
    settings.rerun=1;
    settings.save_en=0;
    % [Greenland,sf_rms]=surf_roughness(Greenland,num_int,repeat_after);
    [Greenland,sf_rms,sf_corr_power,orig_avg_power]=surf_roughness(Greenland,settings);
    
    
    if plots
      figure(3);subplot(3,1,1); plot(lp(orig_avg_power));
      hold on; plot(lp( Greenland.ice_bed_power_avg));
      grid on
      title('Ice Bed Power surface roughness corrected')
      subplot(3,1,2); plot(lp(sf_corr_power));title('Sf corrected power')
      figure(500);plot(sf_rms.rms_height*100); title('Sf Rms height')
    end
    
    
    
    
    %% compensating for bed roughness
    if 0
      settings.type='bed';
      settings.iter=iter;
      if settings.iter==1
        figure(201); subplot(2,1,1); plot(lp(Greenland.ice_bed_power_avg));
      %  Greenland.ice_bed_power_avg=Greenland.ice_bed_power_avg.*(10.^(Attenuation.const_attenuation/10));
       % figure(201); subplot(2,1,1);hold on; plot(lp(Greenland.ice_bed_power_avg));
       % legend('Before const attn','After const attn')
        %figure(201); hold on;subplot(2,1,2); hold on; plot(Attenuation.const_attenuation);
        [Greenland,bed_rms,bed_corr_power]=bed_roughness(Greenland,settings);
      end
   
      
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
  %  Greenland.reference_power = nanmedian((Greenland.ice_bed_power_cgl));
    Greenland.reference_power =median_power;
    %reference_power=-15;
    Greenland.relative_ice_bed_power_G_r_corrected = (Greenland.ice_bed_power_cgl)-Greenland.reference_power;
    
    %relative_ice_bed_power_G_r_corrected= relative_ice_bed_power_G_r_corrected-mean(relative_ice_bed_power_G_r_corrected);
    Greenland.depth_avg = Greenland.depth_avg/1000;
 %  Greenland.relative_depth = nanmean( Greenland.depth_avg);
    Greenland.relative_depth =avg_depth/1000;
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
    %  nanidx2=find(isnan(power_filtered_short));Attenuation.mod_na
    % nanidx=unique([nanidx1 nanidx2]);
    test=zeros(1,length(Greenland.relative_ice_bed_power_G_r_corrected));
    test(nanidx)=nan;
    notnanidx=find(~isnan(test));
    
    %notnanidx=find(power_filtered_long(~nanidx));
    Greenland.relative_ice_bed_power_G_r_corrected_filt= Greenland.relative_ice_bed_power_G_r_corrected(notnanidx);
    % power_filtered_short(nanidx)=[];
    Greenland.depth_avg_filt=Greenland.depth_avg(notnanidx);   %Filtered
    Greenland.Latitude_avg_filt= Greenland.Latitude_avg(notnanidx);
    Greenland.Longitude_avg_filt= Greenland.Longitude_avg(notnanidx);
    Greenland.along_track_avg_filt=Greenland.along_track_avg(notnanidx);
    Greenland.geometric_loss_avg_filt=Greenland.geometric_loss_avg(notnanidx);
    
    
    
    if iter<10
      power_filtered_long=sgolayfilt(Greenland.relative_ice_bed_power_G_r_corrected_filt,2,window1,gausswin(window1));
      power_filtered_short=sgolayfilt(Greenland.relative_ice_bed_power_G_r_corrected_filt,2,window2,gausswin(window2));
    else
      power_filtered_long=Greenland.relative_ice_bed_power_G_r_corrected_filt;
      power_filtered_short=Greenland.relative_ice_bed_power_G_r_corrected_filt;
      
    end
    if plots
      figure;plot(power_filtered_long); hold on; plot(Greenland.relative_ice_bed_power_G_r_corrected_filt)
      legend('filtered','original'); title('Long filter')
      figure;plot(power_filtered_short); hold on; plot(Greenland.relative_ice_bed_power_G_r_corrected_filt)
      legend('filtered','original'); title('short filter')
    end
    
    
    
    
    %    notnanidx=find(power_filtered_long(~nanidx));
    %      nanidx=find(isnan(power_filtered_short));
    %     power_filtered_long(nanidx)=[];
    %     power_filtered_short(nanidx)=[];
    %       Greenland.depth_avg(nanidx)=[];Attenuation.mod_na
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
  %   out.att_method=1;
   %  [Attenuation]=attenuation_calculation_method1(Greenland,power_filtered_long,power_filtered_short);
%     
    %%Method 2 use Na from overall and fit DN for every line
    out.att_method=2;
   [Attenuation]=attenuation_calculation_method2(Greenland,power_filtered_short,Na_bar);
    
     %Method 3 use Na from overall and fit DN for evry 1 km 
 %    out.att_method=3;
  %  [Attenuation]=attenuation_calculation_method3(Greenland,power_filtered_short,Na_bar);
    
    %Method 
    
    if plots
      figure;plot(Attenuation.const_attenuation);title('Const Attenuation')
     hold on;plot(Attenuation.var_attenuation); title('Const and Var attenuation')
      figure;plot(Greenland.depth); title('Depth')
      figure;plot(Attenuation.estimated_dn); title('Estimated DN')
      figure;plot(Attenuation.estimated_na); title('Estimated NA')
      figure;plot(Attenuation.mod_na); title('Modified NA')
      
      if 0
        figure;plot(Attenuation2.const_attenuation);title('Const Attenuation')
        figure;plot(Attenuation2.var_attenuation); title('Var attenuation')
        figure;plot(Attenuation2.mod_na); title('Modified NA')
      end
      
    end
    
    if ~isempty(nanidx)
      const_att=zeros(1,size(Greenland.ice_bed_power_cgl,2));
      power_g_r_corr=zeros(1,size(Greenland.ice_bed_power_cgl,2));
      var_att=zeros(1,size(Greenland.ice_bed_power_cgl,2));
       modi_na=zeros(1,size(Greenland.ice_bed_power_cgl,2));
      const_att(nanidx)=nan;
      var_att(nanidx)=nan;
      power_g_r_corr(nanidx)=nan;
      modi_na(nanidx)=nan;
      const_att(notnanidx)=Attenuation.const_attenuation;
      var_att(notnanidx)=Attenuation.var_attenuation;
      modi_na(notnanidx)=Attenuation.mod_na;
      power_g_r_corr(notnanidx)=Greenland.relative_ice_bed_power_G_r_corrected_filt;
      Attenuation.const_attenuation=const_att;
      Attenuation.var_attenuation=var_att;
      Attenuation.mod_na=modi_na;
      Greenland.relative_ice_bed_power_G_r_corrected=power_g_r_corr;
      
    end
    if iter==1
    constant_attenuation =  cat(2,constant_attenuation, Attenuation.const_attenuation);
    estimated_Na = cat(2,estimated_Na, Attenuation.estimated_na);
    estimated_DN=cat(2,estimated_DN,Attenuation.estimated_dn);
    modified_Na=cat(2,modified_Na,Attenuation.mod_na);
    variable_attenuation=cat(2,variable_attenuation,Attenuation.var_attenuation);
    
    ref=Greenland.relative_ice_bed_power_G_r_corrected+Attenuation.const_attenuation;
    ref2=Greenland.relative_ice_bed_power_G_r_corrected+Attenuation.var_attenuation;
    c_ref=cat(2,c_ref,ref);
    v_ref=cat(2,v_ref,ref2);
    line_num=M*ones(1,length(v_ref));
    line_no=cat(2,line_no,line_num);
     Ice_bed_elevation=cat(2,Ice_bed_elevation, Greenland.ice_bed_elevation);
    end
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
%keyboard
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
  if  strcmp(settings.location,'Peterman')
    xlim([-350 -50]);
    ylim([-1250 -900]);
  else
   xlim([-250 -50]);
    ylim([-2400 -2160]);
  end
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
  [gps.x,gps.y] = projfwd(proj,lat_G_r_corrected,lon_G_r_corrected);
  gps.x = gps.x / 1000;
  gps.y = gps.y / 1000;
  hold on;
  
  scatter(gps.x,gps.y,20,modified_Na,'fill')
  %caxis([-15 15])
  colorbar;
  title('Total Constant Attenuation')
  
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
  [gps.x,gps.y] = projfwd(proj,lat_G_r_corrected,lon_G_r_corrected);
  gps.x = gps.x / 1000;
  gps.y = gps.y / 1000;
  hold on;
  
  scatter(gps.x,gps.y,20,variable_attenuation,'fill')
  %caxis([-15 15])
  colorbar;
  title('Total Variable Attenuation')
  
  %% Plot Value of DN
    figure(7)
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
  
    scatter(gps.x,gps.y,20,modified_Na,'fill')
    %caxis([-15 15])
    colorbar;
    title('Value of mod Na ')
end
keyboard
save_ref_en=0;
if save_ref_en
  
  out.Latitude=lat_G_r_corrected;
  out.Longitude=lon_G_r_corrected;
  out.const_attenuation=constant_attenuation;
  out.var_attenuation=variable_attenuation;
  out.Depth=depth_G_r_corrected;
  out.Refl_const=c_ref;
  out.Refl_var=v_ref;
  out.line_no=line_no;
  out.estimated_Na=estimated_Na;
  out.estimated_DN=estimated_DN;
  out.modified_Na=modified_Na;
  
  if strcmp(settings.location,'Peterman')
    save(['/cresis/snfs1/scratch/manjish/new_peterman/crossline_reflectivity_median_att2_new.mat'],'out')
  else
    save(['/cresis/snfs1/scratch/manjish/new_jacobshavn/cross_line_reflectivity_median_att2_new.mat'],'out')
  end
end
