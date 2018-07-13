clear 
close all
clc
dbstop error

%% Plotting on map
geotiff_fn = '/cresis/snfs1/dataproducts/GIS_data/greenland/Landsat-7/Greenland_natural.tif';
proj = geotiffinfo(geotiff_fn);
%proj = geotiffinfo('X:\GIS_data\antarctica\Landsat-7\Antarctica_LIMA_480m.tif');

[A CMAP R]= geotiffread(geotiff_fn);

%% calculating Na avg
% load(['C:\Users\s343m141\Documents\scripts\matlab\thesis\ice_loss_estimation_paper_data\after_roughness_loss_correction\get heights frames Greenland\Greenland_layerdata_selected_frames_complete_v6.mat'])
if 1
  out_fn=['/cresis/snfs1/scratch/manjish/new_jacobshavn/data/combined_data.mat'];
  load(out_fn);
end
physical_constants
plots =1;

%Roll Criteria
%  Greenland.roll=Greenland.Roll*180/pi;
%  idx= Greenland.roll>5;
%  Greenland.ice_bed_power(idx)=nan;

 %Depth Criteria
%  Greenland.depth = (Greenland.ice_bed_time - Greenland.surface_time)*c/2/sqrt(er_ice);
%  clear idx
%  Greenland.ice_bed_power((Greenland.depth>500))=nan;
 
 %Coherent Integration
numofCohInt=0;
if numofCohInt~=0
      [Greenland]=coh_integration(Greenland,numofCohInt);
end


clear idx
idx = find(isnan(Greenland.ice_bed_power)) ;
Greenland.GPS_time(idx) = [];
Greenland.Latitude(idx) = [];
Greenland.Longitude(idx) = [];
Greenland.Elevation(idx) = [];
Greenland.ice_bed_time(idx) = [];
Greenland.surface_time(idx) = [];
Greenland.ice_bed_power(idx) = [];
Greenland.ice_surface_power(idx) = [];

clear idx
idx = find(isnan(Greenland.ice_surface_power)) ;
Greenland.GPS_time(idx) = [];
Greenland.Latitude(idx) = [];
Greenland.Longitude(idx) = [];
Greenland.Elevation(idx) = [];
Greenland.ice_bed_time(idx) = [];
Greenland.surface_time(idx) = [];
Greenland.ice_bed_power(idx) = [];
Greenland.ice_surface_power(idx) = [];

Greenland.depth = (Greenland.ice_bed_time - Greenland.surface_time)*c/2/sqrt(er_ice);
Greenland.surface_height = (Greenland.surface_time)*c/2;


if plots
  plot(Greenland.depth, lp((Greenland.ice_bed_power)));
  grid on
  title('Depth vs Ice Bed Power')
  % hist(10*log10(Greenland.ice_bed_power),40)
end

[Greenland.depth_sorted Greenland.sortindex]= sort(Greenland.depth);
Greenland.ice_bed_power_sorted = Greenland.ice_bed_power(Greenland.sortindex);
Greenland.surface_height_sorted = Greenland.surface_height(Greenland.sortindex);
Greenland.Latitude_sorted = Greenland.Latitude(Greenland.sortindex);
Greenland.Longitude_sorted = Greenland.Longitude(Greenland.sortindex);

if plots
  figure;plot(Greenland.depth_sorted, 10*log10(abs(Greenland.ice_bed_power_sorted).^2));
  grid on
  title('Depth vs Ice Bed Power after sorting')
  xlabel('Depth in metres'); ylabel('Ice Bed Power in dB')
end
geometric_loss_sorted = (2*(Greenland.surface_height_sorted+Greenland.depth_sorted/sqrt(er_ice))).^2;
Greenland.ice_bed_power_cgl_sorted=(abs(Greenland.ice_bed_power_sorted).^2).*geometric_loss_sorted;
%Greenland.ice_bed_power_cgl_sorted = Greenland.ice_bed_power_cgl(Greenland.sortindex);

if plots
  figure;
  plot(Greenland.depth, 10*log10(abs(Greenland.ice_bed_power).^2));
  grid on
  title('Ice Bed Power vs Depth'); xlabel('Depth')
  
  figure;
  plot(10*log10(abs(Greenland.ice_bed_power).^2));
  grid on
  title('Ice Bed Power vs Along track'); xlabel('Along Track')
  
  figure; plot(Greenland.depth);
  title('Along track vs Depth');
  xlabel('Along Track'); ylabel('Depth')
  
  figure;plot(Greenland.depth_sorted,lp(Greenland.ice_bed_power_cgl_sorted));
  title('Ice Bed Power afer Geometric Loss Corrected')
  xlabel('Depth in m'); ylabel('Power in dB')
end

clear idx
idx = find(isnan(Greenland.ice_bed_power_cgl_sorted)) ;
Greenland.GPS_time(idx) = [];
Greenland.Latitude(idx) = [];
Greenland.Longitude(idx) = [];
Greenland.Elevation(idx) = [];
Greenland.ice_bed_time(idx) = [];
Greenland.surface_time(idx) = [];
Greenland.ice_bed_power_cgl_sorted(idx) = [];
Greenland.depth_sorted (idx)= [];

%Mean -141 for 20110429_01_028
Greenland.ice_bed_power_rgl_sorted = Greenland.ice_bed_power_cgl_sorted./median(Greenland.ice_bed_power_cgl_sorted);
if plots
  figure(6);plot(Greenland.depth_sorted, lp((Greenland.ice_bed_power_rgl_sorted)));
  grid on
  title('Relative Ice Bed Power');  xlabel('Depth in m'); ylabel('Power in dB')
end
median_power = lp(median(Greenland.ice_bed_power_cgl_sorted))  %6.3224 dB
mean_power=lp(mean(Greenland.ice_bed_power_cgl_sorted))
max_power=max(lp(Greenland.ice_bed_power_cgl_sorted))
min_power=min(lp(Greenland.ice_bed_power_cgl_sorted))
dynamic_range=max_power-min_power
avg_depth = mean(Greenland.depth_sorted)     %1.6033 km


[r,m,b]=regression(-2*(Greenland.depth_sorted-mean(Greenland.depth_sorted)),lp((Greenland.ice_bed_power_rgl_sorted)));
val=m*(-2*(Greenland.depth_sorted-mean(Greenland.depth_sorted)))+b;
figure(60); plot(2*(Greenland.depth_sorted-mean(Greenland.depth_sorted)),lp((Greenland.ice_bed_power_rgl_sorted)));
figure(60); hold on; plot(2*(Greenland.depth_sorted-mean(Greenland.depth_sorted)),val);
%Na=10^3*10*log10(exp(1))/m; %
%Na=(val(end)-val(1))/(Greenland.depth_sorted(end)-Greenland.depth_sorted(1))*1000    %One way depth averaged attenuation rate
Na_reg=m*1000

xlabel('Relative Depth in m'); ylabel('Relative Ice Bed Power dB')
title('Regression Fit')

%%  Filtering data

window_size = 10000;
% Greenland.ice_bed_power_frgc = conv(Greenland.ice_bed_power_rgc,gausswin(window_size),'same');
Greenland.ice_bed_power_frgl_sorted  = (sgolayfilt((lp((Greenland.ice_bed_power_rgl_sorted))), 2,window_size+1, gausswin(window_size+1)));
%Greenland.ice_bed_power_frgl_sortedtst  = (sgolayfilt((lp((Greenland.ice_bed_power_rgl_sorted))), 2,window_size+1));

if plots
    figure;plot((1:length(Greenland.ice_bed_power_rgl_sorted)),lp(Greenland.ice_bed_power_rgl_sorted));
    hold on
    plot((1:length(Greenland.ice_bed_power_frgl_sorted)),(Greenland.ice_bed_power_frgl_sorted));
    title('Filtering Data')
end

window_size = 10000;
Greenland.depths_sorted  = sgolayfilt(Greenland.depth_sorted, 2,window_size+1, gausswin(window_size+1));
if plots
    figure;plot((1:length(Greenland.depth_sorted)),Greenland.depth_sorted);
    hold on
    plot((1:length(Greenland.depths_sorted)),(Greenland.depths_sorted));
    title('Filtered Depth ')
    
    
    figure
    subplot(2,1,1)
    plot((1:length(Greenland.depth_sorted)),(Greenland.depths_sorted-mean(Greenland.depths_sorted))/1e3);
    subplot(2,1,2)
    plot((1:length(Greenland.ice_bed_power_frgl_sorted)),(Greenland.ice_bed_power_frgl_sorted));
   
end

%%



Na = 1:0.05:25;    %Least Mean Square Error Method
for  j = 1:length(Na)
  mse(j) = mean((-(Greenland.ice_bed_power_frgl_sorted)- 2*Na(j)*(Greenland.depth_sorted-mean((Greenland.depth_sorted)))/1000).^2);
end
figure;plot(mse);
[~, index1] = min(mse);
Na_bar  = Na(index1)

figure;plot(lp(Greenland.ice_bed_power_rgl_sorted));
hold on; plot(-2*Na_bar*(Greenland.depth_sorted-mean((Greenland.depth_sorted)))/1000);
xlabel('Samples');ylabel('Relative Ice Bed Power dB')
title('Least Mean Square Error Fit Method')

ice_bed_reflectivity=lp(Greenland.ice_bed_power_rgl_sorted)+2*Na_bar*(Greenland.depth_sorted-avg_depth)/1000;
figure;histogram(ice_bed_reflectivity);
xlabel('Relative Reflectivity');ylabel('Frequency')
keyboard


%Ice Bed Reflectivity at 5dB/km att rate
ice_bed_refl_5dB=lp(Greenland.ice_bed_power_rgl_sorted)+2*5*(Greenland.depth_sorted-avg_depth)/1000;
ice_bed_refl_10dB=lp(Greenland.ice_bed_power_rgl_sorted)+2*10*(Greenland.depth_sorted-avg_depth)/1000;
ice_bed_refl_15dB=lp(Greenland.ice_bed_power_rgl_sorted)+2*15*(Greenland.depth_sorted-avg_depth)/1000;
ice_bed_refl_20dB=lp(Greenland.ice_bed_power_rgl_sorted)+2*20*(Greenland.depth_sorted-avg_depth)/1000;


figure;histogram(ice_bed_refl_5dB);
figure;histogram(ice_bed_refl_10dB);
figure;histogram(ice_bed_refl_15dB);
figure;histogram(ice_bed_refl_20dB);


%Decimate
% lat=decimate(Greenland.Latitude_sorted,100);
% lon=decimate(Greenland.Longitude_sorted,100);
% ice_bed_refl=decimate(ice_bed_reflectivity,100);
% ice_bed_refl_5=decimate(ice_bed_refl_5dB,100);
% ice_bed_refl_10=decimate(ice_bed_refl_10dB,100);
% ice_bed_refl_15=decimate(ice_bed_refl_15dB,100);
% ice_bed_refl_20=decimate(ice_bed_refl_20dB,100);
close all
 lat=Greenland.Latitude_sorted;
lon=Greenland.Longitude_sorted;
 ice_bed_refl=ice_bed_reflectivity;
 ice_bed_refl_5=ice_bed_refl_5dB;
ice_bed_refl_10=ice_bed_refl_10dB;
 ice_bed_refl_15=ice_bed_refl_15dB;
 ice_bed_refl_20=ice_bed_refl_20dB;



%% Reflectivity using constant na

figure(101)
mapshow(rgb2gray(A),CMAP/1e3);
xlabel('X (km)');
ylabel('Y (km)');
xlim([-350 -50]);
ylim([-1250 -900]);

hold on
clear gps.x gps.y
[gps.x,gps.y] = projfwd(proj,lat,lon);

gps.x = gps.x / 1000;
gps.y = gps.y / 1000;
hold on;


scatter(gps.x,gps.y,20,ice_bed_refl,'s','fill')
caxis([-15 15])
colorbar;
c=colorbar;
c.Label.String='Reflectivity';
title('Reflectivity using constant na ')

%5dB
figure(102)
mapshow(rgb2gray(A),CMAP/1e3);
xlabel('X (km)');
ylabel('Y (km)');
xlim([-350 -50]);
ylim([-1250 -900]);

hold on
clear gps.x gps.y
[gps.x,gps.y] = projfwd(proj,lat,lon);

gps.x = gps.x / 1000;
gps.y = gps.y / 1000;
hold on;


scatter(gps.x,gps.y,20,ice_bed_refl_5,'s','fill')
caxis([-20 20])
colorbar;
c=colorbar;
c.Label.String='Reflectivity';
title('Reflectivity using 5db/km ')

%10dB
figure(103)
mapshow(rgb2gray(A),CMAP/1e3);
xlabel('X (km)');
ylabel('Y (km)');
xlim([-350 -50]);
ylim([-1250 -900]);

hold on
clear gps.x gps.y
[gps.x,gps.y] = projfwd(proj,lat,lon);

gps.x = gps.x / 1000;
gps.y = gps.y / 1000;
hold on;

scatter(gps.x,gps.y,20,ice_bed_refl_10,'s','fill')
caxis([-20 20])
colorbar;
c=colorbar;
c.Label.String='Reflectivity';
title('Reflectivity using 10dB/km att rate ')

%15db

figure(104)
mapshow(rgb2gray(A),CMAP/1e3);
xlabel('X (km)');
ylabel('Y (km)');
xlim([-350 -50]);
ylim([-1250 -900]);

hold on
clear gps.x gps.y
[gps.x,gps.y] = projfwd(proj,lat,lon);

gps.x = gps.x / 1000;
gps.y = gps.y / 1000;
hold on;


scatter(gps.x,gps.y,20,ice_bed_refl_15,'s','fill')
caxis([-20 20])
colorbar;
c=colorbar;
c.Label.String='Reflectivity';
title('Reflectivity using15 dB/km att rate ')

%20dB
figure(106)
mapshow(rgb2gray(A),CMAP/1e3);
xlabel('X (km)');
ylabel('Y (km)');
xlim([-350 -50]);
ylim([-1250 -900]);

hold on
clear gps.x gps.y
[gps.x,gps.y] = projfwd(proj,lat,lon);

gps.x = gps.x / 1000;
gps.y = gps.y / 1000;
hold on;


scatter(gps.x,gps.y,20,ice_bed_refl_20,'s','fill')
caxis([-20 20])
colorbar;
c=colorbar;
c.Label.String='Reflectivity';
title('Reflectivity using 20 dB/km att rate ')

%Histogram
figure(2), hist(ice_bed_reflectivity,30);
title('Reflectivity using constant na ')
xlabel('Relative Reflectivity');ylabel('Frequency')


%coherence Index Map

idx=find(Greenland.index.coherence>0.2);
coh=Greenland.index.coherence(idx);
lat1=Greenland.index.Latitude_mean(idx);
lon1=Greenland.index.Longitude_mean(idx);


figure(107)
mapshow(rgb2gray(A),CMAP/1e3);
xlabel('X (km)');
ylabel('Y (km)');
xlim([-350 -50]);
ylim([-1250 -900]);
hold on
clear gps.x gps.y
[gps.x,gps.y] = projfwd(proj,lat1,lon1);
gps.x = gps.x / 1000;
gps.y = gps.y / 1000;
hold on;
scatter(gps.x,gps.y,20,coh,'s','fill')
%caxis([min(Greenland.index.coherence) max(Greenland.index.coherence)])
%caxis([0.05 0.08])
colorbar;
c=colorbar;
c.Label.String='Value';
title('Coherence Index ')

%Abruptive Index Map

figure(108)
mapshow(rgb2gray(A),CMAP/1e3);
xlabel('X (km)');
ylabel('Y (km)');
xlim([-350 -50]);
ylim([-1250 -900]);
hold on
clear gps.x gps.y
[gps.x,gps.y] = projfwd(proj,Greenland.index.Latitude_mean,Greenland.index.Longitude_mean);
gps.x = gps.x / 1000;
gps.y = gps.y / 1000;
hold on;
scatter(gps.x,gps.y,20,Greenland.index.abruptness,'s','fill')
%caxis([min(Greenland.index.abruptness) max(Greenland.index.abruptness)])
caxis([0.1 0.3])
colorbar;
c=colorbar;
c.Label.String='value';
title('Abruptive Index ')

%Adjusted Intensity
figure(109)
mapshow(rgb2gray(A),CMAP/1e3);
xlabel('X (km)');
ylabel('Y (km)');
xlim([-350 -50]);
ylim([-1250 -900]);
hold on
clear gps.x gps.y
[gps.x,gps.y] = projfwd(proj,Greenland.index.Latitude_mean,Greenland.index.Longitude_mean);
gps.x = gps.x / 1000;
gps.y = gps.y / 1000;
hold on;
scatter(gps.x,gps.y,20,Greenland.index.Padj,'s','fill')
%caxis([min(Greenland.index.abruptness) max(Greenland.index.abruptness)])
%caxis([0.1 0.3])
colorbar;
c=colorbar;
c.Label.String='value';
title('Adjusted Intensity ')



%Depth
figure(101)
mapshow(rgb2gray(A),CMAP/1e3);
xlabel('X (km)');
ylabel('Y (km)');
xlim([-350 -50]);
ylim([-1250 -900]);

hold on
clear gps.x gps.y
[gps.x,gps.y] = projfwd(proj,lat,lon);

gps.x = gps.x / 1000;
gps.y = gps.y / 1000;
hold on;
scatter(gps.x,gps.y,20,Greenland.depth,'s','fill')
%caxis([-15 15])
colorbar;
c=colorbar;
c.Label.String='Depth';
title('Depth ')

