clear all
close all
clc
dbstop error

%% calculating Na avg
% load(['C:\Users\s343m141\Documents\scripts\matlab\thesis\ice_loss_estimation_paper_data\after_roughness_loss_correction\get heights frames Greenland\Greenland_layerdata_selected_frames_complete_v6.mat'])
if 1
    out_fn=['/cresis/snfs1/scratch/manjish//peterman/completedata_w_idx.mat'];
    load(out_fn);
end
physical_constants
plots =1;

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

[Greenland.depth_sorted Greenland.index]= sort(Greenland.depth);
Greenland.ice_bed_power_sorted = Greenland.ice_bed_power(Greenland.index);
Greenland.surface_height_sorted = Greenland.surface_height(Greenland.index);
Greenland.Latitude_sorted = Greenland.Latitude(Greenland.index);
Greenland.Longitude_sorted = Greenland.Longitude(Greenland.index);

if plots
    figure;plot(Greenland.depth_sorted, 10*log10(abs(Greenland.ice_bed_power_sorted).^2));
    grid on
    title('Depth vd Ice Bed Power after sorting')
end
geometric_loss_sorted = (2*(Greenland.surface_height_sorted+Greenland.depth_sorted/sqrt(er_ice))).^2;
Greenland.ice_bed_power_cgl_sorted=Greenland.ice_bed_power_sorted.*geometric_loss_sorted;
%Greenland.ice_bed_power_cgl_sorted = Greenland.ice_bed_power_cgl(Greenland.index);

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
Greenland.ice_bed_power_rgl_sorted = Greenland.ice_bed_power_cgl_sorted./mean(Greenland.ice_bed_power_cgl_sorted);
if plots
    figure(6);plot(Greenland.depth_sorted, lp((Greenland.ice_bed_power_rgl_sorted)));
    grid on
    title('After Mean Removed')
end
median_power = lp(median(Greenland.ice_bed_power_cgl_sorted))  %6.3224 dB
mean_power=lp(mean(Greenland.ice_bed_power_cgl_sorted))
max_power=max(lp(Greenland.ice_bed_power_cgl_sorted))
avg_depth = mean(Greenland.depth_sorted)     %1.6033 km


[r,m,b]=regression(Greenland.depth_sorted, lp((Greenland.ice_bed_power_rgl_sorted)));
val=m*Greenland.depth+b;
figure(6); hold on; plot(Greenland.depth,val);
%Na=10^3*10*log10(exp(1))/m; %
Na=(val(end)-val(1))/(Greenland.depth(end)-Greenland.depth(1))*1000    %One way depth averaged attenuation rate

ice_bed_reflectivity=lp(Greenland.ice_bed_power_rgl_sorted)+2*Na*(Greenland.depth_sorted-avg_depth)/1000;
figure;histogram(ice_bed_reflectivity);
keyboard

%Decimate
lat=decimate(Greenland.Latitude_sorted,100);
lon=decimate(Greenland.Longitude_sorted,100);
ice_bed_refl=decimate(ice_bed_reflectivity,100);


%% Plotting on map
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
[gps.x,gps.y] = projfwd(proj,lat,lon);

gps.x = gps.x / 1000;
gps.y = gps.y / 1000;
hold on;


scatter(gps.x,gps.y,20,ice_bed_refl,'fill')
caxis([-15 15])
colorbar;
title('Reflectivity using constant na ')

%Histogram
figure(2), hist(ice_bed_reflectivity,30);
title('Reflectivity using constant na ')
