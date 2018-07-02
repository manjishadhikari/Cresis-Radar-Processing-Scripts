%Script: estimation_of_relative_reflectivity_values using englacial atten
%method 2.. FInd Na from MSE and fit DN for every line


clear 
close all
clc
dbstop error

%% calculating Na avg
% load(['C:\Users\s343m141\Documents\scripts\matlab\thesis\ice_loss_estimation_paper_data\after_roughness_loss_correction\get heights frames Greenland\Greenland_layerdata_selected_frames_complete_v6.mat'])
if 1
    in_fn=['/cresis/snfs1/scratch/manjish/new_peterman/combined_data.mat'];
    load(in_fn);
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
Greenland.ice_bed_power_rgl_sorted = Greenland.ice_bed_power_cgl_sorted./median(Greenland.ice_bed_power_cgl_sorted);
if plots
    figure(6);plot(Greenland.depth_sorted, lp((Greenland.ice_bed_power_rgl_sorted)));
    grid on
    title('After Mean Removed')
end
avg_power = lp(median(Greenland.ice_bed_power_cgl_sorted));  %6.3224 dB
max_power=max(lp(Greenland.ice_bed_power_cgl_sorted));
avg_depth = mean(Greenland.depth_sorted);     %1.6033 km


%%
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
    
    %%
    figure
    subplot(2,1,1)
    plot((1:length(Greenland.depth_sorted)),(Greenland.depths_sorted-mean(Greenland.depths_sorted))/1e3);
    subplot(2,1,2)
    plot((1:length(Greenland.ice_bed_power_frgl_sorted)),(Greenland.ice_bed_power_frgl_sorted));
   
end
%%
na = 0:0.01:15;
for N= 1:length(na)
    S(N) = mean((((2*na(N)*(Greenland.depth_sorted-mean(Greenland.depth_sorted))/1e3)+((Greenland.ice_bed_power_frgl_sorted))).^2));
    %       S(N) = mode(round(abs((2*N*(Greenland.depths-mean(Greenland.depths))/1e3)+(10*log10(Greenland.ice_bed_power_frgc)))));
    %        S(N) = sum(abs((2*N*(Greenland.depths-mean(Greenland.depths))/1e3)+(real(Greenland.ice_bed_power_frgc))).^2);
end
[v i] = min(S)
if i==length(S)
    warning('check this')
end
Na = na(i)