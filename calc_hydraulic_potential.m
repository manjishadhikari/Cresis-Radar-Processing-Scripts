%Script: calculation_of_hydraulic_potential.m  
% 
% Purpose: This script is used to calculate the hydraulic potential and equipotential contours of hydraulic potential for the Greenland Glacier. 
% 
% Input data / processed on: Greenland_layerdata_selected_frames_<straight_line/cross_line_number> from the previous step 2 
% 
% Output/data products: Hydraulic potential map and map of equipotential contours of hydraulic potential for the Greenland Glacier 
% 
% see also cross_over_analysis_for_validation, clustering_of_basal_conditions  



clear
close
clc
dbstop error

pw = 1000;
g = 9.81;
pice = 910;
potential_avg = [];
Longitude_avg =[];
Latitude_avg =[];


%geotiff_fn = '/cresis/snfs1/dataproducts/GIS_data/antarctica/Landsat-7/Antarctica_LIMA_480m.tif';
geotiff_fn = '/cresis/snfs1/dataproducts/GIS_data/greenland/Landsat-7/Greenland_natural.tif';
proj = geotiffinfo(geotiff_fn);
%proj = geotiffinfo('X:\GIS_data\antarctica\Landsat-7\Antarctica_LIMA_480m.tif');

[A CMAP R]= geotiffread(geotiff_fn);

%%
settings.location='Jacobshavn';

%%

for M = 1:60
    
    param.radar.fs = 195000000;
%     if cross_lines
%         load(['/cresis/snfs1/scratch/santhosh/thesis/get heights frames Greenland/cross_lines/Greenland_layerdata_selected_frames_v6_' num2str(M,'%03d') '.mat'])
%     else
%         load(['/cresis/snfs1/scratch/santhosh/thesis/get heights frames Greenland/Greenland_layerdata_selected_frames_v6_' num2str(M,'%03d') '.mat'])
%     end
%     

param.radar.fc = 195000000;  %Center Frequency
    if strcmp(settings.location,'Jacobshavn')
      if M<26
        cross_lines = 1;
        M1=0;
        load(['/cresis/snfs1/scratch/manjish/new_jacobshavn/radar_w_idx_new/crossline',num2str(M)]);
      else
        M1=M-25;
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
    
    physical_constants
    %%
    
    
    
    id = ~(isfinite( Greenland.ice_bed_power));
    
    Greenland.ice_bed_time(id) = [];
    Greenland.Latitude(id) = [];
    Greenland.Longitude(id) = [];
    Greenland.surface_time(id) = [];
    Greenland.Elevation(id) = [];
    
    
    
    
    Greenland.depth = (Greenland.ice_bed_time - Greenland.surface_time)*c/2/sqrt(er_ice);
    Greenland.surface_height = (Greenland.surface_time)*c/2;
    surface_elevation = Greenland.Elevation - Greenland.surface_height;
    bed_elevation = surface_elevation - Greenland.depth;
    
    potential = pw*g*bed_elevation + pice*g*Greenland.depth;
    
    k = 1;
    
    for l = 501:500:length(Greenland.ice_bed_time)
        if ((l > 500) && ((l+500) < length(Greenland.ice_bed_time)))
            
            Greenland.Latitude_avg(k) = nanmean(Greenland.Latitude((l-500):(l+499)));
            Greenland.Longitude_avg(k) = nanmean(Greenland.Longitude((l-500):(l+499)));
            Greenland.potential_avg(k) = nanmean(potential((l-500):(l+499)));
            k = k+1;
        end
        
    end
    
    
    
    potential_avg = cat(2,potential_avg, Greenland.potential_avg);
    Latitude_avg = cat(2,Latitude_avg, Greenland.Latitude_avg);
    Longitude_avg = cat(2,Longitude_avg, Greenland.Longitude_avg);
    
    
    % figure(1)
    % hold on
    % clear gps.x gps.y
    % [gps.x,gps.y] = projfwd(proj,Greenland.Latitude_avg,Greenland.Longitude_avg);
    % %         lt = cat(2,lt, Greenland.Latitude_avg ) ;
    % %         ln =cat(2, ln,  Greenland.Longitude_avg) ;
    % gps.x = gps.x / 1000;
    % gps.y = gps.y / 1000;
    % figure(1)
    % hold on;
    % scatter(gps.x,gps.y,20,Greenland.potential_avg,'fill')
    % %scatter(gps.x,gps.y,20,2.*estimated_Na.*depth_G_r_corrected,'fill')
    % ylim([-1000 -655]);
    clear Greenland potential
end
cross_lines = 1;

% %%
% if cross_lines == 1
%     W = 29;
% else
%     W = 19;
% end
% for M = 1:W
%     
%     param.radar.fs = 195000000;
%     if cross_lines
%         load(['/cresis/snfs1/scratch/santhosh/thesis/get heights frames Greenland/cross_lines/Greenland_layerdata_selected_frames_v6_' num2str(M,'%03d') '.mat'])
%     else
%         load(['/cresis/snfs1/scratch/santhosh/thesis/get heights frames Greenland/Greenland_layerdata_selected_frames_v6_' num2str(M,'%03d') '.mat'])
%     end
%     physical_constants
%     %%
%     
%     
%     
%     id = ~(isfinite( Greenland.ice_bed_power));
%     
%     Greenland.ice_bed_time(id) = [];
%     Greenland.Latitude(id) = [];
%     Greenland.Longitude(id) = [];
%     Greenland.surface_time(id) = [];
%     Greenland.Elevation(id) = [];
%     
%     
%     
%     
%     Greenland.depth = (Greenland.ice_bed_time - Greenland.surface_time)*c/2/sqrt(er_ice);
%     Greenland.surface_height = (Greenland.surface_time)*c/2;
%     surface_elevation = Greenland.Elevation - Greenland.surface_height;
%     bed_elevation = surface_elevation - Greenland.depth;
%     
%     potential = pw*g*bed_elevation + pice*g*Greenland.depth;
%     
%     k = 1;
%     
%     for l = 251:500:length(Greenland.ice_bed_time)
%         if ((l > 250) && ((l+1000) < length(Greenland.ice_bed_time)))
%             
%             Greenland.Latitude_avg(k) = nanmean(Greenland.Latitude((l-250):(l+249)));
%             Greenland.Longitude_avg(k) = nanmean(Greenland.Longitude((l-250):(l+249)));
%             Greenland.potential_avg(k) = nanmean(potential((l-250):(l+249)));
%             k = k+1;
%         end
%         
%     end
%     
%     
%     
%     potential_avg = cat(2,potential_avg, Greenland.potential_avg);
%     Latitude_avg = cat(2,Latitude_avg, Greenland.Latitude_avg);
%     Longitude_avg = cat(2,Longitude_avg, Greenland.Longitude_avg);
%     
%     
%     % figure(1)
%     % hold on
%     % clear gps.x gps.y
%     % [gps.x,gps.y] = projfwd(proj,Greenland.Latitude_avg,Greenland.Longitude_avg);
%     % %         lt = cat(2,lt, Greenland.Latitude_avg ) ;
%     % %         ln =cat(2, ln,  Greenland.Longitude_avg) ;
%     % gps.x = gps.x / 1000;
%     % gps.y = gps.y / 1000;
%     % figure(1)
%     % hold on;
%     % scatter(gps.x,gps.y,20,Greenland.potential_avg,'fill')
%     % %scatter(gps.x,gps.y,20,2.*estimated_Na.*depth_G_r_corrected,'fill')
%     % ylim([-1000 -655]);
%     clear Greenland potential
% end

figure(1)
mapshow(rgb2gray(A),CMAP/1e3);
xlabel('X (km)');
ylabel('Y (km)');
% xlim([350 650]);
% ylim([-1000 -655]);

figure(1)
hold on
clear gps.x gps.y
[gps.x,gps.y] = projfwd(proj,Latitude_avg,Longitude_avg);
gps.x = gps.x / 1000;
gps.y = gps.y / 1000;
figure(1)
hold on
scatter(gps.x,gps.y,20,lp(potential_avg),'fill')
% caxis([-2 -0.3])


%%
X = min(gps.x):max(gps.x);
% Z = linspace(0,3800,length(X));
Y = min(gps.y):max(gps.y);
x = repmat(X,length(Y),1);
y = repmat(Y',1,length(X));
% x = repmat(Z,length(Y),1);

% [x, y] = meshgrid(X, Y);


interpolated = griddata(gps.x,gps.y,potential_avg,x,y);
imagesc(X,Y,interpolated)
axis xy
contour(x,y,interpolated)
[fx fy] = gradient(interpolated);
figure(1)
hold on
[cn, h] = contour(x,y,interpolated, 25)
% quiver(x,y,fx,fy)

id = find(fx<0 | -fy<0);
fx(id) = nan;
fy(id) = nan;
figure (1)
hold on
quiver(x,y,fx,-fy,40)
