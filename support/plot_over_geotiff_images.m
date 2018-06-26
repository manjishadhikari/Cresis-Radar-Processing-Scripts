
close all
load(['/cresis/snfs1/scratch/manjish/new_jacobshavn/reflectivity_median_att2.mat'])
%geotiff_fn ='/cresis/snfs1/dataproducts/GIS_data/greenland/Landsat-7/Greenland_natural.tif';
% geotiff_fn='/cresis/snfs1/dataproducts/GIS_data/greenland/plummer_jakobshavn/jak_grid.tif';
% geotiff_fn='/cresis/snfs1/dataproducts/GIS_data/greenland/Jakobshavn/20100709_20100711_mos.tif';
% geotiff_fn='/users/manjish/maps/mog100_2015_grn_v02.tif';
% geotiff_fn='/users/manjish/maps/mog500_2005_gct_v1.1.tif';
%geotiff_fn='/users/manjish/maps/jakob.tif';
% geotiff_fn='/users/manjish/maps/peterman.tif';

  geotiff_fn='/cresis/snfs1/dataproducts/GIS_data/ops/data/geoserver/data/arctic/raster/greenland_bamberV3_bed.tif'; %Bed Elevation Bamber

%   geotiff_fn='/cresis/snfs1/dataproducts/GIS_data/ops/data/geoserver/data/arctic/raster/greenland_bw_45m.tif'; %Natural earth
% geotiff_fn='/cresis/snfs1/dataproducts/GIS_data/ops/data/geoserver/data/arctic/raster/greenland_joughin2011_velocity_magnitude.tif';

proj = geotiffinfo(geotiff_fn);
 [A CMAP R]= geotiffread(geotiff_fn);
figure
%  mapshow(rgb2gray(A),CMAP/1e3);
plot_geotiff(geotiff_fn);
xlim([-220 -50]);
ylim([-2530 -2190]);
xlabel('X (km)');
ylabel('Y (km)');
hold on;
[gps.x,gps.y] = projfwd(proj,out.Latitude,out.Longitude);
gps.x = gps.x/1000 ;
gps.y = gps.y/1000 ;
grid
title('Location of lines');
xlabel('Latitude(deg)');
ylabel('Longitude(deg)');
hold on;
scatter(gps.x,gps.y,20,out.Refl_const,'d','fill')
colorbar
%     if length(bt.coherenceindex.value) >1
%         caxis([min(bt.coherenceindex.value),max(bt.coherenceindex.value)]);
%     elseif length(bt.coherenceindex.value) == 1
%         caxis([bt.coherenceindex.value-50,bt.coherenceindex.value+50]);
%     end
