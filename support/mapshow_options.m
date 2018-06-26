


%geotiff_fn = '/cresis/snfs1/dataproducts/GIS_data/greenland/Landsat-7/Greenland_natural.tif';
geotiff_fn = '/cresis/snfs1/dataproducts/GIS_data/greenland/Greenland1.tif';
  proj = geotiffinfo(geotiff_fn);
  %proj = geotiffinfo('X:\GIS_data\antarctica\Landsat-7\Antarctica_LIMA_480m.tif');
  
  [A CMAP R]= geotiffread(geotiff_fn);
  
  %% Reflectivity using constant na
  
  figure(1)
  mapshow(rgb2gray(A),CMAP/1e3);