source='V:\Thesis\basal_melt_state_oswald\RDBTS4_Greenland_1993_2013_01_basal_thermal_state.nc';
melt_rate= ncread(source,'basal_melt_rate');
melt_state=ncread(source,'likely_basal_thermal_state');
x=ncread(source,'x');
y=ncread(source,'y');

x_l=x(62,:);
y_l=y(62,:);
geotiff_fn ='X:\GIS_data\greenland\Landsat-7\Greenland_natural.tif';
proj = geotiffinfo(geotiff_fn);
[A CMAP R]= geotiffread(geotiff_fn);
figure
mapshow(rgb2gray(A),CMAP/1e3);
[gps.x,gps.y] = projfwd(proj,x_l,y_l);
hold on;
scatter(gps.x,gps.y,melt_state(62,:))
colorbar