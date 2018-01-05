
close all
load('Y:\manjish\peterman\radar_w_index1\verticalline8.mat')




max_length=length(Greenland.Latitude)
start_idx=5000;
stop_idx=280000;

Greenland.segments_length(1)=Greenland.segments_length(1)-start_idx+1;
Greenland.segments_length(end)=Greenland.segments_length(end)-length(Greenland.Latitude)-stop_idx;

Greenland.GPS_time=Greenland.GPS_time(start_idx:stop_idx);
Greenland.Latitude=Greenland.Latitude(start_idx:stop_idx);
Greenland.Longitude=Greenland.Longitude(start_idx:stop_idx);
Greenland.Elevation=Greenland.Elevation(start_idx:stop_idx);
Greenland.ice_bed_time=Greenland.ice_bed_time(start_idx:stop_idx);
Greenland.surface_time=Greenland.surface_time(start_idx:stop_idx);
Greenland.ice_bed_power=Greenland.ice_bed_power(start_idx:stop_idx);
Greenland.ice_surface_power=Greenland.ice_surface_power(start_idx:stop_idx);


geotiff_fn ='X:\GIS_data\greenland\Landsat-7\Greenland_natural.tif';
proj = geotiffinfo(geotiff_fn);
[A CMAP R]= geotiffread(geotiff_fn);
figure
mapshow(rgb2gray(A),CMAP/1e3);
xlim([-300 -50]);
ylim([-1550 -950]);
xlabel('X (km)');
ylabel('Y (km)');
hold on;
[gps.x,gps.y] = projfwd(proj,Greenland.Latitude,Greenland.Longitude);
gps.x = gps.x/1000 ;
gps.y = gps.y/1000 ;
grid
title('Location of lines');
xlabel('Latitude(deg)');
ylabel('Longitude(deg)');
hold on;
scatter(gps.x,gps.y,20,lp(Greenland.ice_bed_power),'fill')
colorbar
%     if length(bt.coherenceindex.value) >1
%         caxis([min(bt.coherenceindex.value),max(bt.coherenceindex.value)]);
%     elseif length(bt.coherenceindex.value) == 1
%         caxis([bt.coherenceindex.value-50,bt.coherenceindex.value+50]);
%     end
figure;plot(lp(Greenland.ice_bed_power));
range(lp(Greenland.ice_bed_power))

Greenland.truncated.bins=[start_idx stop_idx];
Greenland.truncated.original_length=max_length;
if 1
  save('Y:\manjish\peterman\radar_w_index1\verticalline8.mat','Greenland');
end