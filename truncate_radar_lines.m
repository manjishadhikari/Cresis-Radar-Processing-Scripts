
clearvars;
flipped=1;
line_no=20;
close all
if 1
  %load([sprintf('Y:\\manjish\\peterman\\radar_w_index\\verticalline%d.mat',line_no)])
  load([sprintf('/cresis/snfs1/scratch/manjish/peterman/radar_w_index/verticalline%d.mat',line_no)]);
end

lm=length(Greenland.Latitude);
start_idx=1;
stop_idx=lm;

Greenland.segments_length(1)=Greenland.segments_length(1)-start_idx+1;
Greenland.segments_length(end)=Greenland.segments_length(end)-(length(Greenland.Latitude)-stop_idx);

Greenland.GPS_time=Greenland.GPS_time(start_idx:stop_idx);
Greenland.Latitude=Greenland.Latitude(start_idx:stop_idx);
Greenland.Longitude=Greenland.Longitude(start_idx:stop_idx);
Greenland.Elevation=Greenland.Elevation(start_idx:stop_idx);
Greenland.ice_bed_time=Greenland.ice_bed_time(start_idx:stop_idx);
Greenland.surface_time=Greenland.surface_time(start_idx:stop_idx);
Greenland.ice_bed_power=Greenland.ice_bed_power(start_idx:stop_idx);
Greenland.ice_surface_power=Greenland.ice_surface_power(start_idx:stop_idx);

%Plot after truncating

% geotiff_fn ='X:\GIS_data\greenland\Landsat-7\Greenland_natural.tif';
geotiff_fn='/cresis/snfs1/dataproducts/GIS_data/greenland/Landsat-7/Greenland_natural.tif';
proj = geotiffinfo(geotiff_fn);
[A CMAP R]= geotiffread(geotiff_fn);
figure(1)
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
scatter(gps.x(1),gps.y(1),200,'X');
colorbar
%     if length(bt.coherenceindex.value) >1
%         caxis([min(bt.coherenceindex.value),max(bt.coherenceindex.value)]);
%     elseif length(bt.coherenceindex.value) == 1
%         caxis([bt.coherenceindex.value-50,bt.coherenceindex.value+50]);
%     end


Greenland.truncated.bins=[start_idx stop_idx];
Greenland.truncated.original_length=lm;
Greenland.flipped=0;
Greenland
disp(sprintf('1 -->%d',Greenland.ice_bed_power(1)))
disp(sprintf('end -->%d',Greenland.ice_bed_power(end)))

if flipped==1
  Greenland.GPS_time=flip(Greenland.GPS_time);
  Greenland.Latitude=flip(Greenland.Latitude);
  Greenland.Longitude=flip(Greenland.Longitude);
  Greenland.Elevation=flip(Greenland.Elevation);
  Greenland.ice_bed_time=flip(Greenland.ice_bed_time);
  Greenland.surface_time=flip(Greenland.surface_time);
  Greenland.ice_bed_power=flip(Greenland.ice_bed_power);
  Greenland.ice_surface_power=flip(Greenland.ice_surface_power);
  Greenland.segments_length=flip(Greenland.segments_length);
  Greenland.flipped=1;
end

if flipped==1
  figure(2)
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
  title('Location of lines flipped');
  xlabel('Latitude(deg)');
  ylabel('Longitude(deg)');
  hold on;
  scatter(gps.x,gps.y,20,lp(Greenland.ice_bed_power),'fill')
  scatter(gps.x(1),gps.y(1),200,'X');
  colorbar  
  Greenland
  disp(sprintf('1 -->%d',Greenland.ice_bed_power(1)))
disp(sprintf('end -->%d',Greenland.ice_bed_power(end)))
end

figure;plot(lp(Greenland.ice_bed_power));
disp(sprintf('Range of ice bed power= %d',round(range(lp(Greenland.ice_bed_power)))))

keyboard
if 1
  fprintf('Saving ......\n')
%  save([sprintf('Y:\\manjish\\peterman\\radar_w_index_flipped\\verticalline%d.mat',line_no)],'Greenland');
   save([sprintf('/cresis/snfs1/scratch/manjish/peterman/radar_w_index_flipped/verticalline%d.mat',line_no)],'Greenland');
end