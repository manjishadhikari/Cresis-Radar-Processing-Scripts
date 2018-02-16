  
function [success]=laser_radar_roughness_comparision_task(param)


c=3e8;

lat=[];
lon=[];
%laser_elev=[];
laser_rms=[];
radar_rms=[];


debug_flag=0;
cross_line=param.cross_lines_en;
num_int=1000;
coh_int=0;
sf_bin=[0 0]; %[Min Max] [-1 1]
bypass=0;
save_en=0;
if save_en
  warning('Save enabled')
  %keyboard
end
sgolay_filter_length=0;

if ispc
  base_path_sr='Y:';
  base_path_dp='X:';  
  git_path=fullfile('H:\scripts\matlab\total process');

else
  base_path_sr='/cresis/snfs1/scratch';
  base_path_dp='/cresis/snfs1/dataproducts';
  git_path=fullfile('/users/manjish/scripts/matlab/total process');
end

if cross_line==0
  Day_seg={'20100324_01','20110429_01','20110429_02','20110507_02','20130420_02','20140512_01','20140505_01'};
  Day=repelem(Day_seg,[7,2,4,2,3,2,1]);
  direction=[1,1,1,1,1,1,1,-1,1,1,-1,1,-1,1,1,1,1,1,1,1]; %Increasing lat or lon is 1 else 0
  frms={[11 12],[14 15],[17 18],[20 21],[23 24] [30 31],[33 34],[30:32],[33 34],[18:21],[12:15],[2:5],[6:9],[6:8],[17:20],[3 4],[9],[11],[12 13],[17 18],[15 16]};
else
  Day_seg={'20100324_01','20100324_02','20100324_03','20100324_04','20110429_01','20110429_02','20110507_01','20110507_02','20120330_01','20120516_01','20140512_01','20140505_01'};
  Day=repelem(Day_seg,[3,1,1,1,5,1,7,1,1,2,1,5]);
  frms={[36,37],[39,40],[42],[1 2],[1 2],[1 2],[9:12],[13:16],[17:20],[21:24],[25:28],[10:11],[10:14],[15:18],[19:22],[23:26],[27:30],[31:34],[35:37],[1:4],[5 6 7],[13:16],[79:81],[10 11],[12 13],[33 34 35],[38:40],[55 56],[59 60]};
  direction=[1,1,1,1,1,1,1,1,-1,1,-1,-1,1,1,1,1,-1,1,-1,-1,1,1,1,1,1,1,1,1,1];
end



for lno=param.proc_line
  
  
  %load('Y:\manjish\2014_Antarctica_DC8\verticalline2.mat');
  % load(['Y:\manjish\peterman\laser_icessn\verticalline',sprintf('%s',num2str(lno)),'.mat']);   %laser line
  if cross_line==1
    load(fullfile(base_path_sr,'manjish','peterman','laser_icessn', ['crossline',sprintf('%s.mat',num2str(lno))]));
    %load(['Y:\manjish\peterman\laser_icessn\crossline',sprintf('%s',num2str(lno)),'.mat']);   %laser line
  else
     load(fullfile(base_path_sr,'manjish','peterman','laser_icessn', ['verticalline',sprintf('%s.mat',num2str(lno))]));
    % load(['Y:\manjish\peterman\laser_icessn\verticalline',sprintf('%s',num2str(lno)),'.mat']);   %laser line
  end
  % load(['Y:\manjish\peterman\radar\crossline',sprintf('%s',num2str(lno)),'.mat']);          %Radar line
  %   load(['Y:\manjish\peterman\ricefitsurfaceroughness\crossline',sprintf('%s',num2str(lno)),'.mat']);
  % load(['Y:\manjish\peterman\ricefitsurfaceroughness\verticalline',sprintf('%s',num2str(lno)),'.mat']);
  %load(['Y:\manjish\peterman\radarnew\verticalline',sprintf('%s',num2str(lno)),'.mat']);     %new radar line
  
  %    load(['Y:\manjish\2014_Antarctica_DC8\verticalline',sprintf('%s',num2str(lno)),'.mat']);
  %   load(['Y:\manjish\2014_Antarctica_DC8\verticalline',sprintf('%s',num2str(lno)),'.mat']);
  
  if bypass==0
    for p=1:length(frms{lno})
      
      if (cross_line==0 & lno<8 )|(cross_line==1 & lno<7)    %2010
        datapath{p}=(fullfile(base_path_dp,'ct_data','rds','2010_Greenland_DC8','manjish/CSARP_Data',Day{lno},[sprintf('Data_%s_%03d.mat',Day{lno},frms{lno}(p))]));
        layer{p}=(fullfile(base_path_dp,'ct_data','rds','2010_Greenland_DC8','old','CSARP_layerData',Day{lno},[sprintf('Data_%s_%03d.mat',Day{lno},frms{lno}(p))]));
        
       % layer{p}=(['X:\ct_data\rds\2010_Greenland_DC8\CSARP_layerData\',Day{lno},'\',sprintf('Data_%s_%03d.mat',Day{lno},frms{lno}(p))]);
%         prev_path{p}=(['X:\ct_data\rds\2010_Greenland_DC8\manjish\CSARP_Data\',Day{lno},'\',sprintf('Data_%s_%03d.mat',Day{lno},frms{lno}(p)-1)]);
%         prev_layer{p}=(['X:\ct_data\rds\2010_Greenland_DC8\CSARP_layerData\',Day{lno},'\',sprintf('Data_%s_%03d.mat',Day{lno},frms{lno}(p)-1)]);
%         post_path{p}=(['X:\ct_data\rds\2010_Greenland_DC8\manjish\CSARP_Data\',Day{lno},'\',sprintf('Data_%s_%03d.mat',Day{lno},frms{lno}(p)+1)]);
%         post_layer{p}=(['X:\ct_data\rds\2010_Greenland_DC8\CSARP_layerData\',Day{lno},'\',sprintf('Data_%s_%03d.mat',Day{lno},frms{lno}(p)+1)]);
%         
      elseif (cross_line==0 & lno>7 &lno<16 )| (cross_line==1 & lno>6 & lno<21)    %2011
         datapath{p}=(fullfile(base_path_dp,'ct_data','rds','2011_Greenland_P3','CSARP_manjish',Day{lno},[sprintf('Data_%s_%03d.mat',Day{lno},frms{lno}(p))]));
        layer{p}=(fullfile(base_path_dp,'ct_data','rds','2011_Greenland_P3','CSARP_layerData',Day{lno},[sprintf('Data_%s_%03d.mat',Day{lno},frms{lno}(p))]));
       
      elseif (cross_line==1 & lno>20 & lno<24)  %2012
         datapath{p}=(fullfile(base_path_dp,'ct_data','rds','2012_Greenland_P3','CSARP_CSARP_manjish',Day{lno},[sprintf('Data_%s_%03d.mat',Day{lno},frms{lno}(p))]));
        layer{p}=(fullfile(base_path_dp,'ct_data','rds','2012_Greenland_P3','CSARP_post','CSARP_layerData',Day{lno},[sprintf('Data_%s_%03d.mat',Day{lno},frms{lno}(p))]));
        
      elseif (cross_line==0 & lno>15 &lno<19 ) %2013
         datapath{p}=(fullfile(base_path_dp,'ct_data','rds','2013_Greenland_P3','manjish','CSARP_data',Day{lno},[sprintf('Data_%s_%03d.mat',Day{lno},frms{lno}(p))]));
        layer{p}=(fullfile(base_path_dp,'ct_data','rds','2013_Greenland_P3','CSARP_post','CSARP_layerData',Day{lno},[sprintf('Data_%s_%03d.mat',Day{lno},frms{lno}(p))]));
        
      else    %2014
         datapath{p}=(fullfile(base_path_dp,'ct_data','rds','2014_Greenland_P3','CSARP_manjish',Day{lno},[sprintf('Data_%s_%03d.mat',Day{lno},frms{lno}(p))]));
        layer{p}=(fullfile(base_path_dp,'ct_data','rds','2014_Greenland_P3','CSARP_post','CSARP_layerData',Day{lno},[sprintf('Data_%s_%03d.mat',Day{lno},frms{lno}(p))]));
        
        
        
%        datapath{p}=(['X:\ct_data\rds\2011_Greenland_P3\CSARP_manjish\',Day{lno},'\',sprintf('Data_%s_%03d.mat',Day{lno},frms{lno}(p))]);
%        layer{p}=(['X:\ct_data\rds\2011_Greenland_P3\CSARP_layerData\',Day{lno},'\',sprintf('Data_%s_%03d.mat',Day{lno},frms{lno}(p))]);
%       
%         prev_path{p}=(['X:\ct_data\rds\2011_Greenland_P3\CSARP_manjish\',Day{lno},'\',sprintf('Data_%s_%03d.mat',Day{lno},frms{lno}(p)-1)]);
%         prev_layer{p}=(['X:\ct_data\rds\2011_Greenland_P3\CSARP_layerData\',Day{lno},'\',sprintf('Data_%s_%03d.mat',Day{lno},frms{lno}(p)-1)]);
%          post_path{p}=(['X:\ct_data\rds\2011_Greenland_P3\CSARP_manjish\',Day{lno},'\',sprintf('Data_%s_%03d.mat',Day{lno},frms{lno}(p)+1)]);
%         post_layer{p}=(['X:\ct_data\rds\2011_Greenland_P3\CSARP_layerData\',Day{lno},'\',sprintf('Data_%s_%03d.mat',Day{lno},frms{lno}(p)+1)]);
%      
      end
    end
    
    r=roughness_cal_from_radar_v2(lno,datapath,layer,num_int,coh_int,sf_bin,save_en,sgolay_filter_length,cross_line,direction(lno),debug_flag,git_path);
  else
    if cross_line==1
        load(fullfile(base_path_sr,'manjish','peterman','radarnew', ['crossline',sprintf('%s.mat',num2str(lno))]));
        
     % load(['Y:\manjish\peterman\radarnew\crossline',sprintf('%s',num2str(lno)),'.mat']);
    else
        load(fullfile(base_path_sr,'manjish','peterman','radarnew', ['verticalline',sprintf('%s.mat',num2str(lno))]));
  
      % load(['Y:\manjish\peterman\radarnew\verticalline',sprintf('%s',num2str(lno)),'.mat']);
    end
  end
  
  %keyboard
  close all;
 % load([sprintf('Y:\\manjish\\peterman\\sardata\\crossline_roughness%d',lno)])
  if 1
  %latitude_radar_ds=downsample(Greenland.Latitude,1);
  %longitude_radar_ds=downsample(Greenland.Longitude,1);
  %latitude_laser_ds=downsample(latitude_laser,100);
  %longitude_laser_ds=downsample(longitude_laser,100);
  geotiff_fn =fullfile(base_path_dp,'GIS_data','greenland','Landsat-7','Greenland_natural.tif');
  %geotiff_fn='X:\GIS_data\antarctica\Landsat-7\Antarctica_LIMA_480m.tif';
  proj = geotiffinfo(geotiff_fn);
  [A CMAP R]= geotiffread(geotiff_fn);
  
  %% Plot all three
  if debug_flag
    figure(2)
    mapshow(rgb2gray(A),CMAP/1e3);
    % xlim([-300 -50]);
    % ylim([-2350 -2150]);
    xlabel('X (km)');
    ylabel('Y (km)');
    
    
    [gps.x,gps.y] = projfwd(proj,icessn.lat,icessn.lon);
    gps.x = gps.x/1000 ;
    gps.y = gps.y/1000 ;
    figure(2)
    hold on;
    grid
    title('Elevation');
    scatter(gps.x,gps.y,30,[1 0 0],'d')
    scatter(gps.x(1),gps.y(1),50,'c','+')
    
    %Laser line
    %   [gps.xr,gps.yr] = projfwd(proj,latitude_laser_ds,longitude_laser_ds);
    %   gps.xr = gps.xr/1000 ;
    %   gps.yr = gps.yr/1000 ;
    %   figure(2)
    %   hold on
    %   scatter(gps.xr,gps.yr,20,[0 0 1],'o')
    %   scatter(gps.xr(1),gps.yr(1),50,'c','+')
    
    %Radar line
    %     [gps.xrd,gps.yrd] = projfwd(proj,latitude_radar_ds,longitude_radar_ds);
    %     gps.xrd = gps.xrd/1000 ;
    %     gps.yrd = gps.yrd/1000 ;
    %     figure(2)
    %     hold on;
    %     scatter(gps.xrd,gps.yrd,20,[0 1 0],'h')
    %     scatter(gps.xrd(1),gps.yrd(1),50,'c','+')
    
    [gps.xrd,gps.yrd] = projfwd(proj,r.lat,r.lon);
    gps.xrd = gps.xrd/1000 ;
    gps.yrd = gps.yrd/1000 ;
    figure(2)
    hold on;
    scatter(gps.xrd,gps.yrd,40,'k','>')
    scatter(gps.xrd(1),gps.yrd(1),50,'c','+')
    title('All three laser radar and rms ')
  end
  
  %% Selecting TRack 0  from the ICESSN data
  
  idx=find(icessn.track_identifier==0);
  lat_1=icessn.lat(idx);
  lon_1=icessn.lon(idx);
  rms_1=icessn.rms_fit(idx);
  elev_1=icessn.elev(idx);
  
  idx2=find(icessn.track_identifier==1);
  lat_2=icessn.lat(idx2);
  lon_2=icessn.lon(idx2);
  rms_2=icessn.rms_fit(idx2);
  elev_2=icessn.elev(idx2);
  %
  idx4=find(icessn.track_identifier==4);
  lat_4=icessn.lat(idx4);
  lon_4=icessn.lon(idx4);
  rms_4=icessn.rms_fit(idx4);
  elev_4=icessn.elev(idx4);
  
  if debug_flag
    figure(3)
    
    mapshow(rgb2gray(A),CMAP/1e3);
    % xlim([-300 -50]);
    % ylim([-2350 -2150]);
    xlabel('X (km)');
    ylabel('Y (km)');
    
    clear gps.x gps.y
    [gps.x,gps.y] = projfwd(proj,lat_1,lon_1);  %%Track 0
    gps.x = gps.x/1000 ;
    gps.y = gps.y/1000 ;
    figure(3)
    hold on;
    grid
    title('Radar and Track 0 Line');
    scatter(gps.x,gps.y,20,[1 0 0],'d')
    scatter(gps.x(1),gps.y(1),50,'r','+')
    %
    figure(3)
    hold on;
    [gps.xrd,gps.yrd] = projfwd(proj,lat_2,lon_2);   %Track 2
    gps.xrd = gps.xrd/1000 ;
    gps.yrd = gps.yrd/1000 ;
    scatter(gps.xrd,gps.yrd,20,[0 1 0],'*')
    scatter(gps.xrd(1),gps.yrd(1),50,'r','+')
    
    %     figure(3)
    %     hold on;
    %     [gps.xrd,gps.yrd] = projfwd(proj,latitude_radar_ds,longitude_radar_ds); %Radar line
    %     gps.xrd = gps.xrd/1000 ;
    %     gps.yrd = gps.yrd/1000 ;
    %     scatter(gps.xrd,gps.yrd,20,[0 0 1],'h')
    %      scatter(gps.xrd(1),gps.yrd(1),50,'r','+')
    
    figure(3)
    hold on;
    clear gps.xrd gps.yrd;
    [gps.xrd,gps.yrd] = projfwd(proj,lat_4,lon_4);  %Track 4
    gps.xrd = gps.xrd/1000 ;
    gps.yrd = gps.yrd/1000 ;
    scatter(gps.xrd,gps.yrd,20,[1 0 0],'o')
    % scatter(gps.xrd(1),gps.yrd(1),50,'r','+')
    
    figure(3)
    hold on;
    clear gps.xrd gps.yrd;
    [gps.xrd,gps.yrd] = projfwd(proj,r.lat,r.lon); %Radar line
    gps.xrd = gps.xrd/1000 ;
    gps.yrd = gps.yrd/1000 ;
    scatter(gps.xrd,gps.yrd,20,'k','>')
    % scatter(gps.xrd(1),gps.yrd(1),50,'r','+')
    title('Icessn tracks and rms ')
  end
  
  %%
  %only unique lat/lon points
  
  clear idx1 idx2;
  if cross_line==0
    [lat_1, un_idx]=unique(lat_1);
    lon_1=lon_1(un_idx);
    rms_1=rms_1(un_idx);
    elev_1=elev_1(un_idx);
    
    idx1=find(lat_1<=r.lat(1),1,'first');
    if isempty(idx1)
      idx1=1;
    end
    idx2=find(lat_1<=r.lat(end),1,'last');
    if isempty(idx2)
      idx2=length(lat_1);
    end
    
  else
    [lon_1, un_idx]=unique(lon_1);
    lat_1=lat_1(un_idx);
    rms_1=rms_1(un_idx);
    elev_1=elev_1(un_idx);
    if any((r.lon<180))
      r.lon=360+r.lon;
    end
    idx1=find(lon_1<=r.lon(1),1,'first');
     if isempty(idx1)
      idx1=1;
    end
    idx2=find(lon_1<=r.lon(end),1,'last');
    if isempty(idx2)
      idx2=length(lon_1);
    end
    
  end
  
  new_laser_lat=lat_1(idx1:idx2);
  new_laser_lon=lon_1(idx1:idx2);
  new_laser_elevation=elev_1(idx1:idx2);
  new_rms_l=rms_1(idx1:idx2);
  
  %%  if debug_flag
  if debug_flag
    figure(100)
    mapshow(rgb2gray(A),CMAP/1e3);
    % xlim([-300 -50]);
    % ylim([-2350 -2150]);
    xlabel('X (km)');
    ylabel('Y (km)');
    clear gps.x gps.y
    
    [gps.x,gps.y] = projfwd(proj,new_laser_lat,new_laser_lon);
    gps.x = gps.x/1000 ;
    gps.y = gps.y/1000 ;
    figure(100)
    hold on;
    grid
    title('Elevation');
    scatter(gps.x,gps.y,20,new_rms_l)
    
    
    [gps.xrd,gps.yrd] = projfwd(proj,r.lat,r.lon);
    gps.xrd = gps.xrd/1000 ;
    gps.yrd = gps.yrd/1000 ;
    figure(100)
    hold on;
    scatter(gps.xrd,gps.yrd,40,'k','>')
    scatter(gps.xrd(1),gps.yrd(1),50,'c','+')
    title('Truncate for unique rms lines')
  end
  
  %%
  %Remove large gaps
  dist_icessn=geodetic_to_along_track(new_laser_lat,new_laser_lon);
  diff_dist_icessn=diff(dist_icessn);
  if direction(lno)==1
    [d_idx]=find(diff_dist_icessn>10000,1,'first');
    
    if ~isempty(d_idx)
      new_laser_lat=new_laser_lat(1:d_idx-1);
      new_laser_lon=new_laser_lon(1:d_idx-1);
      new_rms_l=new_rms_l(1:d_idx-1);
      new_laser_elevation=new_laser_elevation(1:d_idx-1);
      
    end
  else
    [d_idx]=find(diff_dist_icessn>10000,1,'last');
    
    if ~isempty(d_idx)
      new_laser_lat=new_laser_lat(d_idx+1:end);
      new_laser_lon=new_laser_lon(d_idx+1:end);
      new_rms_l=new_rms_l(d_idx+1:end);
      new_laser_elevation=new_laser_elevation(d_idx+1:end);
      
    end
    
     dist_new=geodetic_to_along_track(new_laser_lat,new_laser_lon);
    diff_dist_new=diff(dist_new);
     clear d_idx;
      [d_idx]=find(diff_dist_new>8000,1,'last');
     if ~isempty(d_idx)
      new_laser_lat=new_laser_lat(1:d_idx-1);
      new_laser_lon=new_laser_lon(1:d_idx-1);
      new_rms_l=new_rms_l(1:d_idx-1);
      new_laser_elevation=new_laser_elevation(1:d_idx-1);
      
    end
    
    
    
    
    
  end
  %%  %%
  if debug_flag
    figure(200)
    mapshow(rgb2gray(A),CMAP/1e3);
    % xlim([-300 -50]);
    % ylim([-2350 -2150]);
    xlabel('X (km)');
    ylabel('Y (km)');
    clear gps.x gps.y
    
    [gps.x,gps.y] = projfwd(proj,new_laser_lat,new_laser_lon);
    gps.x = gps.x/1000 ;
    gps.y = gps.y/1000 ;
    figure(200)
    hold on;
    grid
    title('Elevation');
    scatter(gps.x,gps.y,20,new_rms_l)
    
    
    [gps.xrd,gps.yrd] = projfwd(proj,r.lat,r.lon);
    gps.xrd = gps.xrd/1000 ;
    gps.yrd = gps.yrd/1000 ;
    figure(200)
    hold on;
    scatter(gps.xrd,gps.yrd,40,'k','>')
    scatter(gps.xrd(1),gps.yrd(1),50,'c','+')
    title('Truncate for unique icessn line')
  end
  
  %%
  % load(fullfile(base_path_sr,'manjish','peterman','kuband', 'cl27.mat'));
   
  if cross_line==1
    rms_laser=interp1(new_laser_lon,new_rms_l,r.lon,'linear');  %%Interpolate from laser
    %rms_radar=interp1(r.lon,r.rms_height,new_laser_lon);        %Interpolate from radar
   % rms_kuband=interp1(360+LON,RMS_height,r.lon);                %Kuband 
  else
    rms_laser=interp1(new_laser_lat,new_rms_l,r.lat,'linear');  %%Interpolate from laser
    %rms_radar=interp1(r.lat,r.rms_height,new_laser_lat);        %Interpolate from radar
  end
  
  %
  %Laser
  lat=cat(1,lat,r.lat');
  lon=cat(1,lon,r.lon');
  % laser_elev=cat(1,laser_elev,new_laser_elevation);
  laser_rms=cat(1,laser_rms,rms_laser');
  
  %Radar
  radar_rms=cat(1,radar_rms,r.rms_height'*100);
  
  
  if debug_flag
    
    figure(4)
    mapshow(rgb2gray(A),CMAP/1e3);
    % xlim([-300 -50]);
    % ylim([-2350 -2150]);
    xlabel('X (km)');
    ylabel('Y (km)');
    
    clear gps.x gps.y
    [gps.x,gps.y] = projfwd(proj,r.lat,r.lon);
    gps.x = gps.x/1000 ;
    gps.y = gps.y/1000 ;
    figure(4)
    hold on;
    grid
    title('RMS Fit');
    
    scatter(gps.x,gps.y,50,rms_laser,'fill')
    
    clear gps.x gps.y
    [gps.x,gps.y] = projfwd(proj,r.lat,r.lon);
    gps.x = gps.x/1000 ;
    gps.y = gps.y/1000 ;
    figure(4)
    hold on;
    
    scatter(gps.x,gps.y,20,r.rms_height*100,'fill')
    colorbar;
  end
  

if debug_flag
  
  figure(25)
  mapshow(rgb2gray(A),CMAP/1e3);
  % xlim([-300 -50]);
  % ylim([-2350 -2150]);
  xlabel('X (km)');
  ylabel('Y (km)');
  
  clear gps.x gps.y
  [gps.x,gps.y] = projfwd(proj,lat,lon);
  gps.x = gps.x/1000 ;
  gps.y = gps.y/1000 ;
  figure(25)
  hold on;
  grid
  title('Laser rms(interp)');
  
  scatter(gps.x,gps.y,50,laser_rms,'fill')
  colorbar;
  
  clear gps.x gps.y
  [gps.x,gps.y] = projfwd(proj,lat,lon);
  gps.x = gps.x/1000 ;
  gps.y = gps.y/1000 ;
  figure(26)
  mapshow(rgb2gray(A),CMAP/1e3);
  hold on;
  
  scatter(gps.x,gps.y,20,radar_rms,'fill')
  colorbar;
  title('Radar RMS')
  
  figure(27)
  mapshow(rgb2gray(A),CMAP/1e3);
  % xlim([-300 -50]);
  % ylim([-2350 -2150]);
  xlabel('X (km)');
  ylabel('Y (km)');
  
  clear gps.x gps.y
  [gps.x,gps.y] = projfwd(proj,icessn.lat,icessn.lon);
  gps.x = gps.x/1000 ;
  gps.y = gps.y/1000 ;
  figure(27)
  hold on;
  grid
  title('RMS Fit');
  
  scatter(gps.x,gps.y,50,icessn.rms_fit,'fill')
  colorbar;
end

settings.num_int=num_int;
settings.coh_int=coh_int;
settings.sf_bin=sf_bin;
settings.gitinfo=getGitInfo(git_path);


notnan_idx_laser_rms=find(~isnan(laser_rms));
laser_rms=laser_rms(notnan_idx_laser_rms);
radar_rms=radar_rms(notnan_idx_laser_rms);
%rms_kuband=rms_kuband(notnan_idx_laser_rms);
lat=lat(notnan_idx_laser_rms);
lon=lon(notnan_idx_laser_rms);

  
 
  if 1
    figure;plot(laser_rms-radar_rms);title('Difference in RMS height');
    figure;plot(icessn.lat,icessn.rms_fit); title('RMS Fit');
    figure;plot(icessn.lon,icessn.rms_fit); title('RMS Fit');
  %  figure;plot(r.rms_height*100); title('Radar Rms height');
  %  figure;plot(rms_laser); title('Laser RMS ICESSN');
    % figure;plot(elevation_laser);('Laser Surface Elevation');
    %figure;plot(Greenland.Elevation-Greenland.surface_time*c/2);title('Radar Surface Elevation ');
   % figure;plot(lat,rms_laser,'b');hold on; plot(lat,r.rms_height*100,'r');
    %figure;plot(new_laser_lat,rms_radar*100,'r');hold on ;plot(new_laser_lat,new_rms_l,'b');%hold on; plot(new_laser_lon,rms_kuband_radinterp*100,'g');title('Interp from radar')
    if cross_line==1
      figure;plot(lon,laser_rms,'b','DisplayName','Laser data');hold on; plot(lon,radar_rms,'r','DisplayName','Radar Data');legend('show');
%       figure(7);hold on; plot(lon,rms_kuband*100','g');
    else
      figure;plot(lat,laser_rms,'b','DisplayName','Laser data');hold on; plot(lat,radar_rms,'r','DisplayName','Radar Data');legend('show');
    end
  end

  
  keyboard
  
if save_en                %Save final laser radar roughness
  if cross_line==1
    save_dest=fullfile(base_path_sr,'manjish','peterman','new_process', ['crossline',sprintf('%s.mat',num2str(lno))]);
    save(save_dest,'laser_rms','radar_rms','lat','lon','icessn','settings','r');
    disp(sprintf('Saving crossline_%s\n',num2str(lno)));
  else
      save_dest=fullfile(base_path_sr,'manjish','peterman','new_process', ['verticalline',sprintf('%s.mat',num2str(lno))]);
    save(save_dest,'laser_rms','radar_rms','lat','lon','icessn','settings','r');
   % save(['Y:\manjish\peterman\final\verticalline',sprintf('%s',num2str(lno)),'.mat'],'laser_rms','radar_rms','lat','lon','icessn','settings');
   disp(sprintf('Saving verticalline_%s\n',num2str(lno)));
  end
end

if 0                %Save radar roughness
  if cross_line==1
   out_fn=fullfile(base_path_sr,'manjish','peterman','radarnew', ['crossline',sprintf('%s.mat',num2str(lno))]);
   % out_fn=['Y:\manjish\peterman\radarnew\crossline',num2str(lno)];
  else
     out_fn=fullfile(base_path_sr,'manjish','peterman','radarnew', ['verticalline',sprintf('%s.mat',num2str(lno))]);
    %out_fn=['Y:\manjish\peterman\radarnew\verticalline',num2str(lno)];
  end
  
  out_fn_dir=fileparts(out_fn);
  if ~exist(out_fn_dir,'dir')
    mkdir(out_fn_dir);
  end
  save(out_fn,'r');
end
end
end

if 0
  error=laser_rms-radar_rms;
  figure;plot(error);title('Difference in RMS height');
  figure;plot(radar_rms); title('Radar Rms height');
  figure;plot(laser_rms); title('Laser RMS ICESSN');
  % figure;plot(elevation_laser);('Laser Surface Elevation');
  %  figure;plot(Greenland.Elevation-Greenland.surface_time*c/2);title('Radar Surface Elevation ');
  if cross_line==1
    figure;plot(lon,laser_rms,'b','DisplayName','Laser data');hold on; plot(lon,radar_rms,'r','DisplayName','Radar Data');legend('show');
  else
    figure;plot(lat,laser_rms,'b','DisplayName','Laser data');hold on; plot(lat,radar_rms,'r','DisplayName','Radar Data');legend('show');
    
  end
end



  success=true;
return