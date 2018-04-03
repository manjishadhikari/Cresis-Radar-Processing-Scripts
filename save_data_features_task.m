
function [success]=save_data_features_task(param)



%Script: saving_the_required_data_features
%
% Purpose: This script is used to combine the data of all the corresponding frames in a flight path and save the features required for the estimation of relative reflectivity values for that flight path.
%
% Input data / processed on: 2011_Antarctica_TO pulse compressed data.
%
% Output/data product: data structure consisting of following data fields:
%
% GPS_time, Latitude, Longitude, Elevation, ice_bed_time, surface_time, ice_bed_power, ice_surface_power, ice_bed_time, surface_time, segments_length, ice_bed_echo, and abruptness
%
% Output format: Greenland_layerdata_selected_frames_<straight_line/cross_line_number>
%
% Location saved: Y:\santhosh\thesis\get heights frames Greenland\with abruptness (for data corresponding to straight flight paths) and Y:\santhosh\thesis\get heights frames Greenland\cross_lines\with abruptness (for data corresponding to cross flight paths).
%
% See also verification_of_picked_ice_bed_interface, estimation_of_relative_reflectivity_values

if strcmp(param.location{1},'Peterman')
    Peterman=1;
elseif strcmp(param.location{1},'Jacobshavn')
   Peterman=0;
else
  disp('Location not supported')
  success=true;
  return
end

if Peterman
  disp('Processing DataSet for Peterman')
else
  disp('Processing Dataset for Jacobshavn')
end
dbstop error
cross_lines = param.cross_lines_en;

physical_constants;
debug_flag = 0;


if Peterman
  %Peterman
  if cross_lines
    %   Day_seg={'20100324_01','20120516_01','20120516_01','20140512_01','20140505_01','20140505_01','20140505_01','20140505_01','20140505_01'}
    %   frms={[39 40],[13:16],[79:81],[10 11],[12 13],[33:35],[38:40],[55 56],[59 60]}
  
    Day={'20100324_01','20100324_02','20100324_03','20100324_04','20110429_01','20110429_02','20110507_01','20110507_02','20120330_01','20120516_01','20140512_01','20140505_01'};
    Day_seg=repelem(Day,[3,1,1,1,5,1,7,1,1,2,1,5]);
    frms={[36,37],[39,40],[42],[1 2],[1 2],[1 2],[9:12],[13:16],[17:20],[21:24],[25:28],[10:11],[10:14],[15:18],[19:22],[23:26],[27:30],[31:34],[35:37],[1:4],[5 6 7],[13:16],[79:81],[10 11],[12 13],[33 34 35],[38:40],[55 56],[59 60]};
  
  else
    %      for straight lines
    Day={'20100324_01','20110429_01','20110429_02','20110507_02','20130420_02','20140512_01','20140505_01'};
    Day_seg=repelem(Day,[7,2,4,2,3,2,1]);
    frms={[11 12],[14 15],[17 18],[20 21],[23 24] [30 31],[33 34],[30:32],[33 34],[18:21],[12:15],[2:5],[6:9],[6:8],[17:20],[3 4],[9],[11],[12 13],[17 18],[15 16]};
  %fprintf('%d %d', fsu(m), fsl(m))
  end
  
  
else
  %Jacobshavn
  
  if cross_lines
    %   Day_seg={'20100517_01','20110331_04','20110408_18',}
    %   frms={[6,7],[10,11],[13,14],[16,17],[19],[23,24],[25,26],[29,30],[1],[2],}
    
    Day={'20120502_01','20130410_01','20130415_01','20140419_02','20140419_03'}; 
    Day_seg=repelem(Day,[9,7,1,3,5]);
    frms={[5,6],[8,9],[10,11],[14,15],[16,17],[20,21],[22,23],[28,29],[30],[4,5],[8:10],[11,12],[15],[16,17],[28],[33],[4,5],[4,5],[7,8],[9,10],[2],[3,4],[14],[15,16],[19]};
    
  else
    %      for vertical lines
    
    %Day_seg={'20100514_01','20100514_02'};
    %frms={[5,6],[1,2],[3,4]};
    
    Day={'20110406_01','20110422_01','20120421_01','20130410_01','20130404_02','20140414_02','20140409_01','20140409_02'};
    Day_seg=repelem(Day,[6,1,9,1,8,2,6,4]);
    frms={[4:6],[7:9],[10:12],[13:15],[16:18],[19:21],[2:4],[7:9],[10:12],[13:15],[16:18],[19:21],[22:24],[25:27],[28],[52],[50:52],[1,2],[3:5],[6,7],[8:10],[11:13],[14],[25],[38:40],[5:6],[12],[5:7],[8:10],[11,12],[13;15],[16,17],[18:20],[15],[19:22],[28:29],[34:36]};
  end
end

if param.proc_line>length(Day_seg)
  disp(sprintf('Max line number is: %d',length(Day_seg)))
  success=true;
  return;
end


%%  reading the data

for k =param.proc_line
 
  
  
 %% Peterman
 if Peterman
 if cross_lines
    if k<7
      param1 = read_param_xls(ct_filename_param_v2(param,'rds_param_2010_Greenland_DC8.xls'),Day_seg{k});
      param=mergestruct(param,param1);
      gps_fn = ct_filename_support(param,'','gps',1);
      data_dir = ct_filename_out(param,'','manjish/CSARP_Data');
      layer_dir = ct_filename_out(param,'','old/CSARP_layerData');
    elseif (6<k & k<21)
      param1 = read_param_xls(ct_filename_param_v2(param,'rds_param_2011_Greenland_P3.xls'),Day_seg{k});
      param=mergestruct(param,param1);
      gps_fn = ct_filename_support(param,'','gps',1);
      data_dir = ct_filename_out(param,'','CSARP_manjish');
      layer_dir = ct_filename_out(param,'','CSARP_layerData');
      
    elseif (20<k & k<24)
      param1 = read_param_xls(ct_filename_param_v2(param,'rds_param_2012_Greenland_P3.xls'),Day_seg{k});
      param=mergestruct(param,param1);
      gps_fn = ct_filename_support(param,'','gps',1);
      data_dir = ct_filename_out(param,'','CSARP_CSARP_manjish');
      layer_dir = ct_filename_out(param,'','CSARP_post/CSARP_layerData');
    end
    
  else
    
    if k<8
      param1 = read_param_xls(ct_filename_param_v2(param,'rds_param_2010_Greenland_DC8.xls'),Day_seg{k});
      param=mergestruct(param,param1);
      gps_fn = ct_filename_support(param,'','gps',1);
      data_dir = ct_filename_out(param,'','manjish/CSARP_Data');
      layer_dir = ct_filename_out(param,'','old/CSARP_layerData');
    elseif (7<k & k<16)
      param1 = read_param_xls(ct_filename_param_v2(param,'rds_param_2011_Greenland_P3.xls'),Day_seg{k});
      param=mergestruct(param,param1);
      gps_fn = ct_filename_support(param,'','gps',1);
      data_dir = ct_filename_out(param,'','CSARP_manjish');
      layer_dir = ct_filename_out(param,'','CSARP_layerData');
      
    elseif (15<k & k<19)
      param1 = read_param_xls(ct_filename_param_v2(param,'rds_param_2013_Greenland_P3.xls'),Day_seg{k});
      param=mergestruct(param,param1);
      gps_fn = ct_filename_support(param,'','gps',1);
      data_dir = ct_filename_out(param,'','manjish/CSARP_Data');
      layer_dir = ct_filename_out(param,'','CSARP_post/CSARP_layerData');
    else
      param1 = read_param_xls(ct_filename_param_v2(param,'rds_param_2014_Greenland_P3.xls'),Day_seg{k});
      param=mergestruct(param,param1);
      gps_fn = ct_filename_support(param,'','gps',1);
      data_dir = ct_filename_out(param,'','CSARP_manjish');
      layer_dir = ct_filename_out(param,'','CSARP_post/CSARP_layerData');
    end
end
 else
%% Jacobshavn

 if cross_lines
%     if k<9
%       param1 = read_param_xls(ct_filename_param_v2(param,'rds_param_2010_Greenland_P3.xls'),Day_seg{k});
%       param=mergestruct(param,param1);
%       gps_fn = ct_filename_support(param,'','gps',1);
%       data_dir = ct_filename_out(param,'','CSARP_manjish');
%       layer_dir = ct_filename_out(param,'','CSARP_layerData');
%     elseif 8<k & k<11
%       param1 = read_param_xls(ct_filename_param_v2(param,'rds_param_2011_Greenland_P3.xls'),Day_seg{k});
%       param=mergestruct(param,param1);
%       gps_fn = ct_filename_support(param,'','gps',1);
%       data_dir = ct_filename_out(param,'','CSARP_manjish');
%       layer_dir = ct_filename_out(param,'','CSARP_layerData');
      
    if k<10
      param1 = read_param_xls(ct_filename_param_v2(param,'rds_param_2012_Greenland_P3.xls'),Day_seg{k});
      param=mergestruct(param,param1);
      gps_fn = ct_filename_support(param,'','gps',1);
      data_dir = ct_filename_out(param,'','CSARP_CSARP_manjish');
      layer_dir = ct_filename_out(param,'','CSARP_post/CSARP_layerData');
    elseif (9<k & k<18)
      param1 = read_param_xls(ct_filename_param_v2(param,'rds_param_2013_Greenland_P3.xls'),Day_seg{k});
      param=mergestruct(param,param1);
      gps_fn = ct_filename_support(param,'','gps',1);
      data_dir = ct_filename_out(param,'','manjish/CSARP_Data');
      layer_dir = ct_filename_out(param,'','CSARP_post/CSARP_layerData');
    else 
      param1 = read_param_xls(ct_filename_param_v2(param,'rds_param_2014_Greenland_P3.xls'),Day_seg{k});
      param=mergestruct(param,param1);
      gps_fn = ct_filename_support(param,'','gps',1);
      data_dir = ct_filename_out(param,'','CSARP_manjish');
      layer_dir = ct_filename_out(param,'','CSARP_post/CSARP_layerData');
    end
    
  else
    
 
    if k<8
      param1 = read_param_xls(ct_filename_param_v2(param,'rds_param_2011_Greenland_P3.xls'),Day_seg{k});
      param=mergestruct(param,param1);
      gps_fn = ct_filename_support(param,'','gps',1);
      data_dir = ct_filename_out(param,'','CSARP_manjish');
      layer_dir = ct_filename_out(param,'','CSARP_layerData');
    elseif (7<k & k<17)
      param1 = read_param_xls(ct_filename_param_v2(param,'rds_param_2012_Greenland_P3.xls'),Day_seg{k});
      param=mergestruct(param,param1);
      gps_fn = ct_filename_support(param,'','gps',1);
      data_dir = ct_filename_out(param,'','CSARP_CSARP_manjish');
      layer_dir = ct_filename_out(param,'','CSARP_post/CSARP_layerData');
    elseif (16<k & k<26)
      param1 = read_param_xls(ct_filename_param_v2(param,'rds_param_2013_Greenland_P3.xls'),Day_seg{k});
      param=mergestruct(param,param1);
      gps_fn = ct_filename_support(param,'','gps',1);
      data_dir = ct_filename_out(param,'','manjish/CSARP_Data');
      layer_dir = ct_filename_out(param,'','CSARP_post/CSARP_layerData');
      else
      param1 = read_param_xls(ct_filename_param_v2(param,'rds_param_2014_Greenland_P3.xls'),Day_seg{k});
      param=mergestruct(param,param1);
      gps_fn = ct_filename_support(param,'','gps',1);
      data_dir = ct_filename_out(param,'','CSARP_manjish');
      layer_dir = ct_filename_out(param,'','CSARP_post/CSARP_layerData');
    end
end
  
 end
 %% 
  
  Greenland.GPS_time = [];
  Greenland.Latitude = [];
  Greenland.Longitude = [];
  Greenland.Elevation = [];
  Greenland.Roll=[];
  Greenland.ice_bed_time = [];
  Greenland.surface_time = [];
  Greenland.ice_bed_power = [];
  Greenland.ice_surface_power = [];
  Greenland.ice_bed_time = [];
  Greenland.surface_time = [];
  Greenland.segments_length = [];
  Greenland.index.Latitude_mean=[];
  Greenland.index.Longitude_mean=[];
  Greenland.index.GPS_time_ave=[];
  Greenland.index.coherence=[];
  Greenland.index.abruptness=[];
  Greenland.index.Padj=[];
  Greenland.index.Padj_Na11=[];
  Greenland.index.bt.val=[];
  Greenland.index.bt.idx=[];
  Greenland.index.bt.waveform=[];
  Greenland.index.bt.inc_wf_ave=[];
  
  
  
  
  
  
  for f = 1:length(frms{k})
    frm = frms{k}(f);
    layer_fn = fullfile(layer_dir,sprintf('Data_%s_%03d.mat', param1.day_seg, frm));
    fprintf('Loading data %s\n', layer_fn);
    data_fn = fullfile(data_dir,sprintf('Data_%s_%03d.mat', param1.day_seg, frm));
   
    % Load the file
    data = load(data_fn);
    tmp = load(layer_fn);
    
    surface_twtt=interp1(tmp.GPS_time,tmp.layerData{1}.value{2}.data , data.GPS_time,'linear','extrap');
    surface_twtt(isinf(surface_twtt))=nan;
    
%     extrap_idx=find(data.GPS_time>=tmp.GPS_time(end),1,'first');
%     %extrap_idx=find(isnan(surface_twtt),1,'first');
%     surface_twtt(extrap_idx-50:end)=nanmean(surface_twtt(extrap_idx-50:extrap_idx-1));
     bottom_twtt = interp1(tmp.GPS_time,tmp.layerData{2}.value{2}.data , data.GPS_time,'linear','extrap');
      bottom_twtt(isinf(bottom_twtt))=nan;
%    bottom_twtt(extrap_idx-50:end)=nanmean(bottom_twtt(extrap_idx-50:extrap_idx-1));
     
     if debug_flag
       figure(1);imagesc([],data.Time*1e6,lp(data.Data));
       figure(1);hold on; plot(surface_twtt*1e6);
       figure(1);hold on; plot(bottom_twtt*1e6); title('Before')
     end
  
    %         Elevation = interp1(tmp.GPS_time,tmp.Elevation,data.GPS_time);
    %         surface_twtt = data.Surface + ((data.Elevation-Elevation)/(c/2));
    %         filter_length = 100;
    %         surface_twtt = sgolayfilt(surface_twtt, 2,filter_length+1, hanning(filter_length+1));
    dt= data.Time(2)-data.Time(1);
    dh=dt*c/sqrt(3.14);
    %index_sf = round((surface_twtt-data.Time(1))/dt);
    index_sf=round(interp1(data.Time,1:size(data.Data,1),surface_twtt,'linear','extrap'));
    ice_surface_power  = zeros(1,length(data.Surface));
    
    for i = 1:length(surface_twtt)
      if isnan(surface_twtt(i)) || isinf(surface_twtt(i))
        ice_surface_power(i)= nan;
        continue
      else
        [surface_power, idx] = max(lp(data.Data(index_sf(i)-0:index_sf(i)+0,i)));
      %  [surface_power idx] = max(sqrt(data.Data(index(i)-5:index(i)+5,i).*conj(data.Data(index(i)-5:index(i)+5,i))));
        surface_index = idx + index_sf(i)+0-1;
        if surface_power  == 0
          ice_surface_power(i)= nan;
          surface_twtt(i) = nan;
        else
          ice_surface_power(i) = data.Data(surface_index,i);
          surface_twtt(i) =  interp1([1:length(data.Time)],data.Time,surface_index);
        end
      end
    end
    if any(ice_surface_power ==0 )
      keyboard
      warning('Ice Surface power is zero')
    end
    
    %Bottom Tracking 
   
           % index_bt = round((bottom_twtt-data.Time(1))/dt);
             index_bt=round(interp1(data.Time,1:size(data.Data,1),bottom_twtt,'linear','extrap'));
            ice_bed_power  = zeros(1,length(data.Surface));

            for i = 1:length(bottom_twtt)
             
                if isnan(bottom_twtt(i)) || isinf(bottom_twtt(i))
                    ice_bed_power(i)= nan;
                    continue
                else
                 
                    [bed_power, idx] = max(lp(data.Data(index_bt(i):index_bt(i),i)));
                    bed_index = idx + index_bt(i)-1;
                
                    
                    idx1=bed_index+200;
                    idx2=bed_index+500;
                    if idx1>size(data.Data,1) | idx2>size(data.Data,1)
                      idx2=size(data.Data,1);
                      idx1=idx2-300;
                    end
                    N = mean((data.Data(idx1:idx2,i)));  %Noise floor
                 
                    % Put condition for: if bottom is equal to surface,
                    % reject it as bottom as surface power creates false
                    % bottom analysis results %setting min ice depth of
                    % 500m as bottom criteria
                    depthbt=(index_bt(i)-index_sf(i))*dh;
                    if depthbt<500
                        ice_bed_power(i) = nan;
                        bottom_twtt(i) = nan;
                        continue ;
                    end
                    
                    SNR=bed_power-lp(N);
                    if SNR > 3
                        ice_bed_power(i) = data.Data(bed_index,i);
                        bottom_twtt(i) =  interp1([1:length(data.Time)],data.Time,bed_index);
                        %
                    else
                        ice_bed_power(i) = nan;
                        bottom_twtt(i) = nan;
                        continue ;
                    end
                end
            end
      

    
    if any(ice_bed_power ==0 )
      keyboard
      warning('Ice bed power is zero')
    end
    
    if debug_flag
      figure(2);imagesc([],data.Time*1e6,lp(data.Data));
      figure(2);hold on; plot(surface_twtt*1e6);
      figure(2);hold on; plot(bottom_twtt*1e6);
      title(sprintf('Data-%s-%03d.mat', param1.day_seg, frm))
    end
    
  
    %Save figure as jpg for checking
   
   if param.save_fig_only
     warning('Save figures enabled.. Set it to 0 if not necessary')
      if debug_flag==0
        figure(2);close 
        figure(2);imagesc([],data.Time*1e6,lp(data.Data));
        figure(2);hold on; plot(surface_twtt*1e6);
        figure(2);hold on; plot(bottom_twtt*1e6);
        title(sprintf('Data-%s-%03d.mat', param1.day_seg, frm))
      end
      
      if cross_lines
       if Peterman
        save_path=['/cresis/snfs1/scratch/manjish/peterman/images/',sprintf('crossline%d',k),'/',sprintf('Data_%s_%03d', param1.day_seg, frm)];
       else
         save_path=['/cresis/snfs1/scratch/manjish/jacobshavn/images/',sprintf('crossline%d',k),'/',sprintf('Data_%s_%03d', param1.day_seg, frm)];
       end
        [save_dir] =fileparts(save_path);
        if ~exist(save_dir,'dir')
          
          mkdir(save_dir);
        end
        saveas(figure(2),save_path,'jpg')
       
      else
        if 0
        if Peterman
            save_path=['/cresis/snfs1/scratch/manjish/peterman/images/',sprintf('verticalline%d',k),'/',sprintf('Data_%s_%03d', param1.day_seg, frm)];
        else
           save_path=['/cresis/snfs1/scratch/manjish/jacobshavn/images/',sprintf('verticalline%d',k),'/',sprintf('Data_%s_%03d', param1.day_seg, frm)];
    
        end
        
        [save_dir] =fileparts(save_path);
        if ~exist(save_dir,'dir')
          
          mkdir(save_dir);
        end
        saveas(figure(2),save_path,'jpg')
        end
      end
    end
     
    if ~param.save_fig_only         %Disable to only save figures and not perform save data features
  
      %Coherence Index and Abruptive Index Calculation
    % fc= 195e6;  %radar center frequemcy
    fc=(param1.radar.wfs(1).f0+param1.radar.wfs(1).f1)/2;
    % fs=1.1111e8;  %Sampling frequency
    fs=param1.radar.fs;
    p=4.99;  %Radar Half pulse width in air
    
    %%
    physical_constants;
    
    % L=load(layer_fn);
    
    if debug_flag
      figure(1),imagesc(data.GPS_time,data.Time,10*log10(abs(data.Data)));
      hold on;plot(data.GPS_time,surface_twtt,'--');plot(data.GPS_time,bottom_twtt,'--');
    end
    
    
    %% Coherene Index Calculation with original Datasf(sf==Inf)=NaN;
    
    
    idx=find(isinf(bottom_twtt));
    bottom_twtt(idx)=nan;
    idx2=find(isinf(surface_twtt));
    surface_twtt(idx2)=nan;
    
    MeanDepth=(nanmean(bottom_twtt)-nanmean(surface_twtt))*c/(2*sqrt(er_ice));  %Original Mean Depth
    
    %Mean Ice Depth
    
    if MeanDepth~=0 & MeanDepth>150
      Nx_int_dist=2*sqrt(MeanDepth*p/sqrt(er_ice));                    %Incoherent Integration distance
      distance=geodetic_to_along_track(data.Latitude,data.Longitude,data.Elevation);
      Nx_int = floor(Nx_int_dist/(distance(10)-distance(9)));                %No of lines integrated
      Nx0 = size(data.Data,2);
      Nx = floor(Nx0/Nx_int);                                         %Total lines after integration
      Nx_mod = mod(Nx0,Nx_int);
      if Nx_mod>= Nx_int/2;
        Nx = Nx + 1;
      end
      square_int = zeros(size(data.Data,1),Nx); % incoherent integration, take square of abs first, then sum (phase info lost)
      int_square = zeros(size(data.Data,1),Nx); % coherent integration, sum complex data first, then take square of abs (phase info remains)
      coh_index_ogdata = zeros(1,Nx);
      Abruptiveindex=zeros(1,Nx);
      risingedge=zeros(1,Nx);
      fallingedge=zeros(1,Nx);
      Padj=zeros(1,Nx);
      
      if 1
        for rline =1:Nx
          idx1 = (rline-1)*Nx_int + 1;
          idx2 = rline*Nx_int;
          if Nx0 - idx2 > 0 && Nx0 - idx2 < Nx_int/2;
            idx2 = Nx0;
          else
            idx2 = min(idx2,Nx0);
          end
          square_int(:,rline) = mean(abs(data.Data(:,idx1:idx2)).^2,2);
          square_int_dB = 10*log10(square_int(:,rline));
          int_square(:,rline) = (abs(mean(data.Data(:,idx1:idx2),2))).^2;
          
          meanbt=nansum(bottom_twtt(isfinite(bottom_twtt(idx1:idx2))));
          count=nansum(isfinite(bottom_twtt(idx1:idx2)));
          if meanbt==0 || count==0
            continue;
          end
          meanbt=meanbt/count;
          bt_idx_m = find(data.Time>meanbt,1,'first');
          
          if isempty(bt_idx_m)
            bt_idx_m=round(size(data.Data,1)/2);
          end
          b1=bt_idx_m-5;
          b2=bt_idx_m+5;
          if b2>size(data.Data,1)
            b2=size(data.Data,1);
          end
          
          [bt_val,bt_idx] = max(square_int(b1:b2,rline));
          bt_idx = bt_idx_m-5+bt_idx-1;                 %Peak Index
          bt_pwr = 10*log10(bt_val);                     %Peak Value
          square_int_dB = 10*log10(square_int(:,rline));
          noise_bin1=bt_idx+70;
          noise_bin2=bt_idx+100;
          if noise_bin1>size(data.Data,1)
            noise_bin1=size(data.Data,1)-30;
            noise_bin2=size(data.Data,1);
          end
          
          if noise_bin2>size(data.Data,1)
            noise_bin1=size(data.Data,1)-30;
            noise_bin2=size(data.Data,1);
          end
          noise =mean(square_int_dB(noise_bin1:noise_bin2));
          risingedge(rline) = bt_idx-1;
          SNR = bt_pwr-noise;
          
          % find risingedge(rline) and fallingedge(rline) for integration in range
          while square_int_dB(risingedge(rline))-noise > 0.05*SNR && risingedge(rline)>2
            risingedge(rline) = risingedge(rline) - 1;
          end
          if risingedge(rline)==0
            risingedge(rline)=find(square_int_dB==min(square_int_dB(b1:bt_idx)),1,'first');
          end
          
          fallingedge(rline) = bt_idx+1;
          if fallingedge(rline)>size(data.Data,1)
            fallingedge(rline)=size(data.Data,1) ;
          end
          while square_int_dB(fallingedge(rline))-noise > 0.05*SNR &&  fallingedge(rline)<size(data.Data,1)      %Depth Bins risingedge(rline) and fallingedge(rline)
            fallingedge(rline) = fallingedge(rline) + 1;
          end
          if fallingedge(rline)-bt_idx>50        %Max Allowed bin
            fallingedge(rline)=bt_idx+50;
          end
          if fallingedge(rline)-bt_idx<3           %Guard bin
            fallingedge(rline)=bt_idx+3;
          end
          D1=risingedge(rline);
          D2=fallingedge(rline);
          
          Imeanx=sum((square_int(risingedge(rline):fallingedge(rline),rline)));
          Abruptiveindex(rline)=bt_val/Imeanx;
          coh_index_ogdata(rline) = sum(int_square(risingedge(rline):fallingedge(rline),rline))/sum(square_int(risingedge(rline):fallingedge(rline),rline));
          
          
        end
        distancenx=max(distance).*Nx_int*0.5/(size(distance,2)):max(distance).*Nx_int/size(distance,2):max(distance).*Nx_int/size(distance,2)*Nx;
        if debug_flag
          figure(4);plot(distancenx,coh_index_ogdata,'Displayname','Original Data');
          title('Original Data Coherence Index');   %Coherence Index without Correction
          figure(40);plot(distancenx,Abruptiveindex);
        end
      end
      
      %% Correction for data.Elevation
      
      sf=surface_twtt; %test
      sf_new = zeros(size(sf));
      bt_new = zeros(size(bottom_twtt));
      data.Elevation_new = zeros(size(data.Elevation));
      data.sf_elev_new = zeros(size(data.Elevation));
      data.bt_elev_new = zeros(size(data.Elevation));
      
      % 1)Remove data before zero time
      negative_bins = data.Time < 0;
      data.Time_new = data.Time(~negative_bins);
      data.Data_new = data.Data(~negative_bins,:);
      
      % 2)Create data.Elevation axis to interpolate to
      [max_elev,max_elev_idx] = max(data.Elevation);
      min_elev = min(data.Elevation - sf*c/2 - (data.Time_new(end)-sf)*c/2/sqrt(er_ice));
      dt = data.Time(2)-data.Time(1);
      dr = dt * c/2 / sqrt(er_ice);
      dt_air = dr/(c/2);
      elev_axis = max_elev:-dr:min_elev;
      new_time = zeros(length(elev_axis),length(data.Elevation));
      
      % 3)Zero pad data to create space for interpolated data
      zero_pad_len = length(elev_axis) - length(data.Time_new);
      data.Data_new = cat(1,data.Data_new,zeros(zero_pad_len,size(data.Data_new,2)));
      
      % 4)Determine the corrections to apply to data.Elevation and layers
      dRange = max_elev - data.Elevation;
      dBins = round(dRange / (c/2) / dt);
      dtime = dRange/(c/2);
      
      for rline = 1:size(data.Data_new,2)
        % Determine data.Elevation bins before surface
        sf_elev = data.Elevation(rline) - sf(rline) * c/2;
        time0 = -(max_elev - data.Elevation(rline))/(c/2);
        last_air_idx = find(elev_axis > sf_elev,1,'last');
        new_time_tmp = (time0 + dt_air*(0:last_air_idx-1)).';
        if last_air_idx < length(elev_axis)
          % Determine data.Elevation bins after surface
          dt_ice = dr/(c/2/sqrt(er_ice));
          first_ice_idx = last_air_idx + 1;
          time0 = sf(rline) + (sf_elev - elev_axis(first_ice_idx))/(c/2/sqrt(er_ice));
          new_time(:,rline) = cat(1,new_time_tmp, (time0 + dt_ice*(0:length(elev_axis)-length(new_time_tmp)-1)).');
        end
        data.Data_new(:,rline) = interp1(data.Time_new, data.Data_new(1:length(data.Time_new),rline), new_time(:,rline), 'linear',0);
        data.Elevation_new(rline) = data.Elevation(rline) + dRange(rline);
        sf_new(rline) = sf(rline) + dtime(rline);
        bt_new(rline) = bottom_twtt(rline) + dtime(rline);
        data.sf_elev_new(rline) = data.Elevation_new(rline) - sf_new(rline)*c/2;
        data.bt_elev_new(rline) = data.sf_elev_new(rline) - (bt_new(rline)-sf_new(rline))*c/2/sqrt(er_ice);
        
        
        
      end
      
      if debug_flag
        fh = figure(2);
        figure(fh);imagesc([],elev_axis,10*log10(abs(data.Data_new).^2));title('data.Elevation Correction Data')
        ax = gca;
        ax.YDir = 'normal';
        hold on;plot(data.sf_elev_new,'--');plot(data.bt_elev_new,'--');
        bt_slope = diff(data.bt_elev_new)./diff(distance);
      end
      
      
      
      %% Coherene Index Calculation with data.Elevation Correction
      if 0
        square_int = zeros(size( data.Data_new,1),Nx); % incoherent integration, take square of abs first, then sum
        int_square = zeros(size(data.Data_new,1),Nx); % coherent integration, sum complex data first, then take square of abs
        coh_index_elevcorr = zeros(1,Nx);
        for rline =1:Nx
          idx1 = (rline-1)*Nx_int + 1;
          idx2 = rline*Nx_int;
          if Nx0 - idx2 > 0 & Nx0 - idx2 < Nx_int/2;
            idx2 = Nx0;
          else
            idx2 = min(idx2,Nx0);
          end
          square_int(:,rline) = mean(abs(data.Data_new(:,idx1:idx2)).^2,2);
          int_square(:,rline) = abs(mean(data.Data_new(:,idx1:idx2),2)).^2;
          square_int_dB=lp(square_int(:,rline));
          
          meanbt=nansum(data.bt_elev_new(idx1:idx2));
          count=nansum(isfinite(data.bt_elev_new(idx1:idx2)));
          meanbt=meanbt/count;
          if meanbt==0 || count==0 ||isnan(meanbt)
            continue;
          end
          
          bt_idx_m = find(elev_axis<=meanbt,1,'first');
          
          if isempty(bt_idx_m)
            bt_idx_m=round(size(data.Data_new,1)/2);      %looking for ice bottom
          end
          b1=bt_idx_m-5;
          b2=bt_idx_m+5;
          if b2>size(data.Data_new,1)
            b2=size(data.Data_new,1);
          end
          
          [bt_val,bt_idx]=max(square_int(b1:b2,rline));
          bt_idx = bt_idx_m-5+bt_idx-1;
          bt_pwr = 10*log10(bt_val);
          noise_bin1=bt_idx+70;
          noise_bin2=bt_idx+100;
          if noise_bin1>size(data.Data_new,1)
            noise_bin1=size(data.Data_new,1)-30;
            noise_bin2=size(data.Data_new,1);
          end
          if noise_bin2>size(data.Data_new,1)
            noise_bin2=size(data.Data_new,1);
          end
          noise = mean(square_int_dB(noise_bin1:noise_bin2));
          
          risingedge(rline) = bt_idx-1;
          SNR = bt_pwr-noise;
          
          % find risingedge(rline) and fallingedge(rline) for integration in range
          while square_int_dB(risingedge(rline))-noise > 0.05*SNR & risingedge(rline)>2
            risingedge(rline) = risingedge(rline) - 1;
          end
          if risingedge(rline)==0|square_int_dB(risingedge(rline))==-Inf|square_int_dB(risingedge(rline))==Inf
            risingedge(rline)=find(square_int_dB==min(square_int_dB(b1:bt_idx)),1,'first');
          end
          fallingedge(rline) = bt_idx+1;
          if fallingedge(rline)>size(data.Data_new,1)
            fallingedge(rline)=size(data.Data_new,1) ;
          end
          while square_int_dB(fallingedge(rline))-noise > 0.05*SNR && fallingedge(rline)<size(data.Data_new,1)      %Depth Bins risingedge(rline) and fallingedge(rline)
            fallingedge(rline) = fallingedge(rline) + 1;
          end
          Imeanx=sum((square_int(risingedge(rline):fallingedge(rline),rline)));
          Abruptiveindex(rline)=bt_val/Imeanx;
          coh_index_elevcorr(rline) = sum(int_square(risingedge(rline):fallingedge(rline),rline))/sum(square_int(risingedge(rline):fallingedge(rline),rline));
        end
        
        if debug_flag
          hold on; figure(4); hold on; plot(distancenx,coh_index_elevcorr,'r','Displayname','data.Elevation Corrected');
          figure(41); hold on; plot(distancenx,Abruptiveindex,'r','Displayname','data.Elevation Corrected'); grid;
        end
      end
      
      %% Ice Slope Correction
      
      Nx0 = size(data.Data_new,2);
      Nx = floor(Nx0/Nx_int);
      Nx_mod = mod(Nx0,Nx_int);
      if Nx_mod>= Nx_int/2;
        Nx = Nx + 1;
      end
      slopeerror=zeros(1,size(data.Data_new,1));
      slopeval=zeros(1,size(data.Data_new,1));
      angle=zeros(1,Nx);
      
      for rline =1:Nx
        
        idx1 = (rline-1)*Nx_int + 1;
        idx2 = rline*Nx_int;
        if Nx0 - idx2 > 0 && Nx0 - idx2 < Nx_int/2;
          idx2 = Nx0;
        else
          idx2 = min(idx2,Nx0);
        end
        
        if isinf(mean(data.bt_elev_new(idx1:idx2)))  % If No Ice Bottom, skip
          continue;
        end
        if isnan(mean(data.bt_elev_new(idx1:idx2)))  % If No Ice Bottom, skip
          continue;
        end
        
        p = polyfit(distance(idx1:idx2),data.bt_elev_new(idx1:idx2),1);    % Polygonal fitting
        slopeval(idx1:idx2) = polyval(p,distance(idx1:idx2));           % Ploygonal values after fitting
        base = distance(idx2)-distance(idx1);
        perpendicular = slopeval(idx2)-slopeval(idx1);
        hypotenuse = sqrt(base^2+perpendicular^2);
        angle(rline)=asin(perpendicular/hypotenuse)*180/pi;
        slopeerror = slopeval(idx1:idx2)-slopeval(idx1);         % Error betn original and fitting line
        dtime = 2*slopeerror/c/sqrt(er_ice);
        if debug_flag
          if rline==1
            figure(8);plot([idx1:idx2],data.bt_elev_new(idx1:idx2));
            hold on;plot([idx1:idx2],slopeval(idx1:idx2),'r--')
            figure(9);plot(10*log10(abs(data.Data_new(:,idx2)).^2))
          end
        end
        for idx = idx1:idx2
          data.Data_new(:,idx) = interp1(elev_axis, data.Data_new(:,idx), elev_axis + slopeerror(idx-idx1 +1), 'linear',0);
          data.Elevation_new(idx) = data.Elevation_new(idx) - slopeerror(idx-idx1 +1);
          data.sf_elev_new(idx) = data.sf_elev_new(idx) - slopeerror(idx-idx1 +1);
          data.bt_elev_new(idx) = data.bt_elev_new(idx) - slopeerror(idx-idx1 +1);
          sf_new(idx) = sf_new(idx) + dtime(idx-idx1 +1);
          bt_new(idx) = bt_new(idx) + dtime(idx-idx1 +1);
        end
        if debug_flag
          if rline==1
            figure(8);hold on;plot([idx1:idx2],data.bt_elev_new(idx1:idx2),'g--');
            figure(9);hold on;plot(10*log10(abs(data.Data_new(:,idx2)).^2),'g--');
          end
        end
        end
        
      
      if debug_flag
        fh = figure(3);
        %figure(fh);imagesc([1:size(data.Data_new,2)],elev_axis,10*log10(abs(data.Data_new).^2));title('Slope Correction Data')
        figure(fh);imagesc([],elev_axis,10*log10(abs(data.Data_new).^2));title('Slope Correction Data')
        ax = gca;
        ax.YDir = 'normal';
        hold on;plot(data.sf_elev_new,'--');plot(data.bt_elev_new,'--');
      end
       figure(2);imagesc([],data.Time*1e6,lp(data.Data));
      figure(2);hold on; plot(surface_twtt*1e6);
      figure(2);hold on; plot(bottom_twtt*1e6);
      title(sprintf('Data-%s-%03d.mat', param1.day_seg, frm))
      
      %
      % Truncate data around ice bottom within bt.range_bins
      bt_range_bins =[-50:100];
      bt.val = NaN*ones(1,size(data.Data_new,2));
      bt.idx = NaN*ones(1,size(data.Data_new,2));
      bt.waveform = NaN*ones(length(bt_range_bins),size(data.Data_new,2));
      bt.inc_wf_ave=NaN*ones(size(bt.waveform,1),Nx);
      for rline = 1:size(data.Data_new,2)
        if ~isnan(data.bt_elev_new(rline)) & ~isinf(data.bt_elev_new(rline))
          % bt.idx(rline) = find(elev_axis<=data.bt_elev_new(rline),1,'first');
          bt.idx(rline) = round(interp1(elev_axis,[1:length(elev_axis)],data.bt_elev_new(rline)));
          bt.val(rline) = data.Data_new(bt.idx(rline),rline);
          first_idx = bt.idx(rline) + bt_range_bins(1);
          last_idx = bt.idx(rline) + bt_range_bins(end);
          if first_idx < 1 | last_idx>size(data.Data_new,1)
            bt.idx(rline) = NaN;
            bt.val(rline) = NaN;
            continue
          end
          lower=bt.idx(rline)+bt_range_bins(1);
          upper=bt.idx(rline)+bt_range_bins(end);
          if upper >size(data.Data_new,1)
            upper=size(data.Data_new,1);
          end
          bt.waveform(1:51+upper-bt.idx(rline),rline) = data.Data_new(lower:upper, rline);
        else
          continue
        end
      end
      
      
      if debug_flag
        figure(31);imagesc(lp(data.Data_new));
        hold on;plot(bt.idx,'--');
        figure(33);plot(lp(bt.waveform));
        figure(32);imagesc(lp(bt.waveform));
      end
      
      
      
      %}
      %% Coherene Index Calculation with Slope error Correction
      
      
      %
      
      Nx0 = size(data.Data,2);
      Nx = floor(Nx0/Nx_int);
      Nx_mod = mod(Nx0,Nx_int);
      if Nx_mod>= Nx_int/2;
        Nx = Nx + 1;
      end
      square_int = zeros(size( data.Data_new,1),Nx); % incoherent integration, take square of abs first, then sum
      square_int_dB=zeros(size( data.Data_new,1),Nx);
      int_square = zeros(size(data.Data_new,1),Nx); % coherent integration, sum complex data first, then take square of abs
      coh_index_slopecorr =NaN*ones(1,Nx);
      Abruptiveindex=NaN*ones(1,Nx);
      Padj=NaN*ones(1,Nx);
      Padj2_Na11=NaN*ones(1,Nx);
      Latitude_mean=zeros(1,Nx);
      Longitude_mean=zeros(1,Nx);
      GPS_time_ave=zeros(1,Nx);
      
      risingedge=NaN*ones(1,Nx);
      fallingedge=NaN*ones(1,Nx);
      peakindex=NaN*ones(1,Nx);
      depth=NaN*ones(1,Nx);
      Gmtrc_loss=NaN*ones(1,Nx);
      Power=NaN*ones(1,Nx);
      Dpth=NaN*ones(1,Nx);
      for rline =1:Nx
        idx1 = (rline-1)*Nx_int + 1;
        idx2 = rline*Nx_int;
        if Nx0 - idx2 > 0 & Nx0 - idx2 < Nx_int/2;
          idx2 = Nx0;
        else
          idx2 = min(idx2,Nx0);
        end
        
        square_int(:,rline) = mean(abs(data.Data_new(:,idx1:idx2)).^2,2);
        int_square(:,rline) = abs(mean(data.Data_new(:,idx1:idx2),2)).^2;
        square_int_dB(:,rline) = 10*log10(square_int(:,rline));
        
        Latitude_mean(rline)=mean(data.Latitude(idx1:idx2));  %Mean Latitude
        Longitude_mean(rline)=mean(data.Longitude(idx1:idx2)); %Mean Longitude
        GPS_time_ave(rline)=mean(data.GPS_time(idx1:idx2));   %Mean GPS Time
        
        meanbt=nanmean(data.bt_elev_new(idx1:idx2));
        if meanbt==0 | isnan(meanbt)
          continue;               %skip if no ice bottom
        end
        
        
        meansf=nanmean(data.sf_elev_new(idx1:idx2));  %If bottom=surface; skip
        meanbt1=nanmean(bottom_twtt(idx1:idx2));
        if meanbt1==meansf
          continue;
        end
        depth(rline)=meansf-meanbt;
        if depth(rline)<150
          continue;
        end
        %bt_idx_m = find(elev_axis<=meanbt,1,'first');
        bt_idx_m = round(interp1(elev_axis,[1:length(elev_axis)],meanbt));
        
        if isempty(bt_idx_m)
          continue;  %looking for ice bottom
        end
        b1=bt_idx_m-5;               % Peak search index peak bins
        b2=bt_idx_m+5;
        if b2>size(data.Data_new,1)
          b2=size(data.Data_new,1);
        end
        
        [bt_val,bt_idx] = max(square_int(b1:b2,rline));  %Peak Index and Value
        bt_idx = bt_idx_m-5+bt_idx-1;
        bt_pwr = 10*log10(bt_val);
        
        if bt_pwr==0
          continue;
        end
        
        noise_bin1=bt_idx+470;
        noise_bin2=bt_idx+500;
        if noise_bin1>size(data.Data_new,1) || noise_bin2>size(data.Data_new,1)
          noise_bin1=size(data.Data_new,1)-30;
          noise_bin2=size(data.Data_new,1);
        end
        
        noise = mean(square_int_dB(noise_bin1:noise_bin2,rline));
        if bt_pwr-noise<3
          continue;
        end
        %noise = mean(square_int_dB(noise_bin1:noise_bin2,rline));
        risingedge(rline) = bt_idx-1;
        SNR = bt_pwr-noise;
        
        % find risingedge(rline): risingedge and fallingedge(rline):fallingedge for integration in range
        while square_int_dB(risingedge(rline),rline)-noise > 0.05*SNR  && risingedge(rline)>2
          risingedge(rline) = risingedge(rline) - 1;
        end
        if bt_idx-risingedge(rline)>50
          risingedge(rline)=bt_idx-50;  %Max 50 bins away
        end
        
        if bt_idx-risingedge(rline)<3
          risingedge(rline)=bt_idx-3;  %Guard band of 3
        end
        [r,tmp_idx]=(min(square_int_dB(risingedge(rline):bt_idx-3,rline)));
        risingedge=risingedge+tmp_idx-1;
        
        
        fallingedge(rline) = bt_idx+1;
        if fallingedge(rline)>size(data.Data_new,1)
          fallingedge(rline)=size(data.Data_new,1);
        end
        while square_int_dB(fallingedge(rline),rline)-noise > 0.05*SNR & fallingedge(rline)<size(data.Data_new,1) %Depth Bins risingedge(rline) and fallingedge(rline)
          fallingedge(rline) = fallingedge(rline) + 1;
        end
        if fallingedge(rline)-bt_idx>100
          fallingedge(rline)=bt_idx+100;
        end
        if fallingedge(rline)-bt_idx<3
          fallingedge(rline)=bt_idx+3;
          
        end
        coh_index_slopecorr(rline) = sum(int_square(risingedge(rline):fallingedge(rline),rline))/sum(square_int(risingedge(rline):fallingedge(rline),rline)); %Coherence Index
        
        Imeanx=sum((square_int(risingedge(rline):fallingedge(rline),rline)));  %Power summed from risingedge(rline) to fallingedge(rline)
        Abruptiveindex(rline)=bt_val/Imeanx;    %Abruptive Index
        peakindex(rline)=bt_idx;
        
        D1=risingedge(rline);
        D2=fallingedge(rline);
        Dnoise=noise;
        depth(rline)=meansf-meanbt;
        Dpth(idx1:idx2)=data.sf_elev_new(idx1:idx2)-data.bt_elev_new(idx1:idx2);
        % Power(idx1:idx2)=lp(ice_bed_power(idx1:idx2))+2*lp(2*(480+Dpth(idx1:idx2)/sqrt(3.15)));
        
        
        
        B(rline)=2.3*3000/(depth(rline)+2000);
        Gmtrc_loss(rline)=2*lp(2*((data.Elevation(rline)-data.sf_elev_new(rline))+depth(rline)/sqrt(3.15)));
        Padj(rline)=lp(Imeanx)+Gmtrc_loss(rline)+B(rline)*depth(rline)/100;
        Padj2_Na11(rline)=lp(Imeanx)+Gmtrc_loss(rline)+2*11*(depth(rline)-MeanDepth)/1000; %Using Na=11dB/km
        
        %Saving waveform around the bottom above and below
        bx1=bt_idx-50;
        bx2=bt_idx+100;
        if bx1<=0 |bx2>size(square_int_dB,1)
          continue;
        end
        bt.inc_wf_ave(:,rline) = square_int_dB(bt_idx-50:bt_idx+100,rline);
        
      end
      
      if debug_flag
        %     Padj(Padj<(max(Padj)-40))=max(Padj)-40;
        distancenx=(max(distance).*Nx_int*0.5/(size(distance,2)):max(distance).*Nx_int/size(distance,2):max(distance).*Nx_int/size(distance,2)*Nx);
        hold on; figure(4); hold on; plot(distancenx,coh_index_slopecorr,'g','Displayname','SlopeCorrected'); title('Comparision of Coherence Index')
        legend('show');
        distancenx=distancenx/1000;
        figure(7);plot(distancenx,coh_index_slopecorr,'r','Displayname','Coherence Index'); title('Basal Roughness')
        figure(7); hold on; plot(distancenx,Abruptiveindex,'b','Displayname','Abruptive Index');
        xlabel('Distance(km)');
        legend('show');
        figure(6);plot(distancenx,Padj),title('Adjusted Intesnity')
        figure(11);plot(distancenx,coh_index_slopecorr,'r');
        figure(12); plot(distancenx,Abruptiveindex,'b');
        
        figure(4);plot(distancenx,coh_index_ogdata,'Displayname','Original Data');
        
        figure;subplot(3,1,1);plot(distancenx,Padj)
        subplot(3,1,2);plot(distancenx,coh_index_slopecorr)
        subplot(3,1,3);plot(distancenx,Abruptiveindex)
        figure;plot(depth);
      end
      
      clear index;
      index.Latitude_mean=Latitude_mean;
      index.Longitude_mean=Longitude_mean;
      index.GPS_time_ave=GPS_time_ave;
      index.coherence=coh_index_slopecorr;
      index.abruptness=Abruptiveindex;
      index.Padj=Padj;
      index.Padj_Na11=Padj2_Na11;
      index.bt=bt;
    end
    
    if MeanDepth==0 || isnan(MeanDepth) || MeanDepth<150 || isnan(MeanDepth)
      index.Latitude_mean=[];
      index.Longitude_mean=[];
      index.GPS_time_ave=[];
      index.coherence=[];
      index.abruptness=[];
      index.Padj=[];
      index.Padj_Na11=[];
      index.bt=bt;
    end
    
    
    %%
    Greenland.index.Latitude_mean=cat(2, Greenland.index.Latitude_mean,index.Latitude_mean);
    Greenland.index.Longitude_mean=cat(2, Greenland.index.Longitude_mean,index.Longitude_mean);
    Greenland.index.GPS_time_ave=cat(2, Greenland.index.GPS_time_ave,index.GPS_time_ave);
    Greenland.index.coherence=cat(2, Greenland.index.coherence,index.coherence);
    Greenland.index.abruptness=cat(2, Greenland.index.abruptness,index.abruptness);
    Greenland.index.Padj=cat(2,Greenland.index.Padj,index.Padj);
    Greenland.index.Padj_Na11=cat(2,Greenland.index.Padj_Na11,index.Padj_Na11);
    Greenland.index.bt.val=cat(2,Greenland.index.bt.val,index.bt.val);
    Greenland.index.bt.idx=cat(2,Greenland.index.bt.idx,index.bt.idx);
    Greenland.index.bt.waveform=cat(2,Greenland.index.bt.waveform,index.bt.waveform);
    Greenland.index.bt.inc_wf_ave=cat(2,Greenland.index.bt.inc_wf_ave,index.bt.inc_wf_ave);
    
    end
    Greenland.Latitude = cat(2, Greenland.Latitude, data.Latitude);
    Greenland.Longitude = cat(2, Greenland.Longitude, data.Longitude);
    Greenland.Elevation = cat(2, Greenland.Elevation, data.Elevation);
    Greenland.Roll=cat(2,Greenland.Roll,data.Roll);
    Greenland.ice_bed_time = cat(2, Greenland.ice_bed_time, bottom_twtt);
    Greenland.surface_time = cat(2, Greenland.surface_time, surface_twtt);
    Greenland.ice_bed_power = cat(2, Greenland.ice_bed_power, ice_bed_power);
    Greenland.ice_surface_power = cat(2, Greenland.ice_surface_power, ice_surface_power);
    Greenland.segments_length = cat(2,Greenland.segments_length, length(ice_bed_power));
    Greenland.GPS_time = cat(2, Greenland.GPS_time, data.GPS_time);
    %  Greenland.ice_bed_echo = cat(2, Greenland.ice_bed_echo, ice_bed_echo);
    %Greenland.abruptness =  cat(2,Greenland.abruptness, abruptness)
    
    
    
    
    if ~(length(Greenland.ice_bed_time) ==length(Greenland.ice_bed_power))
      keyboard
      warning('Ice bed time not equal to power')
    end
   
  end
  
  % keyboard
  if ~param.save_fig_only    % ~ changed for saving      %Disable not to save data
  
    %Check if needed to flip to extend from inward towards the coast
    if cross_lines
      % if ismember{param.proc_line,[1,3,5,6])
      if Greenland.Longitude(1)>Greenland.Longitude(end)
         flip_line=0;
       else
         flip_line=1;
       end
    else
     if Peterman
         if Greenland.Latitude(1)<Greenland.Latitude(end)
         flip_line=0;
       else
         flip_line=1;
       end
     else
       if Greenland.Latitude(1)>Greenland.Latitude(end)
         flip_line=0;
       else
         flip_line=1;
       end
     end
    end
%     
     %Flip lines 
   if flip_line
     Greenland.GPS_time=flip(Greenland.GPS_time);
     Greenland.Latitude=flip(Greenland.Latitude);
     Greenland.Longitude=flip(Greenland.Longitude);
     Greenland.Roll=flip(Greenland.Roll);
     Greenland.Elevation=flip(Greenland.Elevation);
     Greenland.ice_bed_time=flip(Greenland.ice_bed_time);
     Greenland.surface_time=flip(Greenland.surface_time);
     Greenland.ice_bed_power=flip(Greenland.ice_bed_power);
     Greenland.ice_surface_power=flip(Greenland.ice_surface_power);
     Greenland.segments_length=flip(Greenland.segments_length);
     Greenland.flipped=1;
   end
    
    
   if 1
      geotiff_fn = '/cresis/snfs1/dataproducts/GIS_data/greenland/Landsat-7/Greenland_natural.tif';
      proj = geotiffinfo(geotiff_fn);
      %proj = geotiffinfo('X:\GIS_data\antarctica\Landsat-7\Antarctica_LIMA_480m.tif');
      
      [A CMAP R]= geotiffread(geotiff_fn);
      
      figure(100)
      mapshow(rgb2gray(A),CMAP/1e3);
      xlabel('X (km)');
      ylabel('Y (km)');
      
      if Peterman
        xlim([-350 -50]);
        ylim([-1250 -900]);
        
      else
        xlim([-350 -50]);
        ylim([-2450 -2150]);
      end
      hold on
      clear gps.x gps.y
      [gps.x,gps.y] = projfwd(proj,Greenland.Latitude,Greenland.Longitude);
      
      gps.x = gps.x / 1000;
      gps.y = gps.y / 1000;
      hold on;
      
      scatter(gps.x,gps.y,20,lp(Greenland.ice_bed_power),'fill')
      scatter(gps.x(1),gps.y(1),200,'X');
      %caxis([-15 15])
      colorbar;
      title('Radar line')
      
      if cross_lines
       if Peterman
        save_path=['/cresis/snfs1/scratch/manjish/peterman/images/',sprintf('crossline%d',k),'/',sprintf('Location_Data_flpd_%s_%03d', param1.day_seg, frm)];
       else
         save_path=['/cresis/snfs1/scratch/manjish/jacobshavn/images/',sprintf('crossline%d',k),'/',sprintf('Location_Data_flpd_%s_%03d', param1.day_seg, frm)];
       end
        [save_dir] =fileparts(save_path);
        if ~exist(save_dir,'dir')
          
          mkdir(save_dir);
        end
        saveas(figure(100),save_path,'jpg')
      else
        if Peterman
        save_path=['/cresis/snfs1/scratch/manjish/peterman/images/',sprintf('verticalline%d',k),'/',sprintf('Location_Data_flpd_%s_%03d', param1.day_seg, frm)];
       else
         save_path=['/cresis/snfs1/scratch/manjish/jacobshavn/images/',sprintf('verticalline%d',k),'/',sprintf('Location_Data_flpd_%s_%03d', param1.day_seg, frm)];
       end
        [save_dir] =fileparts(save_path);
        if ~exist(save_dir,'dir')
          
          mkdir(save_dir);
        end
        saveas(figure(100),save_path,'jpg')
      end
      
      
   end
  
   
  
   
   
  %%Truncate lines to prevent aircraft turns that result power loss
  if 0
   max_length=length(Greenland.Latitude);
  start_idx=1;
 stop_idx=max_length;
  
   Greenland.segments_length(1)=Greenland.segments_length(1)-start_idx+1;
   Greenland.segments_length(end)=Greenland.segments_length(end)-length(Greenland.Latitude)-stop_idx;
   
  Greenland.GPS_time=Greenland.GPS_time(start_idx:stop_idx);
  Greenland.Latitude=Greenland.Latitude(start_idx:stop_idx);
  Greenland.Longitude=Greenland.Longitude(start_idx:stop_idx);
  Greenland.Roll=Greenland.Roll(start_idx:stop_idx);
  Greenland.Elevation=Greenland.Elevation(start_idx:stop_idx);
  Greenland.ice_bed_time=Greenland.ice_bed_time(start_idx:stop_idx);
  Greenland.surface_time=Greenland.surface_time(start_idx:stop_idx);
  Greenland.ice_bed_power=Greenland.ice_bed_power(start_idx:stop_idx);
   Greenland.ice_surface_power=Greenland.ice_surface_power(start_idx:stop_idx);
  end
  
  if 1
  %Save settings
  Greenland.settings.day=Day_seg{k};
  Greenland.settings.frms=frms{k};
  if ispc
    git_path=fullfile('H:\scripts\matlab\total process');
  else
     git_path=fullfile('/users/manjish/scripts/matlab/total process');
  end
  Greenland.settings.gitInfo=getGitInfo(git_path);
  
  %keyboard
  %Peterman
  if cross_lines
    disp(sprintf('Saving Cross Line %d\n',k))
    if Peterman
     save(['/cresis/snfs1/scratch/manjish/peterman/radar_w_idx_new/crossline' num2str(k,'%d') '.mat'],'Greenland');
    else
      save(['/cresis/snfs1/scratch/manjish/jacobshavn/radar_w_index/crossline' num2str(k,'%d') '.mat'],'Greenland');
 
    end
  else
    disp(sprintf('Saving Vertical Line %d\n',k))
    if Peterman
        save(['/cresis/snfs1/scratch/manjish/peterman/radar_w_idx_new/verticalline' num2str(k,'%d') '.mat'],'Greenland');
    else
        save(['/cresis/snfs1/scratch/manjish/jacobshavn/radar_w_index/verticalline' num2str(k,'%d') '.mat'],'Greenland');
 
    end
    end
  
  end
 
  end
end

success=true;

return





