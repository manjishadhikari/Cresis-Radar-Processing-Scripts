
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


clear
close
clc
dbstop error
cross_lines = 0;
global gRadar
physical_constants;
debug_flag = 0;

%lno=1;    %change the line number

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
  
end

%%  reading the data

for k =1
  
  if cross_lines
    if k<7
      param = read_param_xls(ct_filename_param('rds_param_2010_Greenland_DC8.xls'),Day_seg{k});
      gps_fn = ct_filename_support(param,'','gps',1);
      data_dir = ct_filename_out(param,'','manjish/CSARP_Data');
      layer_dir = ct_filename_out(param,'','old/CSARP_layerData');
    elseif (6<k & k<21)
      param = read_param_xls(ct_filename_param('rds_param_2011_Greenland_P3.xls'),Day_seg{k});
      gps_fn = ct_filename_support(param,'','gps',1);
      data_dir = ct_filename_out(param,'','CSARP_manjish');
      layer_dir = ct_filename_out(param,'','CSARP_layerData');
      
    elseif (20<k & k<24)
      param = read_param_xls(ct_filename_param('rds_param_2012_Greenland_P3.xls'),Day_seg{k});
      gps_fn = ct_filename_support(param,'','gps',1);
      data_dir = ct_filename_out(param,'','CSARP_CSARP_manjish');
      layer_dir = ct_filename_out(param,'','CSARP_post/CSARP_layerData');
    end
    
  else
    
    if k<8
      param = read_param_xls(ct_filename_param('rds_param_2010_Greenland_DC8.xls'),Day_seg{k});
      gps_fn = ct_filename_support(param,'','gps',1);
      data_dir = ct_filename_out(param,'','manjish/CSARP_Data');
      layer_dir = ct_filename_out(param,'','old/CSARP_layerData');
    elseif (7<k & k<16)
      param = read_param_xls(ct_filename_param('rds_param_2011_Greenland_P3.xls'),Day_seg{k});
      gps_fn = ct_filename_support(param,'','gps',1);
      data_dir = ct_filename_out(param,'','CSARP_manjish');
      layer_dir = ct_filename_out(param,'','CSARP_layerData');
      
    elseif (15<k & k<19)
      param = read_param_xls(ct_filename_param('rds_param_2013_Greenland_P3.xls'),Day_seg{k});
      gps_fn = ct_filename_support(param,'','gps',1);
      data_dir = ct_filename_out(param,'','manjish/CSARP_Data');
      layer_dir = ct_filename_out(param,'','CSARP_post/CSARP_layerData');
    else
      param = read_param_xls(ct_filename_param('rds_param_2014_Greenland_P3.xls'),Day_seg{k});
      gps_fn = ct_filename_support(param,'','gps',1);
      data_dir = ct_filename_out(param,'','CSARP_manjish');
      layer_dir = ct_filename_out(param,'','CSARP_post/CSARP_layerData');
    end
    
    
  end
  
  frame_fn = ct_filename_support(param,'','frames');
  load(frame_fn);
  records_file = ct_filename_support(param,'','records');
  
  
  
  Greenland.GPS_time = [];
  Greenland.Latitude = [];
  Greenland.Longitude = [];
  Greenland.Elevation = [];
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
    layer_fn = fullfile(layer_dir,sprintf('Data_%s_%03d.mat', param.day_seg, frm));
    fprintf('Loading data %s\n', layer_fn);
    data_fn = fullfile(data_dir,sprintf('Data_%s_%03d.mat', param.day_seg, frm));
    
    % Load the file
    data = load(data_fn);
    tmp = load(layer_fn);
    
    surface_twtt=interp1(tmp.GPS_time,tmp.layerData{1}.value{2}.data , data.GPS_time);
      bottom_twtt = interp1(tmp.GPS_time,tmp.layerData{2}.value{2}.data , data.GPS_time);
    %         Elevation = interp1(tmp.GPS_time,tmp.Elevation,data.GPS_time);
    %         surface_twtt = data.Surface + ((data.Elevation-Elevation)/(c/2));
    %         filter_length = 100;
    %         surface_twtt = sgolayfilt(surface_twtt, 2,filter_length+1, hanning(filter_length+1));
    dt= data.Time(2)-data.Time(1);
    index = round((surface_twtt-data.Time(1))/dt);
    ice_surface_power  = zeros(1,length(data.Surface));
    
    for i = 1:length(surface_twtt)
      if isnan(surface_twtt(i)) || isinf(surface_twtt(i))
        ice_surface_power(i)= nan;
        continue
      else
        [surface_power idx] = max(sqrt(data.Data(index(i):index(i),i).*conj(data.Data(index(i):index(i),i))));
        surface_index = idx + index(i)-1;
        if surface_power  == 0
          ice_surface_power(i)= nan;
          surface_twtt(i) = nan;
        else
          ice_surface_power(i) = data.Data(surface_index);
          surface_twtt(i) =  interp1([1:length(data.Time)],data.Time,surface_index);
        end
      end
    end
    if any(ice_surface_power ==0 )
      keyboard
    end
    
    clear index
            index = round((bottom_twtt-data.Time(1))/dt);
            ice_bed_power  = zeros(1,length(data.Surface));
            bottom_twtt_n = zeros(1,length(data.Surface));
            
            for i = 1:length(bottom_twtt)
                if isnan(bottom_twtt(i)) || isinf(bottom_twtt(i))
                    bottom_twtt_n(i)= nan;
                    continue
                else
                    [bed_power idx] = max(data.Data(index(i)-1:index(i)+1,i));
                    bed_index = idx + index(i)-1-1;
                    bottom_twtt_n(i) =  interp1([1:length(data.Time)],data.Time,bed_index);
                end
            end
            filter_length = 100;
            bottom_twtt_n = sgolayfilt(bottom_twtt_n, 2,filter_length+1, hanning(filter_length+1));
            
            index = round((bottom_twtt_n-data.Time(1))/dt);
            for i = 1:length(bottom_twtt_n)
                
                if isnan(bottom_twtt_n(i)) || isinf(bottom_twtt_n(i))
                    ice_bed_power(i)= nan;
                    continue
                else
                    [bed_power idx] = max(data.Data(index(i):index(i),i));
                    bed_index = idx + index(i)-1;
                    idx1=bed_index+200;
                    idx2=bed_index+500;
                    if idx1>size(data.Data,1) | idx2>size(data.Data,1)
                      idx2=size(data.Data,1);
                      idx1=idx2-300;
                    end
                    N = mean(sqrt(data.Data(idx1:idx2,i)));
                    SNR = 10*log10((sqrt(bed_power))/N);
                    if SNR > 3
                        ice_bed_power(i) = bed_power;
                        %
                    else
                        ice_bed_power(i) = nan;
                        bottom_twtt_n(i) = nan;
                        continue ;
                    end
                end
            end
      
    
%     clear ice_bed_echo
%     if cross_lines
%       load(['/cresis/snfs1/scratch/manjish/test/',sprintf('Data_%s_%03d.mat', param.day_seg, frm)]);
%     else
%       load(['/cresis/snfs1/scratch/manjish/test/',sprintf('Data_%s_%03d.mat', param.day_seg, frm)]);
%     end
    if exist('bottom_twtt_n','var')
      bottom_twtt = bottom_twtt_n;
    end
    
    
    index = round((bottom_twtt-data.Time(1))/dt);
    ice_bed_power  = zeros(1,length(data.Surface));
    abruptness = zeros(1,length(data.Surface));
    
    
    
    for i = 1:length(bottom_twtt)
      if isnan(bottom_twtt(i)) || isinf(bottom_twtt(i))
        ice_bed_power(i)= nan;
        ice_bed_echo{i} = nan;
        abruptness(i) = nan;
        continue
      else
        [bed_power idx] =  max(sqrt(data.Data(index(i)-1:index(i)+1,i).*conj(data.Data(index(i)-1:index(i)+1,i))));
        bed_index = idx + index(i)-1-1;
        ice_bed_power(i)=data.Data(bed_index);
        
      end
      
    end
    
    
    
    if any(ice_bed_power ==0 )
      keyboard
    end
    
    
    if debug_flag
      figure(1);
      clf
      imagesc([],data.Time*1e6,lp(sqrt(data.Data.*conj(data.Data))));
      hold on;plot(surface_twtt*1e6,'--');
      hold on;plot(bottom_twtt*1e6,'--');
      %       keyboard
    end
    
    %Coherence Index and Abruptive Index Calculation
    % fc= 195e6;  %radar center frequemcy
    fc=(param.radar.wfs(1).f1+param.radar.wfs(1).f1)/2;
    % fs=1.1111e8;  %Sampling frequency
    fs=param.radar.fs;
    p=4.99;  %Radar Half pulse width in air
    
    index=indices_calculation(data_fn,surface_twtt,bottom_twtt,fc,fs,p,debug_flag);
    
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
    
    Greenland.Latitude = cat(2, Greenland.Latitude, data.Latitude);
    Greenland.Longitude = cat(2, Greenland.Longitude, data.Longitude);
    Greenland.Elevation = cat(2, Greenland.Elevation, data.Elevation);
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
    end
  end
  disp(sprintf('Saving Line %d\n',k))
  % keyboard
  
  if cross_lines
    save(['/cresis/snfs1/scratch/manjish/peterman/radar_w_index/crossline' num2str(k,'%d') '.mat'],'Greenland');
  else
    save(['/cresis/snfs1/scratch/manjish/peterman/radar_w_index/verticalline' num2str(k,'%d') '.mat'],'Greenland');
  end
end








