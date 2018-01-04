%Script: verification_of_picked_ice_bed_interface
%
% Purpose: This script is used to verify the ice-bed interface if it is picked accurately or if changes can be made to it
%
% Input data / processed on: 2011_Antarctica_TO pulse compressed data
%
%
% Output/data product: two-way travel time to ice-bed interface
%
% Output format: bed_layer_Data_<straight_line/cross_line_number>_<corresponding_frame_number> 
%
% Location saved: Y:\santhosh\byrd\bed_layer (for data corresponding to straight flight paths) and
% Y:\santhosh\byrd\bed_layer\cross_lines (for data corresponding to cross flight paths)
% 
% See also saving_the_required_data_features
%%
clearvars -except p q
close
startup
clc
dbstop error

if 1
    tic
    debug_flag = 0;
    physical_constants;
   
%      Day={'20100324_01','20100324_02','20100324_03','20100324_04','20110429_01','20110429_02','20110507_01','20110507_02','20120330_01','20120516_01','20140512_01','20140505_01'};
%   Day_seg=repelem(Day,[3,1,1,1,5,1,7,1,1,2,1,5]);
%   frms={[36,37],[39,40],[42],[1 2],[1 2],[1 2],[9:12],[13:16],[17:20],[21:24],[25:28],[10:11],[10:14],[15:18],[19:22],[23:26],[27:30],[31:34],[35:37],[1:4],[5 6 7],[13:16],[79:81],[10 11],[12 13],[33 34 35],[38:40],[55 56],[59 60]};
%  
   Day={'20100324_01','20110429_01','20110429_02','20110507_02','20130420_02','20140512_01','20140505_01'};
  Day_seg=repelem(Day,[7,2,4,2,3,2,1]);
   frms={[11 12],[14 15],[17 18],[20 21],[23 24] [30 31],[33 34],[30:32],[33 34],[18:21],[12:15],[2:5],[6:9],[6:8],[17:20],[3 4],[9],[11],[12 13],[17 18],[15 16]};

  
  
  %  Day={'20100324_01','20120516_01','20130420_02','20140512_01','20140505_01'};
   % Day_seg=repelem(Day,[1,2,3,3,6]);
  %  frms={[39 40 ],[13:16],[79:81],[3:4],[8:9],[11],[12,13],[17,18],[10,11],[12,13],[15,16],[33:35],[38:40],[55,56],[59,60]};
    %Day_seg = {'20111213_05','20111201_04','20111214_04','20111213_01','20111213_02','20111212_02','20111212_02','20111212_02','20111205_03','20111205_02','20111206_01','20111205_02','20111206_01','20111205_02','20111213_04','20111201_02','20111206_01','20111228_01','20111228_02','20111206_02'};
    %frms = {[10 11 12 13],[3 4 5 6],[13 14 15 16],[2 3 4],[2 3 4 5],[11 12 13],[7 8 9],[3 4 5],[2 3 4],[10 11 12 13],[7 8 9 10],[6 7 8 9],[12 13 14 15],[1 2 3 4],[5 6 7 8],[5 6 7 8],[1 2 3 4 5],[1 2 3],[1,2],[4 5 6]};
    
    
    
    %       For the cross lines
    %     Day_seg = {'20111213_04','20111213_05','20111213_05','20111201_02','20111201_04','20111206_02','20111214_04','20111214_04','20111212_01','20111214_04','20111206_02','20111216_01','20111206_02','20111222_02','20111212_01','20111209_02','20111212_01','20111220_03','20111220_03','20111219_01','20111219_01','20111219_01','20111218_02','20111213_04','20111212_01','20111219_04','20111222_01','20111222_01','20111214_04'};
    %     frms = {[9 10 11 12],[2 3],[8 9],[9 10 11 12 13],[1 2],[16 17 18 19 20],[11 12],[9 10],[15 16 17 18],[7 8],[12 13 14],[11 12],[9 10],[2 3],[10 11 12],[4 5],[5 6 7],[6 7],[3 4],[13 14 15],[8 9 10 11],[2 3 4 5 6],[11 12 13],[2 3 4],[1 2],[1 2 3 4 5 6 7 8],[7 8],[3 4 5],[2,3,4]};
    
    
    
    global gRadar
    for K =16:21
      if K<8  
      param = read_param_xls(ct_filename_param('rds_param_2010_Greenland_DC8.xls'),Day_seg{K});
        gps_fn = ct_filename_support(param,'','gps',1);
        data_dir = ct_filename_out(param,'','manjish/CSARP_data');
        layer_dir = ct_filename_out(param,'','old/CSARP_layerData');  
      elseif (7<K & K<16)
            param = read_param_xls(ct_filename_param('rds_param_2011_Greenland_P3.xls'),Day_seg{K});
        gps_fn = ct_filename_support(param,'','gps',1);
        data_dir = ct_filename_out(param,'','CSARP_manjish');
        layer_dir = ct_filename_out(param,'','CSARP_layerData');
     
      elseif (15<K & K<19)
           param = read_param_xls(ct_filename_param('rds_param_2013_Greenland_P3.xls'),Day_seg{K});
        gps_fn = ct_filename_support(param,'','gps',1);
        data_dir = ct_filename_out(param,'','manjish/CSARP_Data');
        layer_dir = ct_filename_out(param,'','CSARP_post/CSARP_layerData');
      else
         param = read_param_xls(ct_filename_param('rds_param_2014_Greenland_P3.xls'),Day_seg{K});
        gps_fn = ct_filename_support(param,'','gps',1);
        data_dir = ct_filename_out(param,'','CSARP_manjish');
        layer_dir = ct_filename_out(param,'','CSARP_post/CSARP_layerData');
      end
        frame_fn = ct_filename_support(param,'','frames');
        load(frame_fn);
        records_file = ct_filename_support(param,'','records');
        
        for f =1:length(frms{K})
            
            frm = frms{K}(f);
            layer_fn = fullfile(layer_dir,sprintf('Data_%s_%03d.mat', param.day_seg, frm));
            fprintf('Loading data %s\n', layer_fn);
            data_fn = fullfile(data_dir,sprintf('Data_%s_%03d.mat', param.day_seg, frm));
            
            data = load(data_fn);
            tmp = load(layer_fn);
            
            
            bottom_twtt = interp1(tmp.GPS_time,tmp.layerData{2}.value{2}.data , data.GPS_time);
            surface_twtt=interp1(tmp.GPS_time,tmp.layerData{1}.value{2}.data , data.GPS_time);
            
             idx= find(bottom_twtt>data.Time(end));
              bottom_twtt(idx)=nan;
              
%             Elevation = interp1(tmp.GPS_time,tmp.Elevation,data.GPS_time);
%             bottom_twtt = bed_twtt + ((data.Elevation-Elevation)/(c/2));
%             surface_twtt = data.Surface + ((data.Elevation-Elevation)/(c/2));
%             filter_length = 500;
%             bottom_twtt = sgolayfilt(bottom_twtt, 2,filter_length+1, hanning(filter_length+1)); % points used should be adjusted, more points for inland.
%             surface_twtt = sgolayfilt(surface_twtt, 2,filter_length+1, hanning(filter_length+1));
%             
%             figure;imagesc([],data.Time,lp(data.Data));
%             hold on; plot(surface_twtt);
%             hold on; plot(bottom_twtt);
            dt= data.Time(2)-data.Time(1);
            index = round((surface_twtt-data.Time(1))/dt);
            ice_surface_power  = zeros(1,length(data.Surface));
            
            
            for i = 1:length(surface_twtt)
                if isnan(surface_twtt(i)) || isinf(surface_twtt(i))
                    ice_surface_power(i)= nan;
                    continue
                else
                    [surface_power idx] = max(data.Data(index(i)-1:index(i)+1,i));
                    surface_index = idx + index(i)-1-1;
                    ice_surface_power(i) = surface_power;
                    surface_twtt(i) =  interp1([1:length(data.Time)],data.Time,surface_index);
                end
            end
            
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
            
            
            if debug_flag
                
                figure(1);
                clf
                subplot(2,1,1)
                imagesc([],data.Time*1e6,lp(sqrt(data.Data)));
                                 hold on;plot(surface_twtt*1e6,'--');
%              %   ylim([min(bottom_twtt*1e6)-5 max(bottom_twtt*1e6)+5])
                 subplot(2,1,2)
                 imagesc([],data.Time*1e6,lp(sqrt(data.Data)));
                 hold on;plot(bottom_twtt*1e6,'b');
                 hold on;plot(bottom_twtt_n*1e6,'r');
%            
            end
              save(['/cresis/snfs1/scratch/manjish/test/' ,sprintf('Data_%s_%03d.mat', param.day_seg, frm)],'bottom_twtt_n');
                
%             if choice ==1
%                 save(['/cresis/snfs1/scratch/manjish/test/' ,sprintf('Data_%s_%03d.mat', param.day_seg, frm)],'bottom_twtt_n');
%                 %                            save(['/cresis/snfs1/scratch/santhosh/byrd/bed_layer/cross_lines/bed_layer_' sprintf('Data_%d_%03d.mat',K, f)],'bottom_twtt_n');
%                 
%             else
%                 save(['/cresis/snfs1/scratch/manjish/byrd/bed_layer/bed_layer_' sprintf('Data_%d_%03d.mat',K, f)],'bottom_twtt');
%                 %                                save(['/cresis/snfs1/scratch/santhosh/byrd/bed_layer/cross_lines/bed_layer_' sprintf('Data_%d_%03d.mat',K, f)],'bottom_twtt');
%                 
%             end
         end
    end
end

