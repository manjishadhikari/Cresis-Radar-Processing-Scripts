
clear
close
clc
param.radar.fc = 195e6;
physical_constants;
debug_flag=0;
coh_int=0;


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
Greenland.ice_bed_power_cgl = [];
Greenland.index.Latitude_mean=[];
Greenland.index.Longitude_mean=[];
Greenland.index.GPS_time_ave=[];
Greenland.index.coherence=[];
Greenland.index.abruptness=[];
Greenland.index.Padj=[];
Greenland.Roll=[];
Greenland.index.Padj_Na11=[];
lst=[];
fprintf('Combining all radar lines..\n')
Lines=[67:74 92:103];



for M =1:20
  L=Lines(M);
  if M<9
    tmp=load(['/cresis/snfs1/scratch/manjish/new_jacobshavn/new_lines/crossline', sprintf('%d.mat',L)]);
  else
    % N=M-20;
    tmp= load(['/cresis/snfs1/scratch/manjish/new_jacobshavn/new_lines/verticalline', sprintf('%d.mat',L)]);
  end
  
  if debug_flag
    figure;plot(lp(tmp.Greenland.ice_bed_power));title('Before Coh_int')
  end
  
  clear tmp_data;
  if coh_int~=0
    % tmp.Greenland=coh_integration(tmp.Greenland,coh_int);
    length_data=length(tmp.Greenland.ice_bed_power);
    i=1;
    for m= 1:20: length_data
      %  idx1=(i-1)*numofCohInt+1;
      % idx2=i*numofCohInt;
      idx1=m;
      idx2=m+coh_int;
      if idx2>length_data &  length_data-idx1>=coh_int/2
        idx2= length_data;
      elseif idx2> length_data &  length_data-idx1<coh_int/2
        continue;
      end
      
      tmp_data.GPS_time(i)=nanmean(tmp.Greenland.GPS_time(idx1:idx2));
      tmp_data.Latitude(i)=nanmean(tmp.Greenland.Latitude(idx1:idx2));
      tmp_data.Longitude(i)=nanmean(tmp.Greenland.Longitude(idx1:idx2));
      tmp_data.Elevation(i)=nanmean(tmp.Greenland.Elevation(idx1:idx2));
      tmp_data.Roll(i)=max(tmp.Greenland.Roll(idx1:idx2));
      tmp_data.surface_time(i)=nanmean(tmp.Greenland.surface_time(idx1:idx2));
      tmp_data.ice_bed_time(i)=nanmean(tmp.Greenland.ice_bed_time(idx1:idx2));
      tmp_data.ice_bed_power(i)=nanmean(tmp.Greenland.ice_bed_power(idx1:idx2));
      tmp_data.ice_surface_power(i)=nanmean(tmp.Greenland.ice_surface_power(idx1:idx2));
      i=i+1;
      %  tmp_data.Data(:,i)=nanmean((abs(data.Data(:,idx1:idx2)).^2),2);
      %    i=i+1;
    end
    %    tmp_data.Time=data.Time;
     tmp.Greenland=tmp_data;
  end
 
  
  if debug_flag
    figure;plot(lp(tmp.Greenland.ice_bed_power)); title('After Coh Int')
  end
  
  tmp.Greenland.depth =(-tmp.Greenland.surface_time+tmp.Greenland.ice_bed_time).*3*10^8/(2*sqrt(er_ice));
  tmp.Greenland.surface_height = (tmp.Greenland.surface_time)*c/2;
  tmp.geometric_loss = (2*(tmp.Greenland.surface_height+tmp.Greenland.depth/sqrt(er_ice))).^2;
  tmp.geometric_loss_surface = (2*(tmp.Greenland.surface_height)).^2;
  tmp.Greenland.ice_bed_power_cgl =(tmp.Greenland.ice_bed_power).*tmp.geometric_loss;
  
  if debug_flag
    figure(100);plot(lp(tmp.Greenland.ice_bed_power_cgl))
    power_original=lp(tmp.Greenland.ice_bed_power_cgl);
    figure;plot(tmp.Greenland.depth,lp(tmp.Greenland.ice_bed_power))
    figure;plot( lp(tmp.geometric_loss));
    figure;plot(tmp.Greenland.depth, lp(tmp.geometric_loss))
  end
  
  %   if M<21
  %    load(['/cresis/snfs1/scratch/manjish/new_jacobshavn/surface_roughness/crossline', sprintf('%d.mat',M)]);
  %   else
  %      N=M-20;
  %      load(['/cresis/snfs1/scratch/manjish/new_jacobshavn/surface_roughness/verticalline', sprintf('%d.mat',N)]);
  %   end5
  %   % load(['/cresis/snfs1/scratch/santhosh/thesis/get heights frames Greenland/surface_roughness_v6_' num2str(k,'%03d') '.mat'])
  %   K  = floor(length(tmp.Greenland.ice_bed_power_cgl)/1000);
  %   for l = 1:K
  %     tmp.Greenland.ice_bed_power_cgl(1+(l-1)*1000:(l*1000)) = (tmp.Greenland.ice_bed_power_cgl(1+(l-1)*1000:(l*1000)))./(exp(-(4*pi*r.rms_height(l)/(c/param.radar.fs)*(sqrt(er_ice)-1))^2)*(besseli(0,((4*pi*r.rms_height(l)/(c/param.radar.fs)*(sqrt(er_ice)-1))^2)/2))^2);
  %   end
  %   clear r
  %
  %    if debug_flag
  %   figure(100);hold on;plot(lp( tmp.Greenland.ice_bed_power_cgl))
  %   power_sf_corr=lp( tmp.Greenland.ice_bed_power_cgl);
  %   end
  %
  %  if M<21
  %    load(['/cresis/snfs1/scratch/manjish/new_jacobshavn/bedroughness/crossline', sprintf('%d.mat',M)]);
  %  else
  %     N=M-20;
  %      load(['/cresis/snfs1/scratch/manjish/new_jacobshavn/bedroughness/verticalline', sprintf('%d.mat',N)]);
  %  end
  %  r=rbed;
  %  K  = floor(length(tmp.Greenland.ice_bed_power_cgl)/1000);
  %   for l = 1:K
  %     tmp.Greenland.ice_bed_power_cgl(1+(l-1)*1000:(l*1000)) = (tmp.Greenland.ice_bed_power_cgl(1+(l-1)*1000:(l*1000)))./(exp(-(4*pi*r.rms_height(l)/(c/param.radar.fs*sqrt(er_ice)))^2)*(besseli(0,((4*pi*r.rms_height(l)/(c/param.radar.fs*sqrt(er_ice)))^2)/2))^2);
  %   end
  %   clear r
  
  
  %   if debug_flag
  %   figure(100);hold on; plot(lp( tmp.Greenland.ice_bed_power_cgl))
  %   legend('Original','Sf Roughness Corrected','Bed Roughness Corrected')
  %    power_bed_corr=lp( tmp.Greenland.ice_bed_power_cgl);
  %    figure(101);plot(power_bed_corr-power_sf_corr);title('Bed Roughness Correction values')
  %     figure(102);plot(power_bed_corr-power_original);title('Original minus Roughness Corrected Power')
  %   end
  
  Greenland.GPS_time = cat(2,Greenland.GPS_time, tmp.Greenland.GPS_time);
  Greenland.Latitude = cat(2,Greenland.Latitude, tmp.Greenland.Latitude);
  Greenland.Longitude = cat(2,Greenland.Longitude, tmp.Greenland.Longitude);
  Greenland.Elevation = cat(2,Greenland.Elevation, tmp.Greenland.Elevation);
  Greenland.Roll = cat(2,Greenland.Roll, tmp.Greenland.Roll);
  Greenland.ice_bed_time = cat(2,Greenland.ice_bed_time, tmp.Greenland.ice_bed_time);
  Greenland.surface_time = cat(2,Greenland.surface_time, tmp.Greenland.surface_time);
  Greenland.ice_bed_power = cat(2,Greenland.ice_bed_power, tmp.Greenland.ice_bed_power);
  Greenland.ice_surface_power = cat(2,Greenland.ice_surface_power, tmp.Greenland.ice_surface_power);
  
  Greenland.ice_bed_power_cgl = cat(2,  Greenland.ice_bed_power_cgl, tmp.Greenland.ice_bed_power_cgl);
  %         Greenland.surface_peaks = cat(2,Greenland.surface_peaks, tmp.Greenland.surface_peaks);
  if coh_int==0
    Greenland.segments_length = cat(2,Greenland.segments_length, tmp.Greenland.segments_length);
    Greenland.index.Latitude_mean=cat(2,Greenland.index.Latitude_mean,tmp.Greenland.index.Latitude_mean);
    Greenland.index.Longitude_mean=cat(2,Greenland.index.Longitude_mean,tmp.Greenland.index.Longitude_mean);
    Greenland.index.GPS_time_ave=cat(2,Greenland.index.GPS_time_ave,tmp.Greenland.index.GPS_time_ave);
    Greenland.index.coherence=cat(2,Greenland.index.coherence,tmp.Greenland.index.coherence);
    Greenland.index.abruptness=cat(2,Greenland.index.abruptness,tmp.Greenland.index.abruptness);
    Greenland.index.Padj=cat(2,Greenland.index.Padj,tmp.Greenland.index.Padj);
    if length(tmp.Greenland.index.Padj_Na11)~=length(tmp.Greenland.index.Padj)
      lst(end+1)=M;
    end
    Greenland.index.Padj_Na11=cat(2,Greenland.index.Padj_Na11,tmp.Greenland.index.Padj_Na11);
  end
  
  if (length(Greenland.ice_bed_time)~=length(Greenland.ice_bed_power))
    keyboard
  end
  
end
Greenland.coh_int=coh_int;
Greenland.comb_Lines=Lines;
out_fn=['/cresis/snfs1/scratch/manjish/new_jacobshavn/data/data_2014.mat'];
out_fn_dir=fileparts(out_fn);
if ~exist(out_fn_dir,'dir')
  mkdir(out_fn_dir);
end
save(out_fn,'Greenland','-v7.3');
disp('==Done==')