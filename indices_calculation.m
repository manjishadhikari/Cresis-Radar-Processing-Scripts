

function [index]=indices_calculation(data_fn,surface_twtt,bottom_twtt,fc,fs,p,debug_flag)

%% Determination of coherence Index and abruptive index

if 0
  warning('Test Running ')
  keyboard
  A=load('/cresis/snfs1/dataproducts/ct_data/rds/2011_Greenland_P3/CSARP_manjish/20110507_02/Data_20110507_02_017.mat');
  layer=load('/cresis/snfs1/dataproducts/ct_data/rds/2011_Greenland_P3/CSARP_layerData/20110507_02/Data_20110507_02_017.mat');
   surface_twtt=interp1(layer.GPS_time,layer.layerData{1}.value{2}.data , A.GPS_time);
  load('/cresis/snfs1/scratch/manjish/test/Data_20110507_02_017.mat');
    bottom_twtt = bottom_twtt_n;
    fc=210000000;
    fs= 1.1111e+08;
    p=4.99;
    debug_flag=1;
end

physical_constants;

 A=load(data_fn);
 %L=load(layer_fn);


% distance=geodetic_to_along_track(A.Latitude,A.Longitude,A.Elevation); %Horizontal Distance
%
% bttime= interp1(L.GPS_time, L.layerData{2}.value{2}.data,A.GPS_time,'linear','extrap');
% surface= interp1(L.GPS_time, L.layerData{1}.value{2}.data,A.GPS_time,'linear','extrap');
% Elevation=interp1(L.GPS_time,L.Elevation,A.GPS_time,'linear','extrap');
%
% if debug_flag
% figure(40),imagesc(A.GPS_time,A.Time,10*log10(abs(A.Data)));
% hold on;plot(A.GPS_time,surface,'--');plot(A.GPS_time,bttime,'--');
% end
%
% bottom_twtt=bttime+(A.Elevation-Elevation)/(c/2);
% surface_twtt=surface+(A.Elevation-Elevation)/(c/2);
%
%
% filter_length=100;
% bottom_twtt=sgolayfilt(bottom_twtt,2,filter_length+1,hanning(filter_length+1));
% surface_twtt=sgolayfilt(surface_twtt,2,filter_length+1,hanning(filter_length+1));
%
% if debug_flag
% figure(40),imagesc(A.GPS_time,A.Time,10*log10(abs(A.Data)));
% hold on;plot(A.GPS_time,surface_twtt,'--');plot(A.GPS_time,bottom_twtt,'--');
% end
%
% dt= A.Time(2)-A.Time(1);
% index = round((surface_twtt-A.Time(1))/dt);
%             ice_surface_power  = zeros(1,length(A.Surface));
%

%             for i = 1:length(surface_twtt)
%                   if isnan(surface_twtt(i)) || isinf(surface_twtt(i))
%                     ice_surface_power(i)= nan;
%                     continue
%                   else
%
%
%
%                 [surface_power idx] = max(A.Data(index(i)-1:index(i)+1,i));
%                 surface_index = idx + index(i)-1-1;
%                 ice_surface_power(i) = surface_power;
%                 surface_twtt(i) =  interp1([1:length(A.Time)],A.Time,surface_index);
%
%                   end
%             end
%
%
%             %           dt= data.Time(2)-data.Time(1);
%             index = round((bottom_twtt-A.Time(1))/dt);
%             ice_bed_power  = zeros(1,length(A.Surface));
%
%             for i = 1:length(bottom_twtt)
%                 if isnan(bottom_twtt(i)) || isinf(bottom_twtt(i))
%                     ice_bed_power(i)= nan;
%                     continue
%                 else
%                     [bed_power idx] = max(A.Data(index(i)-5:index(i)+5,i));
%                     bed_index = idx + index(i)-15-1;
%                     %                     ice_bed_echos(:,i)=
%                     %                     data.Data(index(i)-15:index(i)+15,i);
%                     N = mean(sqrt(A.Data(bed_index+200:bed_index+500,i)));
%                     SNR = 10*log10((sqrt(bed_power))/N);
%                     if SNR > 3
%                         ice_bed_power(i) = bed_power;
%                         bottom_twtt(i) =  interp1([1:length(A.Time)],A.Time,bed_index);
%                     else
%                         %                         keyboard
%                         ice_bed_power(i) = nan;
%                         bottom_twtt(i) = nan;
%                         continue ;
%                     end
%                 end
%             end

if debug_flag
  figure(1),imagesc(A.GPS_time,A.Time,10*log10(abs(A.Data)));
  hold on;plot(A.GPS_time,surface_twtt,'--');plot(A.GPS_time,bottom_twtt,'--');
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
  distance=geodetic_to_along_track(A.Latitude,A.Longitude,A.Elevation);
  Nx_int = floor(Nx_int_dist/(distance(10)-distance(9)));                %No of lines integrated
  Nx0 = size(A.Data,2);
  Nx = floor(Nx0/Nx_int);                                         %Total lines after integration
  Nx_mod = mod(Nx0,Nx_int);
  if Nx_mod>= Nx_int/2;
    Nx = Nx + 1;
  end
  square_int = zeros(size(A.Data,1),Nx); % incoherent integration, take square of abs first, then sum (phase info lost)
  int_square = zeros(size(A.Data,1),Nx); % coherent integration, sum complex data first, then take square of abs (phase info remains)
  coh_index_ogdata = zeros(1,Nx);
  Abruptiveindex=zeros(1,Nx);
  risingedge=zeros(1,Nx);
  fallingedge=zeros(1,Nx);
  Padj=zeros(1,Nx);
  for rline =1:Nx
    idx1 = (rline-1)*Nx_int + 1;
    idx2 = rline*Nx_int;
    if Nx0 - idx2 > 0 && Nx0 - idx2 < Nx_int/2;
      idx2 = Nx0;
    else
      idx2 = min(idx2,Nx0);
    end
    square_int(:,rline) = mean(abs(A.Data(:,idx1:idx2)).^2,2);
    square_int_dB = 10*log10(square_int(:,rline));
    int_square(:,rline) = (abs(mean(A.Data(:,idx1:idx2),2))).^2;
    
    meanbt=nansum(bottom_twtt(isfinite(bottom_twtt(idx1:idx2))));
    count=nansum(isfinite(bottom_twtt(idx1:idx2)));
    if meanbt==0 || count==0
      continue;
    end
    meanbt=meanbt/count;
    bt_idx_m = find(A.Time>meanbt,1,'first');
    
    if isempty(bt_idx_m)
      bt_idx_m=round(size(A.Data,1)/2);
    end
    b1=bt_idx_m-5;
    b2=bt_idx_m+5;
    if b2>size(A.Data,1)
      b2=size(A.Data,1);
    end
    
    [bt_val,bt_idx] = max(square_int(b1:b2,rline));
    bt_idx = bt_idx_m-5+bt_idx-1;                 %Peak Index
    bt_pwr = 10*log10(bt_val);                     %Peak Value
    square_int_dB = 10*log10(square_int(:,rline));
    noise_bin1=bt_idx+70;
    noise_bin2=bt_idx+100;
    if noise_bin1>size(A.Data,1)
      noise_bin1=size(A.Data,1)-30;
      noise_bin2=size(A.Data,1);
    end
    
    if noise_bin2>size(A.Data,1)
      noise_bin1=size(A.Data,1)-30;
      noise_bin2=size(A.Data,1);
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
    if fallingedge(rline)>size(A.Data,1)
      fallingedge(rline)=size(A.Data,1) ;
    end
    while square_int_dB(fallingedge(rline))-noise > 0.05*SNR &&  fallingedge(rline)<size(A.Data,1)      %Depth Bins risingedge(rline) and fallingedge(rline)
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
  %% Correction for A.Elevation
  
  sf=surface_twtt; %test
  sf_new = zeros(size(sf));
  bt_new = zeros(size(bottom_twtt));
  A.Elevation_new = zeros(size(A.Elevation));
  A.sf_elev_new = zeros(size(A.Elevation));
  A.bt_elev_new = zeros(size(A.Elevation));
  
  % 1)Remove data before zero time
  negative_bins = A.Time < 0;
  A.Time_new = A.Time(~negative_bins);
  A.Data_new = A.Data(~negative_bins,:);
  
  % 2)Create A.Elevation axis to interpolate to
  [max_elev,max_elev_idx] = max(A.Elevation);
  min_elev = min(A.Elevation - sf*c/2 - (A.Time_new(end)-sf)*c/2/sqrt(er_ice));
  dt = A.Time(2)-A.Time(1);
  dr = dt * c/2 / sqrt(er_ice);
  dt_air = dr/(c/2);
  elev_axis = max_elev:-dr:min_elev;
  new_time = zeros(length(elev_axis),length(A.Elevation));
  
  % 3)Zero pad data to create space for interpolated data
  zero_pad_len = length(elev_axis) - length(A.Time_new);
  A.Data_new = cat(1,A.Data_new,zeros(zero_pad_len,size(A.Data_new,2)));
  
  % 4)Determine the corrections to apply to A.Elevation and layers
  dRange = max_elev - A.Elevation;
  dBins = round(dRange / (c/2) / dt);
  dtime = dRange/(c/2);
  
  for rline = 1:size(A.Data_new,2)
    % Determine A.Elevation bins before surface
    sf_elev = A.Elevation(rline) - sf(rline) * c/2;
    time0 = -(max_elev - A.Elevation(rline))/(c/2);
    last_air_idx = find(elev_axis > sf_elev,1,'last');
    new_time_tmp = (time0 + dt_air*(0:last_air_idx-1)).';
    if last_air_idx < length(elev_axis)
      % Determine A.Elevation bins after surface
      dt_ice = dr/(c/2/sqrt(er_ice));
      first_ice_idx = last_air_idx + 1;
      time0 = sf(rline) + (sf_elev - elev_axis(first_ice_idx))/(c/2/sqrt(er_ice));
      new_time(:,rline) = cat(1,new_time_tmp, (time0 + dt_ice*(0:length(elev_axis)-length(new_time_tmp)-1)).');
    end
    A.Data_new(:,rline) = interp1(A.Time_new, A.Data_new(1:length(A.Time_new),rline), new_time(:,rline), 'linear',0);
    A.Elevation_new(rline) = A.Elevation(rline) + dRange(rline);
    sf_new(rline) = sf(rline) + dtime(rline);
    bt_new(rline) = bottom_twtt(rline) + dtime(rline);
    A.sf_elev_new(rline) = A.Elevation_new(rline) - sf_new(rline)*c/2;
    A.bt_elev_new(rline) = A.sf_elev_new(rline) - (bt_new(rline)-sf_new(rline))*c/2/sqrt(er_ice);
    
    
    
  end
  
  if debug_flag
    fh = figure(2);
    figure(fh);imagesc([],elev_axis,10*log10(abs(A.Data_new).^2));title('A.Elevation Correction Data')
    ax = gca;
    ax.YDir = 'normal';
    hold on;plot(A.sf_elev_new,'--');plot(A.bt_elev_new,'--');
    bt_slope = diff(A.bt_elev_new)./diff(distance);
  end
  
  
  
  %% Coherene Index Calculation with A.Elevation Correction
  
  square_int = zeros(size( A.Data_new,1),Nx); % incoherent integration, take square of abs first, then sum
  int_square = zeros(size(A.Data_new,1),Nx); % coherent integration, sum complex data first, then take square of abs
  coh_index_elevcorr = zeros(1,Nx);
  for rline =1:Nx
    idx1 = (rline-1)*Nx_int + 1;
    idx2 = rline*Nx_int;
    if Nx0 - idx2 > 0 & Nx0 - idx2 < Nx_int/2;
      idx2 = Nx0;
    else
      idx2 = min(idx2,Nx0);
    end
    square_int(:,rline) = mean(abs(A.Data_new(:,idx1:idx2)).^2,2);
    int_square(:,rline) = abs(mean(A.Data_new(:,idx1:idx2),2)).^2;
    square_int_dB=lp(square_int(:,rline));
    
    meanbt=nansum(A.bt_elev_new(idx1:idx2));
    count=nansum(isfinite(A.bt_elev_new(idx1:idx2)));
    meanbt=meanbt/count;
    if meanbt==0 || count==0 ||isnan(meanbt)
      continue;
    end
    
    bt_idx_m = find(elev_axis<=meanbt,1,'first');
    
    if isempty(bt_idx_m)
      bt_idx_m=round(size(A.Data_new,1)/2);      %looking for ice bottom
    end
    b1=bt_idx_m-5;
    b2=bt_idx_m+5;
    if b2>size(A.Data_new,1)
      b2=size(A.Data_new,1);
    end
    
    [bt_val,bt_idx]=max(square_int(b1:b2,rline));
    bt_idx = bt_idx_m-5+bt_idx-1;
    bt_pwr = 10*log10(bt_val);
    noise_bin1=bt_idx+70;
    noise_bin2=bt_idx+100;
    if noise_bin1>size(A.Data_new,1)
      noise_bin1=size(A.Data_new,1)-30;
      noise_bin2=size(A.Data_new,1);
    end
    if noise_bin2>size(A.Data_new,1)
      noise_bin2=size(A.Data_new,1);
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
    if fallingedge(rline)>size(A.Data_new,1)
      fallingedge(rline)=size(A.Data_new,1) ;
    end
    while square_int_dB(fallingedge(rline))-noise > 0.05*SNR && fallingedge(rline)<size(A.Data_new,1)      %Depth Bins risingedge(rline) and fallingedge(rline)
      fallingedge(rline) = fallingedge(rline) + 1;
    end
    Imeanx=sum((square_int(risingedge(rline):fallingedge(rline),rline)));
    Abruptiveindex(rline)=bt_val/Imeanx;
    coh_index_elevcorr(rline) = sum(int_square(risingedge(rline):fallingedge(rline),rline))/sum(square_int(risingedge(rline):fallingedge(rline),rline));
  end
  
  if debug_flag
    hold on; figure(4); hold on; plot(distancenx,coh_index_elevcorr,'r','Displayname','A.Elevation Corrected');
    figure(41); hold on; plot(distancenx,Abruptiveindex,'r','Displayname','A.Elevation Corrected'); grid;
  end
  
  %% Ice Slope Correction
  
  Nx0 = size(A.Data_new,2);
  Nx = floor(Nx0/Nx_int);
  Nx_mod = mod(Nx0,Nx_int);
  if Nx_mod>= Nx_int/2;
    Nx = Nx + 1;
  end
  slopeerror=zeros(1,size(A.Data_new,1));
  slopeval=zeros(1,size(A.Data_new,1));
  angle=zeros(1,Nx);
  
  for rline =1:Nx
    
    idx1 = (rline-1)*Nx_int + 1;
    idx2 = rline*Nx_int;
    if Nx0 - idx2 > 0 && Nx0 - idx2 < Nx_int/2;
      idx2 = Nx0;
    else
      idx2 = min(idx2,Nx0);
    end
    
    if isinf(mean(A.bt_elev_new(idx1:idx2)))  % If No Ice Bottom, skip
      continue;
    end
    if isnan(mean(A.bt_elev_new(idx1:idx2)))  % If No Ice Bottom, skip
      continue;
    end
    
    p = polyfit(distance(idx1:idx2),A.bt_elev_new(idx1:idx2),1);    % Polygonal fitting
    slopeval(idx1:idx2) = polyval(p,distance(idx1:idx2));           % Ploygonal values after fitting
    base = distance(idx2)-distance(idx1);
    perpendicular = slopeval(idx2)-slopeval(idx1);
    hypotenuse = sqrt(base^2+perpendicular^2);
    angle(rline)=asin(perpendicular/hypotenuse)*180/pi;
    slopeerror = slopeval(idx1:idx2)-slopeval(idx1);         % Error betn original and fitting line
    dtime = 2*slopeerror/c/sqrt(er_ice);
    if debug_flag
      if rline==1
        figure(8);plot([idx1:idx2],A.bt_elev_new(idx1:idx2));
        hold on;plot([idx1:idx2],slopeval(idx1:idx2),'r--')
        figure(9);plot(10*log10(abs(A.Data_new(:,idx2)).^2))
      end
    end
    for idx = idx1:idx2
      A.Data_new(:,idx) = interp1(elev_axis, A.Data_new(:,idx), elev_axis + slopeerror(idx-idx1 +1), 'linear',0);
      A.Elevation_new(idx) = A.Elevation_new(idx) - slopeerror(idx-idx1 +1);
      A.sf_elev_new(idx) = A.sf_elev_new(idx) - slopeerror(idx-idx1 +1);
      A.bt_elev_new(idx) = A.bt_elev_new(idx) - slopeerror(idx-idx1 +1);
      sf_new(idx) = sf_new(idx) + dtime(idx-idx1 +1);
      bt_new(idx) = bt_new(idx) + dtime(idx-idx1 +1);
    end
    if debug_flag
      if rline==1
        figure(8);hold on;plot([idx1:idx2],A.bt_elev_new(idx1:idx2),'g--');
        figure(9);hold on;plot(10*log10(abs(A.Data_new(:,idx2)).^2),'g--');
      end
    end
  end
  
  if debug_flag
    fh = figure(3);
    %figure(fh);imagesc([1:size(A.Data_new,2)],elev_axis,10*log10(abs(A.Data_new).^2));title('Slope Correction Data')
    figure(fh);imagesc([],elev_axis,10*log10(abs(A.Data_new).^2));title('Slope Correction Data')
    ax = gca;
    ax.YDir = 'normal';
    hold on;plot(A.sf_elev_new,'--');plot(A.bt_elev_new,'--');
  end
  
  
  %
  % Truncate data around ice bottom within bt.range_bins
  bt_range_bins =[-50:100];
  bt.val = NaN*ones(1,size(A.Data_new,2));
  bt.idx = NaN*ones(1,size(A.Data_new,2));
  bt.waveform = NaN*ones(length(bt_range_bins),size(A.Data_new,2));
  bt.inc_wf_ave=NaN*ones(size(bt.waveform,1),Nx);
  for rline = 1:size(A.Data_new,2)
    if ~isnan(A.bt_elev_new(rline)) & ~isinf(A.bt_elev_new(rline))
      % bt.idx(rline) = find(elev_axis<=A.bt_elev_new(rline),1,'first');
      bt.idx(rline) = round(interp1(elev_axis,[1:length(elev_axis)],A.bt_elev_new(rline)));
      bt.val(rline) = A.Data_new(bt.idx(rline),rline);
      first_idx = bt.idx(rline) + bt_range_bins(1);
      last_idx = bt.idx(rline) + bt_range_bins(end);
      if first_idx < 1 | last_idx>size(A.Data_new,1)
        bt.idx(rline) = NaN;
        bt.val(rline) = NaN;
        continue
      end
      lower=bt.idx(rline)+bt_range_bins(1);
      upper=bt.idx(rline)+bt_range_bins(end);
      if upper >size(A.Data_new,1)
        upper=size(A.Data_new,1);
      end
      bt.waveform(1:51+upper-bt.idx(rline),rline) = A.Data_new(lower:upper, rline);
    else
      continue
    end
  end
  
  
  if debug_flag
    figure(31);imagesc(lp(A.Data_new));
    hold on;plot(bt.idx,'--');
    figure(33);plot(lp(bt.waveform));
    figure(32);imagesc(lp(bt.waveform));
  end
  
  
  
  %}
  %% Coherene Index Calculation with Slope error Correction
  
  
  %
  
  Nx0 = size(A.Data,2);
  Nx = floor(Nx0/Nx_int);
  Nx_mod = mod(Nx0,Nx_int);
  if Nx_mod>= Nx_int/2;
    Nx = Nx + 1;
  end
  square_int = zeros(size( A.Data_new,1),Nx); % incoherent integration, take square of abs first, then sum
  square_int_dB=zeros(size( A.Data_new,1),Nx);
  int_square = zeros(size(A.Data_new,1),Nx); % coherent integration, sum complex data first, then take square of abs
  coh_index_slopecorr =NaN*ones(1,Nx);
  Abruptiveindex=NaN*ones(1,Nx);
  Padj=NaN*ones(1,Nx);
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
   bottom_coh=NaN*ones(1,Nx);
  for rline =1:Nx
    idx1 = (rline-1)*Nx_int + 1;
    idx2 = rline*Nx_int;
    if Nx0 - idx2 > 0 & Nx0 - idx2 < Nx_int/2;
      idx2 = Nx0;
    else
      idx2 = min(idx2,Nx0);
    end
    
    square_int(:,rline) = mean(abs(A.Data_new(:,idx1:idx2)).^2,2);
    int_square(:,rline) = abs(mean(A.Data_new(:,idx1:idx2),2)).^2;
    square_int_dB(:,rline) = 10*log10(square_int(:,rline));
    
    Latitude_mean(rline)=mean(A.Latitude(idx1:idx2));  %Mean Latitude
    Longitude_mean(rline)=mean(A.Longitude(idx1:idx2)); %Mean Longitude
    GPS_time_ave(rline)=mean(A.GPS_time(idx1:idx2));   %Mean GPS Time
    
    meanbt=nanmean(A.bt_elev_new(idx1:idx2));
    if meanbt==0 | isnan(meanbt)
      continue;               %skip if no ice bottom
    end
    
    
    meansf=nanmean(A.sf_elev_new(idx1:idx2));  %If bottom=surface; skip
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
    if b2>size(A.Data_new,1)
      b2=size(A.Data_new,1);
    end
    
    [bt_val,bt_idx] = max(square_int(b1:b2,rline));  %Peak Index and Value
    bt_idx = bt_idx_m-5+bt_idx-1;
    bt_pwr = 10*log10(bt_val);
    
    if bt_pwr==0
      continue;
    end
    
    noise_bin1=bt_idx+470;
    noise_bin2=bt_idx+500;
    if noise_bin1>size(A.Data_new,1) || noise_bin2>size(A.Data_new,1)
      noise_bin1=size(A.Data_new,1)-30;
      noise_bin2=size(A.Data_new,1);
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
    if fallingedge(rline)>size(A.Data_new,1)
      fallingedge(rline)=size(A.Data_new,1);
    end
    while square_int_dB(fallingedge(rline),rline)-noise > 0.05*SNR & fallingedge(rline)<size(A.Data_new,1) %Depth Bins risingedge(rline) and fallingedge(rline)
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
    Dpth(idx1:idx2)=A.sf_elev_new(idx1:idx2)-A.bt_elev_new(idx1:idx2);
    % Power(idx1:idx2)=lp(ice_bed_power(idx1:idx2))+2*lp(2*(480+Dpth(idx1:idx2)/sqrt(3.15)));
    
    
    
    B(rline)=2.3*3000/(depth(rline)+2000);
    Gmtrc_loss(rline)=2*lp(2*((A.Elevation(rline)-A.sf_elev_new(rline))+depth(rline)/sqrt(3.15)));
    Padj(rline)=lp(Imeanx)+Gmtrc_loss(rline)+B(rline)*depth(rline)/100;
    Padj2_Na11(rline)=lp(Imeanx)+Gmtrc_loss(rline)+2*11*(depth(rline)-MeanDepth)/1000; %Using Na=11dB/km
    
    %Saving waveform around the bottom above and below
    bx1=bt_idx-50;
    bx2=bt_idx+100;
    if bx1<=0 |bx2>size(square_int_dB,1)
      continue;
    end
    bt.inc_wf_ave(:,rline) = square_int_dB(bt_idx-50:bt_idx+100,rline);
    bottom_coh(rline)=bt_idx;
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


return