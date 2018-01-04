
%% Determination of coherence Index 

tic
fc = 195e6;
fs = 1.1111e8;  %Sampling Frequency of the radar
c = 3e8;        %Speed of light
er_ice = 3.15;  %Permittivity of ice
p=4.99;   %Pre windowed radar pulse halfwidth in air     

load_flag =1;

if load_flag
     A = load('/cresis/snfs1/dataproducts/ct_data/rds/2010_Greenland_P3/manjish/CSARP_Data/20130404_02/Data_20130404_02_001.mat');
   % L = load('/cresis/snfs1/dataproducts/ct_data/rds/2013_Greenland_P3/CSARP_post/CSARP_layerData/20130404_02/Data_20130404_02_001.mat');
    load('/cresis/snfs1/scratch/manjish/peterman/radar/verticalline1.mat');
    
end
toc

distance=geodetic_to_along_track(Greenland.Latitude,Greenland.Longitude,Greenland.Elevation); %Horizontal Distance

%bttime= interp1(L.GPS_time, L.layerData{2}.value{2}.data,Greenland.GPS_time,'linear','extrap');
%surface= interp1(L.GPS_time, L.layerData{1}.value{2}.data,Greenland.GPS_time,'linear','extrap');
%Elevation=interp1(L.GPS_time,L.Elevation,Greenland.GPS_time,'linear','extrap');
bttime=Greenland.ice_bed_time;
surface=Greenland.surface_time;
Elevation=Greenland.Elevation;


figure(40),imagesc(Greenland.GPS_time,Greenland.Time,10*log10(abs(Greenland.Data)));
hold on;plot(Greenland.GPS_time,surface,'--');plot(Greenland.GPS_time,bttime,'--');
% 
bottom_twtt=bttime+(Greenland.Elevation-Elevation)/(c/2);
surface_twtt=surface+(Greenland.Elevation-Elevation)/(c/2);


filter_length=100;
bottom_twtt=sgolayfilt(bottom_twtt,2,filter_length+1,hanning(filter_length+1));
surface_twtt=sgolayfilt(surface_twtt,2,filter_length+1,hanning(filter_length+1));


figure(40),imagesc(Greenland.GPS_time,Greenland.Time,10*log10(abs(Greenland.Data)));
hold on;plot(Greenland.GPS_time,surface_twtt,'--');plot(Greenland.GPS_time,bottom_twtt,'--');
% 
dt= Greenland.Time(2)-Greenland.Time(1);
index = round((surface_twtt-Greenland.Time(1))/dt);
            ice_surface_power  = zeros(1,length(Greenland.Surface));
            
           
            for i = 1:length(surface_twtt)
                  if isnan(surface_twtt(i)) || isinf(surface_twtt(i))
                    ice_surface_power(i)= nan;
                    continue
                  else
                    
                      
                      
                [surface_power idx] = max(Greenland.Data(index(i)-1:index(i)+1,i));
                surface_index = idx + index(i)-1-1;
                ice_surface_power(i) = surface_power;
                surface_twtt(i) =  interp1([1:length(Greenland.Time)],Greenland.Time,surface_index);
              
                  end
            end
            
            
            %           dt= datGreenland.Time(2)-datGreenland.Time(1);
            index = round((bottom_twtt-Greenland.Time(1))/dt);
            ice_bed_power  = zeros(1,length(Greenland.Surface));
            
            for i = 1:length(bottom_twtt)
                if isnan(bottom_twtt(i)) || isinf(bottom_twtt(i))
                    ice_bed_power(i)= nan;
                    continue
                else
                    [bed_power idx] = max(Greenland.Data(index(i)-5:index(i)+5,i));
                    bed_index = idx + index(i)-15-1;
                    %                     ice_bed_echos(:,i)=
                    %                     datGreenland.Data(index(i)-15:index(i)+15,i);
                    N = mean(sqrt(Greenland.Data(bed_index+200:bed_index+500,i)));
                    SNR = 10*log10((sqrt(bed_power))/N);
                    if SNR > 3
                        ice_bed_power(i) = bed_power;
                        bottom_twtt(i) =  interp1([1:length(Greenland.Time)],Greenland.Time,bed_index);
                    else
                        %                         keyboard
                        ice_bed_power(i) = nan;
                        bottom_twtt(i) = nan;
                        continue ;
                    end
                end
            end


figure(1),imagesc(Greenland.GPS_time,Greenland.Time,10*log10(abs(Greenland.Data)));
hold on;plot(Greenland.GPS_time,surface_twtt,'--');plot(Greenland.GPS_time,bottom_twtt,'--');


%% Coherene Index Calculation with original Datasf(sf==Inf)=NaN;

idx=find(isinf(bottom_twtt));
bottom_twtt(idx)=nan;

MeanDepth=(nanmean(bottom_twtt)-nanmean(surface_twtt))*c/(2*sqrt(er_ice));  %Original Mean Depth
MeanDepth_og=(nanmean(bttime)-nanmean(Greenland.Surface))*c/(2*sqrt(er_ice)); 
 %Mean Ice Depth 

if MeanDepth_og~=0 & MeanDepth_og>150;
    Nx_int_dist=2*sqrt(MeanDepth*p/sqrt(er_ice));                    %Incoherent Integration distance
 
    Nx_int = floor(Nx_int_dist/(distance(10)-distance(9)));                %No of lines integrated
    Nx0 = size(Greenland.Data,2);
    Nx = floor(Nx0/Nx_int);                                         %Total lines after integration
    Nx_mod = mod(Nx0,Nx_int);
    if Nx_mod>= Nx_int/2;
        Nx = Nx + 1;
    end
    square_int = zeros(size(Greenland.Data,1),Nx); % incoherent integration, take square of abs first, then sum (phase info lost)
    int_square = zeros(size(Greenland.Data,1),Nx); % coherent integration, sum complex data first, then take square of abs (phase info remains)
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
        square_int(:,rline) = mean(abs(Greenland.Data(:,idx1:idx2)).^2,2);
        square_int_dB = 10*log10(square_int(:,rline));
        int_square(:,rline) = (abs(mean(Greenland.Data(:,idx1:idx2),2))).^2;
        
        meanbt=nansum(bottom_twtt(isfinite(bottom_twtt(idx1:idx2))));
        count=nansum(isfinite(bottom_twtt(idx1:idx2)));
        if meanbt==0 || count==0
            continue;
        end
        meanbt=meanbt/count;
        bt_idx_m = find(Greenland.Time>meanbt,1,'first');
        
        if isempty(bt_idx_m)
            bt_idx_m=round(size(Greenland.Data,1)/2);
        end
        b1=bt_idx_m-5;
        b2=bt_idx_m+5;
        if b2>size(Greenland.Data,1)
            b2=size(Greenland.Data,1);
        end
        
        [bt_val,bt_idx] = max(square_int(b1:b2,rline));
        bt_idx = bt_idx_m-5+bt_idx-1;                 %Peak Index
        bt_pwr = 10*log10(bt_val);                     %Peak Value
        square_int_dB = 10*log10(square_int(:,rline));
        noise_bin1=bt_idx+70;
        noise_bin2=bt_idx+100;
        if noise_bin1>size(Greenland.Data,1)
            noise_bin1=size(Greenland.Data,1)-30;
            noise_bin2=size(Greenland.Data,1);
        end
        
        if noise_bin2>size(Greenland.Data,1)
            noise_bin1=size(Greenland.Data,1)-30;
            noise_bin2=size(Greenland.Data,1);
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
        if fallingedge(rline)>size(Greenland.Data,1)
            fallingedge(rline)=size(Greenland.Data,1) ;
        end
        while square_int_dB(fallingedge(rline))-noise > 0.05*SNR &&  fallingedge(rline)<size(Greenland.Data,1)      %Depth Bins risingedge(rline) and fallingedge(rline)
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
    figure(4);plot(distancenx,coh_index_ogdata,'Displayname','Original Data');
    title('Original Data Coherence Index');   %Coherence Index without Correction
    figure(40);plot(distancenx,Abruptiveindex);
   
    %% Correction for Greenland.Elevation
 
    sf=surface_twtt; %test
    sf_new = zeros(size(sf));
    bt_new = zeros(size(bottom_twtt));
    Greenland.Elevation_new = zeros(size(Greenland.Elevation));
    Greenland.sf_elev_new = zeros(size(Greenland.Elevation));
    Greenland.bt_elev_new = zeros(size(Greenland.Elevation));
    
    % 1)Remove data before zero time
    negative_bins = Greenland.Time < 0;
    Greenland.Time_new = Greenland.Time(~negative_bins);
    Greenland.Data_new = Greenland.Data(~negative_bins,:);
    
    % 2)Create Greenland.Elevation axis to interpolate to
    [max_elev,max_elev_idx] = max(Greenland.Elevation);
    min_elev = min(Greenland.Elevation - sf*c/2 - (Greenland.Time_new(end)-sf)*c/2/sqrt(er_ice));
    dt = Greenland.Time(2)-Greenland.Time(1);
    dr = dt * c/2 / sqrt(er_ice);
    dt_air = dr/(c/2);
    elev_axis = max_elev:-dr:min_elev;
    new_time = zeros(length(elev_axis),length(Greenland.Elevation));
    
    % 3)Zero pad data to create space for interpolated data
    zero_pad_len = length(elev_axis) - length(Greenland.Time_new);
    Greenland.Data_new = cat(1,Greenland.Data_new,zeros(zero_pad_len,size(Greenland.Data_new,2)));
    
    % 4)Determine the corrections to apply to Greenland.Elevation and layers
    dRange = max_elev - Greenland.Elevation;
    dBins = round(dRange / (c/2) / dt);
    dtime = dRange/(c/2);
    
    for rline = 1:size(Greenland.Data_new,2)
        % Determine Greenland.Elevation bins before surface
        sf_elev = Greenland.Elevation(rline) - sf(rline) * c/2;
        time0 = -(max_elev - Greenland.Elevation(rline))/(c/2);
        last_air_idx = find(elev_axis > sf_elev,1,'last');
        new_time_tmp = (time0 + dt_air*(0:last_air_idx-1)).';
        if last_air_idx < length(elev_axis)
            % Determine Greenland.Elevation bins after surface
            dt_ice = dr/(c/2/sqrt(er_ice));
            first_ice_idx = last_air_idx + 1;
            time0 = sf(rline) + (sf_elev - elev_axis(first_ice_idx))/(c/2/sqrt(er_ice));
            new_time(:,rline) = cat(1,new_time_tmp, (time0 + dt_ice*(0:length(elev_axis)-length(new_time_tmp)-1)).');
        end
        Greenland.Data_new(:,rline) = interp1(Greenland.Time_new, Greenland.Data_new(1:length(Greenland.Time_new),rline), new_time(:,rline), 'linear',0);
        Greenland.Elevation_new(rline) = Greenland.Elevation(rline) + dRange(rline);
        sf_new(rline) = sf(rline) + dtime(rline);
        bt_new(rline) = bottom_twtt(rline) + dtime(rline);
        Greenland.sf_elev_new(rline) = Greenland.Elevation_new(rline) - sf_new(rline)*c/2;
        Greenland.bt_elev_new(rline) = Greenland.sf_elev_new(rline) - (bt_new(rline)-sf_new(rline))*c/2/sqrt(er_ice);
          
       
        
    end
    
    fh = figure(2);
    figure(fh);imagesc([],elev_axis,10*log10(abs(Greenland.Data_new).^2));title('Greenland.Elevation Correction Data')
    ax = gca;
    ax.YDir = 'normal';
    hold on;plot(Greenland.sf_elev_new,'--');plot(Greenland.bt_elev_new,'--');
    bt_slope = diff(Greenland.bt_elev_new)./diff(distance);
   

    
    
    %% Coherene Index Calculation with Greenland.Elevation Correction
    
    square_int = zeros(size( Greenland.Data_new,1),Nx); % incoherent integration, take square of abs first, then sum
    int_square = zeros(size(Greenland.Data_new,1),Nx); % coherent integration, sum complex data first, then take square of abs
    coh_index_elevcorr = zeros(1,Nx);
    for rline =1:Nx
        idx1 = (rline-1)*Nx_int + 1;
        idx2 = rline*Nx_int;
        if Nx0 - idx2 > 0 & Nx0 - idx2 < Nx_int/2;
            idx2 = Nx0;
        else
            idx2 = min(idx2,Nx0);
        end
        square_int(:,rline) = mean(abs(Greenland.Data_new(:,idx1:idx2)).^2,2);
        int_square(:,rline) = abs(mean(Greenland.Data_new(:,idx1:idx2),2)).^2;
        square_int_dB=lp(square_int(:,rline));
        
        meanbt=nansum(Greenland.bt_elev_new(idx1:idx2));
        count=nansum(isfinite(Greenland.bt_elev_new(idx1:idx2)));
        meanbt=meanbt/count;
        if meanbt==0 || count==0 ||isnan(meanbt)
            continue;
        end
        
        bt_idx_m = find(elev_axis<=meanbt,1,'first');
        
        if isempty(bt_idx_m)
            bt_idx_m=round(size(Greenland.Data_new,1)/2);      %looking for ice bottom
        end
        b1=bt_idx_m-5;
        b2=bt_idx_m+5;
        if b2>size(Greenland.Data_new,1)
            b2=size(Greenland.Data_new,1);
        end
        
        [bt_val,bt_idx]=max(square_int(b1:b2,rline));
        bt_idx = bt_idx_m-5+bt_idx-1;
        bt_pwr = 10*log10(bt_val);
        noise_bin1=bt_idx+70;
        noise_bin2=bt_idx+100;
        if noise_bin1>size(Greenland.Data_new,1)
            noise_bin1=size(Greenland.Data_new,1)-30;
            noise_bin2=size(Greenland.Data_new,1);
        end
        if noise_bin2>size(Greenland.Data_new,1)
            noise_bin2=size(Greenland.Data_new,1);
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
        if fallingedge(rline)>size(Greenland.Data_new,1)
            fallingedge(rline)=size(Greenland.Data_new,1) ;
        end
        while square_int_dB(fallingedge(rline))-noise > 0.05*SNR && fallingedge(rline)<size(Greenland.Data_new,1)      %Depth Bins risingedge(rline) and fallingedge(rline)
            fallingedge(rline) = fallingedge(rline) + 1;
        end
        Imeanx=sum((square_int(risingedge(rline):fallingedge(rline),rline)));
        Abruptiveindex(rline)=bt_val/Imeanx;
        coh_index_elevcorr(rline) = sum(int_square(risingedge(rline):fallingedge(rline),rline))/sum(square_int(risingedge(rline):fallingedge(rline),rline));
    end
    hold on; figure(4); hold on; plot(distancenx,coh_index_elevcorr,'r','Displayname','Greenland.Elevation Corrected');
    figure(41); hold on; plot(distancenx,Abruptiveindex,'r','Displayname','Greenland.Elevation Corrected'); grid;
    
    %% Ice Slope Correction
    
    Nx0 = size(Greenland.Data_new,2);
    Nx = floor(Nx0/Nx_int);
    Nx_mod = mod(Nx0,Nx_int);
    if Nx_mod>= Nx_int/2;
        Nx = Nx + 1;
    end
    slopeerror=zeros(1,size(Greenland.Data_new,1));
    slopeval=zeros(1,size(Greenland.Data_new,1));
    angle=zeros(1,Nx);
    
    for rline =1:Nx
        
        idx1 = (rline-1)*Nx_int + 1;
        idx2 = rline*Nx_int;
        if Nx0 - idx2 > 0 && Nx0 - idx2 < Nx_int/2;
            idx2 = Nx0;
        else
            idx2 = min(idx2,Nx0);
        end
        
        if isinf(mean(Greenland.bt_elev_new(idx1:idx2)))  % If No Ice Bottom, skip
            continue;
        end
        if isnan(mean(Greenland.bt_elev_new(idx1:idx2)))  % If No Ice Bottom, skip
            continue;
        end
        
        p = polyfit(distance(idx1:idx2),Greenland.bt_elev_new(idx1:idx2),1);    % Polygonal fitting
        slopeval(idx1:idx2) = polyval(p,distance(idx1:idx2));           % Ploygonal values after fitting
        base = distance(idx2)-distance(idx1);
        perpendicular = slopeval(idx2)-slopeval(idx1);
        hypotenuse = sqrt(base^2+perpendicular^2);
        angle(rline)=asin(perpendicular/hypotenuse)*180/pi;
        slopeerror = slopeval(idx1:idx2)-slopeval(idx1);         % Error betn original and fitting line
        dtime = 2*slopeerror/c/sqrt(er_ice);
        if rline==1
            figure(8);plot([idx1:idx2],Greenland.bt_elev_new(idx1:idx2));
            hold on;plot([idx1:idx2],slopeval(idx1:idx2),'r--')
            figure(9);plot(10*log10(abs(Greenland.Data_new(:,idx2)).^2))
        end
        for idx = idx1:idx2
            Greenland.Data_new(:,idx) = interp1(elev_axis, Greenland.Data_new(:,idx), elev_axis + slopeerror(idx-idx1 +1), 'linear',0);
            Greenland.Elevation_new(idx) = Greenland.Elevation_new(idx) - slopeerror(idx-idx1 +1);
            Greenland.sf_elev_new(idx) = Greenland.sf_elev_new(idx) - slopeerror(idx-idx1 +1);
            Greenland.bt_elev_new(idx) = Greenland.bt_elev_new(idx) - slopeerror(idx-idx1 +1);
            sf_new(idx) = sf_new(idx) + dtime(idx-idx1 +1);
            bt_new(idx) = bt_new(idx) + dtime(idx-idx1 +1);
        end
        if rline==1
            figure(8);hold on;plot([idx1:idx2],Greenland.bt_elev_new(idx1:idx2),'g--');
            figure(9);hold on;plot(10*log10(abs(Greenland.Data_new(:,idx2)).^2),'g--');
        end
    end
    
    
    fh = figure(3);
    %figure(fh);imagesc([1:size(Greenland.Data_new,2)],elev_axis,10*log10(abs(Greenland.Data_new).^2));title('Slope Correction Data')
    figure(fh);imagesc([],elev_axis,10*log10(abs(Greenland.Data_new).^2));title('Slope Correction Data')
    ax = gca;
    ax.YDir = 'normal';
    hold on;plot(Greenland.sf_elev_new,'--');plot(Greenland.bt_elev_new,'--');
    
    
    
    %
     % Truncate data around ice bottom within bt.range_bins
      bt_range_bins =[-50:100];
      bt.val = NaN*ones(1,size(Greenland.Data_new,2));
      bt.idx = NaN*ones(1,size(Greenland.Data_new,2));
      bt.waveform = NaN*ones(length(bt_range_bins),size(Greenland.Data_new,2));
      bt.inc_wf_ave=NaN*ones(size(bt.waveform,1),Nx);
      for rline = 1:size(Greenland.Data_new,2)
        if ~isnan(Greenland.bt_elev_new(rline)) & ~isinf(Greenland.bt_elev_new(rline))
         % bt.idx(rline) = find(elev_axis<=Greenland.bt_elev_new(rline),1,'first');
         bt.idx(rline) = round(interp1(elev_axis,[1:length(elev_axis)],Greenland.bt_elev_new(rline)));
         bt.val(rline) = Greenland.Data_new(bt.idx(rline),rline);
          first_idx = bt.idx(rline) + bt_range_bins(1);
          last_idx = bt.idx(rline) + bt_range_bins(end);
          if first_idx < 1 | last_idx>size(Greenland.Data_new,1)
            bt.idx(rline) = NaN;
            bt.val(rline) = NaN;
            continue
          end
          lower=bt.idx(rline)+bt_range_bins(1);
          upper=bt.idx(rline)+bt_range_bins(end);
          if upper >size(Greenland.Data_new,1)
              upper=size(Greenland.Data_new,1);
          end
          bt.waveform(1:51+upper-bt.idx(rline),rline) = Greenland.Data_new(lower:upper, rline);
        else
          continue
        end
      end
    
      debug_flag=0;
     if debug_flag==1
        figure(31);imagesc(lp(Greenland.Data_new));
        hold on;plot(bt.idx,'--');
        figure(33);plot(lp(bt.waveform));
        figure(32);imagesc(lp(bt.waveform));
     end
    
    
    
    %}
    %% Coherene Index Calculation with Slope error Correction
    

%     

    Nx0 = size(Greenland.Data,2);
    Nx = floor(Nx0/Nx_int);
    Nx_mod = mod(Nx0,Nx_int);
    if Nx_mod>= Nx_int/2;
        Nx = Nx + 1;
    end
    square_int = zeros(size( Greenland.Data_new,1),Nx); % incoherent integration, take square of abs first, then sum
    square_int_dB=zeros(size( Greenland.Data_new,1),Nx);
    int_square = zeros(size(Greenland.Data_new,1),Nx); % coherent integration, sum complex data first, then take square of abs
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
    for rline =1:Nx
        idx1 = (rline-1)*Nx_int + 1;
        idx2 = rline*Nx_int;
        if Nx0 - idx2 > 0 & Nx0 - idx2 < Nx_int/2;
            idx2 = Nx0;
        else
            idx2 = min(idx2,Nx0);
        end
        
        square_int(:,rline) = mean(abs(Greenland.Data_new(:,idx1:idx2)).^2,2);
        int_square(:,rline) = abs(mean(Greenland.Data_new(:,idx1:idx2),2)).^2;
        square_int_dB(:,rline) = 10*log10(square_int(:,rline));
        
        Latitude_mean(rline)=mean(Greenland.Latitude(idx1:idx2));  %Mean Latitude
        Longitude_mean(rline)=mean(Greenland.Longitude(idx1:idx2)); %Mean Longitude
        GPS_time_ave(rline)=mean(Greenland.GPS_time(idx1:idx2));   %Mean GPS Time
        
        meanbt=nanmean(Greenland.bt_elev_new(idx1:idx2));
        if meanbt==0 | isnan(meanbt)
            continue;               %skip if no ice bottom
        end
        
        
        meansf=nanmean(Greenland.sf_elev_new(idx1:idx2));  %If bottom=surface; skip
        meanbt1=nanmean(bttime(idx1:idx2));
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
        if b2>size(Greenland.Data_new,1)
            b2=size(Greenland.Data_new,1);
        end
        
        [bt_val,bt_idx] = max(square_int(b1:b2,rline));  %Peak Index and Value
        bt_idx = bt_idx_m-5+bt_idx-1;
        bt_pwr = 10*log10(bt_val);
        
        if bt_pwr==0
            continue;
        end
        
        noise_bin1=bt_idx+470;
        noise_bin2=bt_idx+500;
        if noise_bin1>size(Greenland.Data_new,1) || noise_bin2>size(Greenland.Data_new,1)
            noise_bin1=size(Greenland.Data_new,1)-30;
            noise_bin2=size(Greenland.Data_new,1);
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
        if fallingedge(rline)>size(Greenland.Data_new,1)
            fallingedge(rline)=size(Greenland.Data_new,1);
        end
        while square_int_dB(fallingedge(rline),rline)-noise > 0.05*SNR & fallingedge(rline)<size(Greenland.Data_new,1) %Depth Bins risingedge(rline) and fallingedge(rline)
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
       % depth(rline)=meansf-meanbt;
       Dpth(idx1:idx2)=Greenland.sf_elev_new(idx1:idx2)-Greenland.bt_elev_new(idx1:idx2); 
       Power(idx1:idx2)=lp(ice_bed_power(idx1:idx2))+2*lp(2*(480+Dpth(idx1:idx2)/sqrt(3.15)));
        
       
       
        B(rline)=2.3*3000/(depth(rline)+2000);
        Gmtrc_loss(rline)=2*lp(2*(480+depth(rline)/sqrt(3.15)));
        Padj(rline)=lp(Imeanx)+Gmtrc_loss(rline)+B(rline)*depth(rline)/100;
        
       
      %Saving waveform around the bottom above and below 
        bx1=bt_idx-50;
        bx2=bt_idx+100;
        if bx1<=0 |bx2>size(square_int_dB,1)
            continue;
        end
        bt.inc_wf_ave(:,rline) = square_int_dB(bt_idx-50:bt_idx+100,rline);
        
    end
    
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
    toc
    figure(4);plot(distancenx,coh_index_ogdata,'Displayname','Original Data');
    
    figure;subplot(3,1,1);plot(distancenx,Padj)
    subplot(3,1,2);plot(distancenx,coh_index_slopecorr)
    subplot(3,1,3);plot(distancenx,Abruptiveindex)
    figure;plot(depth);
end


 Nx_int=find(distancenx>1,1,'first');

  Nx0 = size(distancenx,2);
    Nx = floor(Nx0/Nx_int);
    Nx_mod = mod(Nx0,Nx_int);
    if Nx_mod>= Nx_int/2;
        Nx = Nx + 1;
    end
    

    for rline=1:Nx
           
         idx1 = (rline-1)*Nx_int + 1;
        idx2 = rline*Nx_int;
        if Nx0 - idx2 > 0 & Nx0 - idx2 < Nx_int/2;
            idx2 = Nx0;
        else
            idx2 = min(idx2,Nx0);
        end
        
        
        
         na = 0:0.5:30; 
        S=zeros(1,length(na));
        for N1= 1:length(na)
           S(N1)=sum(abs(((Power(idx1:idx2))+2*na(N1)*((Dpth(idx1:idx2)-mean(Dpth))/1e3))));
        end
        [v i] = min(S);
        if i==length(S)
             warning('check this')
        end
        Na(rline)=na(i);
      
        Ref(idx1:idx2)=Power(idx1:idx2)+ 2*Na(rline)*((Dpth(idx1:idx2)-mean(Dpth)))/1e3;
     


    end

figure;plot(Ref);
figure;plot(Na);
