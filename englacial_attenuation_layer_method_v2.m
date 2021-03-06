
close all
clearvars

debug_flag=1;
frms=[21:24]; %21:24
att_val_bott=[];
lat=[];
lon=[];
gps_time=[];
Na=[];
Count=[];
Depth=[];
for k=1:length(frms)
    
    if 1
      clear data layer
        data=load(sprintf('/cresis/snfs1/dataproducts/ct_data/rds/2011_Greenland_P3/CSARP_standard/20110429_01/Data_20110429_01_%03d.mat',frms(k)));
        layer=load(sprintf('/cresis/snfs1/dataproducts/ct_data/rds/2011_Greenland_P3/CSARP_layerData/20110429_01/Data_20110429_01_%03d.mat',frms(k)));
    end
    
    
    if 0
    figure(100); imagesc(lp(data.Data)); title('Original data')
    
%     roll=data.Roll*180/pi;
%      ridx=find(abs(roll)>5);
%     
%      data.Data(:,ridx)=nan;
    
    inc_B_filter=ones(51,1)/51;
    Incoh_data=(fir_dec(abs(data.Data).^2,inc_B_filter',700));
    figure(101);imagesc(lp(Incoh_data));title('FIR DEC ')
    end
    
    clear tmp_data;
    %Incoherent Integration
    if 0
        
        numofCohInt=700;
        
        Nx=floor(length(data.Elevation)/numofCohInt);
        %        ice_surface_power_tmp=nan*ones(1,Nx);

        %        ice_surface_power_tmp=nan*ones(1,Nx);
        %        Latitude=nan*ones(1,Nx);
        %        Longitude=nan*ones(1,Nx);
        %        Elevation=nan*ones(1,Nx);
        %        surface_twtt_tmp=nan*ones(1,Nx);
%         i=1;
        for i= 1:length(data.Elevation)-numofCohInt
             idx1=(i-1)*numofCohInt+1;
             idx2=i*numofCohInt;
            %idx1=m;
            %idx2=m+numofCohInt;
            if idx2>length(data.GPS_time) & length(data.GPS_time)-idx1>numofCohInt/2
                idx2=length(data.GPS_time);
            elseif idx2>length(data.GPS_time)& length(data.GPS_time)-idx1<numofCohInt/2
                continue;
            end
            tmp_data.GPS_time(i)=mean(data.GPS_time(idx1:idx2));
            tmp_data.Latitude(i)=mean(data.Latitude(idx1:idx2));
            tmp_data.Longitude(i)=mean(data.Longitude(idx1:idx2));
            tmp_data.Elevation(i)=mean(data.Elevation(idx1:idx2));
            tmp_data.Data(:,i)=mean((abs(data.Data(:,idx1:idx2)).^2),2);
           % i=i+1;
        end
        tmp_data.Time=data.Time;
        
        clear data
        data=tmp_data;
        
     data.Data=Incoh_data(:,1:length(data.GPS_time));
    end
    
   dist=geodetic_to_along_track(data.Latitude,data.Longitude); 
   figure;plot(diff(dist));
    
    %layer interpolation
    layer_twtt=[];
    for i=1:length(layer.layerData)
        layer_twtt(i,:)=interp1(layer.GPS_time,layer.layerData{i}.value{2}.data,data.GPS_time,'linear');
        
    end
    
    
    
    layer_power=nan*ones(length(layer.layerData),size(data.Data,2));
    index_layer=nan*ones(length(layer.layerData),size(data.Data,2));
    depth=nan*ones(length(layer.layerData),size(data.Data,2));
    
    if debug_flag
        figure(1);imagesc([],[],lp(data.Data));
    end
    
    for i=1:size(layer_twtt,1)
        %
        if ~isempty(find((~isnan(layer_twtt(i,:))), 1))
            index_layer(i,:)=round(interp1(data.Time,1:size(data.Data,1),layer_twtt(i,:)));
            
            for j=1:length(index_layer(i,:))
                if isnan(index_layer(i,j))
                    layer_power(i,j)=nan;
                    continue;
                end
                %                    if j==1553
                %                        keyboard
                %                    end
                
                %                 if j==36
                %                    keyboard
                %                 end
                %
                [ power_dB,idx] = max(lp(data.Data(index_layer(i,j)-1:index_layer(i,j)+1,j)));
                index_layer(i,j)=index_layer(i,j)-2+idx;
                layer_power(i,j)=data.Data(index_layer(i,j),j);
                
                
                %   layer_power(i,j) = abs(data.Data(index_layer(i,j),i)).^2;
                %             if isfinite(layer_twtt(1,j))
                %                keyboard;
                %             end
                depth(i,j)=(layer_twtt(i,j)-layer_twtt(1,j))*3e8/sqrt(3.2)/2;
            end
            if debug_flag
                hold on; plot(index_layer(i,:));
            end
        end
        
        % hold on; plot(layer.GPS_time,layer.layerData{i}.value{2}.data);
    end
    
    
    %Plot each layer power across track
    if debug_flag
        figure;plot(lp(layer_power(1,:)))
        for i=1:100
            
            hold on; plot(lp(layer_power(i,:)))
        end
          hold on; plot(lp(layer_power(269,:)))
            hold on; plot(lp(layer_power(270,:)))
        title('Power of layers')
    end
    
    
    
    if debug_flag
        figure;plot([1:size(data.Data,1)],lp(data.Data(:,36)));
        hold on; scatter([index_layer(:,36)],lp(layer_power(:,36)))
        
%         figure;plot([1:size(data.Data,1)],lp(data.Data(:,1418)));
%         hold on; scatter([index_layer(:,1418)],lp(layer_power(:,1418)))
        
    end
    
    %Count number of layers in each range line
    count=nan*ones(1,size(layer_power,2));
    for i=1:size(layer_power,2)
        count(i)=length(find(~isnan(layer_power(:,i))));
    end
    
    if debug_flag
        figure;plot(data.GPS_time,count)
        title('No. of internal layers in each range bin')
        figure;plot(data.GPS_time,depth)
        title('Depth of internal layers')
    end
    
    %Fitting to calculate total Ice Attenuation at the ice bottom
    
    att_val_bottom=nan*ones(1,size(layer_power,2));
    na=nan*ones(1,size(layer_power,2));
    [val,idx]=(max(count));
    
    for i=1:size(layer_power,2)-200;
        
        dpth_sc=[];
        power_sc=[];
        ice_bott_dpth=[];
         range_lines_used=0;
        for j=i:i+200
            
            if count(j)>15        %Number of layers used
             
                dpth=depth(:,j);
                nanidx=find(isnan(depth(:,j)));
                dpth(nanidx)=[];
                lpower=layer_power(:,j);
                lpower(nanidx)=[];
                power_loss=lp(lpower(:))-lp(lpower(1));
                %         power_loss=power_loss(3:end);
                %         dpth=dpth(3:end);
                [dpth,idx]=sort(dpth);
                power_loss=power_loss(idx);
                dpth_sc=cat(1,dpth_sc,dpth(1:end-1));
                power_sc=cat(1,power_sc,power_loss(1:end-1));
                ice_bott_dpth=cat(1,ice_bott_dpth,dpth(end));
                  range_lines_used=range_lines_used+1;
                %  dpth_int_layers=dpth(3:end);
                %  power_loss_int_layers=power_loss(3:end);
                
                %  p=polyfit(dpth(1:end-1),power_loss(1:end-1),3);
                %                p=polyfitZero(dpth(1:end-1),power_loss(1:end-1),3);
                %            % fval=polyval(p,dpth(1:end-1));
                %                fval=polyval(p,dpth(1:end-1));
                %             att_val_bottom(i)=polyval(p,dpth(end));
                %             na(i)=att_val_bottom(i)/(2*dpth(end)/1000);
                %             if 1
                %                 if i==36 | i==1418
                %                     figure(5);scatter(dpth(1:end-1),power_loss(1:end-1));
                %                     figure(5);hold on; plot(dpth(1:end-1),fval)
                %                     figure(5);hold on; scatter(dpth(end),att_val_bottom(i),'x');
                %                     title('Fitting for englacial Attenuation')
                %                 end
                %             end
            end
        end
        
%         if i==1212
%             keyboard
%         end
        
        if ~isempty(dpth_sc)
            [dpth_sc,idx]=sort(dpth_sc);
            power_sc=power_sc(idx);
             dpth_sc=dpth_sc(range_lines_used:end);
             power_sc=power_sc(range_lines_used:end);
            p=polyfitZero(dpth_sc,power_sc,3);
            fval=polyval(p,dpth_sc);
            att_val_bottom(i)=polyval(p,nanmean(ice_bott_dpth));
              na(i)=att_val_bottom(i)/(2*nanmean(ice_bott_dpth)/1000);
            if 1
                if i==36 | i==1418
                
                      dpth_plt=reshape(dpth_sc',[1,numel(dpth_sc)]);
                      lpower_plt=reshape(power_sc',[1,numel(power_sc)]);
                      fval_plt=reshape(fval',[1,numel(fval)]);
                          figure(5);scatter(dpth_plt,lpower_plt);
                    figure(5);hold on; plot(dpth_plt,fval_plt)
                    figure(5);hold on; scatter(nanmean(ice_bott_dpth),att_val_bottom(i),'x');
                    title('Fitting for englacial Attenuation')
                end
            end
        end
        
    end
    
    lat=cat(2,lat,data.Latitude);
    lon=cat(2,lon,data.Longitude);
    gps_time=cat(2,gps_time,data.GPS_time);
    att_val_bott=cat(2,att_val_bott,att_val_bottom);
    Na=cat(2,Na,na);
    Count=cat(2,Count,count);
    Depth=cat(2,Depth,depth(2,:));
    
end
% meanpower=nanmean(lp(layer_power),2);
% meandepth=nanmean(depth,2);
% figure;scatter(meandepth,meanpower);
%
% figure;plot(att_val_bottom)
% figure;plot(att_val_bottom-nanmean(att_val_bottom))

% figure;plot(lat,att_val_bott);
% figure;plot(lat,att_val_bott-nanmean(att_val_bott))

%att_val_bott_filt=sgolayfilt(att_val_bott,3,1001);
% figure;plot(lat,att_val_bott);
% figure;plot(lon,att_val_bott_filt);
figure;plot(lat,Na);title('Na')


load('/cresis/snfs1/scratch/manjish/new_peterman/reflectivity/crossline10.mat');

figure;plot(out.Latitude,out.const_attenuation)
nanidx=find(isnan(out.Longitude));
out.const_attenuation(nanidx)=[];
out.Longitude(nanidx)=[];

figure;plot(lat,-att_val_bott-nanmean(-att_val_bott));
hold on; plot(out.Latitude,out.const_attenuation)

%figure;plot(lat,depth)
figure;plot(lat,-att_val_bott);
figure;plot(-att_val_bott)
figure;plot(Depth); title('Depth')
% lastidx=find(out.Longitude<=lon(end),1,'first');
%
% att_val_method2=out.const_attenuation(1:lastidx);
% att_var_method2=out.var_attenuation(1:lastidx);
% LON=out.Longitude(1:lastidx);
% LAT=out.Latitude(1:lastidx);
%
% figure;plot(lon,att_val_bott-nanmean(att_val_bott));
% hold on; plot(LON,att_val_method2)
%
% figure;plot(lon,att_val_bott_filt-nanmean(att_val_bott_filt));
% hold on; plot(LON,att_val_method2)

% figure;plot(out.Longitude,out.const_attenuation);
% att_cal_m1=interp1(out.Longitude,out.const_attenuation,lon);
%
% figure;plot(lat,att_val_bott_filt-nanmean(att_val_bott_filt));
% hold on; plot(lon,att_cal_m1);

% r.rms_height(k) = 0.001;
% options = optimset('PlotFcns',@optimplotfval);
% cost_func=@(r_rms)(abs(rms_fit - exp(-(2*(2*pi/(c/param.radar.fc))*r_rms).^2)/((r_rms).^2)));
% r_rms=[0.00001,1];
% [r_rmsheight,fval]=fminsearch(cost_func,r_rms);
% r.rms_height(k)=mean(abs(r_rmsheight));
% power_loss=lp(lpower(:))-lp(lpower(1));
% [r,m,b]=regression(dpth,power_loss);

% figure;scatter(dpth(2:end),power_loss(2:end))
%%Median Filtering
% figure;
        %        Latitude=nan*ones(1,Nx);
        %        Longitude=nan*ones(1,Nx);
        %        Elevation=nan*ones(1,Nx);
        %        surface_twtt_tmp=nan*ones(1,Nx);
% %         i=1;
%         for i= 1:length(data.Elevation)-numofCohInt
%              idx1=(i-1)*numofCohInt+1;
%              idx2=i*numofCohInt;
%             %idx1=m;
%             %idx2=m+numofCohInt;
%             if idx2>length(data.GPS_time) & length(data.GPS_time)-idx1>numofCohInt/2
%                 idx2=length(data.GPS_time);
%             elseif idx2>length(data.GPS_time)& length(data.GPS_time)-idx1<numofCohInt/2
%                 continue;
%             end
%             tmp_data.GPS_time(i)=mean(data.GPS_time(idx1:idx2));
%             tmp_data.Latitude(i)=mean(data.Latitude(idx1:idx2));
%             tmp_data.Longitude(i)=mean(data.Longitude(idx1:idx2));
%             tmp_data.Elevation(i)=mean(data.Elevation(idx1:idx2));
%             tmp_data.Data(:,i)=mean((abs(data.Data(:,idx1:idx2)).^2),2);
%            % i=i+1;
%         end
%         tmp_data.Time=data.Time;
%         
%         clear data
%         data=tmp_data;
%         
%      data.Data=Incoh_data(:,1:length(data.GPS_time));
%     end
%     
%    dist=geodetic_to_along_track(data.Latitude,data.Longitude); 
%    figure;plot(diff(dist));
%     
%     %layer interpolation
%     layer_twtt=[];
%     for i=1:length(layer.layerData)
%         layer_twtt(i,:)=interp1(layer.GPS_time,layer.layerData{i}.value{2}.data,data.GPS_time,'linear');
%         
%     end
%     
%     
%     
%     layer_power=nan*ones(length(layer.layerData),size(data.Data,2));
%     index_layer=nan*ones(length(layer.layerData),size(data.Data,2));
%     depth=nan*ones(length(layer.layerData),size(data.Data,2));
%     
%     if debug_flag
%         figure(1);imagesc([],[],lp(data.Data));
%     end
%     
%     for i=1:size(layer_twtt,1)
%         %
%         if ~isempty(find((~isnan(layer_twtt(i,:))), 1))
%             index_layer(i,:)=round(interp1(data.Time,1:size(data.Data,1),layer_twtt(i,:)));
%             
%             for j=1:length(index_layer(i,:))
%                 if isnan(index_layer(i,j))
%                     layer_power(i,j)=nan;
%                     continue;
%                 end
%                 %                    if j==1553
%                 %                        keyboard
%                 %                    end
%                 
%                 %                 if j==36
%                 %                    keyboard
%                 %                 end
%                 %
%                 [ power_dB,idx] = max(lp(data.Data(index_layer(i,j)-1:index_layer(i,j)+1,j)));
%                 index_layer(i,j)=index_layer(i,j)-2+idx;
%                 layer_power(i,j)=data.Data(index_layer(i,j),j);
%                 
%                 
%                 %   layer_power(i,j) = abs(data.Data(index_layer(i,j),i)).^2;
%                 %             if isfinite(layer_twtt(1,j))
%                 %                keyboard;
%                 %             end
%                 depth(i,j)=(layer_twtt(i,j)-layer_twtt(1,j))*3e8/sqrt(3.2)/2;
%             end
%             if debug_flag
%                 hold on; plot(index_layer(i,:));
%             end
%         end
%         
%         % hold on; plot(layer.GPS_time,layer.layerData{i}.value{2}.data);
%     end
%     
%     
%     %Plot each layer power across track
% %     if debug_flag
% %         figure;plot(lp(layer_power(1,:)))
% %         for i=1:2
% %             
% %             hold on; plot(lp(layer_power(i,:)))
% %         end
% %           hold on; plot(lp(layer_power(269,:)))
% %             hold on; plot(lp(layer_power(270,:)))
% %         title('Power of layers')
% %     end
%     
%     
%     
%     if debug_flag
%         figure;plot([1:size(data.Data,1)],lp(data.Data(:,36)));
%         hold on; scatter([index_layer(:,36)],lp(layer_power(:,36)))
%         
% %         figure;plot([1:size(data.Data,1)],lp(data.Data(:,1418)));
% %         hold on; scatter([index_layer(:,1418)],lp(layer_power(:,1418)))
%         
%     end
%     
%     %Count number of layers in each range line
%     count=nan*ones(1,size(layer_power,2));
%     for i=1:size(layer_power,2)
%         count(i)=length(find(~isnan(layer_power(:,i))));
%     end
%     
%     if debug_flag
%         figure;plot(data.GPS_time,count)
%         title('No. of internal layers in each range bin')
%         figure;plot(data.GPS_time,depth)
%         title('Depth of internal layers')
%     end
%     
%     %Fitting to calculate total Ice Attenuation at the ice bottom
%     
%     att_val_bottom=nan*ones(1,size(layer_power,2));
%     na=nan*ones(1,size(layer_power,2));
% close all
% clearvars
% 
% debug_flag=1;
% frms=[24]; %21:24
% att_val_bott=[];
% lat=[];
% lon=[];
% gps_time=[];
% Na=[];
% Count=[];
% Depth=[];
% for k=1:length(frms)
%     
%     if 1
%         data=load(sprintf('/cresis/snfs1/dataproducts/ct_data/rds/2011_Greenland_P3/CSARP_standard/20110429_01/Data_20110429_01_%03d.mat',frms(k)));
%         layer=load(sprintf('/cresis/snfs1/dataproducts/ct_data/rds/2011_Greenland_P3/CSARP_layerData/20110429_01/Data_20110429_01_%03d.mat',frms(k)));
%     end
%     
%     
%     if 0
%     figure(100); imagesc(lp(data.Data)); title('Original data')
%     
% %     roll=data.Roll*180/pi;
% %      ridx=find(abs(roll)>5);
% %     
% %      data.Data(:,ridx)=nan;
%     
%     inc_B_filter=ones(51,1)/51;
%     Incoh_data=(fir_dec(abs(data.Data).^2,inc_B_filter',700));
%     figure(101);imagesc(lp(Incoh_data));title('FIR DEC ')
%     end
%     
%     clear tmp_data;
%     %Incoherent Integration
%     if 0
%         
%         numofCohInt=700;
%         
%         Nx=floor(length(data.Elevation)/numofCohInt);
%         %        ice_surface_power_tmp=nan*ones(1,Nx);
%         %        Latitude=nan*ones(1,Nx);
%         %        Longitude=nan*ones(1,Nx);
%         %        Elevation=nan*ones(1,Nx);
%         %        surface_twtt_tmp=nan*ones(1,Nx);
% %         i=1;
%         for i= 1:length(data.Elevation)-numofCohInt
%              idx1=(i-1)*numofCohInt+1;
%              idx2=i*numofCohInt;
%             %idx1=m;
%             %idx2=m+numofCohInt;
%             if idx2>length(data.GPS_time) & length(data.GPS_time)-idx1>numofCohInt/2
%                 idx2=length(data.GPS_time);
%             elseif idx2>length(data.GPS_time)& length(data.GPS_time)-idx1<numofCohInt/2
%                 continue;
%             end
%             tmp_data.GPS_time(i)=mean(data.GPS_time(idx1:idx2));
%             tmp_data.Latitude(i)=mean(data.Latitude(idx1:idx2));
%             tmp_data.Longitude(i)=mean(data.Longitude(idx1:idx2));
%             tmp_data.Elevation(i)=mean(data.Elevation(idx1:idx2));
%             tmp_data.Data(:,i)=mean((abs(data.Data(:,idx1:idx2)).^2),2);
%            % i=i+1;
%         end
%         tmp_data.Time=data.Time;
%         
%         clear data
%         data=tmp_data;
%         
%      data.Data=Incoh_data(:,1:length(data.GPS_time));
%     end
%     
%    dist=geodetic_to_along_track(data.Latitude,data.Longitude); 
%    figure;plot(diff(dist));
%     
%     %layer interpolation
%     layer_twtt=[];
%     for i=1:length(layer.layerData)
%         layer_twtt(i,:)=interp1(layer.GPS_time,layer.layerData{i}.value{2}.data,data.GPS_time,'linear');
%         
%     end
%     
%     
%     
%     layer_power=nan*ones(length(layer.layerData),size(data.Data,2));
%     index_layer=nan*ones(length(layer.layerData),size(data.Data,2));
%     depth=nan*ones(length(layer.layerData),size(data.Data,2));
%     
%     if debug_flag
%         figure(1);imagesc([],[],lp(data.Data));
%     end
%     
%     for i=1:size(layer_twtt,1)
%         %
%         if ~isempty(find((~isnan(layer_twtt(i,:))), 1))
%             index_layer(i,:)=round(interp1(data.Time,1:size(data.Data,1),layer_twtt(i,:)));
%             
%             for j=1:length(index_layer(i,:))
%                 if isnan(index_layer(i,j))
%                     layer_power(i,j)=nan;
%                     continue;
%                 end
%                 %                    if j==1553
%                 %                        keyboard
%                 %                    end
%                 
%                 %                 if j==36
%                 %                    keyboard
%                 %                 end
%                 %
%                 [ power_dB,idx] = max(lp(data.Data(index_layer(i,j)-1:index_layer(i,j)+1,j)));
%                 index_layer(i,j)=index_layer(i,j)-2+idx;
%                 layer_power(i,j)=data.Data(index_layer(i,j),j);
%                 
%                 
%                 %   layer_power(i,j) = abs(data.Data(index_layer(i,j),i)).^2;
%                 %             if isfinite(layer_twtt(1,j))
%                 %                keyboard;
%                 %             end
%                 depth(i,j)=(layer_twtt(i,j)-layer_twtt(1,j))*3e8/sqrt(3.2)/2;
%             end
%             if debug_flag
%                 hold on; plot(index_layer(i,:));
%             end
%         end
%         
%         % hold on; plot(layer.GPS_time,layer.layerData{i}.value{2}.data);
%     end
%     
%     
%     %Plot each layer power across track
%     if debug_flag
%         figure;plot(lp(layer_power(1,:)))
%         for i=1:100
%             
%             hold on; plot(lp(layer_power(i,:)))
%         end
%           hold on; plot(lp(layer_power(269,:)))
%             hold on; plot(lp(layer_power(270,:)))
%         title('Power of layers')
%     end
%     
%     
%     
%     if debug_flag
%         figure;plot([1:size(data.Data,1)],lp(data.Data(:,36)));
%         hold on; scatter([index_layer(:,36)],lp(layer_power(:,36)))
%         
% %         figure;plot([1:size(data.Data,1)],lp(data.Data(:,1418)));
% %         hold on; scatter([index_layer(:,1418)],lp(layer_power(:,1418)))
%         
%     end
%     
%     %Count number of layers in each range line
%     count=nan*ones(1,size(layer_power,2));
%     for i=1:size(layer_power,2)
%         count(i)=length(find(~isnan(layer_power(:,i))));
%     end
%     
%     if debug_flag
%         figure;plot(data.GPS_time,count)
%         title('No. of internal layers in each range bin')
%         figure;plot(data.GPS_time,depth)
%         title('Depth of internal layers')
%     end
%     
%     %Fitting to calculate total Ice Attenuation at the ice bottom
%     
%     att_val_bottom=nan*ones(1,size(layer_power,2));
%     na=nan*ones(1,size(layer_power,2));
%     [val,idx]=(max(count));
%     
%     for i=1:size(layer_power,2)-200;
%         
%         dpth_sc=[];
%         power_sc=[];
%         ice_bott_dpth=[];
%          range_lines_used=0;
%         for j=i:i+200
%             
%             if count(j)>15        %Number of layers used
%              
%                 dpth=depth(:,j);
%                 nanidx=find(isnan(depth(:,j)));
%                 dpth(nanidx)=[];
%                 lpower=layer_power(:,j);
%                 lpower(nanidx)=[];
%                 power_loss=lp(lpower(:))-lp(lpower(1));
%                 %         power_loss=power_loss(3:end);
%                 %         dpth=dpth(3:end);
%                 [dpth,idx]=sort(dpth);
%                 power_loss=power_loss(idx);
%                 dpth_sc=cat(1,dpth_sc,dpth(1:end-1));
%                 power_sc=cat(1,power_sc,power_loss(1:end-1));
%                 ice_bott_dpth=cat(1,ice_bott_dpth,dpth(end));
%                   range_lines_used=range_lines_used+1;
%                 %  dpth_int_layers=dpth(3:end);
%                 %  power_loss_int_layers=power_loss(3:end);
%                 
%                 %  p=polyfit(dpth(1:end-1),power_loss(1:end-1),3);
%                 %                p=polyfitZero(dpth(1:end-1),power_loss(1:end-1),3);
%                 %            % fval=polyval(p,dpth(1:end-1));
%                 %                fval=polyval(p,dpth(1:end-1));
%                 %             att_val_bottom(i)=polyval(p,dpth(end));
%                 %             na(i)=att_val_bottom(i)/(2*dpth(end)/1000);
%                 %             if 1
%                 %                 if i==36 | i==1418
%                 %                     figure(5);scatter(dpth(1:end-1),power_loss(1:end-1));
%                 %                     figure(5);hold on; plot(dpth(1:end-1),fval)
%                 %                     figure(5);hold on; scatter(dpth(end),att_val_bottom(i),'x');
%                 %                     title('Fitting for englacial Attenuation')
%                 %                 end
%                 %             end
%             end
%         end
%         
% %         if i==1212
% %             keyboard
% %         end
%         
%         if ~isempty(dpth_sc)
%             [dpth_sc,idx]=sort(dpth_sc);
%             power_sc=power_sc(idx);
%              dpth_sc=dpth_sc(range_lines_used:end);
%              power_sc=power_sc(range_lines_used:end);
%             p=polyfitZero(dpth_sc,power_sc,3);
%             fval=polyval(p,dpth_sc);
%             att_val_bottom(i)=polyval(p,nanmean(ice_bott_dpth));
%               na(i)=att_val_bottom(i)/(2*nanmean(ice_bott_dpth)/1000);
%             if 1
%                 if i==36 | i==1418
%                 
%                       dpth_plt=reshape(dpth_sc',[1,numel(dpth_sc)]);
%                       lpower_plt=reshape(power_sc',[1,numel(power_sc)]);
%                       fval_plt=reshape(fval',[1,numel(fval)]);
%                           figure(5);scatter(dpth_plt,lpower_plt);
%                     figure(5);hold on; plot(dpth_plt,fval_plt)
%                     figure(5);hold on; scatter(nanmean(ice_bott_dpth),att_val_bottom(i),'x');
%                     title('Fitting for englacial Attenuation')
%                 end
%             end
%         end
%         
%     end
%     
%     lat=cat(2,lat,data.Latitude);
%     lon=cat(2,lon,data.Longitude);
%     gps_time=cat(2,gps_time,data.GPS_time);
%     att_val_bott=cat(2,att_val_bott,att_val_bottom);
%     Na=cat(2,Na,na);
%     Count=cat(2,Count,count);
%     Depth=cat(2,Depth,depth(2,:));
%     
% end
% % meanpower=nanmean(lp(layer_power),2);
% % meandepth=nanmean(depth,2);
% % figure;scatter(meandepth,meanpower);
% %
% % figure;plot(att_val_bottom)
% % figure;plot(att_val_bottom-nanmean(att_val_bottom))
% 
% % figure;plot(lat,att_val_bott);
% % figure;plot(lat,att_val_bott-nanmean(att_val_bott))
% 
% %att_val_bott_filt=sgolayfilt(att_val_bott,3,1001);
% % figure;plot(lat,att_val_bott);
% % figure;plot(lon,att_val_bott_filt);
% figure;plot(lat,Na);title('Na')
% 
% 
% load('/cresis/snfs1/scratch/manjish/new_peterman/reflectivity/crossline10.mat');
% 
% figure;plot(out.Latitude,out.const_attenuation)
% nanidx=find(isnan(out.Longitude));
% out.const_attenuation(nanidx)=[];
% out.Longitude(nanidx)=[];
% 
% figure;plot(lat,-att_val_bott-nanmean(-att_val_bott));
% hold on; plot(out.Latitude,out.const_attenuation)
% 
% %figure;plot(lat,depth)
% figure;plot(lat,-att_val_bott);
% figure;plot(-att_val_bott)
% figure;plot(Depth); title('Depth')
% % lastidx=find(out.Longitude<=lon(end),1,'first');
% %
% % att_val_method2=out.const_attenuation(1:lastidx);
% % att_var_method2=out.var_attenuation(1:lastidx);
% % LON=out.Longitude(1:lastidx);
% % LAT=out.Latitude(1:lastidx);
% %
% % figure;plot(lon,att_val_bott-nanmean(att_val_bott));
% % hold on; plot(LON,att_val_method2)
% %
% % figure;plot(lon,att_val_bott_filt-nanmean(att_val_bott_filt));
% % hold on; plot(LON,att_val_method2)
% 
% % figure;plot(out.Longitude,out.const_attenuation);
% % att_cal_m1=interp1(out.Longitude,out.const_attenuation,lon);
% %
% % figure;plot(lat,att_val_bott_filt-nanmean(att_val_bott_filt));
% % hold on; plot(lon,att_cal_m1);
% 
% % r.rms_height(k) = 0.001;
% % options = optimset('PlotFcns',@optimplotfval);
% % cost_func=@(r_rms)(abs(rms_fit - exp(-(2*(2*pi/(c/param.radar.fc))*r_rms).^2)/((r_rms).^2)));
% % r_rms=[0.00001,1];
% % [r_rmsheight,fval]=fminsearch(cost_func,r_rms);
% % r.rms_height(k)=mean(abs(r_rmsheight));
% % power_loss=lp(lpower(:))-lp(lpower(1));
% % [r,m,b]=regression(dpth,power_loss);
% 
% % figure;scatter(dpth(2:end),power_loss(2:end))
% %%Median Filtering
% % figure;
%     [val,idx]=(max(count));
%     
%     for i=1:size(layer_power,2)-200;
%         
%         dpth_sc=[];
%         power_sc=[];
%         ice_bott_dpth=[];
%          range_lines_used=0;
%         for j=i:i+200
%             
%             if count(j)>15        %Number of layers used
%              
%                 dpth=depth(:,j);
%                 nanidx=find(isnan(depth(:,j)));
%                 dpth(nanidx)=[];
%                 lpower=layer_power(:,j);
%                 lpower(nanidx)=[];
%                 power_loss=lp(lpower(:))-lp(lpower(1));
%                 %         power_loss=power_loss(3:end);
%                 %         dpth=dpth(3:end);
%                 [dpth,idx]=sort(dpth);
%                 power_loss=power_loss(idx);
%                 dpth_sc=cat(1,dpth_sc,dpth(1:end-1));
%                 power_sc=cat(1,power_sc,power_loss(1:end-1));
%                 ice_bott_dpth=cat(1,ice_bott_dpth,dpth(end));
%                   range_lines_used=range_lines_used+1;
%                 %  dpth_int_layers=dpth(3:end);
%                 %  power_loss_int_layers=power_loss(3:end);
%                 
%                 %  p=polyfit(dpth(1:end-1),power_loss(1:end-1),3);
%                 %                p=polyfitZero(dpth(1:end-1),power_loss(1:end-1),3);
%                 %            % fval=polyval(p,dpth(1:end-1));
%                 %                fval=polyval(p,dpth(1:end-1));
%                 %             att_val_bottom(i)=polyval(p,dpth(end));
%                 %             na(i)=att_val_bottom(i)/(2*dpth(end)/1000);
%                 %             if 1
%                 %                 if i==36 | i==1418
%                 %                     figure(5);scatter(dpth(1:end-1),power_loss(1:end-1));
%                 %                     figure(5);hold on; plot(dpth(1:end-1),fval)
%                 %                     figure(5);hold on; scatter(dpth(end),att_val_bottom(i),'x');
%                 %                     title('Fitting for englacial Attenuation')
%                 %                 end
%                 %             end
%             end
%         end
%         
% %         if i==1212
% %             keyboard
% %         end
%         
%         if ~isempty(dpth_sc)
%             [dpth_sc,idx]=sort(dpth_sc);
%             power_sc=power_sc(idx);
%              dpth_sc=dpth_sc(range_lines_used:end);
%              power_sc=power_sc(range_lines_used:end);
%             p=polyfitZero(dpth_sc,power_sc,3);
%             fval=polyval(p,dpth_sc);
%             att_val_bottom(i)=polyval(p,nanmean(ice_bott_dpth));
%               na(i)=att_val_bottom(i)/(2*nanmean(ice_bott_dpth)/1000);
%             if 1
%                 if i==36 | i==1418
%                 
%                       dpth_plt=reshape(dpth_sc',[1,numel(dpth_sc)]);
%                       lpower_plt=reshape(power_sc',[1,numel(power_sc)]);
%                       fval_plt=reshape(fval',[1,numel(fval)]);
%                           figure(5);scatter(dpth_plt,lpower_plt);
%                     figure(5);hold on; plot(dpth_plt,fval_plt)
%                     figure(5);hold on; scatter(nanmean(ice_bott_dpth),att_val_bottom(i),'x');
%                     title('Fitting for englacial Attenuation')
%                 end
%             end
%         end
%         
%     end
%     
%     lat=cat(2,lat,data.Latitude);
%     lon=cat(2,lon,data.Longitude);
%     gps_time=cat(2,gps_time,data.GPS_time);
%     att_val_bott=cat(2,att_val_bott,att_val_bottom);
%     Na=cat(2,Na,na);
%     Count=cat(2,Count,count);
%     Depth=cat(2,Depth,depth(2,:));
%     
% end
% % meanpower=nanmean(lp(layer_power),2);
% % meandepth=nanmean(depth,2);
% % figure;scatter(meandepth,meanpower);
% %
% % figure;plot(att_val_bottom)
% % figure;plot(att_val_bottom-nanmean(att_val_bottom))
% 
% % figure;plot(lat,att_val_bott);
% % figure;plot(lat,att_val_bott-nanmean(att_val_bott))
% 
% %att_val_bott_filt=sgolayfilt(att_val_bott,3,1001);
% % figure;plot(lat,att_val_bott);
% % figure;plot(lon,att_val_bott_filt);
% figure;plot(lat,Na);title('Na')
% 
% 
% load('/cresis/snfs1/scratch/manjish/new_peterman/reflectivity/crossline10.mat');
% 
% figure;plot(out.Latitude,out.const_attenuation)
% nanidx=find(isnan(out.Longitude));
% out.const_attenuation(nanidx)=[];
% out.Longitude(nanidx)=[];
% 
% figure;plot(lat,-att_val_bott-nanmean(-att_val_bott));
% hold on; plot(out.Latitude,out.const_attenuation)
% 
% %figure;plot(lat,depth)
% figure;plot(lat,-att_val_bott);
% figure;plot(-att_val_bott)
% figure;plot(Depth); title('Depth')
% % lastidx=find(out.Longitude<=lon(end),1,'first');
% %
% % att_val_method2=out.const_attenuation(1:lastidx);
% % att_var_method2=out.var_attenuation(1:lastidx);
% % LON=out.Longitude(1:lastidx);
% % LAT=out.Latitude(1:lastidx);
% %
% % figure;plot(lon,att_val_bott-nanmean(att_val_bott));
% % hold on; plot(LON,att_val_method2)
% %
% % figure;plot(lon,att_val_bott_filt-nanmean(att_val_bott_filt));
% % hold on; plot(LON,att_val_method2)
% 
% % figure;plot(out.Longitude,out.const_attenuation);
% % att_cal_m1=interp1(out.Longitude,out.const_attenuation,lon);
% %
% % figure;plot(lat,att_val_bott_filt-nanmean(att_val_bott_filt));
% % hold on; plot(lon,att_cal_m1);
% 
% % r.rms_height(k) = 0.001;
% % options = optimset('PlotFcns',@optimplotfval);
% % cost_func=@(r_rms)(abs(rms_fit - exp(-(2*(2*pi/(c/param.radar.fc))*r_rms).^2)/((r_rms).^2)));
% % r_rms=[0.00001,1];
% % [r_rmsheight,fval]=fminsearch(cost_func,r_rms);
% % r.rms_height(k)=mean(abs(r_rmsheight));
% % power_loss=lp(lpower(:))-lp(lpower(1));
% % [r,m,b]=regression(dpth,power_loss);
% 
% % figure;scatter(dpth(2:end),power_loss(2:end))
% %%Median Filtering
% % figure;