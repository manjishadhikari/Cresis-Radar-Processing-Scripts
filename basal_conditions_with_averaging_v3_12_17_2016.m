clear all
close all
clc
dbstop error
if 1
%% calculating Na avg
% load(['C:\Users\s343m141\Documents\scripts\matlab\thesis\ice_loss_estimation_paper_data\after_roughness_loss_correction\get heights frames Greenland\Greenland_layerdata_selected_frames_complete_v6.mat'])
if 1
    out_fn=['/cresis/snfs1/scratch/manjish//peterman/completedata.mat'];
    load(out_fn);
end
physical_constants
plots =0;

clear idx
 idx = find(isnan(Greenland.ice_bed_power)) ;
        Greenland.GPS_time(idx) = [];
        Greenland.Latitude(idx) = [];
        Greenland.Longitude(idx) = [];
        Greenland.Elevation(idx) = [];
        Greenland.ice_bed_time(idx) = [];
        Greenland.surface_time(idx) = [];
        Greenland.ice_bed_power(idx) = [];
        Greenland.ice_surface_power(idx) = [];

clear idx
 idx = find(isnan(Greenland.ice_surface_power)) ;
        Greenland.GPS_time(idx) = [];
        Greenland.Latitude(idx) = [];
        Greenland.Longitude(idx) = [];
        Greenland.Elevation(idx) = [];
        Greenland.ice_bed_time(idx) = [];
        Greenland.surface_time(idx) = [];
        Greenland.ice_bed_power(idx) = [];
        Greenland.ice_surface_power(idx) = [];

Greenland.depth = (Greenland.ice_bed_time - Greenland.surface_time)*c/2/sqrt(er_ice);
Greenland.surface_height = (Greenland.surface_time)*c/2;
% figure;plot(tmp.Greenland.depth,lp(tmp.Greenland.ice_bed_power_cgl))

if plots
    plot(Greenland.depth, lp((Greenland.ice_bed_power)));
    grid on%verticallines
    title('Depth vs Ice Bed Power')
    % hist(10*log10(Greenland.ice_bed_power),40)
end

[Greenland.depth_sorted Greenland.index]= sort(Greenland.depth);
Greenland.ice_bed_power_sorted = Greenland.ice_bed_power(Greenland.index);
Greenland.surface_height_sorted = Greenland.surface_height(Greenland.index);
Greenland.Latitude_sorted = Greenland.Latitude(Greenland.index);
Greenland.Longitude_sorted = Greenland.Longitude(Greenland.index);

if plots
    figure;plot(Greenland.depth_sorted, 10*log10(abs(Greenland.ice_bed_power_sorted).^2));
    grid on
    title('Depth vd Ice Bed Power after sorting')
end
geometric_loss_sorted = (2*(Greenland.surface_height_sorted+Greenland.depth_sorted)).^2;
Greenland.ice_bed_power_cgl_sorted=Greenland.ice_bed_power_sorted.*geometric_loss_sorted;
%Greenland.ice_bed_power_cgl_sorted = Greenland.ice_bed_power_cgl(Greenland.index);

 if plots
    figure;
    plot(Greenland.depth, 10*log10(abs(Greenland.ice_bed_power).^2));
    grid on
    title('Ice Bed Power vs Depth'); xlabel('Depth')
    
    figure;
    plot(10*log10(abs(Greenland.ice_bed_power).^2));
    grid on
    title('Ice Bed Power vs Along track'); xlabel('Along Track')
    
    figure; plot(Greenland.depth);
    title('Along track vs Depth');
    xlabel('Along Track'); ylabel('Depth')
    
  end

clear idx
idx = find(isnan(Greenland.ice_bed_power_cgl_sorted)) ;
Greenland.GPS_time(idx) = [];
Greenland.Latitude(idx) = [];
Greenland.Longitude(idx) = [];
Greenland.Elevation(idx) = [];
Greenland.ice_bed_time(idx) = []; %figure;plot(tmp.Greenland.depth,lp(tmp.Greenland.ice_bed_power_cgl))
Greenland.surface_time(idx) = [];
Greenland.ice_bed_power_cgl_sorted(idx) = [];
Greenland.depth_sorted (idx)= [];


Greenland.ice_bed_power_rgl_sorted = Greenland.ice_bed_power_cgl_sorted./median(Greenland.ice_bed_power_cgl_sorted);
if plots
    figure;plot(Greenland.depth_sorted, lp((Greenland.ice_bed_power_rgl_sorted)));
    grid on
    title('After Mean Removed')
end
avg_power = median(Greenland.ice_bed_power_cgl_sorted);
avg_depth = mean(Greenland.depth_sorted);

%% assuming constant attenuation rate
% if plots
%     N = [2 5];
%     for k = 1: 2
%         Na  = N(k); % 3 - 20dB attenuation rate
%         Greenland.ice_bed_reflectivity_sorted = 10*log10(Greenland.ice_bed_power_rgl_sorted)+ (2*Na*(Greenland.depth_sorted-mean(Greenland.depth_sorted))/1e3);
%         
%         figure(1)
%         subplot(2,1,k)
%         plot(1:length(Greenland.ice_bed_reflectivity_sorted),Greenland.ice_bed_reflectivity_sorted,'*-');
%         figure(2)
%         subplot(2,1,k)
%         hist(real(Greenland.ice_bed_reflectivity_sorted),40)
%     end
% end

%%
window_size = 10000;
% Greenland.ice_bed_power_frgc = conv(Greenland.ice_bed_power_rgc,gausswin(window_size),'same');
Greenland.ice_bed_power_frgl_sorted  = (sgolayfilt((lp((Greenland.ice_bed_power_rgl_sorted))), 2,window_size+1, gausswin(window_size+1)));
%Greenland.ice_bed_power_frgl_sortedtst  = (sgolayfilt((lp((Greenland.ice_bed_power_rgl_sorted))), 2,window_size+1));

if plots
    figure;plot((1:length(Greenland.ice_bed_power_rgl_sorted)),lp(Greenland.ice_bed_power_rgl_sorted));
    hold on
    plot((1:length(Greenland.ice_bed_power_frgl_sorted)),(Greenland.ice_bed_power_frgl_sorted));
    title('Filtering Data')
end

window_size = 10000;
Greenland.depths_sorted  = sgolayfilt(Greenland.depth_sorted, 2,window_size+1, gausswin(window_size+1));
if plots
    figure;plot((1:length(Greenland.depth_sorted)),Greenland.depth_sorted);
    hold on
    plot((1:length(Greenland.depths_sorted)),(Greenland.depths_sorted));
    title('Filtered Depth ')
    
    %%
    figure
    subplot(2,1,1)
    plot((1:length(Greenland.depth_sorted)),(Greenland.depths_sorted-mean(Greenland.depths_sorted))/1e3);
    subplot(2,1,2)
    plot((1:length(Greenland.ice_bed_power_frgl_sorted)),(Greenland.ice_bed_power_frgl_sorted));
   
end
%%
na = 0:0.01:15;
for N= 1:length(na)
    S(N) = mean((((2*na(N)*(Greenland.depth_sorted-mean(Greenland.depth_sorted))/1e3)+((Greenland.ice_bed_power_frgl_sorted))).^2));
    %       S(N) = mode(round(abs((2*N*(Greenland.depths-mean(Greenland.depths))/1e3)+(10*log10(Greenland.ice_bed_power_frgc)))));
    %        S(N) = sum(abs((2*N*(Greenland.depths-mean(Greenland.depths))/1e3)+(real(Greenland.ice_bed_power_frgc))).^2);
end
[v i] = min(S)
if i==length(S)
    warning('check this')
end
Na = na(i)
%%
if plots
    figure;plot((1:length(Greenland.ice_bed_power_frgl_sorted)),((Greenland.ice_bed_power_frgl_sorted)));
    hold on
    plot((1:length(Greenland.depth_sorted)),2*Na*(Greenland.depths_sorted-mean(Greenland.depths_sorted))/1e3);
    grid on
    title('2Na(d-dx) ')
end
end
%%
close all

Reflectivity_values = [];
attenuation = [];
power_rgc = [];
lt  = []; 
ln = []; 
for M = 1:20
%     keyboard
    %   for M = [5,6,18,19];
    clearvars -except M Na avg_power avg_depth Reflectivity_values attenuation power_rgc param lt ln
    clc
    plots=0;
    param.radar.fs = 1000000000/9;
  load(['/cresis/snfs1/scratch/manjish/peterman/radar/crossline', sprintf('%d.mat',M)]);

    %  load(['/cresis/snfs1/scratch/santhosh/thesis/get heights frames Greenland/Greenland_layerdata_selected_frames_v6_' num2str(M,'%03d') '.mat'])
    physical_constants
    
    
    %     proj = geotiffinfo('X:\GIS_data\antarctica\Landsat-7\Antarctica_LIMA_480m.tif');
    %     [x,y] = projfwd(proj,Greenland.Latitude,Greenland.Longitude);
    %     along_track = [0 cumsum(sqrt(diff(x).^2 + diff(y).^2))]/1e3;
    
    %     along_tack_fn = geodetic_to_along_track(Greenland.Latitude,Greenland.Longitude,Greenland.Elevation);
    
    
    Greenland.depth = (Greenland.ice_bed_time - Greenland.surface_time)*c/2/sqrt(er_ice);
    Greenland.surface_height = (Greenland.surface_time)*c/2;
    
    %  [Greenland.depth_sorted index]= sort(Greenland.depth);
    %  Greenland.ice_bed_power_sorted = Greenland.ice_bed_power(index);
    % Greenland.ice_bed_power = Greenland.ice_bed_power_sorted;
    % Greenland.depth = Greenland.depth_sorted;
    % Greenland.surface_height_sorted = Greenland.surface_height(index);
    % Greenland.surface_height = Greenland.surface_height_sorted ;
    
    geometric_loss = (2*(Greenland.surface_height+Greenland.depth/sqrt(er_ice))).^2;
    geometric_loss_surface = (2*(Greenland.surface_height)).^2;
    Greenland.ice_bed_power_cgl =(Greenland.ice_bed_power).*geometric_loss;
 % Greenland.ice_bed_power_cgl=Greenland.ice_bed_power;
  
  
    if plots
      figure;
        plot(Greenland.depth, 10*log10(abs(Greenland.ice_bed_power).^2));
        grid on
        title('Ice Bed Power vs Depth'); xlabel('Depth')
        
        figure;
         plot(10*log10(abs(Greenland.ice_bed_power).^2));
        grid on
        title('Ice Bed Power vs Along track'); xlabel('Along Track')
        
        figure; plot(Greenland.depth);
        title('Along track vs Depth');
        xlabel('Along Track'); ylabel('Depth')
        
        figure;plot(lp(Greenland.ice_bed_power_cgl))
        title('Geom corrected')
        
    end
    
   
    %% compensating for surface roughness
    
     Greenland.ice_surface_power_cgl = (Greenland.ice_surface_power).*geometric_loss_surface;
    if exist((['/cresis/snfs1/scratch/manjish/peterman/radarnew/',sprintf('crossline%d.mat',M)]),'file')
         load((['/cresis/snfs1/scratch/manjish/peterman/radarnew/',sprintf('crossline%d.mat',M)]));
   
        K  = floor(length(Greenland.ice_bed_power_cgl)/1000);
        for l = 1:K
            
            if any(~isnan(Greenland.ice_bed_power_cgl(1+(l-1)*1000:(l*1000))))
                ice_bed_power_cgl = Greenland.ice_bed_power_cgl(1+(l-1)*1000:(l*1000));
                depth = Greenland.depth(1+(l-1)*1000:(l*1000));
                clear id
                id = find(isnan(Greenland.ice_bed_power_cgl(1+(l-1)*1000:(l*1000))));
                ice_bed_power_cgl(id) = [];
                depth(id) = [];
                Greenland.ice_bed_power_cgl_avg(l) =   mean(ice_bed_power_cgl) ;
                Greenland.depth_avg(l) = mean(depth);
                Greenland.Latitude_avg(l) = mean(Greenland.Latitude(1+(l-1)*1000:(l*1000)));
                Greenland.Longitude_avg(l) = mean(Greenland.Longitude(1+(l-1)*1000:(l*1000)));
            else
                Greenland.ice_bed_power_cgl_avg(l) = nan;
                Greenland.depth_avg(l) = nan;
                Greenland.Latitude_avg(l) = nan;
                Greenland.Longitude_avg(l) = nan;
            end
            
            if isnan(r.rms_height(l))
                continue;
            else
                Greenland.ice_bed_power_cgl_avg(l) = Greenland.ice_bed_power_cgl_avg(l)./(exp(-(4*pi*r.rms_height(l)/(c/param.radar.fs))^2)*(besseli(0,((4*pi*r.rms_height(l)/(c/param.radar.fs))^2)/2))^2);
            end
        end
        
        clearvars r  K
    else
        Greenland.ice_surface_power_cgl = (Greenland.ice_surface_power).*geometric_loss_surface;
        K  = floor(length(Greenland.ice_surface_power_cgl)/1000);
        for l = 1:K
            r.lat(l) = mean(Greenland.Latitude(1+(l-1)*1000:(l*1000)));
            r.lon(l) = mean(Greenland.Longitude(1+(l-1)*1000:(l*1000)));
            s = sqrt(Greenland.ice_surface_power_cgl(1+(l-1)*1000:(l*1000)));
            id = find(isnan(s)|isinf(s)|s==0);
            if length(id) > 500
                r.rms_height(l) = nan;
                r.dielectric_constant(l) = nan;
                r.pn(l) = nan;
                r.pc(l) = nan ;
                continue
            else
                s(id) = [];
            end
            pd = fitdist(double((s)).','Rician')
            % phat = mle(double(abs(s)),'distribution','Rician');
            % x = 0:0.0001:0.2;
            % histogram(abs(s))
            % h = hist((s));
            
            % A = pdf(pd,x);
            a = pd.s;
            % pc = 2*10*log10(a);
            % pn = 10*log10(2*pd.sigma^2);
            S = pd.sigma;
            r.pc(l) = a^2;
            r.pn(l) = 2*2*pd.sigma^2;
            rms_fit = (r.pc(l)/r.pn(l))*4*(2*pi/(c/param.radar.fs))^2;
            
            r.rms_height(l) = 0.0001;
            clear MSE
            for i = 1:5000
                MSE(i) = abs(rms_fit - exp(-(2*(2*pi/(c/param.radar.fs))*r.rms_height(l))^2)/((r.rms_height(l))^2));
                if r.rms_height(l) > 0.40
                    r.rms_height(l) = nan;
                    warning('check this')
                    %                         keyboard
                    break
                else
                    
                    if i>1
                        if MSE(i-1) < MSE(i)
                            break
                        else
                            r.rms_height(l) = r.rms_height(l) + 0.0001;
                            continue  ;
                        end
                    else
                        r.rms_height(l) = r.rms_height(l) + 0.0001;
                    end
                end
            end
            
            
            if isnan(r.rms_height(l))
                r.dielectric_constant(l) = nan;
                continue
            else
                r.dielectric_constant(l) = 1;
                clear mse
                for i = 1:5000
                    mse(i) = abs(r.pc(l) - ((1-sqrt(r.dielectric_constant(l)))/((1+sqrt(r.dielectric_constant(l)))))^2*exp(-(2*(2*pi/(c/param.radar.fs))*r.rms_height(l))^2));
                    if r.dielectric_constant(l) > 4
                        r.dielectric_constant(l) = nan;
                        warning('check this')
                        %                         keyboard
                        break
                    else
                        
                        if i>1
                            if mse(i-1) < mse(i)
                                break
                            else
                                r.dielectric_constant(l) = r.dielectric_constant(l) + 0.01;
                                continue  ;
                            end
                        else
                            r.dielectric_constant(l) = r.dielectric_constant(l) + 0.01;
                        end
                    end
                end
            end
            
            
            if any(~isnan(Greenland.ice_bed_power_cgl(1+(l-1)*1000:(l*1000))))
                ice_bed_power_cgl = Greenland.ice_bed_power_cgl(1+(l-1)*1000:(l*1000));
                depth = Greenland.depth(1+(l-1)*1000:(l*1000));
                clear id
                id = find(isnan(Greenland.ice_bed_power_cgl(1+(l-1)*1000:(l*1000))));
                ice_bed_power_cgl(id) = [];
                depth(id) = [];
                Greenland.ice_bed_power_cgl_avg(l) =   mean(ice_bed_power_cgl) ;
                Greenland.depth_avg(l) = mean(depth);
                Greenland.Latitude_avg(l) = mean(Greenland.Latitude(1+(l-1)*1000:(l*1000)));
                Greenland.Longitude_avg(l) = mean(Greenland.Longitude(1+(l-1)*1000:(l*1000)));
            else
                Greenland.ice_bed_power_cgl_avg(l) = nan;
                Greenland.depth_avg(l) = nan;
                Greenland.Latitude_avg(l) = nan;
                Greenland.Longitude_avg(l) = nan;
            end
            
            if isnan(r.rms_height(l))
                continue;
            else
                Greenland.ice_bed_power_cgl_avg(l) = Greenland.ice_bed_power_cgl_avg(l)./(exp(-(4*pi*r.rms_height(l)/(c/param.radar.fs))^2)*(besseli(0,((4*pi*r.rms_height(l)/(c/param.radar.fs))^2)/2))^2);
            end
        end
        
        save((['/cresis/snfs1/dataproducts/ct_data/rds/2012_Greenland_P3/CSARP_manjish/surfaceroughness_after/',param.day_seg,'/', sprintf('IceBedCoherenceIndex_%s.mat', param.day_seg)]),'r');
        %save(['/cresis/snfs1/scratch/santhosh/thesis/get heights frames Greenland/surface_roughness_v6_' num2str(M,'%03d') '.mat'],'r');
        clearvars r K
    end
    
     if plots
       figure;
       plot(10*log10(abs( Greenland.ice_bed_power_cgl_avg).^2));
        grid on
        title('Ice Bed Power surface roughness corrected')
    end
    
    
    
    %% compensating for bed roughness
    if exist((['/cresis/snfs1/scratch/manjish/peterman/bedroughness/',sprintf('crossline%d.mat',M)]),'file')
         load((['/cresis/snfs1/scratch/manjish/peterman/bedroughness/',sprintf('crossline%d.mat',M)]));
          r=rbed;
        K  = floor(length(Greenland.ice_bed_power_cgl)/1000);
        for l = 1:K
            
            if isnan(r.rms_height(l)) || isnan(Greenland.ice_bed_power_cgl_avg(l))
                continue ;
            else
                Greenland.ice_bed_power_cgl_avg(l) = Greenland.ice_bed_power_cgl_avg(l)./(exp(-(4*pi*r.rms_height(l)/(c/param.radar.fs))^2)*(besseli(0,((4*pi*r.rms_height(l)/(c/param.radar.fs))^2)/2))^2);
            end
        end
        
        %           Greenland.ice_bed_power_cgl(1+(l-1)*1000:(l*1000)) = (Greenland.ice_bed_power_cgl(1+(l-1)*1000:(l*1000)))./(exp(-(4*pi*r.rms_height(l)/(c/param.radar.fs))^2)*(besseli(0,((4*pi*r.rms_height(l)/(c/param.radar.fs))^2)/2))^2);
        
        
        
        clearvars r K
    else
        
        K  = floor(length(Greenland.ice_bed_power_cgl)/1000);
        for l = 1:K
            r.lat(l) = mean(Greenland.Latitude(1+(l-1)*1000:(l*1000)));
            r.lon(l) = mean(Greenland.Longitude(1+(l-1)*1000:(l*1000)));
            s = sqrt(Greenland.ice_bed_power_cgl(1+(l-1)*1000:(l*1000)));
            id = find(isnan(s)|isinf(s)|s==0);
            if length(id) > 500
                r.rms_height(l) = nan;
                r.dielectric_constant(l) = nan;
                r.pn(l) = nan;
                r.pc(l) = nan ;
                continue
            else
                s(id) = [];
            end
            pd = fitdist(double((s)).','Rician')
          
            a = pd.s;
            % pc = 2*10*log10(a);
            % pn = 10*log10(2*pd.sigma^2);
            S = pd.sigma;
            r.pc(l) = a^2;
            r.pn(l) = 2*2*pd.sigma^2;
            rms_fit = (r.pc(l)/r.pn(l))*4*(2*pi/(c/param.radar.fs))^2;
            
            r.rms_height(l) = 0.0001;
            clear MSE
            for i = 1:5000
                MSE(i) = abs(rms_fit - exp(-(2*(2*pi/(c/param.radar.fs))*r.rms_height(l))^2)/((r.rms_height(l))^2));
                if r.rms_height(l) > 0.40
                    r.rms_height(l) = nan;
                    warning('check this')
                    %                         keyboard
                    break
                else
                    
                    if i>1
                        if MSE(i-1) < MSE(i)
                            break
                        else
                            r.rms_height(l) = r.rms_height(l) + 0.0001;
                            continue  ;
                        end
                    else
                        r.rms_height(l) = r.rms_height(l) + 0.0001;
                    end
                end
            end
            
            
            if isnan(r.rms_height(l))
                r.dielectric_constant(l) = nan;
                continue
            else
                r.dielectric_constant(l) = 1;
                clear mse
                for i = 1:5000
                    mse(i) = abs(r.pc(l) - ((1-sqrt(r.dielectric_constant(l)))/((1+sqrt(r.dielectric_constant(l)))))^2*exp(-(2*(2*pi/(c/param.radar.fs))*r.rms_height(l))^2));
                    if r.dielectric_constant(l) > 4
                        r.dielectric_constant(l) = nan;
                        warning('check this')
                        %                         keyboard
                        break
                    else
                        
                        if i>1
                            if mse(i-1) < mse(i)
                                break
                            else
                                r.dielectric_constant(l) = r.dielectric_constant(l) + 0.01;
                                continue  ;
                            end
                        else
                            r.dielectric_constant(l) = r.dielectric_constant(l) + 0.01;
                        end
                    end
                end
            end
            
            
            %        Greenland.ice_bed_power_cgl(1+(l-1)*1000:(l*1000)) = (Greenland.ice_bed_power_cgl(1+(l-1)*1000:(l*1000)))./(exp(-(4*pi*r.rms_height(l)/(c/param.radar.fs))^2)*(besseli(0,((4*pi*r.rms_height(l)/(c/param.radar.fs))^2)/2))^2);
            
            
            if any(~isnan(Greenland.ice_bed_power_cgl(1+(l-1)*1000:(l*1000))))
                ice_bed_power_cgl = Greenland.ice_bed_power_cgl(1+(l-1)*1000:(l*1000));
                depth = Greenland.depth(1+(l-1)*1000:(l*1000));
                clear id
                id = find(isnan(Greenland.ice_bed_power_cgl(1+(l-1)*1000:(l*1000))));
                ice_bed_power_cgl(id) = [];
                depth(id) = [];
                Greenland.ice_bed_power_cgl_avg(l) =   mean(ice_bed_power_cgl) ;
                Greenland.depth_avg(l) = mean(depth);
                Greenland.Latitude_avg(l) = mean(Greenland.Latitude(1+(l-1)*1000:(l*1000)));
                Greenland.Longitude_avg(l) = mean(Greenland.Longitude(1+(l-1)*1000:(l*1000)));
            else
                Greenland.ice_bed_power_cgl_avg(l) = nan;
                Greenland.depth_avg(l) = nan;
                Greenland.Latitude_avg(l) = nan;
                Greenland.Longitude_avg(l) = nan;
                Greenland.Latitude(1+(l-1)*1000:(l*1000)) =[];
                Greenland.Longitude(1+(l-1)*1000:(l*1000))= [];
            end
            
            
            
            
            
            if isnan(r.rms_height(l))
                continue;
            else
                Greenland.ice_bed_power_cgl_avg(l) = Greenland.ice_bed_power_cgl_avg(l)./(exp(-(4*pi*r.rms_height(l)/(c/param.radar.fs*sqrt(er_ice)))^2)*(besseli(0,((4*pi*r.rms_height(l)/(c/param.radar.fs*sqrt(er_ice)))^2)/2))^2);
            end
        end
        
         save((['/cresis/snfs1/dataproducts/ct_data/rds/2012_Greenland_P3/CSARP_manjish/bedroughness_after/',param.day_seg,'/', sprintf('IceBedCoherenceIndex_%s.mat', param.day_seg)]),'r');
       
       % save(['/cresis/snfs1/scratch/santhosh/thesis/get heights frames Greenland/bed_roughness_v6_' num2str(M,'%03d') '.mat'],'r');
        clearvars r K
        
    end
    
      if plots
       hold on;
        plot(10*log10(abs( Greenland.ice_bed_power_cgl_avg).^2));
        grid on
       legend('Sf corrected ',' Bed roughness corrected')
    end
    %%  relative geometrically corrected bed - echo power
    %     for i = 1:length(Greenland.ice_bed_power_cgl)
    %         if i < 3
    %             Greenland.ice_bed_power_rgc(i) = Greenland.ice_bed_power_cgl(i)-mean(Greenland.ice_bed_power_cgl(1:6-i));
    %         elseif i+2 > length(Greenland.ice_bed_power_cgl)
    %             Greenland.ice_bed_power_rgc(i) = Greenland.ice_bed_power_cgl(i)- mean(Greenland.ice_bed_power_cgl(i-5:end));
    %         else
    %             Greenland.ice_bed_power_rgc(i) = Greenland.ice_bed_power_cgl(i)- mean(Greenland.ice_bed_power_cgl(i-2:i+2));
    %         end
    %     end
    
    
   
    
    Greenland.ice_bed_power_rgc = Greenland.ice_bed_power_cgl_avg;%/avg_power;
    if plots
        figure;plot(Greenland.depth_avg,10*log10(Greenland.ice_bed_power_rgc));
        grid on;
    end
    %% assuming constant attenuation rate
    %     if plots
    %         N = [2 10];
    %         for k = 1: 2
    %             Na  = N(k); % 3 - 20dB attenuation rate
    %             Greenland.ice_bed_reflectivity = lp(Greenland.ice_bed_power_rgc)+ (2*Na*(Greenland.depth-mean(Greenland.depth))/1e3);
    %
    %             figure(1)
    %             subplot(2,1,k)
    %             plot(Greenland.ice_bed_reflectivity,'*-');
    %             figure(2)
    %             subplot(2,1,k)
    %             hist(Greenland.ice_bed_reflectivity,40)
    %         end
    %     end
    %%
    if plots
        plot((1:length(Greenland.ice_bed_power_rgc)),10*log10(Greenland.ice_bed_power_rgc));
        figure
        plot(1:length(Greenland.depth_avg),(Greenland.depth_avg/1e3))
    end
    %%
    
    %     clear idx
    idx = find(isnan(Greenland.ice_bed_power_rgc)) ;
    Greenland.ice_bed_power_rgc(idx) = [];
    Greenland.depth_avg(idx) =[];
    Greenland.Latitude_avg(idx) =[];
    Greenland.Longitude_avg(idx) =[];
    %
    %% applying a guassian filter
    %   h = fspecial('gaussian', [1,round(length(Greenland.ice_bed_power_rgc)/5)], 300);
    %  % h = h./max(h); % normalizing
    % %  Greenland.ice_bed_power_frgc = filter(h,1,Greenland.ice_bed_power_rgc);
    %  Greenland.ice_bed_power_frgc = conv(Greenland.ice_bed_power_rgc,h,'same');
    if length(Greenland.ice_bed_power_rgc)< 500
        if mod(length(Greenland.ice_bed_power_rgc),2)==0
            window_size =length(Greenland.ice_bed_power_rgc)-100;
        else
            window_size =length(Greenland.ice_bed_power_rgc)-99;
        end
    else
        window_size = 500;
    end
    if length(Greenland.ice_bed_power_rgc)<50
      window_size=2;
    end
    % Greenland.ice_bed_power_frgc = conv(Greenland.ice_bed_power_rgc,gausswin(window_size),'same');
    Greenland.ice_bed_power_frgc  = sgolayfilt(10*log10(Greenland.ice_bed_power_rgc), 2,window_size+1, gausswin(window_size+1));
    if plots
        plot((1:length(Greenland.ice_bed_power_rgc)),10*log10(Greenland.ice_bed_power_rgc));
        hold on
        plot((1:length(Greenland.ice_bed_power_frgc)),(Greenland.ice_bed_power_frgc));
    end
    if length(Greenland.ice_bed_power_rgc)< 500
        if mod(length(Greenland.ice_bed_power_rgc),2)==0
            window_size =length(Greenland.ice_bed_power_rgc)-100;
        else
            window_size =length(Greenland.ice_bed_power_rgc)-99;
        end
    else
        window_size = 500;
    end
     if length(Greenland.ice_bed_power_rgc)<100
      window_size=2;
    end
    Greenland.depths  = sgolayfilt(Greenland.depth_avg, 2,window_size+1, gausswin(window_size+1));
    if plots
        plot((1:length(Greenland.depth_avg)),Greenland.depth_avg);
        hold on
        plot((1:length(Greenland.depths)),(Greenland.depths));
    end
    %     %
    %     %  for i = 1:length(Greenland.depth)
    %     %      if i < 3
    %     %
    %     %          Greenland.dif_depth(i) = Greenland.depths(i)-mean(Greenland.depths(1:6-i));
    %     %
    %     %     elseif i+2 > length(Greenland.depth)
    %     %     Greenland.dif_depth(i) = Greenland.depths(i)- mean(Greenland.depths(i-5:end));
    %     %     else
    %     %       Greenland.dif_depth(i) = Greenland.depths(i)- mean(Greenland.depths(i-2:i+2));
    %     %      end
    %     %  end
    %     %
    %     % %  Greenland.fdepth  = sgolayfilt(Greenland.dif_depth, 2,window_size+1, gausswin(window_size+1));
    %     %   plot((1:length(Greenland.dif_depth)),(10*log10(Greenland.dif_depth)));
    %     % %    hold on
    %     % % plot((1:length(Greenland.ice_bed_power_rgc)),10*log10(-Greenland.ice_bed_power_rgc)/2);
    %     % % %  hold on
    %     % % %  plot((1:length(Greenland.fdepth)),10*log10(Greenland.fdepth));
    %
    %     %%
    %     %     na =(-(Greenland.ice_bed_power_frgc)- 10*log10(2*Greenland.dif_depth));
    %     %   plot([1:length(na)],na)
    %     % % Na  = mean(na);
    %     %  for na =1:40
    %     %      mse(na) = sum(abs(real(Greenland.ice_bed_power_frgc) + (2*na*((Greenland.depths-mean(Greenland.depths))/1e3))));
    %     %  end
    %     % %  [v i] =min(mse)
    %     % sum(abs((2*12*(Greenland.depths-mean(Greenland.depths))/1e3)+(10*log10(Greenland.ice_bed_power_frgc))));
    %     for N= 1:40
    %          S(N) = mean(abs((2*N*(Greenland.depths-mean(Greenland.depths))/1e3)+(real(Greenland.ice_bed_power_frgc))).^2);
    %  %       S(N) = mode(round(abs((2*N*(Greenland.depths-mean(Greenland.depths))/1e3)+(10*log10(Greenland.ice_bed_power_frgc)))));
    %  %        S(N) = sum(abs((2*N*(Greenland.depths-mean(Greenland.depths))/1e3)+(real(Greenland.ice_bed_power_frgc))).^2);
    %     end
    %     [v i] = min(S)
    %     Na = i;
    %%
    if plots
        figure
        plot((1:length(Greenland.depths)),2*Na*((Greenland.depths-avg_depth)/1e3));
        hold on
        plot((1:length(Greenland.ice_bed_power_frgc)),(-(Greenland.ice_bed_power_frgc)));
    end
    %% to find d(na)/dx
    %     L = cumsum(Greenland.segments_length);
    %     window_size = 36;
    %     debug_flag= 1;
    %     for i = 1: length(L)
    %         if i == 1;
    %             p=1;
    %             q = L(i);
    %         else
    %             p = 1+L(i-1);
    %             q = L(i);
    %         end
    %         along_track_fn{i}= geodetic_to_along_track(Greenland.Latitude(p:q),Greenland.Longitude(p:q));
    %         power{i} = Greenland.ice_bed_power_rgc(p:q);
    %         depth{i} = Greenland.depth(p:q);
    %         depths{i}  = sgolayfilt(depth{i}, 2,window_size+1, gausswin(window_size+1));
    %         power_frgc{i}  = sgolayfilt(10*log10(power{i}), 2,window_size+1, gausswin(window_size+1));
    %         if debug_flag
    %             %             plot((1:length(power{i})),10*log10(power{i}));
    %             %             hold on
    %             %             plot((1:length(power{i})),(power_frgc{i}));
    %             %             keyboardframes
    %             %             clf
    %             %             figure
    %             %             plot((1:length(depth{i})),depth{i});
    %             %             hold on
    %             %             plot((1:length(depths{i})),(depths{i}));
    %             %             keyboard
    %             %             clf
    %             if plots
    %                         subplot(2,1,1)
    %                         plot((1:length(power{i})),(-power_frgc{i}));
    %                         subplot(2,1,2)
    %                         plot((1:length(depths{i})),(depths{i}-mean(depths{i}))/1e3);
    %                         keyboard
    %                         clf
    %             end
    %         end
    %     end
    L = cumsum(Greenland.segments_length);
    debug_flag= 1;
    Greenland.along_track = [];
    for i = 1: length(L)
        if i == 1;
            p=1;
            q = L(i);
            Lat = Greenland.Latitude(p:q);
            Lon = Greenland.Longitude(p:q);
            
            %             idx = find(isnan(Greenland.ice_bed_power(p:q)));
            %             Lat(idx) = [];
            %             Lon(idx) = [];
            along_track_fn{i}= geodetic_to_along_track(Lat,Lon);
            Greenland.along_track = cat(2, Greenland.along_track, along_track_fn{i}) ;
        else
            p = 1+L(i-1);
            q = L(i);
            %             idx = find(isnan(Greenland.ice_bed_power(p:q)));
            Lat = Greenland.Latitude(p:q);
            Lon = Greenland.Longitude(p:q);
            %             Lat(idx) = [];
            %             Lon(idx) = [];
            along_track_fn{i}= geodetic_to_along_track(Lat,Lon);
            Greenland.along_track = cat(2, Greenland.along_track, Greenland.along_track(end)+ along_track_fn{i}) ;
        end
    end
    
    
    K  = floor(length(Greenland.ice_bed_power_cgl)/1000);
    for l = 1:K
        if any(~isnan(Greenland.along_track(1+(l-1)*1000:(l*1000))))
            along_track = Greenland.along_track(1+(l-1)*1000:(l*1000));
            clear id
            id = find(isnan(along_track));
            along_track(id) = [];
            Greenland.along_track_avg(l) =   mean(along_track) ;
        else
            Greenland.along_track_avg(l) = nan;
            continue ;
        end
    end
    
    Greenland.along_track_avg(idx) = [];
    
    if ~(length(Greenland.depth_avg) == length(Greenland.along_track_avg))
        keyboard;
    end
    
    
    
    
    %%
    %     dn = (1:40);
    %     for i = 1: length(power_frgc)
    %         for j = 1:length(dn)
    %
    %             term_1 = 2*dn(j)*((depths{i}-mean(depths{i}))/1e3).*((along_track_fn{i}-mean(along_track_fn{i}))/1e3);
    %             term_2 = 2*Na*((depths{i}-mean(depths{i}))/1e3);
    %             S(i,j) = sum(abs(real(power_frgc{i})+term_1+term_2).^2);
    %             a = S(i,j);
    %             if plots
    %                     plot(real(-power_frgc{i}))
    %                     hold on
    %                     plot(term_1+term_2)
    %                     keyboard
    %                     clf
    %             end
    %             if j > 1
    %                 if S(i,j-1)<S(i,j)
    %                     DN(i) = dn(j-1);
    %                     break
    %                 end
    %             end
    %         end
    %     end
    % end
    % % s = 0;
    % %  for i = 1:length(along_track_fn)
    % % s= s+along_track_fn{i}(end);
    % % end
    
    dn = (-0.75:0.01:0.75);
    for j = 1:length(dn)
        
        term_1 = 2*dn(j).*((Greenland.depths-avg_depth)/1e3).*((Greenland.along_track_avg-mean(Greenland.along_track_avg))/1e3);
        term_2 = 2*Na*((Greenland.depths-avg_depth)/1e3);
        S(j) = sum(abs(real(Greenland.ice_bed_power_frgc)+term_1+term_2).^2);
        if plots
            plot(lp(Greenland.ice_bed_power_frgc))
            hold on
            plot(term_1+term_2)
        %    keyboard
            clf
        end
        if j > 1
            if S(j-1)<S(j)
                DN = dn(j-1);
                break
            end
        end
    end
    
    if j ==length(dn)
        [v id] = min(S);
        DN = dn(id);
    end
    
    
    
    %%
    %     Reflectivity = [];
    %     for i = 1:length(power_frgc)
    %         r  =  real(power_frgc{i})+ 2*DN(i)*((depths{i}-mean(depths{i}))/1e3).*((along_track_fn{i}-mean(along_track_fn{i}))/1e3)+ 2*Na*((depths{i}-mean(depths{i}))/1e3);
    %         Reflectivity = cat(2, Reflectivity, r);
    %     end
    %     if plots
    %         figure
    %         hist(Reflectivity,30)
    %     end
    
    
    Reflectivity  =  real(Greenland.ice_bed_power_frgc)+ 2*DN*((Greenland.depths-avg_depth)/1e3).*((Greenland.along_track_avg-mean(Greenland.along_track_avg))/1e3)+ 2*Na*((Greenland.depths-avg_depth)/1e3);
    Reflectivity_values = cat(2,Reflectivity_values, Reflectivity);
    if plots
        figure
        hist(Reflectivity,30)
       % keyboard
    end
%     L = 2*DN*((Greenland.depths-avg_depth)/1e3).*((Greenland.along_track_avg-mean(Greenland.along_track_avg))/1e3)+ 2*Na*((Greenland.depths-avg_depth)/1e3);
    L =  2*Na*((Greenland.depths-avg_depth)/1e3);

attenuation = cat(2, attenuation, L);
    power_rgc = cat(2,power_rgc, Greenland.ice_bed_power_rgc);
    
    %%
    if plots
        plot((1:length(Greenland.ice_bed_power_rgc)),-10*log10(Greenland.ice_bed_power_rgc));
        hold on
        plot(L);
    end
    %%
    if plots
        figure
        plot(Reflectivity)
    end
    %%
    %     s(:,1) = Reflectivity;
    %     s(:,2) = Greenland.Latitude;
    %     s(:,3) = Greenland.Longitude;
    %     [cidx, cmeans] = kmeans(s,5);
    %
    %     cluster_1 = find(cidx ==1) ;
    %     cluster_2 = find(cidx ==2) ;
    %     cluster_3 = find(cidx ==3) ;
    %     cluster_4 = find(cidx ==4) ;
    %     cluster_5 = find(cidx ==5) ;
    %%
    geotiff_fn = '/cresis/snfs1/dataproducts/GIS_data/greenland/Landsat-7/Greenland_natural.tif';
    proj = geotiffinfo(geotiff_fn);
    %proj = geotiffinfo('X:\GIS_data\antarctica\Landsat-7\Antarctica_LIMA_480m.tif');
    
    clear idx
    idx = find(isnan(Greenland.ice_bed_power)) ;
    Greenland.Latitude(idx) = [];
    Greenland.Longitude(idx) = [];
    
    [A CMAP R]= geotiffread(geotiff_fn);
    if M==1
        figure(1)
        mapshow(rgb2gray(A),CMAP/1e3);
        xlabel('X (km)');
        ylabel('Y (km)');
        xlim([350 650]);
        ylim([-1000 -750]);
        [gps.x,gps.y] = projfwd(proj,Greenland.Latitude_avg,Greenland.Longitude_avg);
        lt = cat(2,lt, Greenland.Latitude_avg ) ;
        ln =cat(2, ln,  Greenland.Longitude_avg) ;
        gps.x = gps.x / 1000;
        gps.y = gps.y / 1000;
        figure(1)
        hold on
        %         coordinates = [146.8025817802702,-80.8425624442354,0 146.859111420776,-80.8605329540687,0 146.9874283189143,-80.87954150883996,0 147.0461921738456,-80.88528971983799,0 147.085089323297,-80.8801308783134,0 147.1482434548151,-80.87775991516887,0 147.1917519299772,-80.87347125778013,0 147.3195670160135,-80.86543377833081,0 147.3645672165925,-80.85787825750721,0 147.4229423576138,-80.84396912935144,0 147.5030744704666,-80.82707270467108,0 147.5070159929441,-80.81896863200822,0 147.5289249162223,-80.80538987794964,0 147.5429653863865,-80.78682274867379,0 147.5569710274941,-80.77885283958237,0 147.5460622232309,-80.76975083868683,0 147.5398376431209,-80.76152454999017,0 147.4901883846538,-80.74866666613673,0 147.463077687949,-80.74179849742697,0 147.4275123289513,-80.7315606505201,0 147.419029972871,-80.72819215478877,0 147.4107047642798,-80.71424368428957,0 147.3945937648576,-80.70589185727766,0 147.3738988490754,-80.69666600228094,0 147.3686825136538,-80.68683355166068,0 147.3580841934058,-80.67774427376065,0 147.3515476296654,-80.66057795038542,0 147.3542032651798,-80.64516483544441,0 147.3665808675984,-80.64045185950934,0 147.3743640833768,-80.63486496040962,0 147.3592598762463,-80.62491050445701,0 147.3487753619517,-80.61583035428667,0 147.2877509501087,-80.58819596650911,0 147.2052854991595,-80.57407513320256,0 147.1304853422634,-80.56491419990755,0 147.0882601171962,-80.56025603767134,0 147.0381745209686,-80.56116475288647,0 146.9747390850922,-80.56837415652973,0 146.9317975728325,-80.57424990621185,0 146.8547588109255,-80.58774212394226,0 146.7999963612884,-80.59667738019473,0 146.7225690973808,-80.61014201919947,0 146.6440153255934,-80.62519963080351,0 146.5869513261213,-80.63732565316801,0 146.5460655379855,-80.64725918255678,0 146.4996423017032,-80.65791600074725,0 146.4640277277719,-80.66711390191517,0 146.4447345040051,-80.67412945820381,0 146.4101560364985,-80.69799150645325,0 146.376820077397,-80.71129504377265,0 146.371877699188,-80.71935771253766,0 146.3379298269023,-80.74161305248747,0 146.350022506786,-80.7548440006107,0 146.371226673235,-80.76985508722662,0 146.3919337788893,-80.77752434878342,0 146.4409640088229,-80.78890589975366,0 146.4966201240507,-80.7979370550982,0 146.5716590267267,-80.80889139304171,0 146.6275978593291,-80.81790800168395,0 146.6698617026529,-80.82426096200722,0 146.7453555743955,-80.83518467260404,0 146.8192317507791,-80.84933482957,0 146.8748444758108,-80.85994717658689,0 146.8025817802702,-80.8425624442354,0];
        %         id = find(coordinates == 0);
        %         coordinates(id)= [];
        %         lon_co = coordinates(1:2:end);
        %         lat_co = coordinates(2:2:end);
        %         [points_x,points_y] = projfwd(proj,lat_co,lon_co);
        %         figure(1)
        %         hold on
        %         points_x = points_x/ 1000;
        %         points_y= points_y / 1000;
        %         scatter(points_x,points_y,10,[0 0 0],'fill')
        
        %         coordinates = [148.4505322453265,-81.17933257473347,0 148.4931742813755,-81.18019774266614,0 148.5217297109509,-81.18047944558919,0 148.5502894307313,-81.18076134675792,0 148.5799056880266,-81.17840938554811,0 148.5829365834471,-81.17799842623899,0 148.6280222387953,-81.17270815019909,0 148.649889034022,-81.16807308629748,0 148.6989177430152,-81.1597303502373,0 148.7237338769781,-81.15468143380964,0 148.7314320951746,-81.14947181216883,0 148.7665542984649,-81.13967797479722,0 148.80062378504,-81.13251116571567,0 148.8212345277378,-81.12301873271431,0 148.8121349296995,-81.11677765994297,0 148.8337675698401,-81.11213448363932,0 148.8459491985473,-81.10212839047021,0 148.8421968768234,-81.09681809805123,0 148.8244592265887,-81.09094297537017,0 148.7984242429184,-81.08455231972782,0 148.7564705524482,-81.07537886416868,0 148.725570271266,-81.06718537330346,0 148.7106430848056,-81.06177583323142,0 148.6666239153455,-81.05126555332161,0 148.6494658016742,-81.04451878137461,0 148.646544367928,-81.03746773528665,0 148.6469405482719,-81.02913256570986,0 148.6472774783445,-81.02825801365908,0 148.6707582305907,-81.01838391479232,0 148.6842966412985,-81.01236591883354,0 148.6918872185219,-81.00717069907741,0 148.7011082385414,-81.0050619691237,0 148.7076928903233,-81.00249009749909,0 148.7086818243151,-80.99986683088764,0 148.7178916136264,-80.99775765274947,0 148.7205177925214,-80.9907626799708,0 148.7060657498585,-80.9844893882953,0 148.6922968278213,-80.97646806870195,0 148.6896740908893,-80.97600532508663,0 148.6778914964372,-80.97019545481462,0 148.6471413807994,-80.9628938636878,0 148.6213334637654,-80.95739231209059,0 148.5851273986036,-80.95004034284733,0 148.5701685675505,-80.94551645798869,0 148.5477623285003,-80.93873078242545,0 148.550849225377,-80.93087051678278,0 148.5366323981828,-80.92460124907186,0 148.5212227745638,-80.91437542116094,0 148.5241605772578,-80.91396536758802,0 148.5246689506209,-80.90564527394662,0 148.5250131937489,-80.90477229083083,0 148.538463540558,-80.89876626874184,0 148.5279253835284,-80.89034284894554,0 148.5031145362246,-80.88309707872655,0 148.4919381201438,-80.8764198251769,0 148.4532741751849,-80.86947708838994,0 148.4129104971972,-80.86689350673683,0 148.387638876286,-80.86795821085867,0 148.3620020669169,-80.86989095579983,0 148.3356281568136,-80.87356291225417,0 148.3110750334782,-80.8794419208873,0 148.2626137474595,-80.88945149933889,0 148.2511944218343,-80.89677702539868,0 148.2331168043176,-80.90666128580418,0 148.2165022391592,-80.91305726047854,0 148.1947151207824,-80.91852242243809,0 148.164039803794,-80.92520749420164,0 148.1447706617504,-80.93113631191804,0 148.1092594700618,-80.93601399524084,0 148.0715253264676,-80.93954856648864,0 148.038895984821,-80.94401006941044,0 148.0181151981174,-80.94684770848069,0 147.9816560789344,-80.95345327608189,0 147.9604303986233,-80.9571596098013,0 147.9383735968227,-80.96260879463416,0 147.9204973059076,-80.96503630789033,0 147.874528476377,-80.97371784106284,0 147.8348717195455,-80.98071452009171,0 147.8173997753759,-80.98796587884296,0 147.8063109888846,-80.99353765483389,0 147.8037765847634,-80.99877095988749,0 147.7845361472137,-81.00380986633186,0 147.7704312073306,-81.00978580359623,0 147.767880530109,-81.01502004619893,0 147.7743300033943,-81.01904422131744,0 147.7735057264172,-81.02649335194612,0 147.7726800181729,-81.03394327401304,0 147.7899017109125,-81.03897288407764,0 147.813638388826,-81.04802869360911,0 147.8442801226546,-81.05453144313663,0 147.8745762336822,-81.0619060824677,0 147.8975547386911,-81.06699932725357,0 147.9409932209669,-81.07671522117538,0 147.9692778657311,-81.08274371218738,0 147.9927856601197,-81.08696016013074,0 148.0046361160123,-81.09192478591433,0 148.0295942144455,-81.0992316852722,0 148.0624694291427,-81.10794327781278,0 148.0709653054862,-81.11418964015348,0 148.0847335887052,-81.12137332424396,0 148.0793005897722,-81.12702935947689,0 148.0673883044252,-81.13437392347107,0 148.0534339361046,-81.14609492334087,0 148.0615584130306,-81.15322189646216,0 148.0867387699101,-81.1605410852018,0 148.1396210976887,-81.16992616022368,0 148.1854636640034,-81.17614838694013,0 148.2268574058771,-81.17967263980717,0 148.288840706019,-81.18208777595341,0 148.339833835722,-81.18349738877274,0 148.3775108661184,-81.18255865393543,0 148.3991602894644,-81.17881447236616,0 148.4223577139819,-81.17816845950685,0 148.4505322453265,-81.17933257473347,0];
        %         id = find(coordinates == 0);
        %         coordinates(id)= [];
        %         lon_co = coordinates(1:2:end);
        %         lat_co = coordinates(2:2:end);
        %         [points_x,points_y] = projfwd(proj,lat_co,lon_co);
        %         figure(1)
        %         hold on
        %         points_x = points_x/ 1000;
        %         points_y= points_y / 1000;
        %         scatter(points_x,points_y,10,[0 0 0],'fill')
        
%         scatter(gps.x,gps.y,20,Reflectivity,'fill')
    scatter(gps.x,gps.y,20,L,'fill')

    else
        
        [gps.x,gps.y] = projfwd(proj,Greenland.Latitude_avg,Greenland.Longitude_avg);
        lt = cat(2,lt, Greenland.Latitude_avg ) ;
        ln =cat(2, ln,  Greenland.Longitude_avg) ;
        gps.x = gps.x / 1000;
        gps.y = gps.y / 1000;
        figure(1)
        hold on;
%         scatter(gps.x,gps.y,20,Reflectivity,'fill')
    scatter(gps.x,gps.y,20,L,'fill')

    end
    % colorbar
end
% figure(1)
% colorbar
for M = 1:20
    % for M = [14,29]
    clearvars -except M Na avg_power avg_depth Reflectivity_values attenuation power_rgc param ln lt 
    clc
    plots =0;
    param.radar.fs = 1000000000/9;
  load(['/cresis/snfs1/scratch/manjish/peterman/radar/verticalline', sprintf('%d.mat',M)]);

    % load(['/cresis/snfs1/scratch/santhosh/thesis/get heights frames Greenland/cross_lines/Greenland_layerdata_selected_frames_v6_' num2str(M,'%03d') '.mat'])
    physical_constants
    
    
    %     proj = geotiffinfo('X:\GIS_data\antarctica\Landsat-7\Antarctica_LIMA_480m.tif');
    %     [x,y] = projfwd(proj,Greenland.Latitude,Greenland.Longitude);
    %     along_track = [0 cumsum(sqrt(diff(x).^2 + diff(y).^2))]/1e3;
    
    %     along_tack_fn = geodetic_to_along_track(Greenland.Latitude,Greenland.Longitude,Greenland.Elevation);
    
    
    Greenland.depth = (Greenland.ice_bed_time - Greenland.surface_time)*c/2/sqrt(er_ice);
    Greenland.surface_height = (Greenland.surface_time)*c/2;
    
    %  [Greenland.depth_sorted index]= sort(Greenland.depth);
    %  Greenland.ice_bed_power_sorted = Greenland.ice_bed_power(index);
    % Greenland.ice_bed_power = Greenland.ice_bed_power_sorted;
    % Greenland.depth = Greenland.depth_sorted;
    % Greenland.surface_height_sorted = Greenland.surface_height(index);
    % Greenland.surface_height = Greenland.surface_height_sorted ;
    
    geometric_loss = (2*(Greenland.surface_height+Greenland.depth)).^2;
    geometric_loss_surface = (2*(Greenland.surface_height)).^2;
    Greenland.ice_bed_power_cgl =(Greenland.ice_bed_power).*geometric_loss;
    
    if plots
        plot(Greenland.depth, 10*log10(Greenland.ice_bed_power_cgl));
        grid on
    end
    %% compensating for surface roughness
    
    Greenland.ice_surface_power_cgl = (Greenland.ice_surface_power).*geometric_loss_surface;
    if exist((['/cresis/snfs1/scratch/manjish/peterman/radarnew/',sprintf('verticalline%d.mat',M)]),'file')
         load((['/cresis/snfs1/scratch/manjish/peterman/radarnew/',sprintf('verticalline%d.mat',M)]));
        K  = floor(length(Greenland.ice_surface_power_cgl)/1000);
        for l = 1:K
            
            if any(~isnan(Greenland.ice_bed_power_cgl(1+(l-1)*1000:(l*1000))))
                ice_bed_power_cgl = Greenland.ice_bed_power_cgl(1+(l-1)*1000:(l*1000));
                depth = Greenland.depth(1+(l-1)*1000:(l*1000));
                clear id
                id = find(isnan(Greenland.ice_bed_power_cgl(1+(l-1)*1000:(l*1000))));
                ice_bed_power_cgl(id) = [];
                depth(id) = [];
                Greenland.ice_bed_power_cgl_avg(l) =   mean(ice_bed_power_cgl) ;
                Greenland.depth_avg(l) = mean(depth);
                Greenland.Latitude_avg(l) = mean(Greenland.Latitude(1+(l-1)*1000:(l*1000)));
                Greenland.Longitude_avg(l) = mean(Greenland.Longitude(1+(l-1)*1000:(l*1000)));
            else
                Greenland.ice_bed_power_cgl_avg(l) = nan;
                Greenland.depth_avg(l) = nan;
                Greenland.Latitude_avg(l) = nan;
                Greenland.Longitude_avg(l) = nan;
            end
            
            if isnan(r.rms_height(l))
                continue;
            else
                Greenland.ice_bed_power_cgl_avg(l) = Greenland.ice_bed_power_cgl_avg(l)./(exp(-(4*pi*r.rms_height(l)/(c/param.radar.fs))^2)*(besseli(0,((4*pi*r.rms_height(l)/(c/param.radar.fs))^2)/2))^2);
            end
        end
        
        clearvars r  K
    else
        Greenland.ice_surface_power_cgl = (Greenland.ice_surface_power).*geometric_loss_surface;
        K  = floor(length(Greenland.ice_surface_power_cgl)/1000);
        for l = 1:K
            r.lat(l) = mean(Greenland.Latitude(1+(l-1)*1000:(l*1000)));
            r.lon(l) = mean(Greenland.Longitude(1+(l-1)*1000:(l*1000)));
            s = sqrt(Greenland.ice_surface_power_cgl(1+(l-1)*1000:(l*1000)));
            id = find(isnan(s)|isinf(s)|s==0);
            if length(id) > 500
                r.rms_height(l) = nan;
                r.dielectric_constant(l) = nan;
                r.pn(l) = nan;
                r.pc(l) = nan ;
                continue
            else
                s(id) = [];
            end
            pd = fitdist(double((s)).','Rician')
            % phat = mle(double(abs(s)),'distribution','Rician');
            % x = 0:0.0001:0.2;
            % histogram(abs(s))
            % h = hist((s));
            
            % A = pdf(pd,x);
            a = pd.s;
            % pc = 2*10*log10(a);
            % pn = 10*log10(2*pd.sigma^2);
            S = pd.sigma;
            r.pc(l) = a^2;
            r.pn(l) = 2*2*pd.sigma^2;
            rms_fit = (r.pc(l)/r.pn(l))*4*(2*pi/(c/param.radar.fs))^2;
            
            r.rms_height(l) = 0.0001;
            clear MSE
            for i = 1:5000
                MSE(i) = abs(rms_fit - exp(-(2*(2*pi/(c/param.radar.fs))*r.rms_height(l))^2)/((r.rms_height(l))^2));
                if r.rms_height(l) > 0.40
                    r.rms_height(l) = nan;
                    warning('check this')
                    %                         keyboard
                    break
                else
                    
                    if i>1
                        if MSE(i-1) < MSE(i)
                            break
                        else
                            r.rms_height(l) = r.rms_height(l) + 0.0001;
                            continue  ;
                        end
                    else
                        r.rms_height(l) = r.rms_height(l) + 0.0001;
                    end
                end
            end
            
            
            if isnan(r.rms_height(l))
                r.dielectric_constant(l) = nan;
                continue
            else
                r.dielectric_constant(l) = 1;
                clear mse
                for i = 1:5000
                    mse(i) = abs(r.pc(l) - ((1-sqrt(r.dielectric_constant(l)))/((1+sqrt(r.dielectric_constant(l)))))^2*exp(-(2*(2*pi/(c/param.radar.fs))*r.rms_height(l))^2));
                    if r.dielectric_constant(l) > 4
                        r.dielectric_constant(l) = nan;
                        warning('check this')
                        %                         keyboard
                        break
                    else
                        
                        if i>1
                            if mse(i-1) < mse(i)
                                break
                            else
                                r.dielectric_constant(l) = r.dielectric_constant(l) + 0.01;
                                continue  ;
                            end
                        else
                            r.dielectric_constant(l) = r.dielectric_constant(l) + 0.01;
                        end
                    end
                end
            end
            
            
            if any(~isnan(Greenland.ice_bed_power_cgl(1+(l-1)*1000:(l*1000))))
                ice_bed_power_cgl = Greenland.ice_bed_power_cgl(1+(l-1)*1000:(l*1000));
                depth = Greenland.depth(1+(l-1)*1000:(l*1000));
                clear id
                id = find(isnan(Greenland.ice_bed_power_cgl(1+(l-1)*1000:(l*1000))));
                ice_bed_power_cgl(id) = [];
                depth(id) = [];
                Greenland.ice_bed_power_cgl_avg(l) =   mean(ice_bed_power_cgl) ;
                Greenland.depth_avg(l) = mean(depth);
                Greenland.Latitude_avg(l) = mean(Greenland.Latitude(1+(l-1)*1000:(l*1000)));
                Greenland.Longitude_avg(l) = mean(Greenland.Longitude(1+(l-1)*1000:(l*1000)));
            else
                Greenland.ice_bed_power_cgl_avg(l) = nan;
                Greenland.depth_avg(l) = nan;
                Greenland.Latitude_avg(l) = nan;
                Greenland.Longitude_avg(l) = nan;
            end
            
            if isnan(r.rms_height(l))
                continue;
            else
                Greenland.ice_bed_power_cgl_avg(l) = Greenland.ice_bed_power_cgl_avg(l)./(exp(-(4*pi*r.rms_height(l)/(c/param.radar.fs))^2)*(besseli(0,((4*pi*r.rms_height(l)/(c/param.radar.fs))^2)/2))^2);
            end
        end
        
        
        save((['/cresis/snfs1/dataproducts/ct_data/rds/2012_Greenland_P3/CSARP_manjish/surfaceroughness_after/',param.day_seg,'/', sprintf('IceBedCoherenceIndex_%s.mat', param.day_seg)]),'r');
        clearvars r K
    end
    
    
    
    
    
    %% compensating for bed roughness
     if exist((['/cresis/snfs1/scratch/manjish/peterman/bedroughness/',sprintf('verticalline%d.mat',M)]),'file')
         load((['/cresis/snfs1/scratch/manjish/peterman/bedroughness/',sprintf('verticalline%d.mat',M)]));
         r=rbed;
        K  = floor(length(Greenland.ice_bed_power_cgl)/1000);
        for l = 1:K
            
            if isnan(r.rms_height(l)) || isnan(Greenland.ice_bed_power_cgl_avg(l))
                continue ;
            else
                Greenland.ice_bed_power_cgl_avg(l) = Greenland.ice_bed_power_cgl_avg(l)./(exp(-(4*pi*r.rms_height(l)/(c/param.radar.fs))^2)*(besseli(0,((4*pi*r.rms_height(l)/(c/param.radar.fs))^2)/2))^2);
            end
        end
        
        %           Greenland.ice_bed_power_cgl(1+(l-1)*1000:(l*1000)) = (Greenland.ice_bed_power_cgl(1+(l-1)*1000:(l*1000)))./(exp(-(4*pi*r.rms_height(l)/(c/param.radar.fs))^2)*(besseli(0,((4*pi*r.rms_height(l)/(c/param.radar.fs))^2)/2))^2);
        
        
        
        clearvars r K
    else
        
        K  = floor(length(Greenland.ice_bed_power_cgl)/1000);
        for l = 1:K
            r.lat(l) = mean(Greenland.Latitude(1+(l-1)*1000:(l*1000)));
            r.lon(l) = mean(Greenland.Longitude(1+(l-1)*1000:(l*1000)));
            s = sqrt(Greenland.ice_bed_power_cgl(1+(l-1)*1000:(l*1000)));
            id = find(isnan(s)|isinf(s)|s==0);
            if length(id) > 500
                r.rms_height(l) = nan;
                r.dielectric_constant(l) = nan;
                r.pn(l) = nan;
                r.pc(l) = nan ;
                continue
            else
                s(id) = [];
            end
            pd = fitdist(double((s)).','Rician')
            % phat = mle(double(abs(s)),'distribution','Rician');
            % x = 0:0.0001:0.2;
            % histogram(abs(s))
            % h = hist((s));
            
            % A = pdf(pd,x);
            a = pd.s;
            % pc = 2*10*log10(a);
            % pn = 10*log10(2*pd.sigma^2);
            S = pd.sigma;
            r.pc(l) = a^2;
            r.pn(l) = 2*2*pd.sigma^2;
            rms_fit = (r.pc(l)/r.pn(l))*4*(2*pi/(c/param.radar.fs))^2;
            
            r.rms_height(l) = 0.0001;
            clear MSE
            for i = 1:5000
                MSE(i) = abs(rms_fit - exp(-(2*(2*pi/(c/param.radar.fs))*r.rms_height(l))^2)/((r.rms_height(l))^2));
                if r.rms_height(l) > 0.40
                    r.rms_height(l) = nan;
                    warning('check this')
                    %                         keyboard
                    break
                else
                    
                    if i>1
                        if MSE(i-1) < MSE(i)
                            break
                        else
                            r.rms_height(l) = r.rms_height(l) + 0.0001;
                            continue  ;
                        end
                    else
                        r.rms_height(l) = r.rms_height(l) + 0.0001;
                    end
                end
            end
            
            
            if isnan(r.rms_height(l))
                r.dielectric_constant(l) = nan;
                continue
            else
                r.dielectric_constant(l) = 1;
                clear mse
                for i = 1:5000
                    mse(i) = abs(r.pc(l) - ((1-sqrt(r.dielectric_constant(l)))/((1+sqrt(r.dielectric_constant(l)))))^2*exp(-(2*(2*pi/(c/param.radar.fs))*r.rms_height(l))^2));
                    if r.dielectric_constant(l) > 4
                        r.dielectric_constant(l) = nan;
                        warning('check this')
                        %                         keyboard
                        break
                    else
                        
                        if i>1
                            if mse(i-1) < mse(i)
                                break
                            else
                                r.dielectric_constant(l) = r.dielectric_constant(l) + 0.01;
                                continue  ;
                            end
                        else
                            r.dielectric_constant(l) = r.dielectric_constant(l) + 0.01;
                        end
                    end
                end
            end
            
            
            %        Greenland.ice_bed_power_cgl(1+(l-1)*1000:(l*1000)) = (Greenland.ice_bed_power_cgl(1+(l-1)*1000:(l*1000)))./(exp(-(4*pi*r.rms_height(l)/(c/param.radar.fs))^2)*(besseli(0,((4*pi*r.rms_height(l)/(c/param.radar.fs))^2)/2))^2);
            
            
            if any(~isnan(Greenland.ice_bed_power_cgl(1+(l-1)*1000:(l*1000))))
                ice_bed_power_cgl = Greenland.ice_bed_power_cgl(1+(l-1)*1000:(l*1000));
                depth = Greenland.depth(1+(l-1)*1000:(l*1000));
                clear id
                id = find(isnan(Greenland.ice_bed_power_cgl(1+(l-1)*1000:(l*1000))));
                ice_bed_power_cgl(id) = [];
                depth(id) = [];
                Greenland.ice_bed_power_cgl_avg(l) =   mean(ice_bed_power_cgl) ;
                Greenland.depth_avg(l) = mean(depth);
                Greenland.Latitude_avg(l) = mean(Greenland.Latitude(1+(l-1)*1000:(l*1000)));
                Greenland.Longitude_avg(l) = mean(Greenland.Longitude(1+(l-1)*1000:(l*1000)));
            else
                Greenland.ice_bed_power_cgl_avg(l) = nan;
                Greenland.depth_avg(l) = nan;
                Greenland.Latitude_avg(l) = nan;
                Greenland.Longitude_avg(l) = nan;
                Greenland.Latitude(1+(l-1)*1000:(l*1000)) =[];
                Greenland.Longitude(1+(l-1)*1000:(l*1000))= [];
            end
            
            
            
            
            
            if isnan(r.rms_height(l))
                continue;
            else
                Greenland.ice_bed_power_cgl_avg(l) = Greenland.ice_bed_power_cgl_avg(l)./(exp(-(4*pi*r.rms_height(l)/(c/param.radar.fs))^2)*(besseli(0,((4*pi*r.rms_height(l)/(c/param.radar.fs))^2)/2))^2);
            end
        end
        
         save((['/cresis/snfs1/dataproducts/ct_data/rds/2012_Greenland_P3/CSARP_manjish/bedroughness_after/',param.day_seg,'/', sprintf('IceBedCoherenceIndex_%s.mat', param.day_seg)]),'r');
       clearvars r K
        
    end
    
    
    %%  relative geometrically corrected bed - echo power
    %     for i = 1:length(Greenland.ice_bed_power_cgl)
    %         if i < 3
    %             Greenland.ice_bed_power_rgc(i) = Greenland.ice_bed_power_cgl(i)-mean(Greenland.ice_bed_power_cgl(1:6-i));
    %         elseif i+2 > length(Greenland.ice_bed_power_cgl)
    %             Greenland.ice_bed_power_rgc(i) = Greenland.ice_bed_power_cgl(i)- mean(Greenland.ice_bed_power_cgl(i-5:end));
    %         else
    %             Greenland.ice_bed_power_rgc(i) = Greenland.ice_bed_power_cgl(i)- mean(Greenland.ice_bed_power_cgl(i-2:i+2));
    %         end
    %     end
    
    
    
    
    Greenland.ice_bed_power_rgc = Greenland.ice_bed_power_cgl_avg;%/avg_power;
    if plots
        plot(Greenland.depth_avg,10*log10(Greenland.ice_bed_power_rgc));
        grid on;
    end
    %% assuming constant attenuation rate
    %     if plots
    %         N = [2 10];
    %         for k = 1: 2
    %             Na  = N(k); % 3 - 20dB attenuation rate
    %             Greenland.ice_bed_reflectivity = lp(Greenland.ice_bed_power_rgc)+ (2*Na*(Greenland.depth-mean(Greenland.depth))/1e3);
    %
    %             figure(1)
    %             subplot(2,1,k)
    %             plot(Greenland.ice_bed_reflectivity,'*-');
    %             figure(2)
    %             subplot(2,1,k)
    %             hist(Greenland.ice_bed_reflectivity,40)
    %         end
    %     end
    %%
    if plots
        plot((1:length(Greenland.ice_bed_power_rgc)),-10*log10(Greenland.ice_bed_power_rgc));
        figure
        plot(1:length(Greenland.depth_avg),(Greenland.depth_avg/1e3))
    end
    %%
    
    %     clear idx
    idx = find(isnan(Greenland.ice_bed_power_rgc)) ;
    Greenland.ice_bed_power_rgc(idx) = [];
    Greenland.depth_avg(idx) =[];
    Greenland.Latitude_avg(idx) =[];
    Greenland.Longitude_avg(idx) =[];
    %
    %% applying a guassian filter
    %   h = fspecial('gaussian', [1,round(length(Greenland.ice_bed_power_rgc)/5)], 300);
    %  % h = h./max(h); % normalizing
    % %  Greenland.ice_bed_power_frgc = filter(h,1,Greenland.ice_bed_power_rgc);
    %  Greenland.ice_bed_power_frgc = conv(Greenland.ice_bed_power_rgc,h,'same');
    if length(Greenland.ice_bed_power_rgc)< 500
        if mod(length(Greenland.ice_bed_power_rgc),2)==0
            window_size =length(Greenland.ice_bed_power_rgc)-100;
        else
            window_size =length(Greenland.ice_bed_power_rgc)-99;
        end
    else
        window_size = 500;
    end
     if length(Greenland.ice_bed_power_rgc)<100
      window_size=2;
    end
    % Greenland.ice_bed_power_frgc = conv(Greenland.ice_bed_power_rgc,gausswin(window_size),'same');
    Greenland.ice_bed_power_frgc  = sgolayfilt(10*log10(Greenland.ice_bed_power_rgc), 2,window_size+1, gausswin(window_size+1));
    if plots
        plot((1:length(Greenland.ice_bed_power_rgc)),10*log10(Greenland.ice_bed_power_rgc));
        hold on
        plot((1:length(Greenland.ice_bed_power_frgc)),(Greenland.ice_bed_power_frgc));
    end
    if length(Greenland.ice_bed_power_rgc)< 500
        if mod(length(Greenland.ice_bed_power_rgc),2)==0
            window_size =length(Greenland.ice_bed_power_rgc)-100;
        else
            window_size =length(Greenland.ice_bed_power_rgc)-99;
        end
    else
        window_size = 500;
    end
     if length(Greenland.ice_bed_power_rgc)<100
      window_size=2;
    end
    Greenland.depths  = sgolayfilt(Greenland.depth_avg, 2,window_size+1, gausswin(window_size+1));
    if plots
        plot((1:length(Greenland.depth_avg)),Greenland.depth_avg);
        hold on
        plot((1:length(Greenland.depths)),(Greenland.depths));
    end
    %     %
    %     %  for i = 1:length(Greenland.depth)
    %     %      if i < 3
    %     %
    %     %          Greenland.dif_depth(i) = Greenland.depths(i)-mean(Greenland.depths(1:6-i));
    %     %
    %     %     elseif i+2 > length(Greenland.depth)
    %     %     Greenland.dif_depth(i) = Greenland.depths(i)- mean(Greenland.depths(i-5:end));
    %     %     else
    %     %       Greenland.dif_depth(i) = Greenland.depths(i)- mean(Greenland.depths(i-2:i+2));
    %     %      end
    %     %  end
    %     %
    %     % %  Greenland.fdepth  = sgolayfilt(Greenland.dif_depth, 2,window_size+1, gausswin(window_size+1));
    %     %   plot((1:length(Greenland.dif_depth)),(10*log10(Greenland.dif_depth)));
    %     % %    hold on
    %     % % plot((1:length(Greenland.ice_bed_power_rgc)),10*log10(-Greenland.ice_bed_power_rgc)/2);
    %     % % %  hold on
    %     % % %  plot((1:length(Greenland.fdepth)),10*log10(Greenland.fdepth));
    %
    %     %%
    %     %     na =(-(Greenland.ice_bed_power_frgc)- 10*log10(2*Greenland.dif_depth));
    %     %   plot([1:length(na)],na)
    %     % % Na  = mean(na);
    %     %  for na =1:40
    %     %      mse(na) = sum(abs(real(Greenland.ice_bed_power_frgc) + (2*na*((Greenland.depths-mean(Greenland.depths))/1e3))));
    %     %  end
    %     % %  [v i] =min(mse)
    %     % sum(abs((2*12*(Greenland.depths-mean(Greenland.depths))/1e3)+(10*log10(Greenland.ice_bed_power_frgc))));
    %     for N= 1:40
    %          S(N) = mean(abs((2*N*(Greenland.depths-mean(Greenland.depths))/1e3)+(real(Greenland.ice_bed_power_frgc))).^2);
    %  %       S(N) = mode(round(abs((2*N*(Greenland.depths-mean(Greenland.depths))/1e3)+(10*log10(Greenland.ice_bed_power_frgc)))));
    %  %        S(N) = sum(abs((2*N*(Greenland.depths-mean(Greenland.depths))/1e3)+(real(Greenland.ice_bed_power_frgc))).^2);
    %     end
    %     [v i] = min(S)
    %     Na = i;
    %%
    if plots
        figure
        plot((1:length(Greenland.depths)),2*Na*((Greenland.depths-avg_depth)/1e3));
        hold on
        plot((1:length(Greenland.ice_bed_power_frgc)),(-(Greenland.ice_bed_power_frgc)));
    end
    %% to find d(na)/dx
    %     L = cumsum(Greenland.segments_length);
    %     window_size = 36;
    %     debug_flag= 1;
    %     for i = 1: length(L)
    %         if i == 1;
    %             p=1;
    %             q = L(i);
    %         else
    %             p = 1+L(i-1);
    %             q = L(i);
    %         end
    %         along_track_fn{i}= geodetic_to_along_track(Greenland.Latitude(p:q),Greenland.Longitude(p:q));
    %         power{i} = Greenland.ice_bed_power_rgc(p:q);
    %         depth{i} = Greenland.depth(p:q);
    %         depths{i}  = sgolayfilt(depth{i}, 2,window_size+1, gausswin(window_size+1));
    %         power_frgc{i}  = sgolayfilt(10*log10(power{i}), 2,window_size+1, gausswin(window_size+1));
    %         if debug_flag
    %             %             plot((1:length(power{i})),10*log10(power{i}));
    %             %             hold on
    %             %             plot((1:length(power{i})),(power_frgc{i}));
    %             %             keyboard
    %             %             clf
    %             %             figure
    %             %             plot((1:length(depth{i})),depth{i});
    %             %             hold on
    %             %             plot((1:length(depths{i})),(depths{i}));
    %             %             keyboard
    %             %             clf
    %             if plots
    %                         subplot(2,1,1)
    %                         plot((1:length(power{i})),(-power_frgc{i}));
    %                         subplot(2,1,2)
    %                         plot((1:length(depths{i})),(depths{i}-mean(depths{i}))/1e3);
    %                         keyboard
    %                         clf
    %             end
    %         end
    %     end
    L = cumsum(Greenland.segments_length);
    debug_flag= 1;
    Greenland.along_track = [];
    for i = 1: length(L)
        if i == 1;
            p=1;
            q = L(i);
            Lat = Greenland.Latitude(p:q);
            Lon = Greenland.Longitude(p:q);
            
            %             idx = find(isnan(Greenland.ice_bed_power(p:q)));
            %             Lat(idx) = [];
            %             Lon(idx) = [];
            along_track_fn{i}= geodetic_to_along_track(Lat,Lon);
            Greenland.along_track = cat(2, Greenland.along_track, along_track_fn{i}) ;
        else
            p = 1+L(i-1);
            q = L(i);
            %             idx = find(isnan(Greenland.ice_bed_power(p:q)));
            Lat = Greenland.Latitude(p:q);
            Lon = Greenland.Longitude(p:q);
            %             Lat(idx) = [];
            %             Lon(idx) = [];
            along_track_fn{i}= geodetic_to_along_track(Lat,Lon);
            Greenland.along_track = cat(2, Greenland.along_track, Greenland.along_track(end)+ along_track_fn{i}) ;
        end
    end
    
    
    K  = floor(length(Greenland.ice_bed_power_cgl)/1000);
    for l = 1:K
        if any(~isnan(Greenland.along_track(1+(l-1)*1000:(l*1000))))
            along_track = Greenland.along_track(1+(l-1)*1000:(l*1000));
            clear id
            id = find(isnan(along_track));
            along_track(id) = [];
            Greenland.along_track_avg(l) =   mean(along_track) ;
        else
            Greenland.along_track_avg(l) = nan;
            continue ;
        end
    end
    
    Greenland.along_track_avg(idx) = [];
    
    if ~(length(Greenland.depth_avg) == length(Greenland.along_track_avg))
        keyboard;
    end
    
    
    
    
    %%
    %     dn = (1:40);
    %     for i = 1: length(power_frgc)
    %         for j = 1:length(dn)
    %
    %             term_1 = 2*dn(j)*((depths{i}-mean(depths{i}))/1e3).*((along_track_fn{i}-mean(along_track_fn{i}))/1e3);
    %             term_2 = 2*Na*((depths{i}-mean(depths{i}))/1e3);
    %             S(i,j) = sum(abs(real(power_frgc{i})+term_1+term_2).^2);
    %             a = S(i,j);
    %             if plots
    %                     plot(real(-power_frgc{i}))
    %                     hold on
    %                     plot(term_1+term_2)
    %                     keyboard
    %                     clf
    %             end
    %             if j > 1
    %                 if S(i,j-1)<S(i,j)
    %                     DN(i) = dn(j-1);
    %                     break
    %                 end
    %             end
    %         end
    %     end
    % end
    % % s = 0;
    % %  for i = 1:length(along_track_fn)
    % % s= s+along_track_fn{i}(end);
    % % end
    
    dn = (-0.75:0.01:0.75);
    for j = 1:length(dn)
        
        term_1 = 2*dn(j).*((Greenland.depths-avg_depth)/1e3).*((Greenland.along_track_avg-mean(Greenland.along_track_avg))/1e3);
        term_2 = 2*Na*((Greenland.depths-avg_depth)/1e3);
        S(j) = sum(abs(real(Greenland.ice_bed_power_frgc)+term_1+term_2).^2);
        if plots
            plot(real(-Greenland.ice_bed_power_frgc))
            hold on
            plot(term_1+term_2)
            keyboard
            clf
        end
        if j > 1
            if S(j-1)<S(j)
                DN = dn(j-1);
                break
            end
        end
    end
    
    if j ==length(dn)
        [v id] = min(S);
        DN = dn(id);
    end
    
    
    
    %%
    %     Reflectivity = [];
    %     for i = 1:length(power_frgc)
    %         r  =  real(power_frgc{i})+ 2*DN(i)*((depths{i}-mean(depths{i}))/1e3).*((along_track_fn{i}-mean(along_track_fn{i}))/1e3)+ 2*Na*((depths{i}-mean(depths{i}))/1e3);
    %         Reflectivity = cat(2, Reflectivity, r);
    %     end
    %     if plots
    %         figure
    %         hist(Reflectivity,30)
    %     end
    
    
    Reflectivity  =  real(Greenland.ice_bed_power_frgc)+ 2*DN*((Greenland.depths-avg_depth)/1e3).*((Greenland.along_track_avg-mean(Greenland.along_track_avg))/1e3)+ 2*Na*((Greenland.depths-avg_depth)/1e3);
    Reflectivity_values = cat(2,Reflectivity_values, Reflectivity);
    if plots
        figure
        hist(Reflectivity,30)
        keyboard
    end
%   L = 2*DN*((Greenland.depths-avg_depth)/1e3).*((Greenland.along_track_avg-mean(Greenland.along_track_avg))/1e3)+ 2*Na*((Greenland.depths-avg_depth)/1e3);
    L =  2*Na*((Greenland.depths-avg_depth)/1e3);
    attenuation = cat(2, attenuation, L);
    power_rgc = cat(2,power_rgc, Greenland.ice_bed_power_rgc);
    
    %%ln
    if plots
        plot((1:length(Greenland.ice_bed_power_rgc)),-10*log10(Greenland.ice_bed_power_rgc));
        hold on
        plot(L);
    end
    %%
    if plots
        figure
        plot(Reflectivity)
    end
    %     S =  [Reflectivity lt ln]; 
    %     [cidx, cmeans] = kmeans(s,5);
    %
    %     cluster_1 = find(cidx ==1) ;
    %     cluster_2 = find(cidx ==2) ;
    %     cluster_3 = find(cidx ==3) ;
    %     cluster_4 = find(cidx ==4) ;
    %     cluster_5 = find(cidx ==5) ;
    %%
    geotiff_fn = '/cresis/snfs1/dataproducts/GIS_data/greenland/Landsat-7/Greenland_natural.tif';
    proj = geotiffinfo(geotiff_fn);
    %proj = geotiffinfo('X:\GIS_data\antarctica\Landsat-7\Antarctica_LIMA_480m.tif');
    [A CMAP R]= geotiffread(geotiff_fn);
    %     if M==1
    %         figure(2)
    %         mapshow(rgb2gray(A),CMAP/1e3);
    %         xlabel('X (km)');
    %         ylabel('Y (km)');
    %         xlim([350 650]);
    %         ylim([-1000 -700]);
    %         [gps.x,gps.y] = projfwd(proj,Greenland.Latitude,Greenland.Longitude);
    %         gps.x = gps.x / 1000;
    %         gps.y = gps.y / 1000;
    %         figure(2)
    %         hold on
    %         scatter(gps.x,gps.y,20,Reflectivity,'fill')
    %     else
    %         figure(1)
    %         hold on
    clear idx
    idx = find(isnan(Greenland.ice_bed_power)) ;
    Greenland.Latitude(idx) = [];
    Greenland.Longitude(idx) = [];
    
    
    [gps.x,gps.y] = projfwd(proj,Greenland.Latitude_avg,Greenland.Longitude_avg);
    lt = cat(2,lt, Greenland.Latitude_avg ) ;
    ln =cat(2, ln,  Greenland.Longitude_avg) ;
    gps.x = gps.x / 1000;
    gps.y = gps.y / 1000;
    figure(1)
    hold on;
%     scatter(gps.x,gps.y,20,Reflectivity,'fill')
    scatter(gps.x,gps.y,20,L,'fill')

end

%  S =  [Reflectivity_values.']; 
%  [cidx, cmeans] = kmeans(S,4);
%         cluster_1 = find(cidx ==1) ;
%         cluster_2 = find(cidx ==2) ;
%         cluster_3 = find(cidx ==3) ;
%         cluster_4 = find(cidx ==4) ;
% % colorbar
% end
%figure(1)
% colorbar
%
geotiff_fn = '/cresis/snfs1/dataproducts/GIS_data/greenland/Landsat-7/Greenland_natural.tif';
% geotiff_fn = 'X:\GIS_data\antarctica\Landsat-7\Antarctica_LIMA_480m.tif';
proj = geotiffinfo(geotiff_fn);
%proj = geotiffinfo('X:\GIS_data\antarctica\Landsat-7\Antarctica_LIMA_480m.tif');
[A CMAP R]= geotiffread(geotiff_fn);
figure(2)
hold on;
mapshow(rgb2gray(A),CMAP/1e3);
xlabel('X (km)');
ylabel('Y (km)');
%  xlim([350 650]);
% ylim([-1000 -750]);

hold on
[gps.x,gps.y] = projfwd(proj,lt,ln);
gps.x = gps.x / 1000;
gps.y = gps.y / 1000;
scatter(gps.x,gps.y,20,Reflectivity_values,'fill')
 colorbar
% caxis([min(lp(brp(cluster_1))), max(lp(brp(cluster_1)))]);
grid
title('reflected bed power');
xlabel('Latitude');
ylabel('Longitude');

%% Attenuation
figure(3)
hold on;
mapshow(rgb2gray(A),CMAP/1e3);
xlabel('X (km)');
ylabel('Y (km)');
%  xlim([350 650]);
% ylim([-1000 -750]);

hold on
[gps.x,gps.y] = projfwd(proj,lt,ln);
gps.x = gps.x / 1000;
gps.y = gps.y / 1000;
scatter(gps.x,gps.y,20,attenuation,'fill')
 colorbar
% caxis([min(lp(brp(cluster_1))), max(lp(brp(cluster_1)))]);
grid
title('Attenuation');
xlabel('Latitude');
ylabel('Longitude');