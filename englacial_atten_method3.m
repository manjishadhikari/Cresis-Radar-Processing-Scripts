%Script: estimation_of_relative_reflectivity_values using englacial atten
%method 2.. Set Na from method 1 and fit DN for every line

%% setup
clear
close
clc
dbstop error

plots =0;
ice_bed_power_G_r_corrected = [];
lat_G_r_corrected = [];
lon_G_r_corrected = [];
depth_G_r_corrected = [];
constant_attenuation = [];
estimated_Na = [];
estimated_Na_mod=[];
var_attn=[];

%% loading the data
cross_lines = 1;

%%

for M =1:19
  
  clearvars -except M plots ice_bed_power_G_r_corrected lat_G_r_corrected lon_G_r_corrected depth_G_r_corrected cross_lines constant_attenuation estimated_Na var_attn estimated_Na_mod
  clc
  
  param.radar.fs = 195000000;
  if M<21
    load(['/cresis/snfs1/scratch/manjish/peterman/radar_w_index/crossline',num2str(M)]);
  else
    M1=M-20;
    load(['/cresis/snfs1/scratch/manjish/peterman/radar_w_index/verticalline',num2str(M1)]);
  end
  physical_constants
  %%
  
  Greenland.depth = (Greenland.ice_bed_time - Greenland.surface_time)*c/2/sqrt(er_ice);
  Greenland.surface_height = (Greenland.surface_time)*c/2;
  
  
  geometric_loss = (2*(Greenland.surface_height+Greenland.depth/sqrt(er_ice))).^2;
  geometric_loss_surface = (2*(Greenland.surface_height)).^2;
  
  if plots
    figure(1);plot(Greenland.depth, lp(Greenland.ice_bed_power));
    grid on; title('Depth vs Power')
    figure(2); plot( lp(Greenland.ice_bed_power));
    grid on; title('Along Track vs Power')
  end
  %% compensating reflected bed power for surface roughness
  file_exist = false;
  if M<21
    if exist((  (['/cresis/snfs1/scratch/manjish/peterman/radarnew/crossline',num2str(M),'.mat'])),'file')
      load(['/cresis/snfs1/scratch/manjish/peterman/radarnew/crossline',num2str(M),'.mat']);
      file_exist = true;
    end
  else
    
    if exist((  (['/cresis/snfs1/scratch/manjish/peterman/radarnew/verticalline',num2str(M1),'.mat'])),'file')
      load(['/cresis/snfs1/scratch/manjish/peterman/radarnew/verticalline',num2str(M1),'.mat']);
      file_exist = true;
    end
  end
  
  %         K  = floor(length(Greenland.ice_surface_power)/1000);
  if file_exist
   k = 1;
    %         for l = 501:250:length(Greenland.ice_surface_power)
    %             if ((l > 500) && ((l+500) < length(Greenland.ice_surface_power)))
    num_int=r.settings.num_int;
    repeat_after=300;
    for l = num_int/2:repeat_after:length(Greenland.ice_surface_power)
      if ((l >= num_int/2) && ((l+num_int/2) < length(Greenland.ice_surface_power)))
        
        if any(~isnan(Greenland.ice_bed_power((l-num_int/2+1):(l+num_int/2))))
          ice_bed_power = Greenland.ice_bed_power((l-num_int/2+1):(l+num_int/2));
          depth = Greenland.depth((l-num_int/2+1):(l+num_int/2));
          clear id
          id = find(isnan(Greenland.ice_bed_power((l-num_int/2+1):(l+num_int/2))));
          ice_bed_power(id) = [];
          depth(id) = [];
          Greenland.ice_bed_power_avg(k) =   nanmean(ice_bed_power) ;
          Greenland.depth_avg(k) = nanmean(depth);
          Greenland.Latitude_avg(k) = nanmean(Greenland.Latitude((l-num_int/2+1):(l+num_int/2)));
          Greenland.Longitude_avg(k) = nanmean(Greenland.Longitude((l-num_int/2+1):(l+num_int/2)));
          Greenland.along_track_avg(k) = nanmedian(geodetic_to_along_track(Greenland.Latitude((l-num_int/2+1):(l+num_int/2)),Greenland.Longitude((l-num_int/2+1):(l+num_int/2))));
          Greenland.geometric_loss_avg(k) = nanmean(geometric_loss((l-num_int/2+1):(l+num_int/2)));
          
          
        else
          Greenland.ice_bed_power_avg(k) = nan;
          Greenland.depth_avg(k) = nan;
          Greenland.Latitude_avg(k) = nan;
          Greenland.Longitude_avg(k) = nan;
          Greenland.along_track_avg(k) = nan;
          Greenland.geometric_loss_avg(k) = nan;
          
        end
        if k>size(r.lat)
          continue
        end
        if isnan(r.rms_height(k))
          k= k+1;
          continue;
        else
          Greenland.ice_bed_power_avg(k) = Greenland.ice_bed_power_avg(k)./(exp(-((4*pi*r.rms_height(k)/(c/param.radar.fs))*(sqrt(r.dielectric_constant(k))-1))^2)*(besseli(0,(((4*pi*r.rms_height(k)/(c/param.radar.fs))*(sqrt(r.dielectric_constant(k))-1))^2)/2))^2);
          k= k+1;
        end
      end
    end
    
    clearvars r  k
  else
    %         K  = floor(length(Greenland.ice_surface_power)/1000);
    k = 1;
    for l = 501:250:length(Greenland.ice_surface_power)
      if ((l > 500) && ((l+500) < length(Greenland.ice_surface_power)))
        
        r.lat(k) = nanmean(Greenland.Latitude((l-500):(l+499)));
        r.lon(k) = nanmean(Greenland.Longitude((l-500):(l+499)));
        s = abs(Greenland.ice_surface_power((l-500):(l+499)));
        id = find(isnan(s)|isinf(s)|s==0);
        if length(id) > 500
          r.rms_height(k) = nan;
          r.dielectric_constant(k) = nan;
          r.pn(k) = nan;
          r.pc(k) = nan ;
          k= k+1;
          continue
        else
          s(id) = [];
        end
        %pd = fitdist(double((s)).','Rician');
        [pd.s, pd.sigma]=ricefit(s);
        % phat = mle(double(abs(s)),'distribution','Rician');
        % x = 0:0.0001:0.2;
        % histogram(abs(s))
        % h = hist((s));
        
        % A = pdf(pd,x);
        a = pd.s;
        % pc = 2*10*log10(a);
        % pn = 10*log10(2*pd.sigma^2);
        S = pd.sigma;
        r.pc(k) = a^2;
        r.pn(k) = 2*2*pd.sigma^2;
        rms_fit = (r.pc(k)/r.pn(k))*4*(2*pi/(c/param.radar.fs))^2;
        
        r.rms_height(k) = 0.001;
        clear MSE
        for i = 1:5000
          MSE(i) = abs(rms_fit - exp(-(2*(2*pi/(c/param.radar.fs))*r.rms_height(k))^2)/((r.rms_height(k))^2));
          if r.rms_height(k) > 0.40
            r.rms_height(k) = nan;
            warning('check this')
            %                         keyboard
            break
          else
            
            if i>1
              if MSE(i-1) < MSE(i)
                break
              else
                r.rms_height(k) = r.rms_height(k) + 0.0001;
                continue  ;
              end
            else
              r.rms_height(k) = r.rms_height(k) + 0.0001;
            end
          end
        end
        
        
        if isnan(r.rms_height(k))
          r.dielectric_constant(k) = nan;
        else
          r.dielectric_constant(k) = 1;
          clear mse
          for i = 1:5000
            mse(i) = abs(r.pc(k) - ((1-sqrt(r.dielectric_constant(k)))/((1+sqrt(r.dielectric_constant(k)))))^2*exp(-(2*(2*pi/(c/param.radar.fs))*r.rms_height(k))^2));
            if r.dielectric_constant(k) > 4
              r.dielectric_constant(k) = nan;
              warning('check this')
              %                         keyboard
              break
            else
              
              if i>1
                if mse(i-1) < mse(i)
                  break
                else
                  r.dielectric_constant(k) = r.dielectric_constant(k) + 0.01;
                  continue  ;
                end
              else
                r.dielectric_constant(k) = r.dielectric_constant(k) + 0.01;
              end
            end
          end
        end
        
        
        if any(~isnan(Greenland.ice_bed_power((l-500):(l+499))))
          ice_bed_power = Greenland.ice_bed_power((l-500):(l+499));
          depth = Greenland.depth((l-500):(l+499));
          clear id
          id = find(isnan(Greenland.ice_bed_power((l-500):(l+499))));
          ice_bed_power(id) = [];
          depth(id) = [];
          Greenland.ice_bed_power_avg(k) =   nanmean(ice_bed_power) ;
          Greenland.depth_avg(k) = nanmean(depth);
          Greenland.Latitude_avg(k) = nanmean(Greenland.Latitude((l-500):(l+499)));
          Greenland.Longitude_avg(k) = nanmean(Greenland.Longitude((l-500):(l+499)));
          Greenland.along_track_avg(k) = nanmedian(geodetic_to_along_track(Greenland.Latitude((l-500):(l+499)),Greenland.Longitude((l-500):(l+499))));
          Greenland.geometric_loss_avg(k) = nanmean(geometric_loss((l-500):(l+499)));
          
        else
          Greenland.ice_bed_power_avg(k) = nan;
          Greenland.depth_avg(k) = nan;
          Greenland.Latitude_avg(k) = nan;
          Greenland.Longitude_avg(k) = nan;
          Greenland.along_track_avg(k) = nan;
          Greenland.geometric_loss_avg(k) = nan;
          
        end
        
        if isnan(r.rms_height(k))
          k = k+1;
          continue;
        else
          Greenland.ice_bed_power_avg(k) = Greenland.ice_bed_power_avg(k)./(exp(-((4*pi*r.rms_height(k)/(c/param.radar.fs))*(sqrt(r.dielectric_constant(k))-1))^2)*(besseli(0,(((4*pi*r.rms_height(k)/(c/param.radar.fs))*(sqrt(r.dielectric_constant(k))-1))^2)/2))^2);
          k= k+1;
        end
        
      end
      if k+1 == length(Greenland.depth_avg)
        keyboard
      end
    end
    
    if cross_lines
      save(['/cresis/snfs1/scratch/manjish/test_r/cross_lines' num2str(M,'%03d') '.mat'],'r');
    else
      save(['/cresis/snfs1/scratch/manjish/test_r/' num2str(M,'%03d') '.mat'],'r');
    end
    clearvars r K
  end
  
   if plots
    figure;
    plot(10*log10(abs( Greenland.ice_bed_power_avg).^2));
    grid on
    title('Ice Bed Power surface roughness corrected')
  end
  
  
  
  %% compensating for bed roughness
  file_exist = false;
  if M<21
    if exist((['/cresis/snfs1/scratch/manjish/peterman/bedroughness/crossline',num2str(M),'.mat']),'file')
      load ((['/cresis/snfs1/scratch/manjish/peterman/bedroughness/crossline',num2str(M),'.mat']))
      file_exist = true;
      r=rbed;
    end
  else
    if exist((  (['/cresis/snfs1/scratch/manjish/peterman/bedroughness/verticalline',num2str(M1),'.mat'])),'file')
      load(['/cresis/snfs1/scratch/manjish/peterman/bedroughness/verticalline',num2str(M1),'.mat']);
      file_exist = true;
    end
    if ~exist('r','var')
        r=rbed;
    end
  end
  
  %         K  = floor(length(Greenland.ice_surface_power)/1000);
  if file_exist
    
    %         K  = floor(length(Greenland.ice_bed_power)/1000);
   k=1;
    num_int=r.settings.num_int;
    repeat_after=300;
    for l = num_int/2:repeat_after:length(Greenland.ice_bed_power)
      if ((l >= num_int/2) && ((l+num_int/2) < length(Greenland.ice_bed_power)))
        if k>size(r.lat)
          continue
        end
        if isnan(r.rms_height(k)) || isnan(Greenland.ice_bed_power_avg(k))
          k= k+1;
          continue ;
        else
          Greenland.ice_bed_power_avg(k) = Greenland.ice_bed_power_avg(k)./(exp(-(4*pi*r.rms_height(k)/(c/(param.radar.fs*sqrt(er_ice))))^2)*(besseli(0,((4*pi*r.rms_height(k)/(c/(param.radar.fs*sqrt(er_ice))))^2)/2))^2);
          k = k+1;
        end
      end
    end
    
    %           Greenland.ice_bed_power(1+(l-1)*1000:(l*1000)) = (Greenland.ice_bed_power(1+(l-1)*1000:(l*1000)))./(exp(-(4*pi*r.rms_height(l)/(c/param.radar.fs))^2)*(besseli(0,((4*pi*r.rms_height(l)/(c/param.radar.fs))^2)/2))^2);
    
    
    
    clearvars r k
  else
    
    %         K  = floor(length(Greenland.ice_bed_power)/1000);
    k = 1;
    for l = 501:250:length(Greenland.ice_bed_power)
      if ((l > 500) && ((l+500) < length(Greenland.ice_surface_power)))
        r.lat(k) = nanmean(Greenland.Latitude((l-500):(l+499)));
        r.lon(k) = nanmean(Greenland.Longitude((l-500):(l+499)));
        s = abs(Greenland.ice_bed_power((l-500):(l+499)));
        id = find(isnan(s)|isinf(s)|s==0);
        if length(id) > 500
          r.rms_height(k) = nan;
          r.dielectric_constant(k) = nan;
          r.pn(k) = nan;
          r.pc(k) = nan ;
          k = k+1;
          continue
        else
          s(id) = [];
        end
        try
          %pd = fitdist(double((s)).','Rician');
          [pd.s,pd.sigma]=ricefit(s);
        catch ME
          warning('unable to fit the distribution')
          r.rms_height(k) = nan;
          r.dielectric_constant(k) = nan;
          r.pn(k) = nan;
          r.pc(k) = nan ;
          k = k+1;
          continue
        end
        
        % phat = mle(double(abs(s)),'distribution','Rician');
        % x = 0:0.0001:0.2;
        % histogram(abs(s))
        % h = hist((s));
        
        % A = pdf(pd,x);
        a = pd.s;
        % pc = 2*10*log10(a);
        % pn = 10*log10(2*pd.sigma^2);
        S = pd.sigma;
        r.pc(k) = a^2;
        r.pn(k) = 2*2*pd.sigma^2;
        rms_fit = (r.pc(k)/r.pn(k))*4*(2*pi/(c/param.radar.fs))^2;
        
        r.rms_height(k) = 0.0001;
        clear MSE
        for i = 1:5000
          MSE(i) = abs(rms_fit - exp(-(2*(2*pi/(c/param.radar.fs))*r.rms_height(k))^2)/((r.rms_height(k))^2));
          if r.rms_height(k) > 0.40B(rline)*depth(rline)/100
            r.rms_height(k) = nan;
            warning('check this')
            %                         keyboard
            break
          else
            
            if i>1
              if MSE(i-1) < MSE(i)
                break
              else
                r.rms_height(k) = r.rms_height(k) + 0.0001;
                continue  ;
              end
            else
              r.rms_height(k) = r.rms_height(k) + 0.0001;
            end
          end
        end
        
        
        if isnan(r.rms_height(k))
          r.dielectric_constant(k) = nan;
        else
          r.dielectric_constant(k) = 1;
          clear mse
          for i = 1:5000
            mse(i) = abs(r.pc(k) - ((1-sqrt(r.dielectric_constant(k)))/((1+sqrt(r.dielectric_constant(k)))))^2*exp(-(2*(2*pi/(c/param.radar.fs))*r.rms_height(k))^2));
            if r.dielectric_constant(k) > 4
              r.dielectric_constant(k) = nan;
              warning('check this')
              %                         keyboard
              break
            else
              
              if i>1
                if mse(i-1) < mse(i)
                  break
                else
                  r.dielectric_constant(k) = r.dielectric_constant(k) + 0.01;
                  continue  ;
                end
              else
                r.dielectric_constant(k) = r.dielectric_constant(k) + 0.01;
              end
            end
          end
        end
        
        
        %        Greenland.ice_bed_power(1+(l-1)*1000:(l*1000)) = (Greenland.ice_bed_power(1+(l-1)*1000:(l*1000)))./(exp(-(4*pi*r.rms_height(l)/(c/param.radar.fs))^2)*(besseli(0,((4*pi*r.rms_height(l)/(c/param.radar.fs))^2)/2))^2);
        
        
        if any(~isnan(Greenland.ice_bed_power((l-500):(l+499))))
          ice_bed_power = Greenland.ice_bed_power((l-500):(l+499));
          depth = Greenland.depth((l-500):(l+499));
          clear id
          id = find(isnan(Greenland.ice_bed_power((l-500):(l+499))));
          ice_bed_power(id) = [];
          depth(id) = [];
          Greenland.ice_bed_power_avg(k) =   nanmean(ice_bed_power) ;
          Greenland.depth_avg(k) = nanmean(depth);
          Greenland.Latitude_avg(k) = nanmean(Greenland.Latitude((l-500):(l+499)));
          Greenland.Longitude_avg(k) = nanmean(Greenland.Longitude((l-500):(l+499)));
          Greenland.along_track_avg(k) = nanmedian(geodetic_to_along_track(Greenland.Latitude((l-500):(l+499)),Greenland.Longitude((l-500):(l+499))));
          Greenland.geometric_loss_avg(k) = nanmean(geometric_loss((l-500):(l+499)));
        else
          Greenland.ice_bed_power_avg(k) = nan;
          Greenland.depth_avg(k) = nan;
          Greenland.Latitude_avg(k) = nan;
          Greenland.Longitude_avg(k) = nan;
          Greenland.geometric_loss_avg(k) = nan;
          %                     Greenland.Latitude((l-500):(l+499)) =[];
          %                     Greenland.Longitude((l-500):(l+499))= [];
          Greenland.along_track_avg(k) = nan;
          
        end
        
        
        
        
        
        if isnan(r.rms_height(k))
          k = k+1;
          continue;
        else
          Greenland.ice_bed_power_avg(k) = Greenland.ice_bed_power_avg(k)./(exp(-(4*pi*r.rms_height(k)/(c/(param.radar.fs*sqrt(er_ice))))^2)*(besseli(0,((4*pi*r.rms_height(k)/(c/(param.radar.fs*sqrt(er_ice))))^2)/2))^2);
          k = k+1;
        end
      end
    end
    
    
    if cross_lines
      save(['/cresis/snfs1/scratch/manjish/cross_lines' num2str(M,'%03d') '.mat'],'r');
    else
      save(['/cresis/snfs1/scratch/manjish/' num2str(M,'%03d') '.mat'],'r');
    end
    clearvars r K
    
  end
 
  if plots
    figure(3);
    hold on;
    plot(10*log10(abs( Greenland.ice_bed_power_avg).^2));
    grid on
    legend('Sf corrected ',' Bed roughness corrected')
  end
  
  %%  relative geometrically corrected bed - echo power
  %     for i = 1:length(Greenland.ice_bed_power)
  %         if i < 3
  %             Greenland.ice_bed_power_rgc(i) = Greenland.ice_bed_power(i)-nanmean(Greenland.ice_bed_power(1:6-i));
  %         elseif i+2 > length(Greenland.ice_bed_power)
  %             Greenland.ice_bed_power_rgc(i) = Greenland.ice_bed_power(i)- nanmean(Greenland.ice_bed_power(i-5:end));
  %         else
  %             Greenland.ice_bed_power_rgc(i) = Greenland.ice_bed_power(i)- nanmean(Greenland.ice_bed_power(i-2:i+2));
  %         end
  %     end
  
  if length(Greenland.ice_bed_power_avg) ~= length(Greenland.geometric_loss_avg)
    keyboard
  end
  
  Greenland.ice_bed_power_cgl =lp(Greenland.ice_bed_power_avg)+lp(Greenland.geometric_loss_avg);
  if plots
    figure(4);subplot(2,1,1);plot(Greenland.ice_bed_power_cgl);title('Bed power after Geom corr')
    subplot(2,1,2);plot(lp(Greenland.geometric_loss_avg)); title('Geom correction')
  end
  
  id = ~(isfinite( Greenland.ice_bed_power_cgl ));
  Greenland.ice_bed_power_cgl(id) = [];
  Greenland.Latitude_avg(id) = [];
  Greenland.Longitude_avg(id) = [];
  Greenland.depth_avg(id) = [];
 
  id = ~(isfinite( Greenland.depth_avg ));
  Greenland.ice_bed_power_cgl(id) = [];
  Greenland.Latitude_avg(id) = [];
  Greenland.Longitude_avg(id) = [];
  Greenland.depth_avg(id) = [];
  
  ice_bed_power_G_r_corrected = cat(2,ice_bed_power_G_r_corrected,Greenland.ice_bed_power_cgl);
  lat_G_r_corrected =  cat(2,lat_G_r_corrected,Greenland.Latitude_avg);
  lon_G_r_corrected =  cat(2,lon_G_r_corrected,Greenland.Longitude_avg);
  depth_G_r_corrected =  cat(2,depth_G_r_corrected,Greenland.depth_avg);
  
  
  %% attenuation_fitting
  %reference_power =  25.1114 ;
  %reference_power = -62.84 ;
  reference_power = median((Greenland.ice_bed_power_cgl));
  relative_ice_bed_power_G_r_corrected = (Greenland.ice_bed_power_cgl)-reference_power;
  
  if plots
    figure;plot(relative_ice_bed_power_G_r_corrected);
    title('Relative Power')
  end
  
  window1=200;
  window2=10;
  power_filtered_long=sgolayfilt(relative_ice_bed_power_G_r_corrected,2,window1+1,gausswin(window1+1));
  power_filtered_short=sgolayfilt(relative_ice_bed_power_G_r_corrected,2,window2+1,gausswin(window2+1));
  
  if plots
    figure;plot(power_filtered_long); hold on; plot(relative_ice_bed_power_G_r_corrected)
    legend('filtered','original'); title('Long filter')
    figure;plot(power_filtered_short); hold on; plot(relative_ice_bed_power_G_r_corrected)
    legend('filtered','original'); title('short filter')
  end
  
   Greenland.depth_avg = Greenland.depth_avg/1000;
   relative_depth = nanmean(isfinite(Greenland.depth_avg));

  %% Attenuation calculation
  along_track = geodetic_to_along_track(Greenland.Latitude_avg,Greenland.Longitude_avg);
  along_track = along_track/1000;
  modified_attenuation = zeros(size(along_track));
  estimated_na = zeros(size(along_track));
  Na=5;
   
      dn = (-10:0.01:10);   %Filter before finding DN????
  for j = 1:length(dn)
    
    term_1 = 2*dn(j).*((Greenland.depth_avg-relative_depth)).*((along_track-mean(along_track)));
    term_2 = 2*Na*((Greenland.depth_avg-relative_depth));
    S(j) = sum(abs((power_filtered_short)+term_1+term_2).^2);
    if 0
      plot(real(-power_filtered_short))
      hold on
      plot(term_1+term_2)
      % keyboard
      %clf
    end
    if j > 1
      if S(j-1)<S(j)
        DN = dn(j-1);
        break
      end
    end
  end
  
  if j ==length(dn)
    [v,id2] = min(S);
    DN = dn(id2);
  end
  
%     if variable_attenuation_mse < constant_attenuation_mse
      %             modified_atteanuation(id) = 2.*variable_attenuation_scale.*Greenland.attenuation_constant(id).*(Greenland.depth_avg(id)-relative_depth);
      modified_attenuation = 2.*Na.*(Greenland.depth_avg-relative_depth);
      
      Na_mod=Na+DN*((along_track-mean(along_track)));
      variable_attenuation=modified_attenuation+ 2*DN.*((Greenland.depth_avg-relative_depth)).*((along_track-mean(along_track)));
      if plots
        figure;plot(modified_attenuation)
        figure;plot(variable_attenuation)
        figure;plot(Greenland.depth)
      end
      
  %%
  
  
  % clf
  % plot(-relative_ice_bed_power_G_r_corrected)
  % hold on
  % plot(2.*variable_attenuation_scale*Greenland.attenuation_constant.*(Greenland.depth_avg-relative_depth))
  % hold on
  % plot(2.*constant_attenuation.*(Greenland.depth_avg-relative_depth))
  % grid
  % keyboard
  
  constant_attenuation =  cat(2,constant_attenuation, modified_attenuation);
  estimated_Na = cat(2,estimated_Na, Na);
  estimated_Na_mod=cat(2,estimated_Na_mod,Na_mod);
  var_attn=cat(2,var_attn,variable_attenuation);
  
  
  close all
  %%
  % along_track = geodetic_to_along_track(Greenland.Latitude_avg,Greenland.Longitude_avg);
  %
  % along_track = along_track/1000;
  % plot( modified_atteanuation+ 0.02*(along_track-median(along_track)))
end
%%
depth_G_r_corrected = depth_G_r_corrected / 1000;

reference_power = median(ice_bed_power_G_r_corrected(isfinite(ice_bed_power_G_r_corrected)));


%relative_ice_bed_power_G_r_corrected = lp(ice_bed_power_G_r_corrected/reference_power);
relative_ice_bed_power_G_r_corrected = (ice_bed_power_G_r_corrected)-(reference_power);

relative_reflectivty_using_constant_attn = relative_ice_bed_power_G_r_corrected + constant_attenuation ;
relative_reflectivty_using_variable_attn = relative_ice_bed_power_G_r_corrected + var_attn ;

geotiff_fn = '/cresis/snfs1/dataproducts/GIS_data/greenland/Landsat-7/Greenland_natural.tif';
proj = geotiffinfo(geotiff_fn);
%proj = geotiffinfo('X:\GIS_data\antarctica\Landsat-7\Antarctica_LIMA_480m.tif');

[A CMAP R]= geotiffread(geotiff_fn);

%% Reflectivity
figure(4)
mapshow(rgb2gray(A),CMAP/1e3);
xlabel('X (km)');
ylabel('Y (km)');
%xlim([350 650]);
%ylim([-1000 -600]);

hold on
clear gps.x gps.y
[gps.x,gps.y] = projfwd(proj,lat_G_r_corrected,lon_G_r_corrected);
%         lt = cat(2,lt, Greenland.Latitude_avg ) ;
%         ln =cat(2, ln,  Greenland.Longitude_avg) ;
gps.x = gps.x / 1000;
gps.y = gps.y / 1000;
hold on;
%scatter(gps.x,gps.y,20,relative_reflectivty_using_modified_na,'fill')
%scatter(gps.x,gps.y,20,2.*estimated_Na.*depth_G_r_corrected,'fill')
%scatter(gps.x,gps.y,20,depth_G_r_corrected,'fill')
scatter(gps.x,gps.y,20,relative_reflectivty_using_variable_attn,'fill')
 caxis([-15 15])
colorbar;
title('Reflectivity')

%ylim([-1000 -655]);
figure, hist(relative_reflectivty_using_variable_attn,30)


%% Attenuation plot
figure(6)
mapshow(rgb2gray(A),CMAP/1e3);
xlabel('X (km)');
ylabel('Y (km)');
%xlim([350 650]);
%ylim([-1000 -600]);

hold on
clear gps.x gps.y
[gps.x,gps.y] = projfwd(proj,lat_G_r_corrected,lon_G_r_corrected);
%         lt = cat(2,lt, Greenland.Latitude_avg ) ;
%         ln =cat(2, ln,  Greenland.Longitude_avg) ;
gps.x = gps.x / 1000;
gps.y = gps.y / 1000;
hold on;
%scatter(gps.x,gps.y,20,relative_reflectivty_using_modified_na,'fill')
%scatter(gps.x,gps.y,20,2.*estimated_Na.*depth_G_r_corrected,'fill')
%scatter(gps.x,gps.y,20,depth_G_r_corrected,'fill')
scatter(gps.x,gps.y,20,estimated_Na_mod,'fill')
%caxis([-15 15])
colorbar;
title('Attenuation')


%% Depth
geotiff_fn = '/cresis/snfs1/dataproducts/GIS_data/greenland/Landsat-7/Greenland_natural.tif';
proj = geotiffinfo(geotiff_fn);
%proj = geotiffinfo('X:\GIS_data\antarctica\Landsat-7\Antarctica_LIMA_480m.tif');

[A CMAP R]= geotiffread(geotiff_fn);
figure(7)
mapshow(rgb2gray(A),CMAP/1e3);
xlabel('X (km)');
ylabel('Y (km)');
%xlim([350 650]);
%ylim([-1000 -600]);

hold on
clear gps.x gps.y
[gps.x,gps.y] = projfwd(proj,lat_G_r_corrected, lon_G_r_corrected);
%         lt = cat(2,lt, Greenland.Latitude_avg ) ;
%         ln =cat(2, ln,  Greenland.Longitude_avg) ;
gps.x = gps.x / 1000;
gps.y = gps.y / 1000;
hold on;
scatter(gps.x,gps.y,20,depth_G_r_corrected,'fill')
%caxis([-15 15])
colorbar;
title('Depth')
