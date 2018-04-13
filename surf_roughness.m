 

function [Greenland,sf_rms,sf_corr_power,orig_avg_power]=surf_roughness(Greenland,settings)
M=21;
M1=1;
physical_constants;
param.radar.fc=195000000;
file_exist = false;
    if M<21
      if exist((  (['/cresis/snfs1/scratch/manjish/new_peterman/surface_roughness/crossline',num2str(M),'.mat'])),'file')
        load(['/cresis/snfs1/scratch/manjish/new_peterman/surface_roughness/crossline',num2str(M),'.mat']);
        file_exist = true;
      end
    else
      
      if exist((  (['/cresis/snfs1/scratch/manjish/new_peterman/surface_roughness/verticalline',num2str(M1),'.mat'])),'file')
        load(['/cresis/snfs1/scratch/manjish/new_peterman/surface_roughness/verticalline',num2str(M1),'.mat']);
        file_exist = true;
      end
    end
    
   
    if ~file_exist
      [r]=roughness_calculation(Greenland,settings);
    end
      k = 1;
      
      num_int=r.num_int;  % ~210 metres
      repeat_after=r.repeat_after;
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
            Greenland.geometric_loss_avg(k) = nanmean(Greenland.geometric_loss((l-num_int/2+1):(l+num_int/2)));
             Greenland.maxroll(k)=max(Greenland.roll((l-num_int/2+1):(l+num_int/2)));
            orig_avg_power(k)=Greenland.ice_bed_power_avg(k); %Original power for comparison
            
          else
            Greenland.ice_bed_power_avg(k) = nan;
            Greenland.depth_avg(k) = nan;
            Greenland.Latitude_avg(k) = nan;
            Greenland.Longitude_avg(k) = nan;
            Greenland.along_track_avg(k) = nan;
            Greenland.geometric_loss_avg(k) = nan;
            Greenland.maxroll(k)=nan;

            orig_avg_power(k)=nan;
          end
          if k>size(r.lat)
            keyboard
            continue
          end
          if isnan(r.rms_height(k))
            k= k+1;
            continue;
          else
            %sf_corr_power(k)=(exp(-((4*pi*r.rms_height(k)/(c/param.radar.fc)))^2)*(besseli(0,(((4*pi*r.rms_height(k)/(c/param.radar.fc)))^2)/2))^2);
           % Greenland.ice_bed_power_avg(k) = Greenland.ice_bed_power_avg(k)*(exp(-((4*pi*r.rms_height(k)/(c/param.radar.fc)))^2)*(besseli(0,(((4*pi*r.rms_height(k)/(c/param.radar.fc)))^2)/2))^2);
           sf_corr_power(k)=(exp(-((4*pi*r.rms_height(k)/(c/param.radar.fc))*(sqrt(er_ice)-1))^2)*(besseli(0,(((4*pi*r.rms_height(k)/(c/param.radar.fc))*(sqrt(er_ice)-1))^2)/2))^2); 
           Greenland.ice_bed_power_avg(k) = Greenland.ice_bed_power_avg(k).*(exp(-((4*pi*r.rms_height(k)/(c/param.radar.fc))*(sqrt(er_ice)-1))^2)*(besseli(0,(((4*pi*r.rms_height(k)/(c/param.radar.fc))*(sqrt(er_ice)-1))^2)/2))^2);
            k= k+1;
          end
        end
      end
      sf_rms=r;
      clearvars r  k
    
end