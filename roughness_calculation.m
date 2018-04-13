
function [r]=roughness_calculation(Greenland,settings)
M=21;
M1=1;
physical_constants;
param.radar.fc=195000000;
num_int=settings.num_int;  % ~210 metres
repeat_after=settings.repeat_after;
cross_lines=0;

if strcmp(settings.type,'surface')
  power=Greenland.ice_surface_power;
else
  power=Greenland.ice_bed_power;
end
k=1;
for l = num_int/2:repeat_after:length(power)
  if ((l >= num_int/2) && ((l+num_int/2) < length(power)))
    if all(abs(Greenland.roll((l-num_int/2+1):(l+num_int/2)))<5)
      r.lat(k) = nanmean(Greenland.Latitude((l-num_int/2+1):(l+num_int/2)));
      r.lon(k) = nanmean(Greenland.Longitude((l-num_int/2+1):(l+num_int/2)));
      s = abs(power((l-num_int/2+1):(l+num_int/2)));
      id = find(isnan(s)|isinf(s)|s==0);
      if length(id) > num_int/2
        r.rms_height(k) = nan;
        r.dielectric_constant(k) = nan;
        r.pn(k) = nan;
        r.pc(k) = nan ;
        k= k+1;
        continue
      else
        s(id) = [];
      end
      
      try
        
        [pd.s, pd.sigma]=ricefit(s);
        
      catch ME
        warning('unable to fit the distribution')
        r.rms_height(k) = nan;
        r.dielectric_constant(k) = nan;
        r.pn(k) = nan;
        r.pc(k) = nan ;
        k = k+1;
        continue
      end
      
      a = pd.s;
      % pc = 2*10*log10(a);
      % pn = 10*log10(2*pd.sigma^2);
      S = pd.sigma;
      r.pc(k) = a^2;
      r.pn(k) = 2*pd.sigma^2;
      rms_fit = (r.pc(k)/r.pn(k))*4*(2*pi/(c/param.radar.fc))^2;
      
      r.rms_height(k) = 0.001;
      
      options = optimset('PlotFcns',@optimplotfval);
      cost_func=@(r_rms)(abs(rms_fit - exp(-(2*(2*pi/(c/param.radar.fc))*r_rms).^2)/((r_rms).^2)));
      r_rms=[0.00001,1];
      [r_rmsheight,fval]=fminsearch(cost_func,r_rms);
      r.rms_height(k)=mean(abs(r_rmsheight));
      %{
            clear MSE
            for i = 1:5000
              MSE(i) = abs(rms_fit - exp(-(2*(2*pi/(c/param.radar.fc))*r.rms_height(k))^2)/((r.rms_height(k))^2));
              if r.rms_height(k) > 0.30
                r.rms_height(k) = nan;
                warning('check this')
                %                         keyboard
                break
              else
                
                if i>1
                  if MSE(i-1) < MSE(i)options = optimset('PlotFcns',@optimplotfval);
                    %          cost_func=@(r_rms)(abs(rms_fit - exp(-(2*(2*pi/(c/param.radar.fc))*r_rms).^2)/((r_rms).^2)));
                    %      r_rms=[0.00001,8];
                    %      [r_rmsheight,fval]=fminsearch(cost_func,r_rms,options);
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
            
      %}
      if isnan(r.rms_height(k))
        r.dielectric_constant(k) = nan;
      else
        r.dielectric_constant(k) = 1;
        clear mse
        for i = 1:5000
          mse(i) = abs(r.pc(k) - ((1-sqrt(r.dielectric_constant(k)))/((1+sqrt(r.dielectric_constant(k)))))^2*exp(-(2*(2*pi/(c/param.radar.fc))*r.rms_height(k))^2));
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
      
      
%       if any(~isnan(Greenland.ice_bed_power((l-num_int/2+1):(l+num_int/2))))
%         ice_bed_power = Greenland.ice_bed_power((l-num_int/2+1):(l+num_int/2));
%         depth = Greenland.depth((l-num_int/2+1):(l+num_int/2));
%         clear id
%         id = find(isnan(Greenland.ice_bed_power((l-num_int/2+1):(l+num_int/2))));
%         ice_bed_power(id) = [];
%         depth(id) = [];
%         Greenland.ice_bed_power_avg(k) =   nanmean(ice_bed_power) ;
%         Greenland.depth_avg(k) = nanmean(depth);
%         Greenland.Latitude_avg(k) = nanmean(Greenland.Latitude((l-num_int/2+1):(l+num_int/2)));
%         Greenland.Longitude_avg(k) = nanmean(Greenland.Longitude((l-num_int/2+1):(l+num_int/2)));
%         Greenland.along_track_avg(k) = nanmedian(geodetic_to_along_track(Greenland.Latitude((l-num_int/2+1):(l+num_int/2)),Greenland.Longitude((l-num_int/2+1):(l+num_int/2))));
%         Greenland.geometric_loss_avg(k) = nanmean(geometric_loss((l-num_int/2+1):(l+num_int/2)));
%         Greenland.maxroll(k)=max(Greenland.roll((l-num_int/2+1):(l+num_int/2)));
%         orig_avg_power(k)=Greenland.ice_bed_power_avg(k); %Original power for comparison
%         
%         
%       else
%         Greenland.ice_bed_power_avg(k) = nan;
%         Greenland.depth_avg(k) = nanmean(depth);
%         Greenland.Latitude_avg(k) = nanmean(Greenland.Latitude((l-num_int/2+1):(l+num_int/2)));
%         Greenland.Longitude_avg(k) = nanmean(Greenland.Longitude((l-num_int/2+1):(l+num_int/2)));
%         Greenland.along_track_avg(k) = nanmedian(geodetic_to_along_track(Greenland.Latitude((l-num_int/2+1):(l+num_int/2)),Greenland.Longitude((l-num_int/2+1):(l+num_int/2))));
%         Greenland.geometric_loss_avg(k) = nanmean(geometric_loss((l-num_int/2+1):(l+num_int/2)));
%         Greenland.maxroll(k)=max(Greenland.roll((l-num_int/2+1):(l+num_int/2)));
%         orig_avg_power(k)=Greenland.ice_bed_power_avg(k); %Original power for comparison
%       end
      
      if isnan(r.rms_height(k))
        k = k+1;
        continue;
      else
%         sf_corr_power(k)=(exp(-((4*pi*r.rms_height(k)/(c/param.radar.fc))*(sqrt(er_ice)-1))^2)*(besseli(0,(((4*pi*r.rms_height(k)/(c/param.radar.fc))*(sqrt(er_ice)-1))^2)/2))^2);
%         
%         Greenland.ice_bed_power_avg(k) = Greenland.ice_bed_power_avg(k).*(exp(-((4*pi*r.rms_height(k)/(c/param.radar.fc))*(sqrt(er_ice)-1))^2)*(besseli(0,(((4*pi*r.rms_height(k)/(c/param.radar.fc))*(sqrt(er_ice)-1))^2)/2))^2);
%         
      end
    else
%       Greenland.ice_bed_power_avg(k) = nan;
%       Greenland.depth_avg(k) = nanmean(Greenland.depth);
%       Greenland.Latitude_avg(k) = nanmean(Greenland.Latitude((l-num_int/2+1):(l+num_int/2)));
%       Greenland.Longitude_avg(k) = nanmean(Greenland.Longitude((l-num_int/2+1):(l+num_int/2)));
%       Greenland.along_track_avg(k) = nanmedian(geodetic_to_along_track(Greenland.Latitude((l-num_int/2+1):(l+num_int/2)),Greenland.Longitude((l-num_int/2+1):(l+num_int/2))));
%       Greenland.geometric_loss_avg(k) = nanmean(geometric_loss((l-num_int/2+1):(l+num_int/2)));
%       Greenland.maxroll(k)=max(Greenland.roll((l-num_int/2+1):(l+num_int/2)));
%       orig_avg_power(k)=Greenland.ice_bed_power_avg(k); %Original power for comparison
%       r.lat(k)= Greenland.Latitude_avg(k);
%       r.lon(k)= Greenland.Longitude_avg(k);
%       r.dielectric_constant(k)=nan;
      r.rms_height(k)=nan;
%       r.pc(k)=nan;
%       r.pn(k)=nan;
%       r.settings.num_int=num_int;
%       r.settings.repeat_after=repeat_after;
    end
    k= k+1;
  end

end

r.num_int=num_int;
r.repeat_after=repeat_after;

disp('Saving surface roughness values')

if strcmp(settings.type,'surface')
  if cross_lines
    save(['/cresis/snfs1/scratch/manjish/new_peterman/surface_roughness/crossline' num2str(M1,'%d') '.mat'],'r');
  else
    save(['/cresis/snfs1/scratch/manjish/new_peterman/surface_roughness/verticalline' num2str(M1,'%d') '.mat'],'r');
  end
elseif strcmp(settings.type,'bed')
  if cross_lines
    save(['/cresis/snfs1/scratch/manjish/new_peterman/bed_roughness/crossline' num2str(M1,'%d') '.mat'],'r');
  else
    save(['/cresis/snfs1/scratch/manjish/new_peterman/bed_roughness/verticalline' num2str(M1,'%d') '.mat'],'r');
  end
end
end
