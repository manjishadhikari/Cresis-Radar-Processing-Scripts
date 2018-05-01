
function [Attenuation]=attenuation_calculation_method2(Greenland,power_filtered_short)

      along_track = geodetic_to_along_track(Greenland.Latitude_avg,Greenland.Longitude_avg);
      along_track = along_track/1000;
      modified_attenuation = zeros(size(along_track));
      estimated_na = zeros(size(along_track));
      Na=11.162;
       
      dn = (-25:0.01:25);   %Filter before finding DN????
      for j = 1:length(dn)
        
        term_1 = 2*dn(j).*((Greenland.depth_avg-Greenland.relative_depth)).*((along_track-mean(along_track)));
        term_2 = 2*Na*((Greenland.depth_avg-Greenland.relative_depth));
        S(j) = sum(abs((power_filtered_short)+term_1+term_2).^2);
        if 0
          figure;plot(real(-power_filtered_short))
          hold on
          plot(term_1+term_2)
          % keyboard
          %clf
        end
%         if j > 1
%           if S(j-1)<S(j)
%             DN = dn(j-1);
%             break
%           end
%         end
      end
      
      if j ==length(dn)
        [v,id2] = min(S);
        DN = dn(id2);
      end
      
      %     if variable_attenuation_mse < constant_attenuation_mse
      %             modified_atteanuation(id) = 2.*varia/N/dcwan/projects/cresis/output/accum/2017_Antarctica_P3/CSARP_qlook/20171103_04/Data_20171103_04_010.matble_attenuation_scale.*Greenland.attenuation_constant(id).*(Greenland.depth_avg(id)-Greenland.relative_depth);
     
      
   
      
      Attenuation.const_attenuation = 2.*Na.*(Greenland.depth_avg-Greenland.relative_depth);
      Attenuation.estimated_na=Na;
      Attenuation.estimated_dn=DN;
      Attenuation.mod_na=Na+DN*((along_track-mean(along_track)));
      Attenuation.var_attenuation=Attenuation.const_attenuation+ 2*DN.*((Greenland.depth_avg-Greenland.relative_depth)).*((along_track-mean(along_track)));
      
      end