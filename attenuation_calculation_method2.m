
function [Attenuation]=attenuation_calculation_method2(Greenland,power_filtered_short,Na_bar)

      along_track = geodetic_to_along_track(Greenland.Latitude_avg_filt,Greenland.Longitude_avg_filt);
      along_track = along_track/1000;
    
      Na= Na_bar;
       
      dn = (-25:0.01:25);   %Filter before finding DN????
      for j = 1:length(dn)
        
        term_1 = 2*dn(j).*((Greenland.depth_avg_filt-Greenland.relative_depth)).*((along_track-mean(along_track)));
        term_2 = 2*Na*((Greenland.depth_avg_filt-Greenland.relative_depth));
        S(j) = mean(((power_filtered_short)+term_1+term_2).^2);
        
%         if j > 1
%           if S(j-1)<S(j)
%             DN = dn(j-1);
%             break
%           end
%         end
      end
      
      if 0
          figure;plot((power_filtered_short))
          hold on
          plot(term_1+term_2)
        figure;plot(S);
      end
      
      
      if j ==length(dn)
        [v,id2] = min(S);
        DN = dn(id2)
      end
      
      %     if variable_attenuation_mse < constant_attenuation_mse
      %             modified_atteanuation(id) = 2.*varia/N/dcwan/projects/cresis/output/accum/2017_Antarctica_P3/CSARP_qlook/20171103_04/Data_20171103_04_010.matble_attenuation_scale.*Greenland.attenuation_constant(id).*(Greenland.depth_avg(id)-Greenland.relative_depth);
     

      
      Attenuation.const_attenuation = 2.*Na.*(Greenland.depth_avg_filt-Greenland.relative_depth);
      Attenuation.estimated_na=Na;
      Attenuation.estimated_dn=DN;
      Attenuation.mod_na=Na+DN*((along_track-mean(along_track)));
      Attenuation.var_attenuation=Attenuation.const_attenuation+ 2*DN.*((Greenland.depth_avg_filt-Greenland.relative_depth)).*((along_track-mean(along_track)));
      
      end