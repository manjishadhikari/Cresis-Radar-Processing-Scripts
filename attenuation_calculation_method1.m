 
function [Attenuation]=attenuation_calculation_method1(Greenland,power_filtered_long,power_filtered_short)

%% Attenuation calculation
    
    %Method 1 fit Na and DN for evry 1 km
    if 1
      along_track = geodetic_to_along_track(Greenland.Latitude_avg,Greenland.Longitude_avg);
      along_track = along_track/1000;
      modified_attenuation = zeros(size(along_track));
      estimated_na = zeros(size(along_track));
      for i = 1:ceil(max(along_track))
        clear id
        if i ==1
          id = find((along_track > 0) & (along_track <= i))  ;
          
        else
          id = find((along_track > i-1) & (along_track <= i))  ; %Every 1 km
        end
        if length(id) < 1
          continue;
        else
          Na = 0:0.05:20;    %filter before finding Na???
          
          for  j = 1:length(Na)
            %         plot(-relative_ice_bed_power_G_r_corrected, '*')% apparent attenuation
            %         hold on
            %         plot(2*Na(j)*(Greenland.depth_avg-relative_depth), '*')
            mse(j) = mean((-power_filtered_long(id)- 2*Na(j)*(Greenland.depth_avg(id)-Greenland.relative_depth)).^2);
            %  keyboard
            %             grid
            %             hold off
          end
          [constant_attenuation_mse, index] = min(mse);
          Na_bar  = Na(index);
        end
        %   Na_bar=11;
        
        dn = (-5:0.01:5);   %Filter before finding DN????
        for j = 1:length(dn)
          
          term_1 = 2*dn(j).*((Greenland.depth_avg(id)-Greenland.relative_depth)).*((along_track(id)-mean(along_track(id))));
          term_2 = 2*Na_bar*((Greenland.depth_avg(id)-Greenland.relative_depth));
          S(j) = sum(abs((power_filtered_short(id))+term_1+term_2).^2);
          if 0
            plot(real(-power_filtered_short(id)))
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
        const_attenuation(id) = 2.*Na_bar.*(Greenland.depth_avg(id)-Greenland.relative_depth);
        estimated_na(id) = Na_bar;
        estimated_dn(id)=DN;
        
        var_attenuation(id)=const_attenuation(id)+ 2*DN.*((Greenland.depth_avg(id)-Greenland.relative_depth)).*((along_track(id)-mean(along_track(id))));
        %     else
        %       modified_atteanuation(id) = 2.*constant_attenuation.*(Greenland.depth_avg(id)-relative_depth);
        %       estimated_na(id) = constant_attenuation;
        %
        %     end
        %     plot(-relative_ice_bed_power_G_r_corrected(id))
        %     hold on
        %     plot(2.*variable_attenuation_scale*Greenland.attenuation_constant(id).*(Greenland.depth_avg(id)-relative_depth))
        %     hold on
        %     plot(2.*constant_attenuation.*(Greenland.depth_avg(id)-relative_depth))
        %     grid
        
      end
    end
    Attenuation.const_attenuation=const_attenuation;
    Attenuation.var_attenuation=var_attenuation;
    Attenuation.estimated_na=estimated_na;
    Attenuation.estimated_dn=estimated_dn;
end