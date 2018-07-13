 
function [Attenuation]=attenuation_calculation_method3(Greenland,power_filtered_short,Na_bar)

%% Attenuation calculation
    
    %Method  DN for evry 10 km and take Na from overall
    if 1
      along_track = geodetic_to_along_track(Greenland.Latitude_avg_filt,Greenland.Longitude_avg_filt);
      along_track = along_track/1000;
  
      dst=floor(along_track(end)/10);
      
      estimated_na = zeros(size(along_track));
      for i = 1:dst
        clear id
        if i ==1
          id = find((along_track > 0) & (along_track <= 10))  ;
        elseif i==dst
 
            id=find((along_track > 10*(i-1)) & (along_track <= along_track(end)));
        else
          id = find((along_track > 10*(i-1)) & (along_track <= 10*i))  ; %Every 1 km
        end
        if length(id) < 1
          continue;
        else
        
          
          %test
%           [r,Na_bar2,b]=regression(-2*(Greenland.depth_avg(id)-Greenland.relative_depth),power_filtered_long(id));
%           if 0
%             figure;plot(-2*(Greenland.depth_avg(id)-Greenland.relative_depth),power_filtered_long(id)); 
%           end
%           
          
        end
        %   Na_bar=11;
      
      
        dn = (-25:0.01:25);   %Filter before finding DN????
        for j = 1:length(dn)
          
          term_1 = 2*dn(j).*((Greenland.depth_avg_filt(id)-Greenland.relative_depth)).*((along_track(id)-mean(along_track(id))));
          term_2 = 2*Na_bar*((Greenland.depth_avg_filt(id)-Greenland.relative_depth));
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
        const_attenuation(id) = 2.*Na_bar.*(Greenland.depth_avg_filt(id)-Greenland.relative_depth);
        estimated_na(id) = Na_bar;
        mod_na(id)=Na_bar+2*DN*((along_track(id)-mean(along_track(id))));
        estimated_dn(id)=DN;
        
        var_attenuation(id)=const_attenuation(id)+ 2*DN.*((Greenland.depth_avg_filt(id)-Greenland.relative_depth)).*((along_track(id)-mean(along_track(id))));
       
      end
    end
    Attenuation.const_attenuation=const_attenuation;
    Attenuation.var_attenuation=var_attenuation;
    Attenuation.estimated_na=estimated_na;
    Attenuation.estimated_dn=estimated_dn;
    Attenuation.mod_na=mod_na;
end