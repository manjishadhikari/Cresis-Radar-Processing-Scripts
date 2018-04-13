 
function [Greenland]=coh_integration(Greenland,numofCohInt)

if numofCohInt~=0
      
      
      Nx=floor(length(Greenland.ice_bed_power)/numofCohInt);
      %        ice_surface_power_tmp=nan*ones(1,Nx);
      %        Latitude=nan*ones(1,Nx);
      %        Longitude=nan*ones(1,Nx);
      %        Elevation=nan*ones(1,Nx);
      %        surface_twtt_tmp=nan*ohttp://www.usagoals.me/c/ncaa/2479253/3/kansas-vs-penn-live-stream/nes(1,Nx);
      
      for i= 1:Nx
        idx1=(i-1)*numofCohInt+1;
        idx2=i*numofCohInt;
        if idx2>length(Greenland.GPS_time) & length(Greenland.GPS_time)-idx1>numofCohInt/2
          idx2=length(Greenland.GPS_time);
        elseif idx2>length(Greenland.GPS_time)& length(Greenland.GPS_time)-idx1<numofCohInt/2
          continue;
        end
        Greenland.GPS_time(i)=mean(Greenland.GPS_time(idx1:idx2));
        Greenland.Latitude(i)=mean(Greenland.Latitude(idx1:idx2));
        Greenland.Longitude(i)=mean(Greenland.Longitude(idx1:idx2));
        Greenland.Elevation(i)=mean(Greenland.Elevation(idx1:idx2));
        Greenland.surface_time(i)=mean(Greenland.surface_time(idx1:idx2));
        Greenland.ice_bed_time(i)=mean(Greenland.ice_bed_time(idx1:idx2));
        Greenland.ice_bed_power(i)=mean(Greenland.ice_bed_power(idx1:idx2));
        Greenland.ice_surface_power(i)=mean(Greenland.ice_surface_power(idx1:idx2));
        
        
      end
      % Greenland.Time=Greenland.Time;
      Greenland.segments_length=round((Greenland.segments_length)/numofCohInt);
end
end