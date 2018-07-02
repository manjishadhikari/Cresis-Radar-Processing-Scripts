
function [Greenlandout]=coh_integration(Greenland,numofCohInt)

if numofCohInt~=0
  
  Nx=floor(length(Greenland.ice_bed_power)/numofCohInt);
  %        ice_surface_power_tmp=nan*ones(1,Nx);
  %        Latitude=nan*ones(1,Nx);
  %        Longitude=nan*ones(1,Nx);
  %        Elevation=nan*ones(1,Nx);
  %        surface_twtt_tmp=nan*ones(1,Nx);
  
  for i= 1:Nx
    idx1=(i-1)*numofCohInt+1;
    idx2=i*numofCohInt;
    if idx2>length(Greenland.GPS_time) & length(Greenland.GPS_time)-idx1>numofCohInt/2
      idx2=length(Greenland.GPS_time);
    elseif idx2>length(Greenland.GPS_time)& length(Greenland.GPS_time)-idx1<numofCohInt/2
      continue;
    end
    
    
    Greenlandout.GPS_time(i)=mean(Greenland.GPS_time(idx1:idx2));
    if isfield(Greenland,'Roll')
      Greenlandout.Roll(i)=max(Greenland.Roll(idx1:idx2));
    end
    if isfield(Greenland,'ice_bed_elevation')
      Greenlandout.ice_bed_elevation(i)=nanmean(Greenland.ice_bed_elevation(idx1:idx2));
    end
    
     if isfield(Greenland,'settings')
      Greenlandout.settings=Greenland.settings;
    end
    
    
    Greenlandout.Latitude(i)=nanmean(Greenland.Latitude(idx1:idx2));
    Greenlandout.Longitude(i)=nanmean(Greenland.Longitude(idx1:idx2));
    Greenlandout.Elevation(i)=nanmean(Greenland.Elevation(idx1:idx2));
    Greenlandout.surface_time(i)=nanmean(Greenland.surface_time(idx1:idx2));
    Greenlandout.ice_bed_time(i)=nanmean(Greenland.ice_bed_time(idx1:idx2));
    if (length(find(isnan(Greenland.ice_bed_power(idx1:idx2))))>numofCohInt/2)
      Greenlandout.ice_bed_power(i)=nan;
    else
      Greenlandout.ice_bed_power(i)=nanmean(Greenland.ice_bed_power(idx1:idx2));
    end
    Greenlandout.ice_surface_power(i)=nanmean(Greenland.ice_surface_power(idx1:idx2));
    
    
  end
  % Greenland.Time=Greenland.Time;
  Greenlandout.segments_length=round((Greenland.segments_length)/numofCohInt);
  Greenlandout.index=Greenland.index;
end
end