
function [Greenland,bed_rms,bed_corr_power]=bed_roughness(Greenland,settings)

bed_corr_power=[];
M=settings.M;
M1=settings.M1;
physical_constants;
param.radar.fc=195000000;
er_bed=6;

if settings.cross_lines==1
  if strcmp(settings.location,'Jacobshavn')
    if exist((  (['/cresis/snfs1/scratch/manjish/new_jacobshavn/bed_roughness/crossline',num2str(M),'.mat'])),'file')
      load(['/cresis/snfs1/scratch/manjish/new_jacobshavn/bed_roughness/crossline',num2str(M),'.mat']);
      file_exist = true;
    end
  elseif strcmp(settings.location,'Peterman')
    if exist((  (['/cresis/snfs1/scratch/manjish/new_peterman/bed_roughness/crossline',num2str(M),'.mat'])),'file')
      load(['/cresis/snfs1/scratch/manjish/new_peterman/bed_roughness/crossline',num2str(M),'.mat']);
      file_exist = true;
    end
  end
  
else
  if strcmp(settings.location,'Jacobshavn')
    if exist((  (['/cresis/snfs1/scratch/manjish/new_jacobshavn/bed_roughness/verticalline',num2str(M1),'.mat'])),'file')
      load(['/cresis/snfs1/scratch/manjish/new_jacobshavn/bed_roughness/verticalline',num2str(M1),'.mat']);
      file_exist = true;
    end
  elseif strcmp(settings.location,'Peterman')
    if exist((  (['/cresis/snfs1/scratch/manjish/new_peterman/bed_roughness/verticalline',num2str(M1),'.mat'])),'file')
      load(['/cresis/snfs1/scratch/manjish/new_peterman/bed_roughness/verticalline',num2str(M1),'.mat']);
      file_exist = true;
    end
  end
end

%Bed roughness calculation
settings.num_int=1000;
settings.repeat_after=10;
settings.type='bed';
settings.save_en=0;
[r]=roughness_calculation(Greenland,settings);

%Bed roughness correction
k=1;
num_int=r.num_int;
repeat_after=r.repeat_after;
for l = num_int/2:repeat_after:length(Greenland.ice_bed_power_avg)
  if ((l >= num_int/2) && ((l+num_int/2) < length(Greenland.ice_bed_power_avg)))
    %         if k>size(r.lat)
    %           continue
    %         end
    if isnan(r.rms_height(k)) || isnan(Greenland.ice_bed_power_avg(k))
      Greenland.ice_bed_power_avg(k)=nan;
      bed_corr_power(k)=nan;
      k= k+1;
      continue ;
    else
      
      bed_corr_power(k) =(exp(-(4*pi*r.rms_height(k)/(c/(param.radar.fc*(sqrt(er_bed)-1))))^2)*(besseli(0,((4*pi*r.rms_height(k)/(c/(param.radar.fc)*(sqrt(er_bed-1)))^2)/2))^2));
      Greenland.ice_bed_power_avg(k) = Greenland.ice_bed_power_avg(k).*(exp(-(4*pi*r.rms_height(k)/(c/(param.radar.fc*(sqrt(er_bed)-1))))^2)*(besseli(0,((4*pi*r.rms_height(k)/(c/(param.radar.fc)*(sqrt(er_bed-1)))^2)/2))^2));
      k = k+1;
    end
  end
end
bed_rms=r;
clearvars r k

end