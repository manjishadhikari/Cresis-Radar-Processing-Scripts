  
function [ice_layer_power]=get_layer_power(data,layer_twtt)

    dt= data.Time(2)-data.Time(1);
%     dh=dt*c/sqrt(3.14);
    %index_sf = round((layer_twtt-data.Time(1))/dt);
    index_sf=round(interp1(data.Time,1:size(data.Data,1),layer_twtt,'linear','extrap'));
    ice_layer_power  = zeros(1,length(layer_twtt));


for i = 1:length(layer_twtt)
      if isnan(layer_twtt(i)) || isinf(layer_twtt(i))
        ice_layer_power(i)= nan;
        continue
      else
        [layer_power, idx] = max(lp(data.Data(index_sf(i)-0:index_sf(i)+0,i)));
      %  [layer_power idx] = max(sqrt(data.data(index(i)-5:index(i)+5,i).*conj(data.data(index(i)-5:index(i)+5,i))));
        layer_index = idx + index_sf(i)+0-1;
        if layer_power  == 0
          ice_layer_power(i)= nan;
          layer_twtt(i) = nan;
        else
          ice_layer_power(i) = data.Data(layer_index,i);
          layer_twtt(i) =  interp1([1:length(data.Time)],data.Time,layer_index);
        end
      end
end
    return 