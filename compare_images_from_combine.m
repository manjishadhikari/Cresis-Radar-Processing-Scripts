if 1
comb=load('X:\ct_data\rds\2011_Greenland_P3\CSARP_manjish\20110429_01\Data_20110429_01_022.mat');
img1=load('X:\ct_data\rds\2011_Greenland_P3\CSARP_manjish\20110429_01\Data_img_01_20110429_01_022.mat');
img2=load('X:\ct_data\rds\2011_Greenland_P3\CSARP_manjish\20110429_01\Data_img_02_20110429_01_022.mat');
layer=load('X:\ct_data\rds\2011_Greenland_P3\CSARP_layerData\20110429_01\Data_20110429_01_022.mat');

end

sf=interp1(layer.GPS_time,layer.layerData{1}.value{2}.data,comb.GPS_time,'linear');   %Surface twtt
bott=interp1(layer.GPS_time,layer.layerData{2}.value{2}.data,comb.GPS_time,'linear');   %Surface twtt

figure(1);imagesc([], comb.Time,lp(comb.Data))
hold on; plot(sf)

sf_power_comb=get_layer_power(comb,sf);
bott_power_comb=get_layer_power(comb,bott);

 figure(2);imagesc([], img1.Time,lp(img1.Data))
hold on; plot(sf)
sf_power_img1=get_layer_power(img1,sf);

figure(3);imagesc([], img2.Time,lp(img2.Data))
hold on; plot(sf)
sf_power_img2=get_layer_power(img2,sf);
bott_power_img2=get_layer_power(img2,bott);

figure;plot(lp(sf_power_img1));
hold on; plot(lp(sf_power_img2));
sf_power_diff=(lp(sf_power_img1)-lp(sf_power_img2));

bott_power_updated=lp(bott_power_img2)+sf_power_diff;
figure;plot(bott_power_updated);

