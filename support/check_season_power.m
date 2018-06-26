

data2010=load(['/cresis/snfs1/scratch/manjish/new_jacobshavn/data_2010.mat']);
data2011=load(['/cresis/snfs1/scratch/manjish/new_jacobshavn/data_2011.mat']);
data2012=load(['/cresis/snfs1/scratch/manjish/new_jacobshavn/data_2012.mat']);

data2013=load(['/cresis/snfs1/scratch/manjish/new_jacobshavn/data_2013.mat']);

data2014=load(['/cresis/snfs1/scratch/manjish/new_jacobshavn/data_2014.mat']);

%   data2010_coh=load(['/cresis/snfs1/scratch/manjish/new_jacobshavn/data_2010_coh_int.mat']);
%     data2011_coh=load(['/cresis/snfs1/scratch/manjish/new_jacobshavn/data_2011_coh_int.mat']);
%
     comb_data=load(['/cresis/snfs1/scratch/manjish/new_jacobshavn/combineddata_w_idx.mat']);
     vert_data=load(['/cresis/snfs1/scratch/manjish/new_jacobshavn/verticalline_w_idx.mat']);
    cross_data=load(['/cresis/snfs1/scratch/manjish/new_jacobshavn/crossline_w_idx.mat']);

figure; histogram(lp(data2010.Greenland.ice_bed_power),100);
hold on;histogram(lp(data2011.Greenland.ice_bed_power),100);
hold on;histogram(lp(data2012.Greenland.ice_bed_power),100);
hold on;histogram(lp(data2013.Greenland.ice_bed_power),100);
hold on;histogram(lp(data2014.Greenland.ice_bed_power),100);
legend('1','2','3','4','5')
%    figure; histogram(lp(data2010_coh.Greenland.ice_bed_power),100);
%     hold on;histogram(lp(data2011_coh.Greenland.ice_bed_power),100);
%
    figure; hold on; histogram(lp(comb_data.Greenland.ice_bed_power),100);
       hold on; histogram(lp(vert_data.Greenland.ice_bed_power),100);
          hold on; histogram(lp(cross_data.Greenland.ice_bed_power),100);
%

mean_2010=nanmean(lp(data2010.Greenland.ice_bed_power))
mean_2011=nanmean(lp(data2011.Greenland.ice_bed_power))
mean_2012=nanmean(lp(data2012.Greenland.ice_bed_power))
mean_2013=nanmean(lp(data2013.Greenland.ice_bed_power))
mean_2014=nanmean(lp(data2014.Greenland.ice_bed_power))



median_2010=nanmedian(lp(data2010.Greenland.ice_bed_power))
median_2011=nanmedian(lp(data2011.Greenland.ice_bed_power))
median_2012=nanmedian(lp(data2012.Greenland.ice_bed_power))
median_2013=nanmedian(lp(data2013.Greenland.ice_bed_power))
median_2014=nanmedian(lp(data2014.Greenland.ice_bed_power))

%    mean_2010_coh=nanmean(lp(data2010_coh.Greenland.ice_bed_power))
%   mean_2011_coh=nanmean(lp(data2011_coh.Greenland.ice_bed_power))
%   median_2010_coh=nanmedian(lp(data2010_coh.Greenland.ice_bed_power))
%   median_2011_coh=nanmedian(lp(data2011_coh.Greenland.ice_bed_power))

mean_comb=nanmean(lp(comb_data.Greenland.ice_bed_power))
median_comb=nanmedian(lp(comb_data.Greenland.ice_bed_power))
mean_vert=nanmean(lp(vert_data.Greenland.ice_bed_power))
mean_cross=nanmean(lp(cross_data.Greenland.ice_bed_power))
median_vert=nanmedian(lp(vert_data.Greenland.ice_bed_power))
median_cross=nanmedian(lp(cross_data.Greenland.ice_bed_power))
