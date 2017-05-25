%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%UNFINISHED

function info=shannon_information(map)
binsize=.2;
map=map(:);
bins=0:binsize:50+binsize/2;
h=hist(map,bins);
h=h/sum(h);
info=nansum(h.*log(h));
stop
h