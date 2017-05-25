function values=MeanFiringRateDistributionofPopulation(files,rat,area,day,prot)

    
    cells=SimultaneousRecordings(area,day);
    spikes_times=SpikesAssembly(files,rat,area,cells,prot);
    
    [n m]=size(spikes_times);
    values=zeros(1,m);
    time=ones(1,m);
    number=zeros(1,m);
    
    for i=1:m
        if spikes_times{i}
        time(i)=max(spikes_times{i});
        number(i)=length(spikes_times{i});
        end
    end
    t_max=max(time);
    values=number/t_max;
end