function matrix=BinarySpikesSignal(spikes_times)

    [n m]=size(spikes_times);
    time=zeros(n,m);
    
    for k=1:n
        for cell=1:m
            if length(spikes_times{k,cell})
                time(k,cell)=max(spikes_times{k,cell}{1});
            end
        end
    end
    total_time=floor(max(max(time))+1);
    matrix=zeros(m,total_time);
    
    for cell=1:m
        matrix(cell,floor(spikes_times{k,cell}{1})+1)=1;
    end
end