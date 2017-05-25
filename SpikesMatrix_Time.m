function  [mean1 var1 mean2 var2]=MeanPhaseCode(spikes_times,pos,totaltime,fast_time,phase)

    [s1 s2]=ReadAlarmTime(pos);
    s1=floor(s1)+1;
    s2=floor(s2)+1;
    e1=s1+totaltime;
    e2=s2+totaltime;
    
    i1=s1+fast_time;
    i2=s2+(totaltime-fast_time);
    
    count1=1;
    count2=1;
    for i=1:length(spike_times)
        for j=1:length(s1)
            if spike_times(i)>s1(j) && spike_times(j)<i1(j)
                aux1(count1)=phase(i);
                count1=count+1;
            end
            if spike_times(i)>i1(j) && spike_times(j)<e1(j)
                aux2(count2)=phase(i);
                count2=count+2;
            end
            if spike_times(i)>s2(j) && spike_times(j)<i2(j)
                aux2(count2)=phase(i);
                count2=count+2;
            end
            if spike_times(i)>i2(j) && spike_times(j)<e2(j)
                aux1(count1)=phase(i);
                count1=count+1;
            end
end