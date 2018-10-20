function  [matrix1,matrix2]=SpikesPhaseMatrix_v2018(phase,stimes,pos,freq,totaltime)
 
   [s1 s2]=ReadAlarmTime(pos); %read alarm time
   s1=floor(s1);s1(1)=s1(1)+1; %add 1 to possible zeros R 
   s2=floor(s2);s2(1)=s2(1)+1; %add 1 to possible zeros L
   e1=s1+totaltime; % end of experiment R
   e2=s2+totaltime; %end of experiment L
   n=floor(floor(totaltime)*freq); % size of matrix
   matrix1=zeros(length(s1),n-1);
   matrix2=zeros(length(s2),n-1);
   time1=zeros(length(s1),n);
   time2=zeros(length(s2),n);
   matrix1=zeros(length(s1),n-1,150);
   matrix2=zeros(length(s2),n-1,150);
   phase_spike=phase(floor(stimes{1}*4800));
   
   for i=1:length(s1)
       time1(i,:)=s1(i):1/freq:s1(i)+n/freq-1/freq;
   end
   for i=1:length(s1)
       time2(i,:)=s2(i):1/freq:s2(i)+n/freq-1/freq;
   end
   
   for i=1:length(s1)
       for j=1:n-1
           index=find(stimes{1}>time1(i,j) & stimes{1}<time1(i,j+1));
           spike_phase=phase_spike(index);
           matrix1(i,j,1:length(spike_phase))=spike_phase;
       end
   end
    
   for i=1:length(s2)
       for j=1:n-1
           index=find(stimes{1}>time2(i,j) & stimes{1}<time2(i,j+1));
           spike_phase=phase_spike(index);
           matrix2(i,j,1:length(spike_phase))=spike_phase;
       end
   end

end