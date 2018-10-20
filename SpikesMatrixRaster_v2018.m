%This function gives the spikes separated in trials repetitions.

function  [matrix1 matrix2]=SpikesMatrixRaster_v2018(stimes,pos,totaltime)

   freq=4800;
   [s1 s2]=ReadAlarmTime(pos); %read alarm time
   s1=floor(s1);s1(1)=s1(1)+1; %add 1 to possible zeros R 
   s2=floor(s2);s2(1)=s2(1)+1; %add 1 to possible zeros L
   e1=s1+totaltime; % end of experiment R
   e2=s2+totaltime; %end of experiment L
   n=floor(floor(totaltime)*freq); % size of matrix
   matrix1=[];matrix2=[];
         
   for i=1:length(s1)
       spikes=stimes{1}(find(stimes{1}>s1(i)&stimes{1}<e1(i)));
       matrix1{i}=spikes;
   end
   for i=1:length(s2)
       spikes=stimes{1}(find(stimes{1}>s2(i)&stimes{1}<e2(i)));
       matrix2{i}=spikes;
   end
   
end