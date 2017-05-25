function  [matrix_r,matrix_l]=SpikesPhaseMatrix(phase_spike,spikes_times,register,session_numbers,freq)
    
    [fil sessions]=size(session_numbers);
    [index nprot]=size(spikes_times);
    matrix_r=[];
    matrix_l=[];
    first=0;
    for number=1:sessions
         pos=Loading_Pos(register,session_numbers(number));
         tf = isfield(pos.data, 'log');
         if tf==1
            [s1 s2]=ReadAlarmTime(pos);
            %[sr sl]=ReadStartTime(pos);
            s1=floor(s1)+1;
            s2=floor(s2)+1;
            %sr=floor(sr)+1;
            %sl=floor(sl)+1;
            %[p1 p2]=ReadEndTime(pos);
            %e1=sr+min(p1-sr);
            %e2=sl+min(p2-sl);
            %totaltime1=min(p1-sr);
            %totaltime2=min(p2-sl);
            totaltime1=48;
            totaltime2=48;
            e1=s1+totaltime1;
            e2=s2+totaltime2;
            
            %aux1=min(sr-s1);
            %aux2=min(sl-s2);
            %n1=floor(floor(totaltime1)*freq+floor(aux1)*freq);
            %n2=floor(floor(totaltime2)*freq+floor(aux2)*freq);
            %n=min(n1,n2);
            n=floor(floor(totaltime1)*freq);
            matrix1=zeros(length(s1),n,150);
            matrix2=zeros(length(s2),n,150);
            time1=zeros(length(s1),n);
            time2=zeros(length(s2),n);

    
            for i=1:length(s1)
                time1(i,:)=s1(i):1/freq:s1(i)+n/freq-1/freq;
            end

            for i=1:length(s1)
                time2(i,:)=s2(i):1/freq:s2(i)+n/freq-1/freq;
            end

            for i=1:length(s1)
                for j=1:n-1
                    contador=1;
                    for k=1:length(spikes_times{number})
                        if spikes_times{number}(k)>time1(i,j) && spikes_times{number}(k)<time1(i,j+1)
                            matrix1(i,j,contador)=phase_spike{number}(k);
                            contador=contador+1;
                        end
                    end
                end
            end
    
            for i=1:length(s2)
                for j=1:n-1
                    contador=1;
                    for k=1:length(spikes_times{number})
                        if spikes_times{number}(k)>time2(i,j) && spikes_times{number}(k)<time2(i,j+1)
                            matrix2(i,j,contador)=phase_spike{number}(k);
                            contador=contador+1;
                        end
                    end
                end        
            end
            if number==first 
             matrix_r=matrix1;
             matrix_l=matrix2;
            else
              matrix_r=[matrix_r;matrix1];
              matrix_l=[matrix_l;matrix2];
            end
         else
             first=first+1;
         end
    end
end