

function TemporalHistograms(directory,files,rat,area,list,freq,totaltime)
    
    n=size(list);
    for index=1:n 
        [spikes_times session_number]=SpikesAssembly(files,rat,area,list(index,:),'C');
        [ncells nse]=size(session_number);
        if length(spikes_times)
            for se=1:nse
            order = Wished_Register_Order(area,list(index,:));
            register=Rat_Register(files,area,rat,order);
            st= register.sessions{se}.spike_times;
                        
            [matrix_1 matrix_2]=SpikesMatrix(spikes_times{1},register,session_number,freq,totaltime);
            [nr sigmar nl sigmal]=GaussianAnalysisStimators(matrix_1,matrix_2);

            if (sum(sum(matrix_1))+sum(sum(matrix_2))>300)
            signal=[nr,nl];
            error=[sigmar,sigmal];
            time1=0:1/freq:length(nr)/freq-1/freq;
            time2=0:1/freq:length(nl)/freq-1/freq;

            figure(index)
            subplot(2,1,1)
            plot([time1,time1(end)+1/freq:1/freq:time1(end)+6*freq],[nr,zeros(1,6*freq)],'b','Linewidth',2)
            hold on
            plot([6 6],[0 max(signal)],'k','Linewidth',2')
            legend('-->','Start R')
            axis([0 totaltime+6 0 max(signal)])
            title(strcat('SPIKE COUNT : ',int2str(length(st)),'-name:',list(index,:)));
            xlabel('Time [s]')
            ylabel('<n>')
            subplot(2,1,2)
            plot([-6*freq:1/freq:0-1/freq,time2],[zeros(1,6*freq),wrev(nl(6*freq+1:end)),wrev(nl(1:6*freq))],'r','Linewidth',2);
            axis([-6 totaltime 0 max(signal)])
            hold on
            plot([42 42],[0 max(signal)],'k','Linewidth',2')
            xlabel('Time [s]')
            ylabel('<n>')
            legend('<--','Start L')
            file_name=strcat(directory,list(index,:));
            saveas(index,file_name,'png');
            close(index);   
            end
            end
        end
    end
end