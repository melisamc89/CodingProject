function PhaseHistograms(directory,files,rat,list,area,n_fs,low,high,freq)

    n=size(area);
    for index= 1:n 
        order=Wished_Register_Order(list,area(index,:));
        register=Rat_Register(files,list,rat,order);
        [fil col]=size(register.sessions);
        for se=1:col
            protocol=register.sessions{se}.sess_type;
            if protocol == 'C'
                spikes_times= register.sessions{se}.spike_times;
                pos=Loading_Pos(register,se);
                tf = isfield(pos.data, 'log');
                if (tf==1 && length(spikes_times)>50)
                    [eeg_vector,fs]=ReadEEG(register,se);
                    eeg_vector=resample(eeg_vector,1,floor(fs/n_fs));
                    n_fs=fs/floor(fs/n_fs);
                    [phase energy]=Hilbert_phases_energy(eeg_vector,n_fs,low,high);
                    phase_spike=Phase_Calculator(spikes_times,phase,n_fs);
                    [matrix1 matrix2]=SpikesPhaseMatrix(phase_spike,spikes_times,pos,42,freq);
                    
                    [theta variance]=CircularAnalysisStimators(matrixa,matrixb)

                    [nr sigmar nl sigmal]=FiringTimeData(spikes_times,pos,41,freq);
                    
                    time1=0:1/freq:82-1/freq;
                    signal=[nr,nl];
                    error=[sigmar,sigmal];
                    figure(index)
                    subplot(2,1,1)
                    errorbar(time1,signal,error)
                    title(strcat('SPIKE COUNT  ',int2str(length(spikes_times)),area(index,:)));
                    xlabel('Time [s]')
                    ylabel('<n>')
                    time2=0:1/freq:82-2/freq;
                    subplot(2,1,2)
                    errorbar(time2,phase,var,'r');
                    title(strcat('PHASE',area(index,:)));
                    xlabel('Time [s]')
                    ylabel('Phase [rad]')
                    file_name=strcat(directory,area(index,:),'-',int2str(low),'-',int2str(high));
                    saveas(index,file_name,'png');
                    close(figure(index));
                    
                    name1=strcat(directory,area(index,:),'-',int2str(low),'-',int2str(high));
                    save(name1,'phase','var','-mat');
                    save(name1,'phase','var','-ASCII');
                end
            end
        end
    end
end