function PhaseAndTemporalHistograms(directory,files,rat,list,area,n_fs,low,high,freq)

    n=size(list);
    for index= 1:n 
        [spikes_times session_number]=SpikesAssembly(files,rat,area,list(index,:),'C');
        [ncells nse]=size(session_number);
        matrix=[];
        phasematrix=[];
        for se=1:nse
            order = Wished_Register_Order(area,list(index,:));
            register=Rat_Register(files,area,rat,order);
            [eeg fs]=ReadEEG(register,session_number(se));
            eeg_vector=resample(eeg,1,floor(fs/n_fs));
            n_fs=fs/floor(fs/n_fs);
            [fase energy]=Hilbert_phases_energy(eeg_vector,n_fs,low,high);
            phase{se}=fase;
        end
        phase_spike=Phase_Calculator(spikes_times{1},phase,1);
        [matrix1,matrix2]=SpikesPhaseMatrix(phase_spike,spikes_times{1},register,session_number(1,:),freq);     
        [theta variance]=CircularAnalysisStimators(matrix1,matrix2);
        [matrix_1 matrix_2]=SpikesMatrix(spikes_times{1},register,session_number(1,:),freq);
        [nr sigmar nl sigmal]=GaussianAnalysisStimators(matrix_1,matrix_2);
        
        %matrix=[matrix_1,matrix_2];
        %[pn pnt pt dn_1]=GaussianProbability(matrix);
        %[p_theta p_thetat pt2]=VonMisesProbability2(matrix1,matrix2,20);
                    
        signal=[nr,nl];
        error=[sigmar,sigmal];
        time1=0:1/freq:length(signal)/freq-1/freq;
        figure(index)
        subplot(2,1,1)
        plot(time1,signal,'b','Linewidth',2)
        %imagesc(pnt)
        axis([0 85 0 max(signal)])
        title(strcat('SPIKE COUNT  ',int2str(length(spikes_times)),area(index,:)));
        xlabel('Time [s]')
        ylabel('<n>')
        time2=0:1/freq:length(theta)/freq-1/freq;
        subplot(2,1,2)
        %theta=theta*360/(2*pi);
        %theta=mod(theta,180);
        errorbar(time2,theta,variance,'g','Linewidth',2);
        %plot(time2,theta,'g','Linewidth',2);
        axis([0 85 -pi pi])
        %imagesc(p_thetat)
        title(strcat('PHASE',list(index,:)));
        xlabel('Time [s]')
        ylabel('Phase [rad]')
        file_name=strcat(directory,list(index,:),'-',int2str(low),'-',int2str(high));
        %file_name=strcat(directory,list(index,:),'-',int2str(low),'-',int2str(high),'Probability');
        saveas(index,file_name,'png');
        close(figure(index));                    
        mat=[time1' signal' error' theta' variance'];
        name1=strcat(directory,list(index,:),'-',int2str(low),'-',int2str(high),'-dat');
        %save(name1,'theta','variance','-mat');
        save(name1,'mat','-ASCII');
    end
end
