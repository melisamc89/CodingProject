function PhaseAndTemporalHistograms3(directory,files,rat,list,area,n_fs,low1,high1,low2,high2,freq,legends)

    n=size(list);
    for index=1:n 
        [spikes_times session_number]=SpikesAssembly(files,rat,area,list(index,:),'C');
        [ncells nse]=size(session_number);
        matrix=[];
        phasematrix=[];
        if length(spikes_times)
            for se=1:nse
            order = Wished_Register_Order(area,list(index,:));
            register=Rat_Register(files,area,rat,order);
            [eeg fs]=ReadEEG(register,session_number(se));
            eeg_vector=resample(eeg,1,floor(fs/n_fs));
            n_fs=fs/floor(fs/n_fs);
            [fase1 energy]=Hilbert_phases_energy(eeg_vector,n_fs,low1,high1);
            [fase2 energy]=Hilbert_phases_energy(eeg_vector,n_fs,low2,high2);
            phase1{se}=fase1;
            phase2{se}=fase2;
            end
        phase_spike1=Phase_Calculator(spikes_times,phase1,1);
        phase_spike2=Phase_Calculator(spikes_times,phase2,1);
        
        [matrix1,matrix2]=SpikesPhaseMatrix(phase_spike1{1},spikes_times{1},register,session_number,freq);     
        [matrix1_2,matrix2_2]=SpikesPhaseMatrix(phase_spike2{1},spikes_times{1},register,session_number,freq);     
        [theta1 variance1]=CircularAnalysisStimators(matrix1,matrix2);
        [theta2 variance2]=CircularAnalysisStimators(matrix1_2,matrix2_2);
        [matrix_1 matrix_2]=SpikesMatrix(spikes_times{1},register,session_number,freq);
        [nr sigmar nl sigmal]=GaussianAnalysisStimators(matrix_1,matrix_2);
 
        signal=[nr,nl];
        error=[sigmar,sigmal];
        time1=0:1/freq:length(signal)/freq-1/freq;
       
        figure(index)
        subplot(3,1,1)
        plot(time1,signal,'k','Linewidth',2)
        legend(num2str(legends(1)))
        %imagesc(pnt)
        axis([0 100 0 max(signal)])
        title(strcat('SPIKE COUNT  ',int2str(length(spikes_times)),list(index,:)));
        xlabel('Time [s]')
        ylabel('<n>')
        time2=0:1/freq:length(theta1)/freq-1/freq;
        subplot(3,1,2)
        errorbar(time2,theta1,variance1,'r','Linewidth',2);
        axis([0 100 -pi pi])
        xlabel('Time [s]')
        ylabel('Phase_d_e_l_t_a[rad]')
        legend(num2str(legends(3)))
        subplot(3,1,3)
        errorbar(time2,theta2,variance2,'b','Linewidth',2);
        axis([0 100 -pi pi])
        xlabel('Time [s]')
        ylabel('Phase_t_h_e_t_a [rad]')
        legend(num2str(legends(2)))
        file_name=strcat(directory,list(index,:));
        saveas(index,file_name,'png');
        close(index);                    
        %mat=[time1' signal' error' theta1' variance1' theta2' variance2'];
        %name1=strcat(directory,list(index,:),'-dat');
        %save(name1,'theta','variance','-mat');
        %save(name1,'mat','-ASCII');
        end
    end
end