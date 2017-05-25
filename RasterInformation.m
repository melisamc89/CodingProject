function RasterInformation(directory,files,rat,list,area,n_fs,low1,high1,low2,high2,freq,totaltime)

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
            %[fase1 energy]=Hilbert_phases_energy(eeg_vector,n_fs,low1,high1);
            %[fase2 energy]=Hilbert_phases_energy(eeg_vector,n_fs,low2,high2);
            %phase1{se}=fase1;
            %phase2{se}=fase2;
            end
        %phase_spike1=Phase_Calculator(spikes_times,phase1,1);
        %phase_spike2=Phase_Calculator(spikes_times,phase2,1);
        
        %[matrix1,matrix2]=SpikesPhaseMatrix(phase_spike1{1},spikes_times{1},register,session_number,freq);     
        %[matrix1_2,matrix2_2]=SpikesPhaseMatrix(phase_spike2{1},spikes_times{1},register,session_number,freq);     
        [matrix_1 matrix_2]=SpikesMatrix(spikes_times{1},register,session_number,freq,totaltime);
        %[n m l]=size([matrix1 matrix2]);
        %for i=1:n
        %mat2(i,:)=CircularAnalysisStimators(matrix1(i,:,:),matrix2(i,:,:));
        %end
        %[n m l]=size([matrix1_2 matrix2_2]);
        %for i=1:n
        %mat3(i,:)=CircularAnalysisStimators(matrix1_2(i,:,:),matrix2_2(i,:,:));
        %end
        
        
        mat=[matrix_1 matrix_2];
        [n m]=size(mat);
        time1=0:1/freq:m/freq-1/freq;
        figure(1)
        %imagesc(time1,1:n,mat)
        plotSpikeRaster(logical(mat),'PlotType','vertline');
        %cmap=[[1,1,1];[0,0,0]];
        %colormap(cmap)
        title(strcat('SPIKE COUNT :',int2str(length(spikes_times)),' _',list(index,:)));
        xlabel('Time (revisar unidades)')
        ylabel('trials')
        %[n m]=size(mat2);
        %time2=0:1/freq:m/freq-1/freq;
        %figure(2)
        %imagesc(time2,1:n,mat2);
        %colormap(jet)
        %xlabel('Time [s]')
        %ylabel('trials')
        %caxis([-pi pi])
        %figure(3)
        %imagesc(time2,1:n,mat3);
        %colormap(jet)
        %xlabel('Time [s]')
        %ylabel('trials')
        %caxis([-pi pi])
        
        file_name1=strcat(directory,'firing',list(index,:));
        saveas(1,file_name1,'png');              
        a=[time1' mat'];
        name1=strcat(directory,'firing',list(index,:),'.dat');
        save(name1,'a','-ASCII')
        
        %file_name2=strcat(directory,'delta_',list(index,:));
        %saveas(2,file_name2,'png');              
        %b=[time2' mat2'];
        %name2=strcat(directory,'delta_',list(index,:),'.dat');
        %save(name2,'b','-ASCII')

        %file_name3=strcat(directory,'theta_',list(index,:));
        %saveas(3,file_name3,'png');              
        %c=[time2' mat3'];
        %name3=strcat(directory,'theta_',list(index,:),'.dat');
        %save(name3,'c','-ASCII')
        
        close(1)
        %close(2)
        %close(3)
        end
    end
end