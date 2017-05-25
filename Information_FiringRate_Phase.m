%This function computes the information between firing rate and phase of
%the spikes.
low=6;
high=12;
n_fs=4800;
load('3CA1','-mat');
areaname='CA1';
area=CA1_3;
list=area;
files=Data_Listing();
rat=3;
%load('CA1_complete','-mat');
freq=1;

    lfp=Days(list);
    for index= 8:length(lfp)
        cells=SimultaneousRecordings(area,lfp(index,1:7));
        [spikes_times session_number]=SpikesAssembly(files,rat,area,cells,'C');
        [ncells nse]=size(session_number);
        index
        ncells
        if ncells
            matrix=[];
            phasematrix=[];
            for se=1:nse
                order = Wished_Register_Order(area,lfp(index,:));
                register=Rat_Register(files,area,rat,order);
                [eeg fs]=ReadEEG(register,session_number(1,se));
                eeg_vector=resample(eeg,1,floor(fs/n_fs));
                n_fs=fs/floor(fs/n_fs);
                [fase energy]=Hilbert_phases_energy(eeg_vector,n_fs,low,high);
                phase{se}=fase;
            end
            pos=Loading_Pos(register,session_number(1,1));
            for i=1:ncells
                recording_data{index}{i}.name=cells(i,:);
                order = Wished_Register_Order(area,cells(i,:));
                register=Rat_Register(files,area,rat,order);
                cell_spikes_times=spikes_times{i};
                phase_spike=Phase_Calculator(cell_spikes_times,phase,4800);
                [matrix1 matrix2]=SpikesMatrix(cell_spikes_times,register,session_number(i,:),freq);
                [pmatrix1 pmatrix2]=SpikesPhaseMatrix(phase_spike,cell_spikes_times,register,session_number(i,:),freq);               
                matrix(:,:,i)=[matrix1,matrix2];
                [p_theta p_theta_n ptheta_nt pntheta pntheta_t pnthetat pn pt]=GaussianVonMisesDistribution(matrix(:,:,i),pmatrix1,pmatrix2);
                [I hy hyx hx]=Information(p_theta,pn,p_theta_n);        
                recording_data{index}{i}.I_firing_phase=I;
                recording_data{index}{i}.h_phase=hy;
                recording_data{index}{i}.ht_phase_firingrate=hyx;
                recording_data{index}{i}.h_firingrate=hx;
                 %for aspect=1:6
                  %    [I0 h0 h0t hj]=BehaviouralPairInformationSquareProtocol(pntheta,pntheta_t,pt,1,pos,freq,aspect-1);
                   %   recording_data{index}{i}.I_np_time(aspect)=I0;
                    %  recording_data{index}{i}.h_np(aspect)=h0;
                     % recording_data{index}{i}.h_np_t(aspect)=h0t;
                      %recording_data{index}{i}.h_t(aspect)=hj;
                 %end
            end
        end
    end
    name=strcat(directory,areaname,'_firingrate_phase');
    save(name,'recording_data','-mat')
%% 