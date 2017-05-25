%This function computes the information of pair of cells about one behavioural aspect of
%the rat with the square protocol in it's velocity.
%The variable type is the one which determinates the behavioural aspect.
% (0) Time information
% (1) Space Information
% (2) Directional Information
% (3) Speed Information
% (4) Space Directional Information
% (5) Speed Directional Information
low=6;
high=12;
n_fs=4800;
load('3CA1','-mat');
areaname='CA1';
area=EC1_3;
list=area;
files=Data_Listing();
rat=3;
load('CA1_complete','-mat');
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
            order = Wished_Register_Order(area,cells(i,:));
            register=Rat_Register(files,area,rat,order);
            cell_spikes_times=spikes_times{i};
            phase_spike=Phase_Calculator(cell_spikes_times,phase,4800);
            [matrix1 matrix2]=SpikesMatrix(cell_spikes_times,register,session_number(i,:),freq);
            [pmatrix1 pmatrix2]=SpikesPhaseMatrix(phase_spike,cell_spikes_times,register,session_number(i,:),freq);                
            recording_data{index}{i}.estimators.time_firing_rate=mean([matrix1,matrix2]);
            recording_data{index}{i}.estimators.time_firing_rate_variance=var([matrix1,matrix2]);
            [circular_mean,circular_variance]=CircularAnalysis([pmatrix1,pmatrix2]);
            recording_data{index}{i}.estimators.time_firing_phase=circular_mean;
            recording_data{index}{i}.estimators.time_firing_phase_variance=circular_variance;
        end
        end
    end
    name=strcat(directory,areaname,'_complete');
    save(name,'recording_data','-mat')