%This function computes the information between firing rate and phase for
%each cell.

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
    for index=8:length(lfp)
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
                %recording_data{index}{i}.name=cells(i,:);
                order = Wished_Register_Order(area,cells(i,:));
                register=Rat_Register(files,area,rat,order);
                cell_spikes_times=spikes_times{i};
                phase_spike=Phase_Calculator(cell_spikes_times,phase,4800);
                [matrix1 matrix2]=SpikesMatrix(cell_spikes_times,register,session_number(i,:),freq);
                [pmatrix1 pmatrix2]=SpikesPhaseMatrix(phase_spike,cell_spikes_times,register,session_number(i,:),freq);                
                matrix=[matrix1,matrix2];
                [p_theta p_theta_n pn vec phase_mean kappa]=GaussianVonMisesDistribution2(matrix,pmatrix1,pmatrix2);
                [I hy hyx hx]=Information(p_theta,pn,p_theta_n);
                I
                data{index}{i}.I_np=I;
                data{index}{i}.hn=hx;
                data{index}{i}.hp=hy;
                data{index}{i}.hpn=hyx;
                data{index}{i}.n=vec;
                data{index}{i}.phase=phase_mean;
                data{index}{i}.phase_error=kappa;
            end
        end
    end
    name=strcat(directory,'CA1','_firing_phase');
    save(name,'data','-mat')