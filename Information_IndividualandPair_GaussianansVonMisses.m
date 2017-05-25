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
    for index= 21:length(lfp)
        cells=SimultaneousRecordings(area,lfp(index,1:7));
        [spikes_times session_number]=SpikesAssembly(files,rat,area,cells,'C');
        %binary_data=BinarySpikesSignal(spikes_times);
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
            %[matrix1 matrix2]=SpikesMatrix(cell_spikes_times,register,session_number(i,:),freq);
            [pmatrix1 pmatrix2]=SpikesPhaseMatrix(phase_spike,cell_spikes_times,register,session_number(i,:),freq);                
            %recording_data{index}{i}.estimators.time_firing_rate=mean([matrix1,matrix2]);
            %recording_data{index}{i}.estimators.time_firing_rate_variance=var([matrix1,matrix2]);
            %[circular_mean,circular_variance]=CircularAnalysis([pmatrix1,pmatrix2]);
            %recording_data{index}{i}.estimators.time_firing_phase=circular_mean;
            %recording_data{index}{i}.estimators.time_firing_phase_variance=circular_variance;
            %matrix(:,:,i)=[matrix1,matrix2];
            %[pn pnt pt dn_1]=GaussianProbability(matrix(:,:,i));
            [p_theta p_thetat pt2]=VonMisesProbability2(pmatrix1,pmatrix2,20);
            %[p_theta p_theta_n ptheta_nt pntheta pntheta_t pnthetat]=GaussianVonMisesDistribution(matrix(:,:,i),pmatrix1,pmatrix2);
            for aspect=1:6
                   %[I h ht hvar]=BehaviouralInformationSquareProtocol(pn,pnt,pt,1,pos,freq,aspect-1);
                   [Ip hp htp hvarp]=BehaviouralInformationSquareProtocol(p_theta,p_thetat,pt2,1,pos,freq,aspect-1);        
                   %recording_data{index}{i}.I(aspect)=I;
                   %recording_data{index}{i}.h(aspect)=h;
                   %recording_data{index}{i}.ht(aspect)=ht;
                   %recording_data{index}{i}.h_var(aspect)=hvar;
                   recording_data{index}{i}.I_theta_10(aspect)=Ip;
                   recording_data{index}{i}.h_theta_10(aspect)=hp;
                   recording_data{index}{i}.ht_theta_10(aspect)=htp;
                   recording_data{index}{i}.h_var_theta_10(aspect)=hvarp;
            end
        end
        end
    end
        for i=1:ncells
            for j=1:i-1
                aux=corrcoef(binary_data(i,:),binary_data(j,:));
                    if size(matrix(:,:,i))==size(matrix(:,:,j))
                        [pn1n2 pn1n2t pt corr_matrix dn1 dn2]=BidimentionalGaussianProbability(matrix(:,:,i),matrix(:,:,j));
                        [pn1 pn1t]=MarginalTwoDimentional(pn1n2,pn1n2t,pt,1);
                        recording_data{index}{i}.corr_in_time(j,:)=corr_matrix;
                        for aspect=1:6
                              [I0 h0 h0t hj]=BehaviouralPairInformationSquareProtocol(pn1n2,pn1n2t,pt,1,pos,freq,aspect-1);
                              [I1 h1 h1t h1j]=BehaviouralInformationSquareProtocol(pn1,pn1t,pt,1,pos,freq,aspect-1);
                              recording_data{index}{i}.I_marginal(j,aspect)=I1;
                              recording_data{index}{i}.h_marginal(j,aspect)=h1;
                              recording_data{index}{i}.ht_marginal(j,aspect)=h1t;
                              recording_data{index}{i}.h_var_marginal(j,aspect)=h1j;
                              recording_data{index}{i}.I_pair(j,aspect)=I0;
                              recording_data{index}{i}.h_pair(j,aspect)=h0;
                              recording_data{index}{i}.ht_pair(j,aspect)=h0t;
                              recording_data{index}{i}.h_var_pair(j,aspect)=hj;
                        end
                    end
                    recording_data{index}{i}.corr(j)=aux(1,2);
            end
        end
    end
    name=strcat(directory,areaname,'_complete');
    save(name,'recording_data','-mat')
%% 