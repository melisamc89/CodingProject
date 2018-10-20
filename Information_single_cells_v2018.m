%This function computes the information of pair of cells about one behavioural aspect of
%the rat with the square protocol in it's velocity.
%The variable type is the one which determinates the behavioural aspect.
% (0) Time information
% (1) Space Information
% (2) Directional Information
% (3) Speed Information
% (4) Space Directional Information
% (5) Speed Directional Information

clear all
directory='/home/melisa/Escritorio/Melisa/Doctorado/Information/';
%Rat for the experiment
rat=3;
%setting parameters for filters
low1=1.5;
low2=4;
high1=5;
high2=12;
n_fs=4800;
%loading avaiable data and setting names
%for CA1
load('3CA1','-mat');
areaname='CA1';
area=CA1_3;
%for EC1
load('3EC1','-mat');
areaname='EC1';
area=EC1_3;
list=area;
%Loading the avaible files
files=Data_Listing();

%SET FREQUENCY FOR THE IN TIME HISTOGRAMS
freq=1;
%creating the list of different LFP of the recordings. It is not the same
%as the number of recodings and cells.
lfp=Days(list);

for index= 1:length(lfp) %goes for every recorded session
        count=0; 
        order = Wished_Register_Order(area,lfp(index,:)); %choose the corresponding day of register
        register=Rat_Register(files,area,rat,order); %load the day register
        cells=SimultaneousRecordings(area,lfp(index,1:7)); %create a list of cells that where measured the same day
        [spikes_times session_number]=SpikesAssembly(files,rat,area,cells,'C'); % load spike times of all cells measured
        %in protocol 'C', measured the same day.
        [ncells nse]=size(session_number); %number of cells and repetitions of protocol C
        if ncells
            for se=1:nse
                pos=Loading_Pos(register,session_number(1,se));
                tf = isfield(pos.data, 'log');
                if tf==1
                    [eeg fs]=ReadEEG(register,session_number(1,se));
                    eeg_vector=resample(eeg,1,floor(fs/n_fs));
                    n_fs=fs/floor(fs/n_fs);
                    [qual4,qual,qual2]=egf_theta_quality_2016(eeg);
                    if qual>0.5
                        [fase1 energy1]=Hilbert_phases_energy(eeg_vector,n_fs,low1,high1);
                        phase1=fase1;
                        [fase2 energy2]=Hilbert_phases_energy(eeg_vector,n_fs,low2,high2);
                        phase2=fase2;
                        for i=1:ncells
                            count=count+1;
                            data{index}{count}.name=cells(i,:);
                            cell_spikes_times=spikes_times{i};
                            %phase_spike=Phase_Calculator(cell_spikes_times,phase,4800);
                            [matrix1 matrix2]=SpikesMatrix_v2018(cell_spikes_times,pos,freq,42);
                            matrix=[matrix1,matrix2];
                            %[matrix1_raster matrix2_raster]=SpikesMatrixRaster_v2018(stimes,pos,42);
                            %filename=strcat(directory,'Raster/',areaname,'_',cells(i,:));
                            %save(filename,'matrix_raster','-ASCII');
                            
                            [pmatrix1_d pmatrix2_d]=SpikesPhaseMatrix_v2018(phase1,cell_spikes_times,pos,freq,42);
                            [pmatrix1_t pmatrix2_t]=SpikesPhaseMatrix_v2018(phase2,cell_spikes_times,pos,freq,42);
                            
                            %compute mean firing rate parameters
                            data{index}{count}.estimators.mean_firing_rate=sum(sum([matrix1 matrix2]'))/(size(matrix1,1)*42);
                            data{index}{count}.estimators.time_firing_rate=mean([matrix1,matrix2]);
                            data{index}{count}.estimators.time_firing_rate_variance=var([matrix1,matrix2]);
                            
                            %compute mean_circular_variance_parameters
                            %delta
                            [circular_mean_d,circular_variance_d]=CircularAnalysis_v2018([pmatrix1_d,pmatrix2_d]);
                            pmatrix_d=[pmatrix1_d pmatrix2_d];
                            aux=reshape(pmatrix_d,[size(pmatrix_d,1)*size(pmatrix_d,2)*size(pmatrix_d,3) 1]);
                            aux2=aux(find(aux));   
                            data{index}{count}.estimators.mean_phase_delta=circ_mean(aux2);
                            data{index}{count}.estimators.var_phase_delta=circ_var(aux2);
                            data{index}{count}.estimators.time_phase_delta=circular_mean_d;
                            data{index}{count}.estimators.time_phase_variance_delta=circular_variance_d;
                            
                            %compute mean_circular_variance_parameters
                            %theta
                            [circular_mean_t,circular_variance_t]=CircularAnalysis_v2018([pmatrix1_t,pmatrix2_t]);
                            pmatrix_t=[pmatrix1_t pmatrix2_t];
                            aux=reshape(pmatrix_t,[size(pmatrix_t,1)*size(pmatrix_t,2)*size(pmatrix_t,3) 1]);
                            aux2=aux(find(aux));   
                            data{index}{count}.estimators.mean_phase_theta=circ_mean(aux2);
                            data{index}{count}.estimators.var_phase_theta=circ_var(aux2);
                            data{index}{count}.estimators.time_phase_theta=circular_mean_t;
                            data{index}{count}.estimators.time_phase_variance_theta=circular_variance_t;
                            
                            
                            [pn pnt pt]=GaussianProbability_v2018(matrix);
                            [p_delta p_deltat pt2]=VonMisesProbability_v2018(matrix,pmatrix_d,20);
                            [p_theta p_thetat pt2]=VonMisesProbability_v2018(matrix,pmatrix_t,20);

                            %[p_theta p_theta_n ptheta_nt pntheta pntheta_t pnthetat]=GaussianVonMisesDistribution(matrix(:,:,i),pmatrix1,pmatrix2);
                            for aspect=1:6
                                   %firing_rate
                                   [I h ht hvar]=BehaviouralInformationSquareProtocol(pn,pnt,pt,1,pos,freq,aspect-1);
                                   data{index}{count}.I(aspect)=I;
                                   data{index}{count}.h(aspect)=h;
                                   data{index}{count}.ht(aspect)=ht;
                                   data{index}{count}.h_var(aspect)=hvar;
                                   data{index}{count}.I_norm(aspect)=I/min(h,hvar); 
                                   
                                   %delta phase
                                   [I h ht hvar]=BehaviouralInformationSquareProtocol(p_delta,p_deltat,pt2,1,pos,freq,aspect-1);        
                                   data{index}{count}.I_delta(aspect)=I;
                                   data{index}{count}.h_delta(aspect)=h;
                                   data{index}{count}.ht_delta(aspect)=ht;
                                   data{index}{count}.h_var_delta(aspect)=hvar;
                                   data{index}{count}.I_delta_norm(aspect)=I/min(h,hvar);
                                   
                                   %theta phase
                                   [I h ht hvar]=BehaviouralInformationSquareProtocol(p_theta,p_thetat,pt2,1,pos,freq,aspect-1);        
                                   data{index}{count}.I_theta(aspect)=I;
                                   data{index}{count}.h_theta(aspect)=h;
                                   data{index}{count}.ht_theta(aspect)=ht;
                                   data{index}{count}.h_var_theta(aspect)=hvar;
                                   data{index}{count}.I_theta_norm(aspect)=I/min(h,hvar);
                            end
                            for test=1:100
                                new_matrix=RandomSpikesTimes_v2018(matrix);
                                [pn pnt pt]=GaussianProbability_v2018(new_matrix);
                                
                                new_pmatrix_d=RandomSpikesPhase_v2018(pmatrix_d);
                                [p_delta p_deltat pt2]=VonMisesProbability_v2018(matrix,new_pmatrix_d,20);
                                
                                new_pmatrix_t=RandomSpikesPhase_v2018(pmatrix_t);
                                [p_theta p_thetat pt2]=VonMisesProbability_v2018(matrix,new_pmatrix_t,20);
                                for aspects=1:6
                                       [I h ht hvar]=BehaviouralInformationSquareProtocol(pn,pnt,pt,1,pos,freq,aspect-1);
                                       data{index}{count}.I_test(aspect,test)=I;
                                       data{index}{count}.h_test(aspect,test)=h;
                                       data{index}{count}.ht_test(aspect,test)=ht;
                                       data{index}{count}.h_var_test(aspect,test)=hvar;
                                       data{index}{count}.I_norm_test(aspect,test)=I/min(h,hvar); 
                                       
                                       
                                       [I h ht hvar]=BehaviouralInformationSquareProtocol(p_delta,p_deltat,pt2,1,pos,freq,aspect-1);        
                                       data{index}{count}.I_test(aspect,test)=I;
                                       data{index}{count}.h_test(aspect,test)=h;
                                       data{index}{count}.ht_test(aspect,test)=ht;
                                       data{index}{count}.h_var_test(aspect,test)=hvar;
                                       data{index}{count}.I_norm_test(aspect,test)=I/min(h,hvar); 
                                       
                                       [I h ht hvar]=BehaviouralInformationSquareProtocol(p_theta,p_thetat,pt2,1,pos,freq,aspect-1);        
                                       data{index}{count}.I_test(aspect,test)=I;
                                       data{index}{count}.h_test(aspect,test)=h;
                                       data{index}{count}.ht_test(aspect,test)=ht;
                                       data{index}{count}.h_var_test(aspect,test)=hvar;
                                       data{index}{count}.I_norm_test(aspect,test)=I/min(h,hvar);    
                                   end
                            end
                            save_matrix(:,:,count)=matrix;
                        end
                    end
                end
            end
        end
end
                        
name=strcat(directory,areaname,'_complete');
save(name,'recording_data','-mat')