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
                        for f=1:80
                        [fase1 energy1]=Hilbert_phases_energy(eeg_vector,n_fs,f,f+1);
                        phase1=fase1;
                        for i=24:24
                             aspect=1;

                            data{index}{i}.name=cells(i,:);
                            cell_spikes_times=spikes_times{i};
                            [matrix1 matrix2]=SpikesMatrix_v2018(cell_spikes_times,pos,freq,42);
                            matrix=[matrix1,matrix2];
                            [pmatrix1 pmatrix2]=SpikesPhaseMatrix_v2018(phase1,cell_spikes_times,pos,freq,42);
                            
                            [pn pnt pt]=GaussianProbability_v2018(matrix);
                            pmatrix=[pmatrix1 pmatrix2];
                            [pf pft pt2]=VonMisesProbability_v2018(matrix,pmatrix,20);

                            %firing_rate
                            [I h ht hvar]=BehaviouralInformationSquareProtocol(pn,pnt,pt,1,pos,freq,0);
                            data{index}{i}.I(aspect)=I;
                            data{index}{i}.h(aspect)=h;
                            data{index}{i}.ht(aspect)=ht;
                            data{index}{i}.h_var(aspect)=hvar;
                            data{index}{i}.I_norm(aspect)=I/min(h,hvar); 
                                   
                            % phase
                            [I h ht hvar]=BehaviouralInformationSquareProtocol(pf,pft,pt2,1,pos,freq,0);        
                            data{index}{i}.I_phase(f)=I;
                            data{index}{i}.h_phase(f)=h;
                            data{index}{i}.ht_phase(f)=ht;
                            data{index}{i}.h_var_phase(f)=hvar;
                            data{index}{i}.I_phase_norm(f)=I/min(h,hvar);
                            
                            for test=1:100
                                new_matrix=RandomSpikesTimes_v2018(matrix);
                                [pn pnt pt]=GaussianProbability_v2018(new_matrix);
                                
                                new_pmatrix=RandomSpikesPhase_v2018(pmatrix);
                                [p_delta p_deltat pt2]=VonMisesProbability_v2018(matrix,new_pmatrix,20);
                                
                                       [I h ht hvar]=BehaviouralInformationSquareProtocol(pn,pnt,pt,1,pos,freq,aspect-1);
                                       data{index}{i}.TEST.I(test)=I;
                                       data{index}{i}.TEST.h(test)=h;
                                       data{index}{i}.TEST.ht(test)=ht;
                                       data{index}{i}.TEST.h_var(test)=hvar;
                                       data{index}{i}.TEST.I_norm(test)=I/min(h,hvar); 
                                       
                                       
                                       [I h ht hvar]=BehaviouralInformationSquareProtocol(p_delta,p_deltat,pt2,1,pos,freq,aspect-1);        
                                       data{index}{i}.TEST.I_phase(f,test)=I;
                                       data{index}{i}.TEST.h_phase(f,test)=h;
                                       data{index}{i}.TEST.ht_phase(f,test)=ht;
                                       data{index}{i}.TEST.h_var_phase(f,test)=hvar;
                                       data{index}{i}.TEST.I_norm_phase(f,test)=I/min(h,hvar);      
                            end
                        end
                    end
                end
            end
            end
        end
end

name=strcat(directory,areaname);
save(name,'data','-mat')
name2=strcat(directory,areaname,'_matrix');
save(name2,'save_matrix','-mat')
