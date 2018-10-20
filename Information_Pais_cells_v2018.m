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
                               for j=i+1:ncells
                                   if size(matrix(:,:,i))==size(matrix(:,:,j))
                                      [pn1n2 pn1n2t pt corr_matrix]=BidimentionalGaussianProbability(matrix(:,:,i),matrix(:,:,j));
                                      [pn1 pnt1 pt1]=GaussianProbability_v2018(matrix(:,:,i));
                                      [pn2 pnt2 pt2]=GaussianProbability_v2018(matrix(:,:,j));                                        
                                      for aspect=1:6
                                           [I0 h0 h0t hj]=BehaviouralPairInformationSquareProtocol(pn1n2,pn1n2t,pt,1,pos,freq,aspect-1);
                                           [I1 h1 h1t hj1]=BehaviouralInformationSquareProtocol(pn1,pnt1,pt1,1,pos,freq,aspect-1);
                                           [I2 h2 h2t hj2]=BehaviouralInformationSquareProtocol(pn2,pnt2,pt2,1,pos,freq,aspect-1);
                                           data{index}{i}.Pair.I_pair(j-1,aspect)=I0;
                                           data{index}{i}.Pair.h_pair(j-1,aspect)=h0;
                                           data{index}{i}.Pair.ht_pair(j-1,aspect)=h0t;
                                           data{index}{i}.Pair.h_var_pair(j-1,aspect)=hj;
                                           data{index}{i}.Pair.I_norm(j-1,aspect)=I0/min(h0,hj);
                                           
                                           data{index}{i}.Pair.I_1(j-1,aspect)=I1;
                                           data{index}{i}.Pair.h_1(j-1,aspect)=h1;
                                           data{index}{i}.Pair.ht_1(j-1,aspect)=ht1;
                                           data{index}{i}.Pair.h_var_1(j-1,aspect)=hj1;
                                           data{index}{i}.Pair.I_norm_1(j-1,aspect)=I1/min(h1,hj1);
                                           
                                           data{index}{i}.Pair.I_2(j-1,aspect)=I2;
                                           data{index}{i}.Pair.h_2(j-1,aspect)=h2;
                                           data{index}{i}.Pair.ht_2(j-1,aspect)=ht2;
                                           data{index}{i}.Pair.h_var_2(j-1,aspect)=hj2;
                                           data{index}{i}.Pair.I_norm_2(j-1,aspect)=I2/min(h2,hj2);
                                      end
                                      for test=1:100
                                        [new_matrix1 new_matrix2]=RandomPairSpikesTimes_v2018(matrix(:,:,i),matrix(:,:,j));
                                        [pn1n2 pn1n2t pt corr_matrix]=BidimentionalGaussianProbability(new_matrix1,mew_matrix2);
                                        [pn1 pnt1 pt1]=GaussianProbability_v2018(new_matrix1);
                                        [pn2 pnt2 pt2]=GaussianProbability_v2018(new_matrix2);      
                                        for aspect=1:6
                                            [I0 h0 h0t hj]=BehaviouralPairInformationSquareProtocol(pn1n2,pn1n2t,pt,1,pos,freq,aspect-1);
                                            [I1 h1 h1t hj1]=BehaviouralInformationSquareProtocol(pn1,pnt1,pt1,1,pos,freq,aspect-1);
                                            [I2 h2 h2t hj2]=BehaviouralInformationSquareProtocol(pn2,pnt2,pt2,1,pos,freq,aspect-1);
                                              
                                            data{index}{i}.TEST.Pair.I_pair(j-1,aspect,test)=I0;
                                            data{index}{i}.TEST.Pair.h_pair(j-1,aspect)=h0;
                                            data{index}{i}.TEST.Pair.ht_pair(j-1,aspect)=h0t;
                                            data{index}{i}.TEST.Pair.h_var_pair(j-1,aspect)=hj;
                                            data{index}{i}.TEST.Pair.I_norm(j-1,aspect)=I0/min(h0,hj);
                                           
                                            data{index}{i}.TEST.Pair.I_1(j-1,aspect)=I1;
                                            data{index}{i}.TEST.Pair.h_1(j-1,aspect)=h1;
                                            data{index}{i}.TEST.Pair.ht_1(j-1,aspect)=ht1;
                                            data{index}{i}.TEST.Pair.h_var_1(j-1,aspect)=hj1;
                                            data{index}{i}.TEST.Pair.I_norm_1(j-1,aspect)=I1/min(h1,hj1);
                                           
                                            data{index}{i}.TEST.Pair.I_2(j-1,aspect)=I2;
                                            data{index}{i}.TEST.Pair.h_2(j-1,aspect)=h2;
                                            data{index}{i}.TEST.Pair.ht_2(j-1,aspect)=ht2;
                                            data{index}{i}.TEST.Pair.h_var_2(j-1,aspect)=hj2;
                                            data{index}{i}.TEST.Pair.I_norm_2(j-1,aspect)=I2/min(h2,hj2);  
                                        end
                                      end
                                   end
                               end
                        end
                    end
                end
            end
        end
end