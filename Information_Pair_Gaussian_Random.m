%This function computes the information of pair of cells about one behavioural aspect of
%the rat with the square protocol in it's velocity.
%The variable type is the one which determinates the behavioural aspect.
% (0) Time information
% (1) Space Information
% (2) Directional Information
% (3) Speed Information
% (4) Space Directional Information
% (5) Speed Directional Information

load('CA1_complete','-mat');
low=6;
high=12;
n_fs=4800;
load('3CA1','-mat');
areaname='CA1';
area=CA1_3;
list=area;
files=Data_Listing();
rat=3;
freq=1;

    lfp=Days(list);
    for index= 8:length(lfp)
        cells=SimultaneousRecordings(area,lfp(index,1:7));
        [spikes_times session_number]=SpikesAssembly(files,rat,area,cells,'C');
        binary_data=BinarySpikesSignal(spikes_times);
        [ncells nse]=size(session_number);
        index
        for aleat=1:100
            matrix=[];
            for i=1:ncells
                order = Wished_Register_Order(area,cells(i,:));
                register=Rat_Register(files,area,rat,order);
                pos=Loading_Pos(register,session_number(1,1));
                times=RandomSpikesTimes(spikes_times{i},session_number);
                [matrix1 matrix2]=SpikesMatrix(times,register,session_number(i,:),freq);
                matrix(:,:,i)=[matrix1,matrix2];
            end
            for i=1:ncells
                for j=1:ncells
                    if i~=j
                    if size(matrix(:,:,i))==size(matrix(:,:,j))
                        [pn1n2 pn1n2t pt corr_matrix dn1 dn2]=BidimentionalGaussianProbability(matrix(:,:,i),matrix(:,:,j));
                        for aspect=1:1
                            [I0 h0 h0t hj]=BehaviouralPairInformationSquareProtocol(pn1n2,pn1n2t,pt,1,pos,freq,aspect-1);
                            recording_data{index}{i}.random.I_pair(aleat,j,aspect)=I0;
                            recording_data{index}{i}.random.h_pair(aleat,j,aspect)=h0;
                            recording_data{index}{i}.random.ht_pair(aleat,j,aspect)=h0t;
                            recording_data{index}{i}.random.h_var_pair(aleat,j,aspect)=hj;
                        end    
                    end
                    end
                end
            end
        end
    end 
    name=strcat(directory,areaname);
    save(name,'recording_data','-mat')