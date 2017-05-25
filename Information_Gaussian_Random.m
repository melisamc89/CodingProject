% This code computes de information about behavioural aspects of cells with
% the same firing rate as the real ones, but with random spike_times

low=6;
high=12;
n_fs=4800;
load('3CA1','-mat');
areaname='CA1';
area=CA1_3;
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
        if ncells
            for i=1:ncells
                matrix=[];
                order = Wished_Register_Order(area,cells(i,:));
                register=Rat_Register(files,area,rat,order);
                pos=Loading_Pos(register,session_number(1,1));
                for j=1:100
                    times=RandomSpikesTimes(spikes_times{i},session_number);
                    [matrix1 matrix2]=SpikesMatrix(times,register,session_number(i,:),freq);
                    matrix(:,:,i)=[matrix1,matrix2];
                    [pn pnt pt dn_1]=GaussianProbability(matrix(:,:,i));
                    for aspect=1:6
                           [a b c d]=BehaviouralInformationSquareProtocol(pn,pnt,pt,1,pos,freq,aspect-1);
                           I(j,aspect)=a;
                           h(j,aspect)=b;
                           ht(j,aspect)=c;
                           h_var(j,aspect)=d;
                    end
                end
                for aspect=1:6
                    I(:,aspect)=sort(I(:,aspect));
                end
                recording_data{index}{i}.random.I_mean=mean(I);
                recording_data{index}{i}.random.I_std=std(I);
                recording_data{index}{i}.random.I=I;
            end
        end
    end
    name=strcat(directory,areaname,'_complete');
    save(name,'recording_data','-mat')