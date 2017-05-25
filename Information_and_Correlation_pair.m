%This function computes the mutual information between pair of neurons and
%their correlation

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
    for index= [60,62:74]
        index
        cells=SimultaneousRecordings(area,lfp(index,1:7));
        [spikes_times session_number]=SpikesAssembly(files,rat,area,cells,'C');
        size(spikes_times)
        binary_data=BinarySpikesSignal(spikes_times);
        [ncells nse]=size(session_number);
        corr_and_info{index}.corr_2=corrcoef(binary_data');
    end
    
    for index= [8:60,62:length(lfp)]
        cells=SimultaneousRecordings(area,lfp(index,1:7));
        [spikes_times session_number]=SpikesAssembly(files,rat,area,cells,'C');
        [ncells nse]=size(session_number);
        index
        ncells
        if ncells
            matrix=[];
            for i=1:ncells
                order = Wished_Register_Order(area,cells(i,:));
                register=Rat_Register(files,area,rat,order);
                cell_spikes_times=spikes_times{i};
                [matrix1 matrix2]=SpikesMatrix(cell_spikes_times,register,session_number(i,:),freq);
                matrix(:,:,i)=[matrix1,matrix2];
            end
        end
        for i=1:ncells
            for j=1:ncells
                if size(matrix(:,:,i))==size(matrix(:,:,j)) 
                    if i~=j
                    pn1n2=BidimentionalGaussianProbability2(matrix(:,:,i),matrix(:,:,j));
                    I=Information2(pn1n2)
                    corr_and_info{index}.info(i,j)=I;
                    end
                end
            end
        end
    end
    name=strcat(directory,areaname,'_corr_and_info');
    save(name,'corr_and_info','-mat')