%This function computes the information of pair of cells about one behavioural aspect of
%the rat with the square protocol in it's velocity.
%The variable type is the one which determinates the behavioural aspect.
% (0) Time information
% (1) Space Information
% (2) Directional Information
% (3) Speed Information
% (4) Space Directional Information
% (5) Speed Directional Information

function recording_data=InformationPairBinarySquareProtocol(directory,files,rat,areaname,area,list,freq)

    lfp=Days(list);
    for index= 9:9
        cells=SimultaneousRecordings(area,lfp(index,1:7));
        [spikes_times session_number]=SpikesAssembly(files,rat,area,cells,'C');
        binary_data=BinarySpikesSignal(spikes_times);
        [ncells nse]=size(session_number);
        matrix=[];
        for i=1:ncells
            order = Wished_Register_Order(area,cells(i,:));
            register=Rat_Register(files,area,rat,order);
            for se=1:nse
                cell_spikes_times{se}=spikes_times{i,se};
            end
               [matrix1 matrix2]=SpikesMatrix(cell_spikes_times,register,session_number(i,:),freq);
               matrix(:,:,i)=[matrix1,matrix2];
               [p pt]=BinaryProbability(matrix(:,:,i));
               tam=size(matrix(:,:,i));
               p_t=ones(1,tam(3))/tam(3);
               for aspect=1:1
                   [I h ht hvar]=BehaviouralInformationSquareProtocol(p,pt,p_t,1,pos,freq,aspect-1);
                   recording_data{index}{i}.I(aspect)=I;
                   recording_data{index}{i}.h(aspect)=h;
                   recording_data{index}{i}.ht(aspect)=ht;
                   recording_data{index}{i}.h_var(aspect)=hvar;
               end
        end
        for i=1:ncells
            for j=1:ncells
                aux=corrcoef(binary_data(i,:),binary_data(j,:));
                if i~=j
                    if size(matrix(:,:,i))==size(matrix(:,:,j))
                        pt=zeros(1,86);
                        pt(1,1:10)=0.1;
                        [p12 p12t]=BidimentionalBinaryProbability(matrix(:,:,i),matrix(:,:,j));
                        [p1 p1t]=MarginalTwoDimentional(p12,p12t,pt,1);
                        recording_data{index}{i}.corr_in_time(j,:)=corr_matrix;
                        for aspect=1:1
                            [I0 h0 h0t hj]=BehaviouralPairInformationSquareProtocol(p12,p12t,pt,1,pos,freq,aspect-1);
                            [I1 h1 h1t h1j]=BehaviouralInformationSquareProtocol(p1,p1t,pt,1,pos,freq,aspect-1);
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
                end
            end
            recording_data{index}{i}.corr(j)=aux(1,2);
        end
    end 
    name=strcat(directory,areaname);
    save(name,'recording_data','-mat')
end