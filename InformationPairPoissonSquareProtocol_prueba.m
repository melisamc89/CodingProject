%This function computes the information of pair of cells about one behavioural aspect of
%the rat with the square protocol in it's velocity.
%The variable type is the one which determinates the behavioural aspect.
% (0) Time information
% (1) Space Information
% (2) Directional Information
% (3) Speed Information
% (4) Space Directional Information
% (5) Speed Directional Information

function recording_data=InformationPairPoissonSquareProtocol_prueba(directory,files,rat,areaname,area,list,freq)

    lfp=Days(list);
    for index= 9:9
        index
        cells=SimultaneousRecordings(area,lfp(index,1:7));
        spikes_times=SpikesAssembly(files,rat,area,cells,'C');
        binary_data=BinarySpikesSignal(spikes_times);
        for i=1:size(cells)
            order = Wished_Register_Order(area,cells(i,:));
            register=Rat_Register(files,area,rat,order);        
            [fil col]=size(register.sessions);
            for se=1:col
                protocol=register.sessions{se}.sess_type;
                if protocol == 'C'
                    spikes_times= register.sessions{se}.spike_times;
                    pos=Loading_Pos(register,se);
                    tf = isfield(pos.data, 'log');
                    if (tf==1 && length(spikes_times)>50)
                        [matrix1 matrix2]=SpikesMatrix(spikes_times,pos,freq);
                        matrix{:,:,i}=[matrix1,matrix2];
                        %[pn pnt pt dn_1]=GaussianProbability(matrix1,matrix2,times);
                        for j=1:i-1
                            aux=corrcoef(binary_data(i,:),binary_data(j,:));
                            if size(matrix{:,:,i})==size(matrix{:,:,j})
                                [pn1n2 pn1n2t pt corr_matrix dn1 dn2 times]=BidimentionalGaussianProbability_prueba(matrix{:,:,i},matrix{:,:,j});
                                [pn1 pn1t]=MarginalTwoDimentional(pn1n2,pn1n2t,pt,1);
                                [pn pnt pt dn_1]=GaussianProbability_prueba(matrix1,matrix2,times);
                                for aspect=1:1
                                    [I h ht hvar]=BehaviouralInformationSquareProtocol(pn,pnt,pt,1,pos,freq,aspect-1);
                                     recording_data{index}{i}.I(j,aspect)=I;
                                     recording_data{index}{i}.h(j,aspect)=h;
                                     recording_data{index}{i}.ht(j,aspect)=ht;
                                     recording_data{index}{i}.h_var(j,aspect)=hvar;
                                end   
                                %name=strcat('C:\Users\meli__000\Desktop\Melisa\Doctorado\CellsCorrelation\Pn-Pn_marginal2\',int2str(i),'_',int2str(j));
                                %figure(j)
                                %subplot(2,10,[1:9])
                                %imagesc(pnt)
                                %subplot(2,10,10)
                                %imagesc(pn)
                                %subplot(2,10,[11:19])
                                %imagesc(pn1t)
                                %subplot(2,10,20)
                                %imagesc(pn1')
                                %saveas(j,name,'png')
                                %close(j)
                                %prueba=prueba+pn1t;
                                recording_data{index}{i}.corr_in_time(j,:)=corr_matrix;
                                for aspect=1:1
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
                        %figure(100)
                        %name=strcat('C:\Users\meli__000\Desktop\Melisa\Doctorado\CellsCorrelation\Pn-Pn_marginal\suma_',int2str(i));
                        %subplot(2,1,1)
                        %imagesc(pnt)
                        %subplot(2,1,2)
                        %imagesc(prueba)
                        %saveas(100,name,'png')
                        %close(100)
                    end
                end
            end
        end
    end
    %name=strcat(directory,areaname);
    %save(name,'recording_data','-mat')
end