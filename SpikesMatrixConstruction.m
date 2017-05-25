
function SpikesMatrixConstruction(directory,files,rat,area,freq)

    for index= 1:size(area)
        register=Rat_Register(files,area,rat,index);        
        [fil col]=size(register.sessions);
        for se=1:col
            protocol=register.sessions{se}.sess_type;
            if protocol == 'C'
                spikes_times= register.sessions{se}.spike_times;
                pos=Loading_Pos(register,se);
                tf = isfield(pos.data, 'log');
                if (tf==1 && length(spikes_times)>50)
                    [matrixr matrixl]=SpikesMatrix(spikes_times,pos,40,freq);
                    name=strcat(directory,area(index,:),'C',int2str(se));
                    save(name,'matrixr','matrixl','-mat');
                end
            end
        end
    end
end
