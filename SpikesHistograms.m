function SpikesHistograms(directory,files,rat,area,freq)

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
                    name=strcat(directory,area(index,:),protocol,int2str(se));
                    matrix=load(name,'-mat');
                    matrix=[matrix.matrixr,matrix.matrixl];
                    h=sum(matrix);
                    %h=h/sum(h);
                    figure(index)
                    plot(h)
                    hold on
                    plot(h,'lineStyle','none','Marker','o','MarkerFaceColor','b','MarkerSize',5,'MarkerEdgeColor','b');
                    xlabel ('Time [s]');
                    ylabel ('Spike Count')
                    title(area(index,:))
                    saveas(index,name,'png');
                    close(figure(index));
                end
            end
        end
    end
end