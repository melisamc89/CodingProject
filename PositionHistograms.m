
function PositionHistograms(directory1,directory2,files,rat,area)

    for index= 1:size(area)
        register=Rat_Register(files,area,rat,index);        
        [fil col]=size(register.sessions);
        for se=1:col
            protocol=register.sessions{se}.sess_type;
            if protocol == 'C'
                spikes_times= register.sessions{se}.spike_times;
                position=Loading_Pos(register,se);
                tf = isfield(position.data, 'log');
                if (tf==1 && length(spikes_times)>50)
                    disc=80;
                    registername=strcat(area(index,:),'C',int2str(se));
                    [xr sigmaxr xl sigmaxl]=FiringPositionData(directory2,registername,'square');
                    name=strcat(directory1,area(index,:),'C',int2str(se));
                    save(name,'xr','sigmaxr','xl','sigmaxl','-ASCII');
                    meanr([1:length(xr)])=mean(xr);
                    stdr([1:length(xr)])=std(xr);
                    meanl([1:length(xl)])=mean(xl);
                    stdl([1:length(xl)])=std(xl);
                    
                    pos=0:18:length(xr)*18-18;
                    figure(index)
                    subplot(2,1,1)
                    errorbar(pos,xr,sigmaxr)
                    xlim([0 400])
                    hold on
                    errorbar(pos,meanr,stdr,'r')
                    title(strcat(area(index,:),' Right'))
                    xlabel('X [cm]')
                    ylabel('<n>')
                    subplot(2,1,2)
                    errorbar(pos,xl,sigmaxl)
                    hold on
                    errorbar(pos,meanl,stdl,'r')
                    title(strcat(area(index,:),' Left'))
                    xlim([0 400])
                    xlabel('X [cm]')
                    ylabel('<n>')
                    saveas(index,name,'png');
                    close(index)
                end
            end
        end
    end
end