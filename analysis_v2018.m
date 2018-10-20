%Program for analysis of information of CA1 and EC1 neurons
clear all
load('EC1','-mat')
areaname='EC1';

load('CA1','-mat')
areaname='CA1';
%% COMPARE DIFFERENT CODES

count=ones(6,3);
percent=99;
total_cells=0;
for index=1:length(data)
    if isempty(data{index})~=1
        for j=1:length(data{index})
            total_cells=total_cells+1;
            for aspect=1:6
                test=sort(data{index}{j}.TEST.I(aspect,:));
                if data{index}{j}.I(aspect)>test(percent)
                    I{aspect}(count(aspect,1))=data{index}{j}.I(aspect);
                    count(aspect,1)=count(aspect,1)+1;
                end
                test=sort(data{index}{j}.TEST.I_delta(aspect,:));
                if data{index}{j}.I_delta(aspect)>test(99)
                    I_delta{aspect}(count(aspect,2))=data{index}{j}.I_delta(aspect);
                    count(aspect,2)=count(aspect,2)+1;
                end
                test=sort(data{index}{j}.TEST.I_theta(aspect,:));
                if data{index}{j}.I_theta(aspect)>test(99)
                    I_theta{aspect}(count(aspect,3))=data{index}{j}.I_theta(aspect);
                    count(aspect,3)=count(aspect,3)+1;
                end
            end
        end
    end
end

directory='/home/melisa/Escritorio/Datos_Ines/';
beh=['tim';'pos';'dir';'vel';'dip';'dis'];
for aspect=1:6
    name=strcat(directory,areaname,'_information_firing_',beh(aspect,:));
    I_save=I{aspect};
    save(name,'I_save','-ASCII');
    name=strcat(directory,areaname,'_information_delta_',beh(aspect,:));
    I_save=I_delta{aspect};
    save(name,'I_save','-ASCII');
    name=strcat(directory,areaname,'_information_theta_',beh(aspect,:));
    I_save=I_theta{aspect};
    save(name,'I_save','-ASCII');
end

%fraction of significative informative cells in each code
for aspect=1:6
    number(1,aspect)=length(I{aspect})/total_cells;
    number(2,aspect)=length(I_delta{aspect})/total_cells;
    number(3,aspect)=length(I_theta{aspect})/total_cells;
end
bar(number')
name=strcat(directory,areaname,'_information_fraction');
save(name,'number','-ASCII')

%using normalize information
count=ones(6,3);
percent=99;
total_cells=0;
for index=1:length(data)
    if isempty(data{index})~=1
        for j=1:length(data{index})
            total_cells=total_cells+1;
            for aspect=1:6
                test=sort(data{index}{j}.TEST.I(aspect,:));
                if data{index}{j}.I(aspect)>test(percent)
                    I{aspect}(count(aspect,1))=data{index}{j}.I_norm(aspect);
                    count(aspect,1)=count(aspect,1)+1;
                end
                test=sort(data{index}{j}.TEST.I_delta(aspect,:));
                if data{index}{j}.I_delta(aspect)>test(99)
                    I_delta{aspect}(count(aspect,2))=data{index}{j}.I_delta_norm(aspect);
                    count(aspect,2)=count(aspect,2)+1;
                end
                test=sort(data{index}{j}.TEST.I_theta(aspect,:));
                if data{index}{j}.I_theta(aspect)>test(99)
                    I_theta{aspect}(count(aspect,3))=data{index}{j}.I_theta_norm(aspect);
                    count(aspect,3)=count(aspect,3)+1;
                end
            end
        end
    end
end


beh=['tim';'pos';'dir';'vel';'dip';'dis'];
%histograms of information
for aspect=1:6
    subplot(3,2,aspect)
    a=histogram(I{aspect},'Binwidth',0.005,'Normalization','probability');
    hold on
    c=histogram(I_theta{aspect},'Binwidth',0.005,'Normalization','probability');
    b=histogram(I_delta{aspect},'Binwidth',0.005,'Normalization','probability');
    title(beh(aspect,:))
    xlabel(strcat('I(rta;',beh(aspect,:),')/H'))
    xlim([0 0.3])
    ylim([0 1])
    legend('Firing rate','Theta','Delta')
end

for aspect=1:6
    info(1,aspect)=mean(I{aspect});
    info(2,aspect)=mean(I_delta{aspect});
    info(3,aspect)=mean(I_theta{aspect});
end




