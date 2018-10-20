
clear all
load('EC1','-mat')
load('CA1','-mat')
%%sinergy analysis


%%compare informative and non-informative
count=ones(6,3);
count2=ones(6);

percent=99;
total_cells=0;
for index=1:length(data)
    if isempty(data{index})~=1
        for j=1:length(data{index})
            total_cells=total_cells+1;
            for aspect=1:6
                I_all{aspect}(count2(aspect))=data{index}{j}.I(aspect);
                I_delta_all{aspect}(count2(aspect))=data{index}{j}.I_delta(aspect);
                I_theta_all{aspect}(count2(aspect))=data{index}{j}.I_theta(aspect);
                count2(aspect)=count2(aspect)+1;
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

count=ones(6,3);
count2=ones(6);

percent=99;
total_cells=0;
for index=1:length(data)
    if isempty(data{index})~=1
        for j=1:length(data{index})
            total_cells=total_cells+1;
            for aspect=1:6
                I_all{aspect}(count2(aspect))=data{index}{j}.I(aspect);
                I_delta_all{aspect}(count2(aspect))=data{index}{j}.I_delta(aspect);
                I_theta_all{aspect}(count2(aspect))=data{index}{j}.I_theta(aspect);
                count2(aspect)=count2(aspect)+1;
                test=sort(data{index}{j}.TEST.I(aspect,:));
                h=min(data{index}{j}.h_var(aspect),data{index}{j}.h(aspect));
                if data{index}{j}.I(aspect)>test(percent) && h
                    I{aspect}(count(aspect,1))=data{index}{j}.I(aspect)/h;
                    count(aspect,1)=count(aspect,1)+1;
                end
                test=sort(data{index}{j}.TEST.I_delta(aspect,:));
                h=min(data{index}{j}.h_var_delta(aspect),data{index}{j}.h_delta(aspect));
                if data{index}{j}.I_delta(aspect)>test(99) && h
                    h=min(data{index}{j}.h_var_delta(aspect),data{index}{j}.h_delta(aspect));
                    I_delta{aspect}(count(aspect,2))=data{index}{j}.I_delta(aspect)/h;
                    count(aspect,2)=count(aspect,2)+1;
                end
                test=sort(data{index}{j}.TEST.I_theta(aspect,:));
                h=min(data{index}{j}.h_var_theta(aspect),data{index}{j}.h_theta(aspect));
                if data{index}{j}.I_theta(aspect)>test(99) && h
                    I_theta{aspect}(count(aspect,3))=data{index}{j}.I_theta(aspect)/h;
                    h_prueba{aspect}(count(aspect,3))=h;
                    count(aspect,3)=count(aspect,3)+1;
                end
            end
        end
    end
end


figure(1)
beh=['tim';'pos';'dir';'vel';'dip';'dis'];
%histograms of information firing
for aspect=1:6
    subplot(3,2,aspect)
    a=histogram(I_all{aspect},'Binwidth',0.005,'Normalization','cdf');
    hold on
    a=histogram(I{aspect},'Binwidth',0.005,'Normalization','cdf');    
    title(beh(aspect,:))
    xlabel(strcat('I(rta;',beh(aspect,:),')'))
    xlim([0 1])
    ylim([0 1])
    legend('Non Inf','Inf')
end

figure(2)
beh=['tim';'pos';'dir';'vel';'dip';'dis'];
%histograms of information firing
for aspect=1:6
    subplot(3,2,aspect)
    a=histogram(I_delta_all{aspect},'Binwidth',0.001,'Normalization','cdf');
    hold on
    a=histogram(I_delta{aspect},'Binwidth',0.001,'Normalization','cdf');    
    title(beh(aspect,:))
    xlabel(strcat('I(rta;',beh(aspect,:),')'))
    xlim([0 0.15])
    ylim([0 1])
    legend('Non Inf','Inf')
end

figure(3)
beh=['tim';'pos';'dir';'vel';'dip';'dis'];
%histograms of information firing
for aspect=1:6
    subplot(3,2,aspect)
    a=histogram(I_theta_all{aspect},'Binwidth',0.001,'Normalization','cdf');
    hold on
    a=histogram(I_theta{aspect},'Binwidth',0.001,'Normalization','cdf');    
    title(beh(aspect,:))
    xlabel(strcat('I(rta;',beh(aspect,:),')'))
    xlim([0 0.15])
    ylim([0 1])
    legend('Non Inf','Inf')
end