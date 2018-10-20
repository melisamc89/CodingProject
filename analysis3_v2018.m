%Program for analysis of information of CA1 and EC1 neurons
clear all
load('EC1','-mat')
load('CA1','-mat')
%% COMPARE DIFFERENT CODES
%using normalize information
count=ones(6,3);
percent=99;
total_cells=0;
for index=1:length(data)
    if isempty(data{index})~=1
        for j=1:length(data{index})
            total_cells=total_cells+1;
            fr_all(total_cells)=data{index}{j}.estimators.mean_firing_rate;
            for aspect=1:6
                test=sort(data{index}{j}.TEST.I(aspect,:));
                if data{index}{j}.I(aspect)>test(percent)
                    I{aspect}(count(aspect,1))=data{index}{j}.I_norm(aspect);
                    fr{aspect}(count(aspect,1))=data{index}{j}.estimators.mean_firing_rate;
                    count(aspect,1)=count(aspect,1)+1;
                end
                test=sort(data{index}{j}.TEST.I_delta(aspect,:));
                if data{index}{j}.I_delta(aspect)>test(99)
                    I_delta{aspect}(count(aspect,2))=data{index}{j}.I_delta_norm(aspect);
                    fr_delta{aspect}(count(aspect,2))=data{index}{j}.estimators.mean_firing_rate;
                    var_delta{aspect}(count(aspect,2))=data{index}{j}.estimators.var_phase_delta;
                    count(aspect,2)=count(aspect,2)+1;
                end
                test=sort(data{index}{j}.TEST.I_theta(aspect,:));
                if data{index}{j}.I_theta(aspect)>test(99)
                    I_theta{aspect}(count(aspect,3))=data{index}{j}.I_theta_norm(aspect);
                    fr_theta{aspect}(count(aspect,3))=data{index}{j}.estimators.mean_firing_rate;
                    var_theta{aspect}(count(aspect,3))=data{index}{j}.estimators.var_phase_theta;
                    count(aspect,3)=count(aspect,3)+1;
                end
            end
        end
    end
end

beh=['tim';'pos';'dir';'vel';'dip';'dis'];
%firing rate vs information
for aspect=1:6
    subplot(3,2,aspect)
    scatter(I{aspect},fr{aspect},'filled')
    hold on
    scatter(I_theta{aspect},fr_theta{aspect},'filled')
    scatter(I_delta{aspect},fr_delta{aspect},'filled')
    box on
    title(beh(aspect,:))
    ylabel('Firing Rate [Hz]')
    xlabel('Inf')
    legend('Firing','Theta','Delta')
end
    
%circ variance vs information
for aspect=1:6
    subplot(3,2,aspect)
    scatter(I_theta{aspect},var_theta{aspect},'filled')
    hold on
    scatter(I_delta{aspect},var_delta{aspect},'filled')
    box on
    title(beh(aspect,:))
    ylabel('Cvar')
    xlabel('Inf')
    legend('Theta','Delta')
end

%circ variance vs firing_rate
for aspect=1:6
    subplot(3,2,aspect)
    scatter(fr_theta{aspect},var_theta{aspect},'filled')
    hold on
    scatter(fr_delta{aspect},var_delta{aspect},'filled')
    box on
    title(beh(aspect,:))
    ylabel('Cvar')
    xlabel('Frate')
    legend('Theta','Delta')
end
