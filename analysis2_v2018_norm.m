
clear all
load('EC1','-mat')
%load('CA1','-mat')
%%sinergy analysis

count=ones(3,2);
count2=0;
total_cells=0;
inf_cells=ones(3,1);

percent=99;
total_cells=0;
for index=1:length(data)
    if isempty(data{index})~=1
        for j=1:length(data{index})
            total_cells=total_cells+1;
            count2=count2+1;
            
            Syn(count2,1)=data{index}{j}.I(5)-data{index}{j}.I(2)-data{index}{j}.I(3);
            Syn(count2,2)=data{index}{j}.I(6)-data{index}{j}.I(4)-data{index}{j}.I(3);
            
            Syn(count2,3)=data{index}{j}.I_delta(5)-data{index}{j}.I_delta(2)-data{index}{j}.I_delta(3);
            Syn(count2,4)=data{index}{j}.I_delta(6)-data{index}{j}.I_delta(4)-data{index}{j}.I_delta(3);
            
            Syn(count2,5)=data{index}{j}.I_theta(5)-data{index}{j}.I_theta(2)-data{index}{j}.I_theta(3);
            Syn(count2,6)=data{index}{j}.I_theta(6)-data{index}{j}.I_theta(4)-data{index}{j}.I_theta(3);
            
            Syn_all(total_cells,1)=data{index}{j}.I(1)-data{index}{j}.I(5)-data{index}{j}.I(6);
            Syn_all(total_cells,2)=data{index}{j}.I_delta(1)-data{index}{j}.I_delta(5)-data{index}{j}.I_delta(6);
            Syn_all(total_cells,3)=data{index}{j}.I_theta(1)-data{index}{j}.I_theta(5)-data{index}{j}.I_theta(6);

            
            test1=sort(data{index}{j}.TEST.I(5,:)-data{index}{j}.TEST.I(2,:)-data{index}{j}.TEST.I(3,:));
            test2=sort(data{index}{j}.TEST.I(6,:)-data{index}{j}.TEST.I(4,:)-data{index}{j}.TEST.I(3,:));
            test3=sort(data{index}{j}.TEST.I_delta(5,:)-data{index}{j}.TEST.I_delta(2,:)-data{index}{j}.TEST.I_delta(3,:));
            test4=sort(data{index}{j}.TEST.I_delta(6,:)-data{index}{j}.TEST.I_delta(4,:)-data{index}{j}.TEST.I_delta(3,:));
            test5=sort(data{index}{j}.TEST.I_theta(5,:)-data{index}{j}.TEST.I_theta(2,:)-data{index}{j}.TEST.I_theta(3,:));
            test6=sort(data{index}{j}.TEST.I_theta(6,:)-data{index}{j}.TEST.I_theta(4,:)-data{index}{j}.TEST.I_theta(3,:));
            test7=sort(data{index}{j}.TEST.I(1,:)-data{index}{j}.TEST.I(5,:)-data{index}{j}.TEST.I(6,:));
            test8=sort(data{index}{j}.TEST.I_delta(1,:)-data{index}{j}.TEST.I_delta(5,:)-data{index}{j}.TEST.I_delta(6,:));
            test9=sort(data{index}{j}.TEST.I_theta(1,:)-data{index}{j}.TEST.I_theta(5,:)-data{index}{j}.TEST.I_theta(6,:));
            
            %firing
            if Syn(count2,1)>test1(99)
                Syn_firing{1}(count(1,1))=Syn(count2,1)/data{index}{j}.I(5);
                count(1,1)=count(1,1)+1;
            end
            if Syn(count2,2)>test2(99)
                Syn_firing{2}(count(1,2))=Syn(count2,2)/data{index}{j}.I(6);
                count(1,2)=count(1,2)+1;
            end   
            %delta
            if Syn(count2,3)>test3(99)
                Syn_delta{1}(count(2,1))=Syn(count2,3)/data{index}{j}.I_delta(5);
                count(2,1)=count(2,1)+1;
            end
            if Syn(count2,4)>test4(99)
                Syn_delta{2}(count(2,2))=Syn(count2,4)/data{index}{j}.I_delta(6);
                count(2,2)=count(2,2)+1;
            end   
            %theta
            if Syn(count2,5)>test5(99)
                Syn_theta{1}(count(3,1))=Syn(count2,5)/data{index}{j}.I_theta(5);
                count(3,1)=count(3,1)+1;
            end
            if Syn(count2,6)>test2(99)
                Syn_theta{2}(count(3,2))=Syn(count2,6)/data{index}{j}.I_theta(6);
                count(3,2)=count(3,2)+1;
            end
            
            if Syn_all(total_cells,1)>test7(99)
                Syn_firing{3}(inf_cells(1))=Syn_all(total_cells,1)/data{index}{j}.I(1);
                inf_cells(1)=inf_cells(1)+1;
            end
            if Syn_all(total_cells,2)>test8(99)
                Syn_delta{3}(inf_cells(2))=Syn_all(total_cells,2)/data{index}{j}.I_delta(1);
                inf_cells(2)=inf_cells(2)+1;
            end
            if Syn_all(total_cells,3)>test9(99)
                Syn_theta{3}(inf_cells(3))=Syn_all(total_cells,3)/data{index}{j}.I_theta(1);
                inf_cells(3)=inf_cells(3)+1;
            end
            
            Syn(count2,1)=Syn(count2,1)/data{index}{j}.I(5);
            Syn(count2,2)=Syn(count2,2)/data{index}{j}.I(6);
            
            Syn(count2,3)=Syn(count2,3)/data{index}{j}.I_delta(5);
            Syn(count2,4)=Syn(count2,4)/data{index}{j}.I_delta(6);
            
            Syn(count2,5)=Syn(count2,5)/data{index}{j}.I_theta(5);
            Syn(count2,6)=Syn(count2,6)/data{index}{j}.I_theta(6);
            
            Syn_all(total_cells,1)=Syn_all(total_cells,1)/data{index}{j}.I(1);
            Syn_all(total_cells,2)=Syn_all(total_cells,2)/data{index}{j}.I_delta(1);
            Syn_all(total_cells,3)=Syn_all(total_cells,3)/data{index}{j}.I_theta(1);
        end
    end
end



figure(1)
%beh=['tim';'pos';'dir';'vel';'dip';'dis'];
%histograms of information firing
subplot(3,3,1)
a=histogram(Syn(:,1),'Binwidth',0.1);%,'Normalization','probability');
hold on
a=histogram(Syn_firing{1},'Binwidth',0.01);%,'Normalization','probability');    
title('Synergy Firing PosDir')
legend('Non Inf','Inf')
xlim([-1 1])


subplot(3,3,2)
a=histogram(Syn(:,2),'Binwidth',0.1);%,'Normalization','probability');
hold on
a=histogram(Syn_firing{2},'Binwidth',0.1);%,'Normalization','probability');    
title('Synergy Firing VelDir')
legend('Non Inf','Inf')
xlim([-1 1])


subplot(3,3,3)
a=histogram(Syn_all(:,1),'Binwidth',0.1);%,'Normalization','probability');
hold on
a=histogram(Syn_firing{3},'Binwidth',0.1);%,'Normalization','probability');    
title('Synergy Firing All')
legend('Non Inf','Inf')
xlim([-1 1])



subplot(3,3,4)
a=histogram(Syn(:,3),'Binwidth',0.1);%,'Normalization','probability');
hold on
a=histogram(Syn_delta{1},'Binwidth',0.1);%,'Normalization','probability');    
title('Synergy Delta PosDir')
legend('Non Inf','Inf')
xlim([-1 1])


subplot(3,3,5)
a=histogram(Syn(:,4),'Binwidth',0.1);%,'Normalization','probability');
hold on
a=histogram(Syn_firing{2},'Binwidth',0.1);%,'Normalization','probability');    
title('Synergy Delta VelDir')
legend('Non Inf','Inf')
xlim([-1 1])

           
subplot(3,3,6)
a=histogram(Syn_all(:,2),'Binwidth',0.1);%,'Normalization','probability');
hold on
a=histogram(Syn_delta{3},'Binwidth',0.1);%,'Normalization','probability');    
title('Synergy Firing All')
legend('Non Inf','Inf')
xlim([-1 1])


subplot(3,3,7)
a=histogram(Syn(:,5),'Binwidth',0.1);%,'Normalization','probability');
hold on
a=histogram(Syn_theta{1},'Binwidth',0.1);%,'Normalization','probability');    
title('Synergy Theta PosDir')
legend('Non Inf','Inf')
xlim([-1 1])


subplot(3,3,8)
a=histogram(Syn(:,6),'Binwidth',0.1);%,'Normalization','probability');
hold on
a=histogram(Syn_theta{2},'Binwidth',0.1);%,'Normalization','probability');    
title('Synergy Theta VelDir')
legend('Non Inf','Inf')
xlim([-1 1])


subplot(3,3,9)
a=histogram(Syn_all(:,3),'Binwidth',0.1);%,'Normalization','probability');
hold on
a=histogram(Syn_theta{3},'Binwidth',0.1);%,'Normalization','probability');    
title('Synergy Firing All')
legend('Non Inf','Inf')
xlim([-1 1])


%%fraction of synergistic cells

for i=1:3
    number(1,1)=length(Syn_firing{1})/total_cells;
    number(1,2)=length(Syn_firing{2})/total_cells;
    number(1,3)=length(Syn_firing{3})/total_cells;
    
    number(2,1)=length(Syn_delta{1})/total_cells;
    number(2,2)=length(Syn_delta{2})/total_cells;
    number(2,3)=length(Syn_delta{3})/total_cells;
    
    number(3,1)=length(Syn_theta{1})/total_cells;
    number(3,2)=length(Syn_theta{2})/total_cells;
    number(3,3)=length(Syn_theta{3})/total_cells; 
end
bar(number')

%% Comparison of synergies
clear all
load('EC1','-mat')
%load('CA1','-mat')
%%sinergy analysis


percent=99;
total_cells=0;
count2=0;
for index=1:length(data)
    if isempty(data{index})~=1
        for j=1:length(data{index})
            total_cells=total_cells+1;
            count2=count2+1;
            Syn(count2,1)=data{index}{j}.I(5)-data{index}{j}.I(2)-data{index}{j}.I(3);
            Syn(count2,2)=data{index}{j}.I(6)-data{index}{j}.I(4)-data{index}{j}.I(3);
            
            Syn(count2,3)=data{index}{j}.I_delta(5)-data{index}{j}.I_delta(2)-data{index}{j}.I_delta(3);
            Syn(count2,4)=data{index}{j}.I_delta(6)-data{index}{j}.I_delta(4)-data{index}{j}.I_delta(3);
            
            Syn(count2,5)=data{index}{j}.I_theta(5)-data{index}{j}.I_theta(2)-data{index}{j}.I_theta(3);
            Syn(count2,6)=data{index}{j}.I_theta(6)-data{index}{j}.I_theta(4)-data{index}{j}.I_theta(3);
            
            test1=sort(data{index}{j}.TEST.I(5,:)-data{index}{j}.TEST.I(2,:)-data{index}{j}.TEST.I(3,:));
            test2=sort(data{index}{j}.TEST.I(6,:)-data{index}{j}.TEST.I(4,:)-data{index}{j}.TEST.I(3,:));
            test3=sort(data{index}{j}.TEST.I_delta(5,:)-data{index}{j}.TEST.I_delta(2,:)-data{index}{j}.TEST.I_delta(3,:));
            test4=sort(data{index}{j}.TEST.I_delta(6,:)-data{index}{j}.TEST.I_delta(4,:)-data{index}{j}.TEST.I_delta(3,:));
            test5=sort(data{index}{j}.TEST.I_theta(5,:)-data{index}{j}.TEST.I_theta(2,:)-data{index}{j}.TEST.I_theta(3,:));
            test6=sort(data{index}{j}.TEST.I_theta(6,:)-data{index}{j}.TEST.I_theta(4,:)-data{index}{j}.TEST.I_theta(3,:));
            vector(count2,1:6)=zeros(1,6);
            %firing
            if Syn(count2,1)>test1(99)
                vector(count2,1)=1;
            end
            if Syn(count2,2)>test2(99)
                vector(count2,2)=1;
            end   
            %delta
            if Syn(count2,3)>test3(99)
                vector(count2,3)=1;
            end
            if Syn(count2,4)>test4(99)
                vector(count2,4)=1;
            end   
            %theta
            if Syn(count2,5)>test5(99)
                vector(count2,5)=1;
            end
            if Syn(count2,6)>test2(99)
                vector(count2,6)=1;
            end
            
            Syn(count2,1)=Syn(count2,1)/data{index}{j}.I(5);
            Syn(count2,2)=Syn(count2,2)/data{index}{j}.I(6);
            
            Syn(count2,3)=Syn(count2,3)/data{index}{j}.I_delta(5);
            Syn(count2,4)=Syn(count2,4)/data{index}{j}.I_delta(6);
            
            Syn(count2,5)=Syn(count2,5)/data{index}{j}.I_theta(5);
            Syn(count2,6)=Syn(count2,6)/data{index}{j}.I_theta(6);
        end
    end
end

Syn2=Syn.*vector;

figure(1)
%beh=['tim';'pos';'dir';'vel';'dip';'dis'];
%histograms of information firing
subplot(1,3,1)
scatter(Syn(:,1),Syn(:,2))
hold on
scatter(Syn2(:,1),Syn2(:,2),'*')
title('Syn VD vs PD firing')
xlabel('Syn(P,D)')
ylabel('Syn(V,D)')
xlim([-1 1])
ylim([-1 1])
box on

subplot(1,3,2)
scatter(Syn(:,3),Syn(:,4))
hold on
scatter(Syn2(:,3),Syn2(:,4),'*')
title('Syn VD vs PD delta')
xlim([-1 1])
ylim([-1 1])
xlabel('Syn(P,D)')
ylabel('Syn(V,D)')
box on

subplot(1,3,3)
scatter(Syn(:,5),Syn(:,6))   
hold on
scatter(Syn2(:,5),Syn2(:,6),'*')
title('Syn VD vs PD theta')
xlim([-1 1])
ylim([-1 1])
xlabel('Syn(P,D)')
ylabel('Syn(V,D)')
box on

%% another attempt


count=ones(3,2);
count2=0;
total_cells=0;
inf_cells=ones(3,1);

percent=99;
total_cells=0;
for index=1:length(data)
    if isempty(data{index})~=1
        for j=1:length(data{index})
            total_cells=total_cells+1;
            count2=count2+1;
            
            Syn(count2,1)=data{index}{j}.I(5)-data{index}{j}.I(2)-data{index}{j}.I(3);
            Syn(count2,2)=data{index}{j}.I(6)-data{index}{j}.I(4)-data{index}{j}.I(3);
            I(count2,1)=data{index}{j}.h_var(2)+data{index}{j}.h_var(3)-data{index}{j}.h_var(5);
            I(count2,2)=data{index}{j}.h_var(4)+data{index}{j}.h_var(3)-data{index}{j}.h_var(6);            
            
            Syn(count2,3)=data{index}{j}.I_delta(5)-data{index}{j}.I_delta(2)-data{index}{j}.I_delta(3);
            Syn(count2,4)=data{index}{j}.I_delta(6)-data{index}{j}.I_delta(4)-data{index}{j}.I_delta(3);
            I(count2,3)=data{index}{j}.h_var_deta(2)+data{index}{j}.h_var_delta(3)-data{index}{j}.h_var_delta(5);
            I(count2,4)=data{index}{j}.h_var_delta(4)+data{index}{j}.h_var_delta(3)-data{index}{j}.h_var_delta(6);            
            
            Syn(count2,5)=data{index}{j}.I_theta(5)-data{index}{j}.I_theta(2)-data{index}{j}.I_theta(3);
            Syn(count2,6)=data{index}{j}.I_theta(6)-data{index}{j}.I_theta(4)-data{index}{j}.I_theta(3);
            I(count2,5)=data{index}{j}.h_var_theta(2)+data{index}{j}.h_var_theta(3)-data{index}{j}.h_var_theta(5);
            I(count2,6)=data{index}{j}.h_var_theta(4)+data{index}{j}.h_var_theta(3)-data{index}{j}.h_var_theta(6);                    
            
            Syn_all(total_cells,1)=data{index}{j}.I(1)-data{index}{j}.I(5)-data{index}{j}.I(6);
            Syn_all(total_cells,2)=data{index}{j}.I_delta(1)-data{index}{j}.I_delta(5)-data{index}{j}.I_delta(6);
            Syn_all(total_cells,3)=data{index}{j}.I_theta(1)-data{index}{j}.I_theta(5)-data{index}{j}.I_theta(6);
            I_all(total_cells,1)=data{index}{j}.h_var(5)+data{index}{j}.h_var(6)-data{index}{j}.h_var(1);
            I_all(total_cells,2)=data{index}{j}.h_var_delta(5)+data{index}{j}.h_var_delta(6)-data{index}{j}.h_var_delta(1);
            I_all(total_cells,3)=data{index}{j}.h_var_theta(5)+data{index}{j}.h_var_theta(6)-data{index}{j}.h_var_theta(1);

            %%%%% vertificar a partir de aca!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
            
            test1=sort(data{index}{j}.TEST.I(5,:)-data{index}{j}.TEST.I(2,:)-data{index}{j}.TEST.I(3,:));
            test2=sort(data{index}{j}.TEST.I(6,:)-data{index}{j}.TEST.I(4,:)-data{index}{j}.TEST.I(3,:));
            test3=sort(data{index}{j}.TEST.I_delta(5,:)-data{index}{j}.TEST.I_delta(2,:)-data{index}{j}.TEST.I_delta(3,:));
            test4=sort(data{index}{j}.TEST.I_delta(6,:)-data{index}{j}.TEST.I_delta(4,:)-data{index}{j}.TEST.I_delta(3,:));
            test5=sort(data{index}{j}.TEST.I_theta(5,:)-data{index}{j}.TEST.I_theta(2,:)-data{index}{j}.TEST.I_theta(3,:));
            test6=sort(data{index}{j}.TEST.I_theta(6,:)-data{index}{j}.TEST.I_theta(4,:)-data{index}{j}.TEST.I_theta(3,:));
            test7=sort(data{index}{j}.TEST.I(1,:)-data{index}{j}.TEST.I(5,:)-data{index}{j}.TEST.I(6,:));
            test8=sort(data{index}{j}.TEST.I_delta(1,:)-data{index}{j}.TEST.I_delta(5,:)-data{index}{j}.TEST.I_delta(6,:));
            test9=sort(data{index}{j}.TEST.I_theta(1,:)-data{index}{j}.TEST.I_theta(5,:)-data{index}{j}.TEST.I_theta(6,:));
            
            %firing
            if Syn(count2,1)>test1(99)
                Syn_firing{1}(count(1,1))=Syn(count2,1);
                I_firing{1}(count(1,1))=I(count2,1);                
                count(1,1)=count(1,1)+1;
            end
            if Syn(count2,2)>test2(99)
                Syn_firing{2}(count(1,2))=Syn(count2,2)/data{index}{j}.I(6);
                count(1,2)=count(1,2)+1;
            end   
            %delta
            if Syn(count2,3)>test3(99)
                Syn_delta{1}(count(2,1))=Syn(count2,3)/data{index}{j}.I_delta(5);
                count(2,1)=count(2,1)+1;
            end
            if Syn(count2,4)>test4(99)
                Syn_delta{2}(count(2,2))=Syn(count2,4)/data{index}{j}.I_delta(6);
                count(2,2)=count(2,2)+1;
            end   
            %theta
            if Syn(count2,5)>test5(99)
                Syn_theta{1}(count(3,1))=Syn(count2,5)/data{index}{j}.I_theta(5);
                count(3,1)=count(3,1)+1;
            end
            if Syn(count2,6)>test2(99)
                Syn_theta{2}(count(3,2))=Syn(count2,6)/data{index}{j}.I_theta(6);
                count(3,2)=count(3,2)+1;
            end
            
            if Syn_all(total_cells,1)>test7(99)
                Syn_firing{3}(inf_cells(1))=Syn_all(total_cells,1)/data{index}{j}.I(1);
                inf_cells(1)=inf_cells(1)+1;
            end
            if Syn_all(total_cells,2)>test8(99)
                Syn_delta{3}(inf_cells(2))=Syn_all(total_cells,2)/data{index}{j}.I_delta(1);
                inf_cells(2)=inf_cells(2)+1;
            end
            if Syn_all(total_cells,3)>test9(99)
                Syn_theta{3}(inf_cells(3))=Syn_all(total_cells,3)/data{index}{j}.I_theta(1);
                inf_cells(3)=inf_cells(3)+1;
            end
            
            Syn(count2,1)=Syn(count2,1)/data{index}{j}.I(5);
            Syn(count2,2)=Syn(count2,2)/data{index}{j}.I(6);
            
            Syn(count2,3)=Syn(count2,3)/data{index}{j}.I_delta(5);
            Syn(count2,4)=Syn(count2,4)/data{index}{j}.I_delta(6);
            
            Syn(count2,5)=Syn(count2,5)/data{index}{j}.I_theta(5);
            Syn(count2,6)=Syn(count2,6)/data{index}{j}.I_theta(6);
            
            Syn_all(total_cells,1)=Syn_all(total_cells,1)/data{index}{j}.I(1);
            Syn_all(total_cells,2)=Syn_all(total_cells,2)/data{index}{j}.I_delta(1);
            Syn_all(total_cells,3)=Syn_all(total_cells,3)/data{index}{j}.I_theta(1);
        end
    end
end


%% Comparison of synergies



