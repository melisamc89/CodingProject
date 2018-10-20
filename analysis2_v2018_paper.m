
clear all
load('EC1','-mat')
load('CA1','-mat')
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
            I1(count2,1)=-data{index}{j}.ht(5)+data{index}{j}.ht(2)+data{index}{j}.ht(3);
            I1(count2,2)=-data{index}{j}.ht(6)+data{index}{j}.ht(4)+data{index}{j}.ht(3);
            
            Syn(count2,3)=data{index}{j}.I_delta(5)-data{index}{j}.I_delta(2)-data{index}{j}.I_delta(3);
            Syn(count2,4)=data{index}{j}.I_delta(6)-data{index}{j}.I_delta(4)-data{index}{j}.I_delta(3);
            I1(count2,3)=-data{index}{j}.ht_delta(5)+data{index}{j}.ht_delta(2)+data{index}{j}.ht_delta(3);
            I1(count2,4)=-data{index}{j}.ht_delta(6)+data{index}{j}.ht_delta(4)+data{index}{j}.ht_delta(3);
            
            Syn(count2,5)=data{index}{j}.I_theta(5)-data{index}{j}.I_theta(2)-data{index}{j}.I_theta(3);
            Syn(count2,6)=data{index}{j}.I_theta(6)-data{index}{j}.I_theta(4)-data{index}{j}.I_theta(3);
            I1(count2,5)=-data{index}{j}.ht_theta(5)+data{index}{j}.ht_theta(2)+data{index}{j}.ht_theta(3);
            I1(count2,6)=-data{index}{j}.ht_theta(6)+data{index}{j}.ht_theta(4)+data{index}{j}.ht_theta(3);
            
            test1=sort(data{index}{j}.TEST.I(5,:)-data{index}{j}.TEST.I(2,:)-data{index}{j}.TEST.I(3,:));
            test2=sort(data{index}{j}.TEST.I(6,:)-data{index}{j}.TEST.I(4,:)-data{index}{j}.TEST.I(3,:));
            test3=sort(data{index}{j}.TEST.I_delta(5,:)-data{index}{j}.TEST.I_delta(2,:)-data{index}{j}.TEST.I_delta(3,:));
            test4=sort(data{index}{j}.TEST.I_delta(6,:)-data{index}{j}.TEST.I_delta(4,:)-data{index}{j}.TEST.I_delta(3,:));
            test5=sort(data{index}{j}.TEST.I_theta(5,:)-data{index}{j}.TEST.I_theta(2,:)-data{index}{j}.TEST.I_theta(3,:));
            test6=sort(data{index}{j}.TEST.I_theta(6,:)-data{index}{j}.TEST.I_theta(4,:)-data{index}{j}.TEST.I_theta(3,:));
            
            %firing
            if Syn(count2,1)>test1(end)
                Syn_firing{1}(count(1,1))=Syn(count2,1);
                I_firing{1}(count(1,1))=I1(count2,1);
                count(1,1)=count(1,1)+1;
            end
            if Syn(count2,2)>test2(end)
                Syn_firing{2}(count(1,2))=Syn(count2,2);
                I_firing{2}(count(1,2))=I1(count2,2);
                count(1,2)=count(1,2)+1;
            end   
            %delta
            if Syn(count2,3)>test3(end)
                Syn_delta{1}(count(2,1))=Syn(count2,3);
                I_delta{1}(count(2,1))=I1(count2,3);
                count(2,1)=count(2,1)+1;
            end
            if Syn(count2,4)>test4(end)
                Syn_delta{2}(count(2,2))=Syn(count2,4);
                  I_delta{2}(count(2,2))=I1(count2,4);
                count(2,2)=count(2,2)+1;
            end   
            %theta
            if Syn(count2,5)>test5(end)
                Syn_theta{1}(count(3,1))=Syn(count2,5);
                I_theta{1}(count(3,1))=I1(count2,5);
                count(3,1)=count(3,1)+1;
            end
            if Syn(count2,6)>test2(end)
                Syn_theta{2}(count(3,2))=Syn(count2,6);
                I_theta{2}(count(3,2))=I1(count2,6);
                count(3,2)=count(3,2)+1;
            end       
        end
    end
end



figure(1)
x=[0:0.1:6];
subplot(2,3,1)
scatter(I1(:,1),Syn(:,1)+I1(:,1))
hold on
plot(x,x,'Linewidth',2)
xlim([0 6])
ylim([0 6])
box on
xlabel('<I(P;D|n1)>')
ylabel('I(P;D)')
title('Firing')

subplot(2,3,4)
scatter(I1(:,2),Syn(:,2)+I1(:,2))
hold on
plot(x,x,'Linewidth',2)
xlim([0 6])
ylim([0 6])
box on
xlabel('<I(V;D|n1)>')
ylabel('I(V;D)')

subplot(2,3,2)
scatter(I1(:,3),Syn(:,3)+I1(:,3))
hold on
plot(x,x,'Linewidth',2)
xlim([0 6])
ylim([0 6])
box on
xlabel('<I(P;D|n1)>')
ylabel('I(P;D)')
title('Delta')

subplot(2,3,5)
scatter(I1(:,4),Syn(:,4)+I1(:,4))
hold on
plot(x,x,'Linewidth',2)
box on
xlabel('<I(V;D|n1)>')
ylabel('I(V;D)')

subplot(2,3,3)
scatter(I1(:,5),Syn(:,5)+I1(:,5))
hold on
plot(x,x,'Linewidth',2)
xlim([0 6])
ylim([0 6])
box on
xlabel('<I(P;D|n1)>')
ylabel('I(P;D)')
title('Theta')

subplot(2,3,6)
scatter(I1(:,6),Syn(:,6)+I1(:,6))
hold on
plot(x,x,'Linewidth',2)
xlim([0 6])
ylim([0 6])
box on
xlabel('<I(V;D|n1)>')
ylabel('I(V;D)')

%%

figure(2)
x=[0:0.1:6];
subplot(2,3,1)
scatter(I1(:,1),Syn(:,1)+I1(:,1))
hold on
scatter(I_firing{1},Syn_firing{1}+I_firing{1},'*')
plot(x,x,'Linewidth',2)
xlim([0 6])
ylim([0 6])
box on
xlabel('<I(P;D|n1)>')
ylabel('I(P;D)')
title('Firing')

subplot(2,3,4)
scatter(I1(:,2),Syn(:,2)+I1(:,2))
hold on
scatter(I_firing{2},Syn_firing{2}+I_firing{2},'*')
plot(x,x,'Linewidth',2)
xlim([0 6])
ylim([0 6])
box on
xlabel('<I(V;D|n1)>')
ylabel('I(V;D)')

subplot(2,3,2)
scatter(I1(:,3),Syn(:,3)+I1(:,3))
hold on
scatter(I_delta{1},Syn_delta{1}+I_delta{1},'*')
plot(x,x,'Linewidth',2)
xlim([0 6])
ylim([0 6])
box on
xlabel('<I(P;D|n1)>')
ylabel('I(P;D)')
title('Delta')

subplot(2,3,5)
scatter(I1(:,4),Syn(:,4)+I1(:,4))
hold on
scatter(I_delta{2},Syn_delta{2}+I_delta{2},'*')
plot(x,x,'Linewidth',2)
box on
xlabel('<I(V;D|n1)>')
ylabel('I(V;D)')

subplot(2,3,3)
scatter(I1(:,5),Syn(:,5)+I1(:,5))
hold on
scatter(I_theta{1},Syn_theta{1}+I_theta{1},'*')
plot(x,x,'Linewidth',2)
xlim([0 6])
ylim([0 6])
box on
xlabel('<I(P;D|n1)>')
ylabel('I(P;D)')
title('Theta')

subplot(2,3,6)
scatter(I1(:,6),Syn(:,6)+I1(:,6))
hold on
scatter(I_theta{2},Syn_theta{2}+I_theta{2},'*')
plot(x,x,'Linewidth',2)
xlim([0 6])
ylim([0 6])
box on
xlabel('<I(V;D|n1)>')
ylabel('I(V;D)')

