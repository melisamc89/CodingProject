
clear all
%load('EC1','-mat')
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
            
            Syn(count2,3)=data{index}{j}.I_delta(5)-data{index}{j}.I_delta(2)-data{index}{j}.I_delta(3);
            Syn(count2,4)=data{index}{j}.I_delta(6)-data{index}{j}.I_delta(4)-data{index}{j}.I_delta(3);
            
            Syn(count2,5)=data{index}{j}.I_theta(5)-data{index}{j}.I_theta(2)-data{index}{j}.I_theta(3);
            Syn(count2,6)=data{index}{j}.I_theta(6)-data{index}{j}.I_theta(4)-data{index}{j}.I_theta(3);
            
            Syn_all(total_cells,1)=data{index}{j}.I(1)-data{index}{j}.I(5)-data{index}{j}.I(6);
            Syn_all(total_cells,2)=data{index}{j}.I_delta(1)-data{index}{j}.I_delta(5)-data{index}{j}.I_delta(6);
            Syn_all(total_cells,3)=data{index}{j}.I_theta(1)-data{index}{j}.I_theta(5)-data{index}{j}.I_theta(6);

            
            test1=sort(abs(data{index}{j}.TEST.I(5,:)-data{index}{j}.TEST.I(2,:)-data{index}{j}.TEST.I(3,:)));
            test2=sort(abs(data{index}{j}.TEST.I(6,:)-data{index}{j}.TEST.I(4,:)-data{index}{j}.TEST.I(3,:)));
            test3=sort(abs(data{index}{j}.TEST.I_delta(5,:)-data{index}{j}.TEST.I_delta(2,:)-data{index}{j}.TEST.I_delta(3,:)));
            test4=sort(abs(data{index}{j}.TEST.I_delta(6,:)-data{index}{j}.TEST.I_delta(4,:)-data{index}{j}.TEST.I_delta(3,:)));
            test5=sort(abs(data{index}{j}.TEST.I_theta(5,:)-data{index}{j}.TEST.I_theta(2,:)-data{index}{j}.TEST.I_theta(3,:)));
            test6=sort(abs(data{index}{j}.TEST.I_theta(6,:)-data{index}{j}.TEST.I_theta(4,:)-data{index}{j}.TEST.I_theta(3,:)));
            test7=sort(data{index}{j}.TEST.I(1,:)-data{index}{j}.TEST.I(5,:)-data{index}{j}.TEST.I(6,:));
            test8=sort(data{index}{j}.TEST.I_delta(1,:)-data{index}{j}.TEST.I_delta(5,:)-data{index}{j}.TEST.I_delta(6,:));
            test9=sort(data{index}{j}.TEST.I_theta(1,:)-data{index}{j}.TEST.I_theta(5,:)-data{index}{j}.TEST.I_theta(6,:));
            
            %firing
            if abs(Syn(count2,1))>test1(99)
                Syn_firing{1}(count(1,1))=Syn(count2,1);
                count(1,1)=count(1,1)+1;
            end
            if abs(Syn(count2,2))>test2(99)
                Syn_firing{2}(count(1,2))=Syn(count2,2);
                count(1,2)=count(1,2)+1;
            end   
            %delta
            if abs(Syn(count2,3))>test3(99)
                Syn_delta{1}(count(2,1))=Syn(count2,3);
                count(2,1)=count(2,1)+1;
            end
            if abs(Syn(count2,4))>test4(99)
                Syn_delta{2}(count(2,2))=Syn(count2,4);
                count(2,2)=count(2,2)+1;
            end   
            %theta
            if abs(Syn(count2,5))>test5(99)
                Syn_theta{1}(count(3,1))=Syn(count2,5);
                count(3,1)=count(3,1)+1;
            end
            if abs(Syn(count2,6))>test2(99)
                Syn_theta{2}(count(3,2))=Syn(count2,6);
                count(3,2)=count(3,2)+1;
            end
            
            if Syn_all(total_cells,1)>test7(99)
                Syn_firing{3}(inf_cells(1))=Syn_all(total_cells,1);
                inf_cells(1)=inf_cells(1)+1;
            end
            if Syn_all(total_cells,2)>test8(99)
                Syn_delta{3}(inf_cells(2))=Syn_all(total_cells,2);
                inf_cells(2)=inf_cells(2)+1;
            end
            if Syn_all(total_cells,3)>test9(99)
                Syn_theta{3}(inf_cells(3))=Syn_all(total_cells,3);
                inf_cells(3)=inf_cells(3)+1;
            end
        end
    end
end


for i=1:3
    x1(i)=length(find(Syn_firing{i}<0))/220;
    x2(i)=length(find(Syn_delta{i}<0))/220;
    x3(i)=length(find(Syn_theta{i}<0))/220;
end

data=Syn_firing{1}';
save('EC_syn_firing_posdir','data','-ASCII')
data=Syn_firing{2}';
save('EC_syn_firing_veldir','data','-ASCII')
data=Syn_delta{1}';
save('EC_syn_delta_posdir','data','-ASCII')
data=Syn_delta{2}';
save('EC_syn_delta_veldir','data','-ASCII')
data=Syn_theta{1}';
save('EC_syn_theta_posdir','data','-ASCII')
data=Syn_theta{2}';
save('EC_syn_theta_veldir','data','-ASCII')

figure(1)
%beh=['tim';'pos';'dir';'vel';'dip';'dis'];
%histograms of information firing
subplot(3,3,1)
a=histogram(Syn(:,1),'Binwidth',0.005);%,'Normalization','probability');
hold on
a=histogram(Syn_firing{1},'Binwidth',0.005);%,'Normalization','probability');    
title('Synergy Firing PosDir')
legend('Non Inf','Inf')
subplot(3,3,2)
a=histogram(Syn(:,2),'Binwidth',0.005);%,'Normalization','probability');
hold on
a=histogram(Syn_firing{2},'Binwidth',0.005);%,'Normalization','probability');    
title('Synergy Firing VelDir')
legend('Non Inf','Inf')
subplot(3,3,3)
a=histogram(Syn_all(:,1),'Binwidth',0.005);%,'Normalization','probability');
hold on
a=histogram(Syn_firing{3},'Binwidth',0.005);%,'Normalization','probability');    
title('Synergy Firing All')
legend('Non Inf','Inf')


subplot(3,3,4)
a=histogram(Syn(:,3),'Binwidth',0.005);%,'Normalization','probability');
hold on
a=histogram(Syn_delta{1},'Binwidth',0.005);%,'Normalization','probability');    
title('Synergy Delta PosDir')
legend('Non Inf','Inf')

subplot(3,3,5)
a=histogram(Syn(:,4),'Binwidth',0.005);%,'Normalization','probability');
hold on
a=histogram(Syn_firing{2},'Binwidth',0.005);%,'Normalization','probability');    
title('Synergy Delta VelDir')
legend('Non Inf','Inf')
           
subplot(3,3,6)
a=histogram(Syn_all(:,2),'Binwidth',0.005);%,'Normalization','probability');
hold on
a=histogram(Syn_delta{3},'Binwidth',0.005);%,'Normalization','probability');    
title('Synergy Firing All')
legend('Non Inf','Inf')

subplot(3,3,7)
a=histogram(Syn(:,5),'Binwidth',0.005);%,'Normalization','probability');
hold on
a=histogram(Syn_theta{1},'Binwidth',0.005);%,'Normalization','probability');    
title('Synergy Theta PosDir')
legend('Non Inf','Inf')

subplot(3,3,8)
a=histogram(Syn(:,6),'Binwidth',0.005);%,'Normalization','probability');
hold on
a=histogram(Syn_theta{2},'Binwidth',0.005);%,'Normalization','probability');    
title('Synergy Theta VelDir')
legend('Non Inf','Inf')

subplot(3,3,9)
a=histogram(Syn_all(:,3),'Binwidth',0.005);%,'Normalization','probability');
hold on
a=histogram(Syn_theta{3},'Binwidth',0.005);%,'Normalization','probability');    
title('Synergy Firing All')
legend('Non Inf','Inf')

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

areaname='CA1';
directory='/home/melisa/Escritorio/Datos_Ines/';
name=strcat(directory,areaname,'_syn_firing_posdir')
save_data=Syn_firing{1};
save(name,'save_data','-ASCII')
name=strcat(directory,areaname,'_syn_firing_veldir')
save_data=Syn_firing{2};
save(name,'save_data','-ASCII')
name=strcat(directory,areaname,'_syn_delta_posdir')
save_data=Syn_delta{1};
save(name,'save_data','-ASCII')
name=strcat(directory,areaname,'_syn_delta_veldir')
save_data=Syn_delta{2};
save(name,'save_data','-ASCII')
name=strcat(directory,areaname,'_syn_theta_posdir')
save_data=Syn_theta{1};
save(name,'save_data','-ASCII')
name=strcat(directory,areaname,'_syn_theta_veldir')
save_data=Syn_theta{2};
save(name,'save_data','-ASCII')



areaname='CA1';
directory='/home/melisa/Escritorio/Datos_Ines/';
name=strcat(directory,areaname,'_syn_firing_posdir_all')
save_data=Syn(:,1);
save(name,'save_data','-ASCII')
name=strcat(directory,areaname,'_syn_firing_veldir_all')
save_data=Syn(:,2);
save(name,'save_data','-ASCII')
name=strcat(directory,areaname,'_syn_delta_posdir_all')
save_data=Syn(:,3);
save(name,'save_data','-ASCII')
name=strcat(directory,areaname,'_syn_delta_veldir_all')
save_data=Syn(:,4);
save(name,'save_data','-ASCII')
name=strcat(directory,areaname,'_syn_theta_posdir_all')
save_data=Syn(:,5);
save(name,'save_data','-ASCII')
name=strcat(directory,areaname,'_syn_theta_veldir_all')
save_data=Syn(:,6);
save(name,'save_data','-ASCII')