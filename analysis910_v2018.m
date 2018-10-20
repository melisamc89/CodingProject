
clear all
load('EC1','-mat')
%load('CA1','-mat')

count=ones(6,3);
percent=99;
total_cells=0;
for index=1:length(data)
    if isempty(data{index})~=1
        for j=1:length(data{index})
            total_cells=total_cells+1;
            for aspect=1:6
                test=sort(data{index}{j}.TEST.I(aspect,:));
                information(aspect,total_cells)=data{index}{j}.I(aspect);
                information_delta(aspect,total_cells)=data{index}{j}.I_delta(aspect);
                information_theta(aspect,total_cells)=data{index}{j}.I_theta(aspect);                
                I(aspect,total_cells)=0;
                I_delta(aspect,total_cells)=0;
                I_theta(aspect,total_cells)=0;
                if data{index}{j}.I(aspect)>test(percent)
                    I(aspect,total_cells)=1;
                end
                test=sort(data{index}{j}.TEST.I_delta(aspect,:));
                if data{index}{j}.I_delta(aspect)>test(99)
                    I_delta(aspect,total_cells)=1;
                end
                test=sort(data{index}{j}.TEST.I_theta(aspect,:));
                if data{index}{j}.I_theta(aspect)>test(99)
                    I_theta(aspect,total_cells)=1;
                end
            end
        end
    end
end

number=sum(I);
for i=0:6
    list=I(:,find(number==i));
    prob(i+1,:)=sum(list,2)/length(find(number==i));
    clear list
end


number=sum(I_delta);
for i=0:6
    list=I_delta(:,find(number==i));
    prob(i+1,:)=sum(list,2)/length(find(number==i));
    clear list
end


number=sum(I_theta);
for i=0:6
    list=I_theta(:,find(number==i));
    prob(i+1,:)=sum(list,2)/length(find(number==i));
    clear list
end