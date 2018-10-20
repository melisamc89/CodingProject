
clear all
%load('EC1','-mat')
load('CA1','-mat')

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

matrix1=information.*I;
matrix2=information_delta.*I_delta;
matrix3=information_theta.*I_theta;

%al possible combinations

vector=[1,1,1,1,1,1,0,0,0,0,0,0];
matrix=nchoosek(vector,6);
number=2;
element(1,:)=matrix(1,:);
for i=1:length(matrix)
    index=0;
    for j=1:size(element,1)
        if ~isequal(matrix(i,:),element(j,:))
            index=index+1;
        end
    end
    if index==size(element,1)
        element(number,:)=matrix(i,:);
        number=number+1;
    end
end

new_element(1,:)=element(1,:);
number2=2;
for k=1:size(element)
    permutations=perms(element(k,:));
    for i=1:length(permutations)
        index=0;
        for j=1:size(new_element,1)
            if ~isequal(permutations(i,:),new_element(j,:))
                index=index+1;
            end
        end
        if index==size(new_element,1)
            new_element(number2,:)=permutations(i,:);
            number2=number2+1;
        end
    end
end

% generate a 6d random vector from a uniform distribution
random_vector=rand(6,total_cells);
random_vector(find(random_vector>0.5))=1;
random_vector(find(random_vector<=0.5))=0;
%generate random data
probability_random=zeros(size(new_element,1),1);
for i=1:length(matrix1)
    for j=1:size(new_element,1)
        if isequal(random_vector(:,i),new_element(j,:)')
            probability_random(j)=probability_random(j)+1;
        end
    end
end
probability_random=probability_random/length(I);
h=-sum((probability_random+eps).*log((probability_random+eps)));
matrix_random=zeros(7,20);
matrix_random(1,1)=probability_random(1)/sum(probability_random(1));
matrix_random(2,1:6)=probability_random(2:7)/sum(probability_random(2:7));
matrix_random(3,1:15)=probability_random(8:22)/sum(probability_random(8:22));
matrix_random(4,1:20)=probability_random(23:42)/sum(probability_random(23:42));
matrix_random(5,1:15)=probability_random(43:57)/sum(probability_random(43:57));
matrix_random(6,1:6)=probability_random(58:63)/sum(probability_random(58:63));
matrix_random(7,1)=probability_random(64)/sum(probability_random(64));
bar3(matrix_random')
h_random=-sum((matrix_random+eps)'.*log(matrix_random+eps)');


%compare with I
probability1=zeros(size(new_element,1),1);
for i=1:length(matrix1)
    for j=1:size(new_element,1)
        if isequal(I(:,i),new_element(j,:)')
            probability1(j)=probability1(j)+1;
        end
    end
end
probability1=probability1/length(I);
dkl_atributes(1)=KLDivergence_v2018(probability1,probability_random)/h;
%compare with I_delta
probability=zeros(size(new_element,1),1);
for i=1:length(matrix1)
    for j=1:size(new_element,1)
        if isequal(I_delta(:,i),new_element(j,:)')
            probability(j)=probability(j)+1;
        end
    end
end
probability2=probability/length(I_delta);
dkl_atributes(2)=KLDivergence_v2018(probability2,probability_random)/h;
%compare with I_theta
probability=zeros(size(new_element,1),1);
for i=1:length(matrix1)
    for j=1:size(new_element,1)
        if isequal(I_theta(:,i),new_element(j,:)')
            probability(j)=probability(j)+1;
        end
    end
end
probability3=probability/length(I_theta);
dkl_atributes(3)=KLDivergence_v2018(probability3,probability_random)/h;

subplot(1,4,1)
bar3(matrix_random')

subplot(1,4,2)
matrix_z1=zeros(7,20);
matrix_z1(1,1)=probability1(1)/sum(probability1(1));
matrix_z1(2,1:6)=probability1(2:7)/sum(probability1(2:7));
matrix_z1(3,1:15)=probability1(8:22)/sum(probability1(8:22));
matrix_z1(4,1:20)=probability1(23:42)/sum(probability1(23:42));
matrix_z1(5,1:15)=probability1(43:57)/sum(probability1(43:57));
matrix_z1(6,1:6)=probability1(58:63)/sum(probability1(58:63));
matrix_z1(7,1)=probability1(64)/sum(probability1(64));
bar3(matrix_z1')

subplot(1,4,3)
matrix_z2=zeros(7,20);
matrix_z2(1,1)=probability2(1)/sum(probability2(1));
matrix_z2(2,1:6)=probability2(2:7)/sum(probability2(2:7));
matrix_z2(3,1:15)=probability2(8:22)/sum(probability2(8:22));
matrix_z2(4,1:20)=probability2(23:42)/sum(probability2(23:42));
matrix_z2(5,1:15)=probability2(43:57)/sum(probability2(43:57));
matrix_z2(6,1:6)=probability2(58:63)/sum(probability2(58:63));
matrix_z2(7,1)=probability2(64)/sum(probability2(64));
bar3(matrix_z2')

subplot(1,4,4)
matrix_z3=zeros(7,20);
matrix_z3(1,1)=probability3(1)/sum(probability3(1));
matrix_z3(2,1:6)=probability3(2:7)/sum(probability3(2:7));
matrix_z3(3,1:15)=probability3(8:22)/sum(probability3(8:22));
matrix_z3(4,1:20)=probability3(23:42)/sum(probability3(23:42));
matrix_z3(5,1:15)=probability3(43:57)/sum(probability3(43:57));
matrix_z3(6,1:6)=probability3(58:63)/sum(probability3(58:63));
matrix_z3(7,1)=probability3(64)/sum(probability3(64));
bar3(matrix_z3')

for condition=1:7
    dkl(condition,1)=KLDivergence_v2018(matrix_z1(condition,:),matrix_random(condition,:))/(h_random(condition)+eps);
    dkl(condition,2)=KLDivergence_v2018(matrix_z2(condition,:),matrix_random(condition,:))/(h_random(condition)+eps);
    dkl(condition,3)=KLDivergence_v2018(matrix_z3(condition,:),matrix_random(condition,:))/(h_random(condition)+eps);
end
bar(dkl)    
