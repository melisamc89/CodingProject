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

matrix1=information.*I;
matrix2=information_delta.*I_delta;
matrix3=information_theta.*I_theta;

%comparacion de codigos usando la idea de probabilidad

vector=[1,1,1,0,0,0];
matrix=nchoosek(vector,3);
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

% generate a 3d random vector from a uniform distribution
random_vector=rand(3,total_cells);
random_vector(find(random_vector>0.5))=1;
random_vector(find(random_vector<=0.5))=0;
%compare with I
probability1=zeros(size(new_element,1),1);
for i=1:length(random_vector)
    for j=1:size(new_element,1)
        if isequal(random_vector(:,i),new_element(j,:)')
            probability1(j)=probability1(j)+1;
        end
    end
end
probability1=probability1/length(random_vector);
h=-sum((probability1+eps).*log((probability1+eps)));
matrix_random=zeros(4,3);
matrix_random(1,1)=probability1(1)/sum(probability1(1));
matrix_random(2,1:3)=probability1(2:4)/sum(probability1(2:4));
matrix_random(3,1:3)=probability1(5:7)/sum(probability1(5:7));
matrix_random(4,1)=probability1(8)/sum(probability1(8));
h_random=-sum((matrix_random+eps)'.*log(matrix_random+eps)');


for aspect=1:6
    matrix_a(1,:)=I(aspect,:);
    matrix_a(2,:)=I_delta(aspect,:);
    matrix_a(3,:)=I_theta(aspect,:);
    %compare with I
    probability=zeros(size(new_element,1),1);
    for i=1:length(matrix_a)
        for j=1:size(new_element,1)
            if isequal(matrix_a(:,i),new_element(j,:)')
                probability(j)=probability(j)+1;
            end
        end
    end
    probability=probability/length(I);
    %dkl(aspect)=KLDivergence_v2018(probability,probability1)/h;
    %subplot(2,3,aspect)
    matrix_z=zeros(4,3);
    matrix_z(1,1)=probability(1)/sum(probability(1));
    matrix_z(2,1:3)=probability(2:4)/sum(probability(2:4));
    matrix_z(3,1:3)=probability(5:7)/sum(probability(5:7));
    matrix_z(4,1)=probability(8)/sum(probability(8));
    for condition=1:4
        dkl_codes(aspect,condition)=KLDivergence_v2018(matrix_z(condition,:),matrix_random(condition,:))/(h_random(condition)+eps);
    end
    %bar3(matrix_z')
end
bar3(dkl_codes')
