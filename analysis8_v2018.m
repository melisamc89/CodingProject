
clear all
load('EC1','-mat')
load('CA1','-mat')

count=ones(6,3);
percent=99;
total_cells=0;
for index=1:length(data)
    if isempty(data{index})~=1
        for j=1:length(data{index})
            total_cells=total_cells+1;
            cell_name(total_cells,:)=data{index}{j}.name;
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


%compute the data matrix
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

h1=-sum((probability1+eps).*log((probability1+eps)));
h2=-sum((probability2+eps).*log((probability2+eps)));
h3=-sum((probability3+eps).*log((probability3+eps)));


for k=1:200
    % generate a 6d random vector from a uniform distribution
    random_vector1=zeros(6,length(I));
    random_vector2=zeros(6,length(I));
    random_vector3=zeros(6,length(I));
    for l=1:6
        x=randperm(length(I));
        random_vector1(l,:)=I(l,x);
        random_vector2(l,:)=I_delta(l,x);
        random_vector3(l,:)=I_theta(l,x);
    end
    %generate random data
    probability_random1=zeros(size(new_element,1),1);
    probability_random2=zeros(size(new_element,1),1);
    probability_random3=zeros(size(new_element,1),1);

    for i=1:length(I)
        for j=1:size(new_element,1)
            if isequal(random_vector1(:,i),new_element(j,:)')
                probability_random1(j)=probability_random1(j)+1;
            end
            if isequal(random_vector2(:,i),new_element(j,:)')
                probability_random2(j)=probability_random2(j)+1;
            end
            if isequal(random_vector3(:,i),new_element(j,:)')
                probability_random3(j)=probability_random3(j)+1;
            end
        end
    end
    probability_random1=probability_random1/length(I);
    probability_random2=probability_random2/length(I_delta);
    probability_random3=probability_random3/length(I_theta);

    h1_r=-sum((probability_random1+eps).*log((probability_random1+eps)));
    h2_r=-sum((probability_random2+eps).*log((probability_random2+eps)));
    h3_r=-sum((probability_random3+eps).*log((probability_random3+eps)));

    h_random=-sum((matrix_random+eps)'.*log(matrix_random+eps)');
    dkl_atributes(1,k)=KLDivergence_v2018(probability1,probability_random1);
    dkl_atributes(2,k)=KLDivergence_v2018(probability2,probability_random2);
    dkl_atributes(3,k)=KLDivergence_v2018(probability3,probability_random3);
end
figure(1)
histogram(dkl_atributes(1,:),'Binwidth',0.1,'Normalization','probability')
hold on
histogram(dkl_atributes(2,:),'Binwidth',0.1,'Normalization','probability')
histogram(dkl_atributes(3,:),'Binwidth',0.1,'Normalization','probability')
legend('Firing','Delta','Theta')

number=1;
for i=1:3
    for j=i+1:3
       [h p]=ttest(dkl_atributes(1,:),dkl_atributes(2,:));
       prob(number)=p;
       number=number+1;
    end
end

%example

matrix_random1=zeros(7,20);
matrix_random1(1,1)=probability_random1(1);%/sum(probability_random(1));
matrix_random1(2,1:6)=probability_random1(2:7);%/sum(probability_random(2:7));
matrix_random1(3,1:15)=probability_random1(8:22);%/sum(probability_random(8:22));
matrix_random1(4,1:20)=probability_random1(23:42);%/sum(probability_random(23:42));
matrix_random1(5,1:15)=probability_random1(43:57);%/sum(probability_random(43:57));
matrix_random1(6,1:6)=probability_random1(58:63);%/sum(probability_random(58:63));
matrix_random1(7,1)=probability_random1(64);%/sum(probability_random(64));

matrix_random2=zeros(7,20);
matrix_random2(1,1)=probability_random2(1);%/sum(probability_random(1));
matrix_random2(2,1:6)=probability_random2(2:7);%/sum(probability_random(2:7));
matrix_random2(3,1:15)=probability_random2(8:22);%/sum(probability_random(8:22));
matrix_random2(4,1:20)=probability_random2(23:42);%/sum(probability_random(23:42));
matrix_random2(5,1:15)=probability_random2(43:57);%/sum(probability_random(43:57));
matrix_random2(6,1:6)=probability_random2(58:63);%/sum(probability_random(58:63));
matrix_random2(7,1)=probability_random2(64);%/sum(probability_random(64));

matrix_random3=zeros(7,20);
matrix_random3(1,1)=probability_random3(1);%/sum(probability_random(1));
matrix_random3(2,1:6)=probability_random3(2:7);%/sum(probability_random(2:7));
matrix_random3(3,1:15)=probability_random3(8:22);%/sum(probability_random(8:22));
matrix_random3(4,1:20)=probability_random3(23:42);%/sum(probability_random(23:42));
matrix_random3(5,1:15)=probability_random3(43:57);%/sum(probability_random(43:57));
matrix_random3(6,1:6)=probability_random3(58:63);%/sum(probability_random(58:63));
matrix_random3(7,1)=probability_random3(64);%/sum(probability_random(64));

matrix_z1=zeros(7,20);
matrix_z1(1,1)=(probability1(1)+eps);
matrix_z1(2,1:6)=(probability1(2:7)+eps);
matrix_z1(3,1:15)=(probability1(8:22)+eps);
matrix_z1(4,1:20)=(probability1(23:42)+eps);
matrix_z1(5,1:15)=(probability1(43:57)+eps);
matrix_z1(6,1:1:6)=(probability1(58:63)+eps);
matrix_z1(7,1)=(probability1(64)+eps);

matrix_z2=zeros(7,20);
matrix_z2(1,1)=(probability2(1)+eps);
matrix_z2(2,1:6)=(probability2(2:7)+eps);
matrix_z2(3,1:15)=(probability2(8:22)+eps);
matrix_z2(4,1:20)=(probability2(23:42)+eps);
matrix_z2(5,1:15)=(probability2(43:57)+eps);
matrix_z2(6,1:6)=(probability2(58:63)+eps);
matrix_z2(7,1)=(probability2(64)+eps);

matrix_z3=zeros(7,20);
matrix_z3(1,1)=(probability3(1)+eps);
matrix_z3(2,1:6)=(probability3(2:7)+eps);
matrix_z3(3,1:15)=(probability3(8:22)+eps);
matrix_z3(4,1:20)=(probability3(23:42)+eps);
matrix_z3(5,1:15)=(probability3(43:57)+eps);
matrix_z3(6,1:6)=(probability3(58:63)+eps);
matrix_z3(7,1)=(probability3(64)+eps);

subplot(2,3,1)
bar3(matrix_random1')
ylim([0 20])

subplot(2,3,4)
bar3(matrix_z1')
ylim([0 20])

subplot(2,3,2)
bar3(matrix_random1')
ylim([0 20])

subplot(2,3,5)
bar3(matrix_z2')
ylim([0 20])

subplot(2,3,3)
bar3(matrix_random3')
ylim([0 20])

subplot(2,3,6)
bar3(matrix_z3')
ylim([0 20])

%% the same but for conditioned encoded atributes


p_marg1=sum(matrix_z1,2)/sum(sum(matrix_z1,2));
p_marg2=sum(matrix_z2,2)/sum(sum(matrix_z2,2));
p_marg3=sum(matrix_z3,2)/sum(sum(matrix_z3,2));

%plot marginls! code this!

matrix_z1=zeros(7,20);
matrix_z1(1,1)=(probability1(1)+eps)/sum(probability1(1)+eps);
matrix_z1(2,1:6)=(probability1(2:7)+eps)/sum(probability1(2:7)+eps);
matrix_z1(3,1:15)=(probability1(8:22)+eps)/sum(probability1(8:22)+eps);
matrix_z1(4,1:20)=(probability1(23:42)+eps)/sum(probability1(23:42)+eps);
matrix_z1(5,1:15)=(probability1(43:57)+eps)/sum(probability1(43:57)+eps);
matrix_z1(6,1:6)=(probability1(58:63)+eps)/sum(probability1(58:63)+eps);
matrix_z1(7,1)=(probability1(64)+eps)/sum(probability1(64)+eps);

matrix_z2=zeros(7,20);
matrix_z2(1,1)=(probability2(1)+eps)/sum(probability2(1)+eps);
matrix_z2(2,1:6)=(probability2(2:7)+eps)/sum(probability2(2:7)+eps);
matrix_z2(3,1:15)=(probability2(8:22)+eps)/sum(probability2(8:22)+eps);
matrix_z2(4,1:20)=(probability2(23:42)+eps)/sum(probability2(23:42)+eps);
matrix_z2(5,1:15)=(probability2(43:57)+eps)/sum(probability2(43:57)+eps);
matrix_z2(6,1:6)=(probability2(58:63)+eps)/sum(probability2(58:63)+eps);
matrix_z2(7,1)=(probability2(64)+eps)/sum(probability2(64)+eps);

matrix_z3=zeros(7,20);
matrix_z3(1,1)=(probability3(1)+eps)/sum(probability3(1)+eps);
matrix_z3(2,1:6)=(probability3(2:7)+eps)/sum(probability3(2:7)+eps);
matrix_z3(3,1:15)=(probability3(8:22)+eps)/sum(probability3(8:22)+eps);
matrix_z3(4,1:20)=(probability3(23:42)+eps)/sum(probability3(23:42)+eps);
matrix_z3(5,1:15)=(probability3(43:57)+eps)/sum(probability3(43:57)+eps);
matrix_z3(6,1:6)=(probability3(58:63)+eps)/sum(probability3(58:63)+eps);
matrix_z3(7,1)=(probability3(64)+eps)/sum(probability3(64)+eps);

for k=1:200
    % generate a 6d random vector from a uniform distribution
    random_vector1=zeros(6,length(I));
    random_vector2=zeros(6,length(I));
    random_vector3=zeros(6,length(I));
    for l=1:6
        x=randperm(length(I));
        random_vector1(l,:)=I(l,x);
        random_vector2(l,:)=I_delta(l,x);
        random_vector3(l,:)=I_theta(l,x);
    end
    %generate random data
    probability_random1=zeros(size(new_element,1),1);
    probability_random2=zeros(size(new_element,1),1);
    probability_random3=zeros(size(new_element,1),1);

    for i=1:length(I)
        for j=1:size(new_element,1)
            if isequal(random_vector1(:,i),new_element(j,:)')
                probability_random1(j)=probability_random1(j)+1;
            end
            if isequal(random_vector2(:,i),new_element(j,:)')
                probability_random2(j)=probability_random2(j)+1;
            end
            if isequal(random_vector3(:,i),new_element(j,:)')
                probability_random3(j)=probability_random3(j)+1;
            end
        end
    end
    probability_random1=probability_random1/length(I);
    probability_random2=probability_random2/length(I_delta);
    probability_random3=probability_random3/length(I_theta);
    
    matrix_random1=zeros(7,20);
    matrix_random1(1,1)=probability_random1(1)/sum(probability_random1(1)+eps);
    matrix_random1(2,1:6)=probability_random1(2:7)/sum(probability_random1(2:7)+eps);
    matrix_random1(3,1:15)=probability_random1(8:22)/sum(probability_random1(8:22)+eps);
    matrix_random1(4,1:20)=probability_random1(23:42)/sum(probability_random1(23:42)+eps);
    matrix_random1(5,1:15)=probability_random1(43:57)/sum(probability_random1(43:57)+eps);
    matrix_random1(6,1:6)=probability_random1(58:63)/sum(probability_random1(58:63)+eps);
    matrix_random1(7,1)=probability_random1(64)/sum(probability_random1(64)+eps);

    matrix_random2=zeros(7,20);
    matrix_random2(1,1)=probability_random2(1)/sum(probability_random2(1)+eps);
    matrix_random2(2,1:6)=probability_random2(2:7)/sum(probability_random2(2:7)+eps);
    matrix_random2(3,1:15)=probability_random2(8:22)/sum(probability_random2(8:22)+eps);
    matrix_random2(4,1:20)=probability_random2(23:42)/sum(probability_random2(23:42)+eps);
    matrix_random2(5,1:15)=probability_random2(43:57)/sum(probability_random2(43:57)+eps);
    matrix_random2(6,1:6)=probability_random2(58:63)/sum(probability_random2(58:63)+eps);
    matrix_random2(7,1)=probability_random2(64)/sum(probability_random2(64)+eps);

    matrix_random3=zeros(7,20);
    matrix_random3(1,1)=probability_random3(1)/sum(probability_random3(1)+eps);
    matrix_random3(2,1:6)=probability_random3(2:7)/sum(probability_random3(2:7)+eps);
    matrix_random3(3,1:15)=probability_random3(8:22)/sum(probability_random3(8:22)+eps);
    matrix_random3(4,1:20)=probability_random3(23:42)/sum(probability_random3(23:42)+eps);
    matrix_random3(5,1:15)=probability_random3(43:57)/sum(probability_random3(43:57)+eps);
    matrix_random3(6,1:6)=probability_random3(58:63)/sum(probability_random3(58:63)+eps);
    matrix_random3(7,1)=probability_random3(64)/sum(probability_random3(64)+eps);
    
    for condition=1:7
        dkl(condition,1,k)=p_marg1(condition)*KLDivergence_v2018(matrix_z1(condition,:),matrix_random1(condition,:));
        dkl(condition,2,k)=p_marg2(condition)*KLDivergence_v2018(matrix_z2(condition,:),matrix_random2(condition,:));
        dkl(condition,3,k)=p_marg3(condition)*KLDivergence_v2018(matrix_z3(condition,:),matrix_random3(condition,:));
    end
end

dkl_mean=mean(dkl,3)
bar3(dkl_mean')

%examples plot

subplot(2,3,1)
bar3(matrix_random1')
ylim([0 20])

subplot(2,3,4)
bar3(matrix_z1')
ylim([0 20])

subplot(2,3,2)
bar3(matrix_random1')
ylim([0 20])

subplot(2,3,5)
bar3(matrix_z2')
ylim([0 20])

subplot(2,3,3)
bar3(matrix_random3')
ylim([0 20])

subplot(2,3,6)
bar3(matrix_z3')
ylim([0 20])


figure(2)
number=1;
for i=1:3
    for condition=1:7
        subplot(3,7,number)
        number=number+1;
        histogram(dkl(condition,i,:),20,'Normalization','probability')
    end
end
    
