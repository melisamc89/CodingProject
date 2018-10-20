function [circular_mean circular_variance]=CircularAnalysis_v2018(matrix)

    [n m l]=size(matrix);    
    for i=1:n
        for j=1:m
            index=find(matrix(i,j,:));
            if isempty(index)~=1
                c_mean=circ_mean(matrix(i,j,index));
                c_var=circ_var(reshape(matrix(i,j,index),[length(index) 1]));
                circular_mean(i,j)=c_mean;
                circular_variance(i,j)=c_var;
            else
                circular_mean(i,j)=0;
                circular_variance(i,j)=1;
            end
        end
    end
end