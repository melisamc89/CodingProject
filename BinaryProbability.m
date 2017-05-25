function [p pt]=BinaryProbability(matrix)

    [n m]=size(matrix);
    r=mean(reshape(matrix,n*m,1));
    newmatrix=zeros(n,m);
    
    for i=1:n
        for j=1:m
            if matrix(i,j)>=r
                newmatrix(i,j)=1;
            end
        end
    end     
    pt(1,:)=1.0*sum(newmatrix)/n;
    pt(2,:)=1-pt(1,:);
    p(1)=1.0*sum(pt(1,:))/m;
    p(2)=1-p(1);
end