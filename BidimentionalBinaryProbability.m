function [p pt]=BidimentionalBinaryProbability(matrix1,matrix2)

    [n m]=size(matrix1);
    r1=mean(reshape(matrix1,n*m,1));
    r2=mean(reshape(matrix2,n*m,1));
    newmatrix=zeros(n,m);

    for i=1:n
        for t=1:m
            if matrix1(i,t)>=r1 && matrix2(i,t)>=r2
                newmatrix(i,t)=1;
            end
            if matrix2(i,t)>=r1 && matrix2(i,t)<r2
                newmatrix(i,t)=3;
            end
            if matrix2(i,t)<r1 && matrix2(i,t)>=r2
                newmatrix(i,t)=0;
            end
            if matrix2(i,t)<r1 && matrix2(i,t)<r2
                newmatrix(i,t)=2;
            end        
        end
    end
    for i=1:2
        for j=1:2
            for t=1:m
                pt(i,j,t)=1.0*length(find(newmatrix(:,t)==i*2-j));
            end
        end
    end
    
    for t=1:m
        pt(:,:,t)=pt(:,:,t)/sum(sum(pt(:,:,t)));
    end
    p=sum(pt,3)/m;
end