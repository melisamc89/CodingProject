function [matrix1 matrix2]=BinaryCountDistribution(matrixa,matrixb)

    [n1 m1]=size(matrixa);
    [n2 m2]=size(matrixb);
    totalmatrix=[matrixa,matrixb];
    [n m]=size(totalmatrix);
    
    for i=1:n
        for j=1:m
            if totalmatrix(i,j)==0
                newmatrix(i,j)=totalmatrix(i,j);
            else
                newmatrix(i,j)=1;
            end

        end
    end
    
    for i=1:n1
        for j=1:m1
            matrix1(i,j)=newmatrix(i,j);
        end
    end
    
    for i=1:n2
        for j=1:m2
            matrix2(i,j)=newmatrix(i,m1+j);
        end
    end

end