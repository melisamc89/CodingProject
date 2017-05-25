function [circular_mean circular_variance]=CircularAnalysis(matrix)

    [n m]=size(matrix);
    
    for j=1:m
        sum1=0;
        sum2=0;
        count=0;
        for i=1:n
               if matrix(i,j)~=0
                    sum1=sum1+cos(matrix(i,j));
                    sum2=sum2+sin(matrix(i,j));
                    count=count+1;
               end
        end
        phase(j)=atan2(sum2,sum1);
        r(j)=sqrt((sum1/count)^2+(sum2/count)^2);
        
        if find(matrix(:,j))<10
        phase(j)=0;
        r(j)=1;
        end
    end

    circular_mean=phase;
    circular_variance=1-r;
end