%computes information from join probability

function I=Information2(pyx)

    px=sum(pyx,2);
    py=sum(pyx);
    
    [n m]=size(pyx);
    I=0;
    for i=1:n
        for j=1:m
           if py(i)~=0 && px(j)~=0 && pyx(i,j)~=0
            I=I+pyx(i,j)*log((pyx(i,j))/(px(j)*py(i)));
           end
        end
    end
    
    %I=hy-hyx;
    
end