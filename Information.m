%computes information from conditional probability

function [I hy hyx hx]=Information(py,px,py_x)

    hy=0;
    for i=1:length(py)
        hy=hy-py(i)*log(py(i)+eps);
    end
    hx=0;
    for i=1:length(px)
        hx=hx-px(i)*log(px(i)+eps);
    end
    
    [n m]=size(py_x);
    hyx=0;
    for i=1:n
        for j=1:m
            hyx=hyx-px(j)*py_x(i,j)*log((py_x(i,j)+eps));
        end
    end

    I=hy-hyx;
end