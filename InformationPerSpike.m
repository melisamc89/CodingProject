function I0=InformationPerSpike(matrix,nt)

    [n m]=size(matrix);
    [r var]=MeanRate(matrix,nt);

    media=mean(r);
    
    I0=0;
    for i=1:m
        if r(i)>0
            I0=I0+r(i)*log(r(i)/media)/media;
        end
    end

end