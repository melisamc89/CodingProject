% This function computes information by the suposition that the p(n) has an
% gaussian distribution and that p(n|t) also is gaussian. Where the n0 ans
% de sigma0 are taken from the data

function I0=GaussianInformation(matrix,nt)

    [l m]=size(matrix);
    [r var]=MeanRate(matrix,nt);
    
    r0=mean(r);
    sd0=std(r);
    
    I0=0;
    if sd0 ~= 0 && r0~= 0
        for i=1:m
            nmin=r(i)-2*sqrt(var(i));
            nmax=r(i)+2*sqrt(var(i));
            dn=abs(nmax-nmin)/100;
            n=nmin;
            hn=0;
            if var(i) ~= 0 && r(i) ~= 0
                while n < nmax
                    aux_nt=exp(-(n-r(i))*(n-r(i))/(2*var(i)))/sqrt(2*pi*var(i));
                    aux_n=0;
                    for j=1:m
                        if var(j) ~= 0 && r(j) ~= 0
                        aux_n=aux_n+exp(-(n-r(j))*(n-r(j))/(2*var(j)))/sqrt(2*pi*var(j));
                        end
                    end
                    aux_n=aux_n/m;
                    hn=hn+aux_nt*log(aux_nt/aux_n+eps)*dn;
                    n=n+dn;
                end
            I0=I0+hn/m;
            end
        end
    end
    
end