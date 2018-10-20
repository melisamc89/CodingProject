function [pn1 pn1t pt n1min n1max]=GaussianProbability_v2018(matrix)
 
    nbins=100;
    [ntrials ns]=size(matrix);
    r_1=mean(matrix,1);
    var_1=var(matrix,1);
    auxvar=min(var_1(find(var_1)));
    
    r=mean(r_1);
    variance=var(reshape(matrix,ntrials*ns,1));
    n1min=r-3*sqrt(variance);
    n1max=r+3*sqrt(variance);
    dn1=1.0*(n1max-n1min)/nbins;
    
    ns=ns-1;
    pn1=zeros(nbins,1);
    pn1t=zeros(nbins,ns);
    pt=ones(1,ns);
    for t=1:ns
        x1=find(matrix(:,t));
        total=sum(matrix(:,t));
        r1=mean(matrix(:,t));
        var1=var(matrix(:,t));
        if length(x1)<=5 || r1== 0 || total <= 20
            var1=auxvar;
        end
        n1=n1min;
        i=1;
        while n1 < n1max && i<=nbins
            pn1t(i,t)=exp(-(n1-r1)*(n1-r1)/(2*var1))/sqrt(2*pi*var1);
            n1=n1+dn1;
            i=i+1;
        end
    end
    pt=pt/ns;
    suma=sum(pn1t,1);
    for t=1:ns
        pn1t(:,t)=pn1t(:,t)/suma(t);
    end
    pn1=sum(pn1t,2);
    suma=sum(pn1);
    pn1=pn1/suma;
end