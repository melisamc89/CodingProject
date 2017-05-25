function pn1n2=BidimentionalGaussianProbability2(matrix1,matrix2)
 
    nbins=50;
    [ntrials ns]=size(matrix1);
    
    x=reshape(matrix1,ntrials*ns,1);
    r_1=mean(x);
    variance_1=var(x);
    n1min=r_1-3*sqrt(variance_1);
    n1max=r_1+3*sqrt(variance_1);
    x=reshape(matrix2,ntrials*ns,1);
    r_2=mean(x);
    variance_2=var(x);
    n2min=r_2-3*sqrt(variance_2);
    n2max=r_2+3*sqrt(variance_2);
    dn1=1.0*(n1max-n1min)/nbins;
    dn2=1.0*(n2max-n2min)/nbins;
    
    pn1n2=zeros(nbins,nbins);        
    
    cov_matrix=cov(reshape(matrix1,ntrials*ns,1),reshape(matrix2,ntrials*ns,1)); 
    mu=[r_1;r_2];
    determinant=det(cov_matrix);
    n1=n1min;
    i=1;
    while n1<n1max && i<nbins
        j=1;
        n2=n2min;
        while n2<n2max && j<nbins
             n=[n1;n2];
             pn1n2(i,j)=exp(-(n-mu)'*inv(cov_matrix)*(n-mu)/2)/(2*pi*sqrt(determinant));
             n2=n2+dn2;
             j=j+1;
        end
        n1=n1+dn1;
        i=i+1;
    end
    suma=sum(sum(pn1n2));
    pn1n2=pn1n2/suma;
end