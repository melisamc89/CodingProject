function [pn1n2 pn1n2t pt corr_matrix_t]=BidimentionalGaussianProbability_v2018(matrix1,matrix2)
 
    nbins=50;
    [ntrials ns]=size(matrix1);
    variance1=var(matrix1,1);
    variance2=var(matrix2,1);
    auxvar1=min(variance1(find(variance1)));
    auxvar2=min(variance2(find(variance2)));

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
    
    ns=ns-1;
    pn1n2=zeros(nbins,nbins);
    pn1n2t=zeros(nbins,nbins,ns);
    pt=ones(1,ns);
    
    for t=1:ns
        corr_t=corrcoef(matrix1(:,t),matrix2(:,t));
        corr_matrix_t(t)=corr_t(1,2);
        cov_matrix_t=cov(matrix1(:,t),matrix2(:,t)); 
        r1=mean(matrix1(:,t));
        r2=mean(matrix2(:,t));
        mu=[r1;r2];
        x1=find(matrix1(:,t));
        x2=find(matrix2(:,t));
        total1=sum(matrix1(:,t));
        total2=sum(matrix2(:,t));
        if length(x1)<=5 || r1== 0 || total1<=20
            cov_matrix_t(1,1)=auxvar1;
            cov_matrix_t(1,2)=0;
            cov_matrix_t(2,1)=0;
        end
        if length(x2)<=5 || r2== 0 || total2<=20
            cov_matrix_t(2,2)=auxvar2; 
            cov_matrix_t(1,2)=0;
            cov_matrix_t(2,1)=0;
        end
        if (length(x1)<=5 && length(x2)<=5) || (r1== 0 && r2== 0) || (total1<=20 && total2<=20)
            cov_matrix_t(1,1)=auxvar1; 
            cov_matrix_t(2,2)=auxvar2;
            cov_matrix_t(1,2)=0;
            cov_matrix_t(2,1)=0;
        end
        %cov_matrix_t(1,1)
        %var(matrix1(:,t))
        %cov_matrix_t*inv(cov_matrix_t)
        determinant=det(cov_matrix_t);
            n1=n1min;
            i=1;
            while n1<n1max && i<nbins
                j=1;
                n2=n2min;
                while n2<n2max && j<nbins
                    n=[n1;n2];
                    pn1n2t(i,j,t)=exp(-(n-mu)'*inv(cov_matrix_t)*(n-mu)/2)/(2*pi*sqrt(determinant));
                    n2=n2+dn2;
                    j=j+1;
                end
                n1=n1+dn1;
                i=i+1;
            end
    end
    pt=pt/ns;      
    suma=sum(sum(pn1n2t,2),1);
    for t=1:ns
        pn1n2t(:,:,t)=pn1n2t(:,:,t)/suma(t);
    end
    pn1n2=sum(pn1n2t,3);
    suma=sum(sum(pn1n2));
    pn1n2=pn1n2/suma;
end