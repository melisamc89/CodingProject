function [ptheta_n pn p_theta vec phase_mean kappa]=ConditionalVonMisesProbability2(n_matrix,matrixa,matrixb)

    [n m]=size(n_matrix);
    nmax=max(max(n_matrix));
    phase=reshape([matrixa,matrixb],n*m,150);

    phase_mean=zeros(nmax,1);
    phase_std=zeros(nmax,1);
    pn=zeros(1,nmax);
    counter=0;
    
    for i=1:nmax
        pos=find(n_matrix==i);
        if length(pos)
            pn(i)=length(pos);
            counter=counter+length(pos);
            aux=phase(pos(1),1);
            for k=2:150
                if phase(pos(1),k)~=0
                    aux=cat(2,aux,phase(pos(1),k));
                end
            end
            for j=2:length(pos)
                for k=1:length(phase(pos(j),:))
                    if phase(pos(j),k)~=0
                        aux=cat(2,aux,phase(pos(j),k));
                    end
                end
            end
            phase_mean(i)=circ_mean(aux');
            phase_std(i)=circ_var(aux');
            clear aux
        end
    end
    pn=pn/counter;
   
    kappa=zeros(1,length(phase_std));
    for i=1:length(phase_std)
        if phase_std(i)~=0
            kappa(i)=(-1/(2*log(1-phase_std(i))));
        end
    end
    
    phase_max=150;
    ptheta_n=zeros(phase_max,nmax);
    for n=1:nmax
        for i=1:phase_max
            o=i*2*pi/phase_max;
            ptheta_n(i,n)=exp(cos(o-phase_mean(n))*kappa(n))/(2*pi*besseli(0,kappa(n)));
            if phase_mean(n)==0
                ptheta_n(i,n)=0;
            end
        end
    end
    [n m]=size(ptheta_n);
    suma=sum(ptheta_n);
    for i=1:m
        ptheta_n(:,i)=ptheta_n(:,i)/suma(i);
    end
    p_theta=ptheta_n*pn';
    sum(p_theta)
    vec=1:nmax;
end