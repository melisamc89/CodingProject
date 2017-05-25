%This calculation of information satisfices the processing information
%theorem
%This code computes the information between two neurons and the different
%behabioural aspects. For that, it took into account the correlation of two
%neurons activity (that inforamtion is in the matrix pn, pnt)
function [I hn hnt ht]=BehaviouralTriInformationSquareProtocol(pn_a,pnt_a,pt,dn,pos,freq,aspect)
    
    time1=5;
    h1=time1*freq;
    time2=30;
    h2=time2*freq;
    I=0;
    hn=0;
    hnt=0;
    ht=0;
    [a b c ns]=size(pnt_a);
    pnt=reshape(pnt_a,a*b*c,ns);
    sum(pnt);
    %imagesc(pnt)
    size(pnt_a);
    size(pnt);
    [a b c]=size(pn_a);
    pn=reshape(pn_a,a*b*c,1);
    sum(pn);
    nmax=a*b*c;
    
    [s1 s2]=ReadAlarmTime(pos);
    [sr sl]=ReadStartTime(pos);
    s1=floor(s1)+1;
    s2=floor(s2)+1;
    sr=floor(sr)+1;
    sl=floor(sl)+1;
    cero1=min(sr-s1);
    cero2=min(sl-s2);
    cant1=floor(cero1*freq);
    cant2=floor(cero2*freq);
    cant=min(cant1,cant2)-2;
    l=2*cant+2*h1+2*h2;

    
    %if ns>=2*cant+2*h1+2*h2
    
        switch aspect

            case 0
                for i=1:ns-1
                    aux(i)=0;
                    size(pnt);
                    for n=0:nmax-1
                        if pnt(n+1,i)~=0    
                            aux(i)=aux(i)-dn*pnt(n+1,i)*log(pnt(n+1,i));
                        end
                    end
                    aux(i)=aux(i)*pt(i);
                end
                hn=Entropy(pnt,pn,dn);
                hnt=sum(aux);
                I=hn-hnt;
                ht=0;
                for i=1:length(pt)
                    if pt(i)~=0
                    ht=ht-pt(i)*log(pt(i));
                    end
                end

            case 1

                [s1 s2]=ReadAlarmTime(pos);
                [sr sl]=ReadStartTime(pos);
                s1=floor(s1)+1;
                s2=floor(s2)+1;
                sr=floor(sr)+1;
                sl=floor(sl)+1;
                cero1=min(sr-s1);
                cero2=min(sl-s2);
                cant1=floor(cero1*freq);
                cant2=floor(cero2*freq);
                cant=min(cant1,cant2);

                m=2+h1+h2/6;

                px=zeros(1,m);
                    px(1)=sum(pt(1:cant));
                    for i=1:h1
                        px(i+1)=(pt(cant+i)+pt(ns-i));
                    end
                    for i=1:h2/6
                       aux1=sum(pt(cant+h1+1+(i-1)*6:cant+h1+i*6));
                       aux2=sum(pt(ns-h1-i*6:ns-h1-(i-1)*6-1));
                       px(h1+i+1)=(aux1+aux2);
                    end
                    px(m)=sum(pt(cant+h1+h2+1:2*cant+h1+h2));

                sum2=sum(px);
                if sum2~=0
                    px=px/sum2;
                end

                pnx=zeros(nmax,m);
                for n=0:nmax-1
                    a=find(pnt(n+1,1:cant));
                    if length(a)~=0
                        pnx(n+1,1)=sum(pnt(n+1,1:cant))/length(a);
                    end
                    for i=1:h1
                        if pnt(n+1,cant+i)~=0 && pnt(n+1,ns-i)~=0
                            pnx(n+1,i+1)=0.5*(pnt(n+1,cant+i)+pnt(n+1,ns-i));
                        else
                            if pnt(n+1,cant+i)~=0 && pnt(n+1,ns-i)==0
                                pnx(n+1,i+1)=pnt(n+1,cant+i);
                            else
                                pnx(n+1,i+1)=pnt(n+1,ns-i);
                            end
                        end
                    end
                    for i=1:h2/6
                        b=find(pnt(n+1,cant+h1+1+(i-1)*6:cant+h1+i*6));
                        c=find(pnt(n+1,ns-h1-i*6:ns-h1-(i-1)*6-1));
                        aux1=0;
                        aux2=0;
                        if length(b)~=0
                            aux1=sum(pnt(n+1,cant+h1+1+(i-1)*6:cant+h1+i*6))/length(b);
                        end
                        if length(c)~=0
                            aux2=sum(pnt(n+1,ns-h1-i*6:ns-h1-(i-1)*6-1))/length(c);
                        end
                        pnx(n+1,h1+i+1)=0.5*(aux1+aux2);
                    end
                    d=find(pnt(n+1,cant+h1+h2+1:2*cant+h1+h2));
                    if length(d)~=0
                        pnx(n+1,m)=sum(pnt(n+1,cant+h1+h2+1:2*cant+h1+h2))/length(d);
                    end
                end
                
                for i=1:m
                    sum1=sum(pnx(:,i));
                    if sum1~=0
                        pnx(:,i)=pnx(:,i)/sum1;
                    end
                end

                pn=zeros(1,nmax);
                for n=0:nmax-1
                    pn(n+1)=sum(pnx(n+1,:).*px);
                end
                if sum(pn)~=0
                    pn=pn/sum(pn);
                end

                hn=Entropy(pnx,pn,dn);

                hnt=0;
                for i=1:m
                    hnx=0;
                    for n=0:nmax-1
                        if pnx(n+1,i)~=0
                            hnx=hnx-dn*pnx(n+1,i)*log(pnx(n+1,i));
                        end
                    end
                    hnt=hnt+hnx*px(i);
                end
                I=hn-hnt;

                ht=0;
                for i=1:length(px)
                    if px(i)~=0
                    ht=ht-px(i)*log(px(i));
                    end
                end

            case 2

                pd(1)=sum(pt(1:floor(ns/2)));
                pd(2)=sum(pt(floor(ns/2)+1:ns));
                sum2=sum(pd);
                if sum2~=0
                    pd=pd/sum2;
                end

                pndir=zeros(nmax,2);
                for n=0:nmax-1
                    a=length(find(pnt(n+1,1:floor(ns/2))));
                    if a~=0
                        pndir(n+1,1)=sum(pnt(n+1,1:floor(ns/2)))/a;
                    end
                    b=length(find(pnt(n+1,floor(ns/2)+1:ns)));
                    if b~=0
                        pndir(n+1,2)=sum(pnt(n+1,floor(ns/2)+1:ns))/b;
                    end
                end
                
                for i=1:2
                    sum1=sum(pndir(:,i));
                    if sum1~=0
                        pndir(:,i)=pndir(:,i)/sum1;
                    end
                end

  
                pn=zeros(1,nmax);
                for n=0:nmax-1
                    pn(n+1)=sum(pndir(n+1,:).*pd);
                end

                if sum(pn)~=0
                    pn=pn/sum(pn);
                end
                
               hndir1=0;
               hndir2=0;
                for n=0:nmax-1
                    if pndir(n+1,1)~=0
                        hndir1=hndir1-dn*pndir(n+1,1)*log(pndir(n+1,1));
                    end
                    if pndir(n+1,2)~=0                
                        hndir2=hndir2-dn*pndir(n+1,2)*log(pndir(n+1,2));
                    end
                end
                
                hndir=pd(1)*hndir1+pd(2)*hndir2;
                hn=Entropy(pndir,pn,dn);

                hnt=hndir;
                I=hn-hnt;
                ht=0;
                for i=1:length(pd)
                    if pd(i)~=0
                    ht=ht-pd(i)*log(pd(i));
                    end
                end

            case 3

                [s1 s2]=ReadAlarmTime(pos);
                [sr sl]=ReadStartTime(pos);
                s1=floor(s1)+1;
                s2=floor(s2)+1;
                sr=floor(sr)+1;
                sl=floor(sl)+1;
                cero1=min(sr-s1);
                cero2=min(sl-s2);
                cant1=floor(cero1*freq);
                cant2=floor(cero2*freq);
                cant=min(cant1,cant2);

                m=3;
                ps=zeros(1,3);
                ps(1)=sum(pt(1:cant))+sum(pt(1+cant+h1+h2:2*cant+h1+h2));
                ps(2)=sum(pt(cant+1:cant+h1))+sum(pt(2*cant+h1+2*h2+1:ns));
                ps(3)=sum(pt(cant+h1+1:cant+h1+h2))+sum(pt(2*cant+h1+h2+1:ns-h1));
                sum2=sum(ps);

                if sum2~=0
                    ps=ps/sum2;
                end

                pns=zeros(nmax,3);
                aux=zeros(1,6);
                
                for n=0:nmax-1
                    a=length(find(pnt(n+1,1:cant)));
                    b=length(find(pnt(n+1,1+cant+h1+h2:2*cant+h1+h2)));
                    if a~=0
                        aux(1)=sum(pnt(n+1,1:cant))/a;
                    end
                    if b~=0
                        aux(2)=sum(pnt(n+1,1+cant+h1+h2:2*cant+h1+h2))/b;
                    end
                    pns(n+1,1)=0.5*(aux(1)+aux(2));
                    c=length(find(pnt(n+1,cant+1:cant+h1)));
                    d=length(find(pnt(n+1,2*cant+h1+2*h2+1:ns)));
                    if c~=0
                        aux(3)=sum(pnt(n+1,cant+1:cant+h1))/c;
                    end
                    if d~=0
                        aux(4)=sum(pnt(n+1,2*cant+h1+2*h2+1:ns))/d;
                    end
                    pns(n+1,2)=0.5*(aux(3)+aux(4));
                    e=length(find(pnt(n+1,cant+h1+1:cant+h1+h2)));
                    f=length(find(pnt(n+1,2*cant+h1+h2+1:ns-h1)));
                    if e~=0
                        aux(5)=sum(pnt(n+1,cant+h1+1:cant+h1+h2))/e;
                    end
                    if f~=0
                        aux(6)=sum(pnt(n+1,2*cant+h1+h2+1:ns-h1))/f;
                    end
                    pns(n+1,3)=0.5*(aux(5)+aux(6));
                end
                
                               
                for i=1:3
                    sum1=sum(pns(:,i));
                    if sum1~=0
                        pns(:,i)=pns(:,i)/sum1;
                    end
                end
                
                pn=zeros(1,nmax);
                for n=0:nmax-1
                    pn(n+1)=sum(pns(n+1,:).*ps);
                end
                if sum(pn)~=0
                    pn=pn/sum(pn);
                end

                hnt=0;
                for i=1:3
                    hns=0;
                    for n=0:nmax-1
                        if pns(n+1,i)~=0
                            hns=hns-dn*pns(n+1,i)*log(pns(n+1,i));
                        end
                    end
                    hnt=hnt+hns*ps(i);
                end

                hn=Entropy(pns,pn,dn);
                I=hn-hnt;
                ht=0;
                for i=1:length(ps)
                    if ps(i)~=0
                    ht=ht-ps(i)*log(ps(i));
                    end
                end
            case 4

                [s1 s2]=ReadAlarmTime(pos);
                [sr sl]=ReadStartTime(pos);
                s1=floor(s1)+1;
                s2=floor(s2)+1;
                sr=floor(sr)+1;
                sl=floor(sl)+1;
                cero1=min(sr-s1);
                cero2=min(sl-s2);
                cant1=floor(cero1*freq);
                cant2=floor(cero2*freq);
                cant=min(cant1,cant2);

                v=[cant,h1,h2,h2,h1,cant];
                m=2+2*h1+2*h2/6;

                px=zeros(1,m);
                    px(1)=sum(pt(1:cant));
                    for i=1:h1
                        px(i+1)=pnt(cant+i);
                    end
                    for i=1:h2/6
                        aux1=sum(pt(cant+h1+1+(i-1)*6:cant+h1+i*6));
                        px(h1+i+1)=aux1;
                    end
                    for i=1:h2/6
                        aux2=sum(pt(ns-h1-i*6:ns-h1-(i-1)*6-1));
                        px(1+h1+i+h2/6)=aux2;
                    end
                    for i=1:h1
                        px(1+h1+2*h2/6+i)=pt(ns-h1+i);
                    end
                    px(m)=sum(pt(cant+h1+h2+1:2*cant+h1+h2));

                sum2=sum(px);
                if sum2~=0
                    px=px/sum2;
                end

 
                pnx=zeros(nmax,m);
                for n=0:nmax-1
                    a=find(pnt(n+1,1:cant));
                    if length(a)~=0
                        pnx(n+1,1)=sum(pnt(n+1,1:cant))/length(a);
                    end
                    for i=1:h1
                            pnx(n+1,i+1)=pnt(n+1,cant+i)+pnt(n+1,ns-i);
                    end
                    for i=1:h2/6
                        b=find(pnt(n+1,cant+h1+1+(i-1)*6:cant+h1+i*6));
                        aux1=0;
                        if length(b)~=0
                            aux1=sum(pnt(n+1,cant+h1+1+(i-1)*6:cant+h1+i*6))/length(b);
                        end
                        pnx(n+1,h1+i+1)=aux1;
                    end
                    for i=1:h2/6
                        c=find(pnt(n+1,ns-h1-i*6:ns-h1-(i-1)*6-1));
                        aux2=0;
                        if length(c)~=0
                            aux2=sum(pnt(n+1,ns-h1-i*6:ns-h1-(i-1)*6-1))/length(c);
                        end
                        pnx(n+1,1+h1+i+h2/6)=aux2;
                    end
                    for i=1:h1
                        pnx(n+1,h1+2*h2/6+i)=pnt(n+1,ns-h1+i);
                    end
                    d=find(pnt(n+1,cant+h1+h2+1:2*cant+h1+h2));
                    if length(d)~=0
                        pnx(n+1,m)=sum(pnt(n+1,cant+h1+h2+1:2*cant+h1+h2))/length(d);
                    end
                end
                
                         
                for i=1:m
                    sum1=sum(pnx(:,i));
                    if sum1~=0
                        pnx(:,i)=pnx(:,i)/sum1;
                    end
                end

                pn=zeros(1,nmax);
                for n=0:nmax-1
                    pn(n+1)=sum(pnx(n+1,:).*px);
                end
                if sum(pn)~=0
                    pn=pn/sum(pn);
                end

                hnt=0;
                for i=1:m
                    hnx=0;
                    for n=0:nmax-1
                        if pnx(n+1,i)~=0
                            hnx=hnx-dn*pnx(n+1,i)*log(pnx(n+1,i));
                        end
                    end
                    hnt=hnt+hnx*px(i);
                end

                hn=Entropy(pnx,pn,dn);

                I=hn-hnt;
                ht=0;
                for i=1:length(pt)
                    if pt(i)~=0
                    ht=ht-pt(i)*log(pt(i));
                    end
                end

            case 5

                [s1 s2]=ReadAlarmTime(pos);
                [sr sl]=ReadStartTime(pos);
                s1=floor(s1)+1;
                s2=floor(s2)+1;
                sr=floor(sr)+1;
                sl=floor(sl)+1;
                cero1=min(sr-s1);
                cero2=min(sl-s2);
                cant1=floor(cero1*freq);
                cant2=floor(cero2*freq);
                cant=min(cant1,cant2);

                m=6;

                ps=zeros(1,6);

                    ps(1)=sum(pt(1:cant));
                    ps(2)=sum(pt(1+cant+h1+h2:2*cant+h1+h2));
                    ps(3)=sum(pt(cant+1:cant+h1));
                    ps(4)=sum(pt(2*cant+h1+2*h2+1:ns));
                    ps(5)=sum(pt(cant+h1+1:cant+h1+h2));
                    ps(6)=sum(pt(2*cant+h1+h2+1:ns-h1));

                sum2=sum(ps);
                if sum2~=0
                    ps=ps/sum2;
                end
                pns=zeros(nmax,6);
                aux=zeros(1,6);

                for n=0:nmax-1
                    a=length(find(pnt(n+1,1:cant)));
                    b=length(find(pnt(n+1,1+cant+h1+h2:2*cant+h1+h2)));
                    if a~=0
                        aux(1)=sum(pnt(n+1,1:cant))/a;
                    end
                    if b~=0
                        aux(2)=sum(pnt(n+1,1+cant+h1+h2:2*cant+h1+h2))/b;
                    end
                    c=length(find(pnt(n+1,cant+1:cant+h1)));
                    d=length(find(pnt(n+1,2*cant+h1+2*h2+1:ns)));
                    if c~=0
                        aux(3)=sum(pnt(n+1,cant+1:cant+h1))/c;
                    end
                    if d~=0
                        aux(4)=sum(pnt(n+1,2*cant+h1+2*h2+1:ns))/d;
                    end
                    e=length(find(pnt(n+1,cant+h1+1:cant+h1+h2)));
                    f=length(find(pnt(n+1,2*cant+h1+h2+1:ns-h1)));
                    if e~=0
                        aux(5)=sum(pnt(n+1,cant+1+h1:cant+h1+h2))/e;
                    end
                    if f~=0
                        aux(6)=sum(pnt(n+1,2*cant+h1+h2+1:ns-h1))/f;
                    end
                    pns(n+1,:)=aux;
                end
                
                 for i=1:6
                    sum1=sum(pns(:,i));
                    if sum1~=0
                        pns(:,i)=pns(:,i)/sum1;
                    end
                end
                
                pn=zeros(1,nmax);
                for n=0:nmax-1
                    pn(n+1)=sum(pns(n+1,:).*ps);
                end
                if sum(pn)~=0
                    pn=pn/sum(pn);
                end


                hnt=0;
                for i=1:m
                    hns=0;
                    for n=0:nmax-1
                        if pns(n+1,i)~=0
                            hns=hns-dn*pns(n+1,i)*log(pns(n+1,i));
                        end
                    end
                    hnt=hnt+hns*ps(i);
                end

                hn=Entropy(pns,pn,dn);

                I=hn-hnt;
                ht=0;
                for i=1:length(ps)
                    if ps(i)~=0
                    ht=ht-ps(i)*log(ps(i));
                    end
                end
        end
    %end
end