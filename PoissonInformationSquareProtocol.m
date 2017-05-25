%This calculation of information satisfices the processing information
%theorem

function [I hn hnt ht]=PoissonInformationSquareProtocol(pn,pnt,pos,freq,aspect)
    
    time1=6;
    h1=time1*freq;
    time2=30;
    h2=time2*freq;
    
    [nmax ns]=size(pnt);
    hn=0;
    for n=0:nmax-1
        if pn(n+1)~=0
           hn=hn-pn(n+1)*log(pn(n+1));
        end
    end
              
    switch aspect
 
        case 0
            
            for i=1:ns-1
                aux(i)=0;
                for n=0:nmax-1
                    if pnt(n+1,i)~=0
                        aux(i)=aux(i)-pnt(n+1,i)*log(pnt(n+1,i));
                    end
                end
            end
            hnt=sum(aux)/ns;
            I=hn-hnt;
            ht=log(ns);
            
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
            
            v=[cant,h1,h2,h2,h1,cant];
            m=2+h1+h2/6;
            
            pnx=zeros(nmax,m);
            
            for n=0:nmax-1
                pnx(n+1,1)=mean(pnt(n+1,1:cant));
                for i=1:h1
                    pnx(n+1,i+1)=(pnt(n+1,cant+i)+pnt(n+1,ns-(cant+i)))*0.5;
                end
                for i=1:h2/6
                    aux1=mean(pnt(n+1,cant+h1+1+(i-1)*6:cant+h1+1+i*6));
                    aux2=mean(pnt(n+1,ns-(cant+h1+1+i*6):ns-(cant+h1+1+(i-1)*6)));
                    pnx(n+1,h1+i+1)=0.5*(aux1+aux2);
                end
                pnx(n+1,m)=mean(pnt(n+1,ns-cant:ns));

            end

            hnt=0;
            for i=1:m
                hnx=0;
                for n=0:nmax-1
                    if pnx(n+1,i)~=0
                        hnx=hnx-pnx(n+1,i)*log(pnx(n+1,i));
                    end
                end
                hnt=hnt+hnx/m;
            end
            I=hn-hnt;
            ht=log(m);
        case 2
            
            for n=0:nmax-1
                pndir(n+1,1)=mean(pnt(n+1,1:floor(ns/2)));
                pndir(n+1,2)=mean(pnt(n+1,floor(ns/2):ns));
            end
            hndir=0;
            for n=0:nmax-1
                if pndir(n+1,1)~=0
                    hndir=hndir-pndir(n+1,1)*log(pndir(n+1,1));
                end
                if pndir(n+1,2)~=0                
                    hndir=hndir-pndir(n+1,2)*log(pndir(n+1,2));
                end
            end
            hnt=hndir/2;
            I=hn-hnt;
            ht=log(2);
            
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
            
            pns=zeros(nmax,3);
            
            for n=0:nmax-1
                pns(n+1,1)=0.5*mean(pnt(n+1,1:cant))+0.5*mean(pnt(n+1,ns-cant:ns));
                pns(n+1,2)=(mean(pnt(n+1,cant:cant+h1))+mean(pnt(n+1,ns-(cant+h1):ns-cant)))*0.5;
                pns(n+1,3)=mean(pnt(n+1,cant+h1+1:ns-cant));
            end
            
            hnt=0;
            for i=1:3
                hns=0;
                for n=0:nmax-1
                    if pns(n+1,i)~=0
                        hns=hns-pns(n+1,i)*log(pns(n+1,i));
                    end
                end
                hnt=hnt+hns/m;
            end
            I=hn-hnt;
            ht=log(3);
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
            
            pnx=zeros(nmax,m);
            
            for n=0:nmax-1
                pnx(n+1,1)=mean(pnt(n+1,1:cant));
                for i=1:h1
                    pnx(n+1,i+1)=pnt(n+1,cant+i);
                end
                for i=1:h2/6
                    aux1=mean(pnt(n+1,cant+h1+1+(i-1)*6:cant+h1+1+i*6));
                    pnx(n+1,h1+i+1)=aux1;
                end
                for i=1:h2/6
                    aux2=mean(pnt(n+1,ns-(cant+h1+1+i*6):ns-(cant+h1+1+(i-1)*6)));
                    pnx(n+1,h1+i+1+h2/6)=aux2;
                end
                for i=1:h1
                    pnx(n+1,h1+2*h2/6+i+1)=pnt(n+1,ns-(cant+i));
                end
                pnx(n+1,m)=mean(pnt(n+1,ns-cant:ns));
            end
            
            hnt=0;
            for i=1:m
                hnx=0;
                for n=0:nmax-1
                    if pnx(n+1,i)~=0
                        hnx=hnx-pnx(n+1,i)*log(pnx(n+1,i));
                    end
                end
                hnt=hnt+hnx/m;
            end
            I=hn-hnt;
            ht=log(m);
            
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
            
            pns=zeros(nmax,6);
            
            for n=0:nmax-1
                pns(n+1,1)=mean(pnt(n+1,1:cant));
                pns(n+1,2)=mean(pnt(n+1,ns-cant:ns));
                pns(n+1,3)=mean(pnt(n+1,cant+1:cant+h1));
                pns(n+1,4)=mean(pnt(n+1,ns-(cant+h1):ns-cant));
                pns(n+1,5)=mean(pnt(n+1,cant+h1+1:cant+h1+h2));
                pns(n+1,6)=mean(pnt(n+1,cant+h1+h2+1:ns-cant-h1));
            end
            
            hnt=0;
            for i=1:m
                hns=0;
                for n=0:nmax-1
                    if pns(n+1,i)~=0
                        hns=hns-pns(n+1,i)*log(pns(n+1,i));
                    end
                end
                hnt=hnt+hns/m;
            end
            I=hn-hnt;
            ht=log(m);
    end
      
        
end