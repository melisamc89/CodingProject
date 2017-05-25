function [I0 I1]=InformationCalculusMontemurro(spikes,nt)

        [ntrials, nbins]=size(spikes);
        ns=floor(nbins);
        spk=zeros(1,1,ntrials,ns);
        for i=1:ntrials
            spk(1,:,i,:)=reshape(spikes(i,1:ns),1,[]);
        end

        [ncells,L,ntrials,ns]=size(spk);
        
        I0(1)=hr(spk,nt,0)-hrs(spk,nt,0);
        I1(1)=hr(spk,nt,2)-hrs(spk,nt,2);
        
        
        I0(2)=hr(spk,nt,0)-(hrs(spk,nt,0)+hrsind(spk,nt,0)-hrs_shuff(spk,nt,0));
        I1(2)=hr(spk,nt,2)-(hrs(spk,nt,2)+hrsind(spk,nt,2)-hrs_shuff(spk,nt,2));
        
        I0(3)=hrn(spk,nt,0)-hrsn(spk,nt,0);
        I1(3)=hrn(spk,nt,2)-hrsn(spk,nt,2);

        %I0=hr(spk,nt,0)-hrs(spk,nt,0);
        %I1=hr(spk,nt,4)-hrs(spk,nt,4);
        %Hn0=hrn(spk,nt,0);
        %Hnt0=hrsn(spk,nt,0);
        %Hn1=hrn(spk,nt,4);
        %Hnt1=hrsn(spk,nt,4);
        %I0=Hn0-Hnt0;
        %I1=Hn1-Hnt1;

end