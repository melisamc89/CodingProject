function [H0 H1]=EntropyNCalculusMontemurro(spikes,nt)

        [ntrials, nbins]=size(spikes);
        ns=floor(nbins);
        spk=zeros(1,1,ntrials,ns);
        for i=1:ntrials
            spk(1,:,i,:)=reshape(spikes(i,1:ns),1,[]);
        end

        [ncells,L,ntrials,ns]=size(spk);

        H0=hrn(spk,nt,0);
        H1=hrsn(spk,nt,0);
       

end