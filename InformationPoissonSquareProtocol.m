%This function computes the information about one behavioural aspect of
%the rat with the square protocol in it's velocity.
%The variable type is the one which determinates the behavioural aspect.
% (0) Time information
% (1) Space Information
% (2) Directional Information
% (3) Speed Information
% (4) Space Directional Information
% (5) Speed Directional Information

function InformationPoissonSquareProtocol(directory,files,rat,areaname,area,list,freq)
    count=1;
    for index= 1:size(list)
        order = Wished_Register_Order(area,list(index,:));
        register=Rat_Register(files,area,rat,order);        
        [fil col]=size(register.sessions);
        for se=1:col
            protocol=register.sessions{se}.sess_type;
            if protocol == 'C'
                spikes_times= register.sessions{se}.spike_times;
                pos=Loading_Pos(register,se);
                tf = isfield(pos.data, 'log');
                if (tf==1 && length(spikes_times)>50)
                    [matrix1 matrix2]=SpikesMatrix(spikes_times,pos,freq);
                    [pn pnt pt]=PoissonProbability(matrix1,matrix2);
                    
                    for j=0:5
                        [I0 h0 h0t hj]=BehaviouralInformationSquareProtocol(pn,pnt,pt,pos,freq,j);
                        I(count,j+1)=I0;
                        hn(count,j+1)=h0;
                        hnt(count,j+1)=h0t;
                        h_var(count,j+1)=hj;
                    end
                    data(count,:)=area(index,:);
                    count=count+1;
                end
            end
        end
    end
    name=strcat(directory,areaname);
    save(name,'data','I','hn','hnt','h_var','-mat')

end