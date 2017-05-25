%This function computes the information about one behavioural aspect of
%the rat with the square protocol in it's velocity.
%The variable type is the one which determinates the behavioural aspect.
%The mutal information between eeg phase with any behaviural aspect is
%computed.
% (0) Time information
% (1) Space Information
% (2) Directional Information
% (3) Speed Information
% (4) Space Directional Information
% (5) Speed Directional Information

function InformationVonMisesPhaseSquareProtocol(directory,files,rat,areaname,area,list,freq,n_fs,low,high,num)

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
                    [eeg_vector,fs]=ReadEEG(register,se);
                    eeg_vector=resample(eeg_vector,1,floor(fs/n_fs));
                    n_fs=fs/floor(fs/n_fs);
                    [phase energy]=Hilbert_phases_energy(eeg_vector,n_fs,low,high);
                    phase_spike=Phase_Calculator(spikes_times,phase,n_fs);
                    [matrix1 matrix2]=SpikesPhaseMatrix(phase_spike,spikes_times,pos,42,freq);
                    [p_theta p_thetat pt]=VonMisesProbability(matrix1,matrix2,num);
                    
                    for j=0:5
                        [I0 h0 h0t hj]=BehaviouralInformationSquareProtocol(p_theta,p_thetat,pt,pos,freq,j);
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