%This function computes the information about one behavioural aspect of
%the rat with the square protocol in it's velocity.
%The variable type is the one which determinates the behavioural aspect.
% (0) Time information
% (1) Space Information
% (2) Directional Information
% (3) Speed Information
% (4) Space Directional Information
% (5) Speed Directional Information

function InformationSquareProtocol(directory,files,rat,areaname,area,list,freq,aspect)

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
                if (tf==1 && length(spikes_times)>1500)
                    [matrix1 matrix2]=SpikesMatrix(spikes_times,pos,42,freq);
                    
                    [newmatrix nt]=SquareProtocolMatrix(matrix1,matrix2,pos,aspect,freq);
                    [n m]=size(newmatrix);
                    x=reshape(newmatrix,n*m,1);
                    y=1:m;
                    for i=1:n-1
                         y=[y 1:m];
                    end
                    I0=Information2(x,y,max(y));
               
                    I(count,:)=I0;
                    data(count,:)=area(index,:);
                    count=count+1;
                end
            end
        end
    end
    switch aspect
        case 0
            name=strcat(directory,areaname);
        case 1
            name=strcat(directory,areaname,'_Space');
        case 2
            name=strcat(directory,areaname,'_Directional');
        case 3
            name=strcat(directory,areaname,'_Speed');
        case 4
            name=strcat(directory,areaname,'_SpaceDirectional');
        case 5
            name=strcat(directory,areaname,'_SpeedDirectional');
    end
    save(name,'data','I','-mat')

end