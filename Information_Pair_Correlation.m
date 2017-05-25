
    counter=1;
    for index=[8:60,62:71]
        index;
        data=corr_and_info{index};
        [n m]=size(data.corr);
        for i=1:n
            for j=i+1:m
                corr_val(counter)=data.corr_2(i,j);
                info(counter)=data.info(i,j);
                if data.corr(i,j)>0.2
                    index
                    i
                    j
                end
                counter=counter+1;
            end
        end
    end
    
    
low=6;
high=12;
n_fs=4800;
load('3CA1','-mat');
areaname='CA1';
area=CA1_3;
list=area;
files=Data_Listing();
rat=3;
%load('CA1_complete','-mat');
freq=1;
    lfp=Days(list);
    for index=9:60
        data=corr_and_info{index};
        cells=SimultaneousRecordings(area,lfp(index,1:7));
        [spikes_times session_number]=SpikesAssembly(files,rat,area,cells,'C');
        [ncells nse]=size(session_number);
            matrix=[];
            for i=1:ncells
                order = Wished_Register_Order(area,cells(i,:));
                register=Rat_Register(files,area,rat,order);
                cell_spikes_times=spikes_times{i};
                [matrix1 matrix2]=SpikesMatrix(cell_spikes_times,register,session_number(i,:),freq);
                matrix(:,:,i)=[matrix1,matrix2];
            end
            for i=1:ncells
                for j=i+1:ncells
                    pn1n2=BidimentionalGaussianProbability2(matrix(:,:,i),matrix(:,:,j));
                    I=data.info(i,j);
                    c=data.corr(i,j);
                    figure(1)
                    imagesc(pn1n2)
                    n1=sum(sum(matrix(:,:,i)));
                    n2=sum(sum(matrix(:,:,j)));
                    title(strcat('I=',num2str(I),' _ Corr=',num2str(c),' _ n1=',int2str(n1),' _ n2=',int2str(n2)));
                    dir='C:\Users\meli__000\Desktop\Melisa\Doctorado\CellsCorrelation\pn1n2_2\';
                    nombre=strcat(dir,num2str(index),'_',int2str(i),'_',int2str(j));
                    saveas(1,nombre,'png')
                    close(1)
                end
            end
    end