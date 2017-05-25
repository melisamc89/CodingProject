%This function computes the information between firing rate and phase for
%each cell.

load('CA1_complete','-mat');
low=6;
high=12;
n_fs=4800;
load('3CA1','-mat');
areaname='CA1';
area=CA1_3;
list=area;
files=Data_Listing();
rat=3;
freq=1;

    lfp=Days(list);
    for index=8:length(lfp)
        cells=SimultaneousRecordings(area,lfp(index,1:7));
        [spikes_times session_number]=SpikesAssembly(files,rat,area,cells,'C');
        [ncells nse]=size(session_number);
        index
        ncells
        if ncells
            matrix=[];
            for i=1:ncells
                i
                %recording_data{index}{i}.name=cells(i,:);
                order = Wished_Register_Order(area,cells(i,:));
                register=Rat_Register(files,area,rat,order);
                cell_spikes_times=spikes_times{i};
                if length(cell_spikes_times)
                    for j=1:100
                        phase_spike=Random_Phase(cell_spikes_times);
                        [matrix1 matrix2]=SpikesMatrix(cell_spikes_times,register,session_number(i,:),freq);
                        [pmatrix1 pmatrix2]=SpikesPhaseMatrix(phase_spike,cell_spikes_times,register,session_number(i,:),freq);                
                        matrix=[matrix1,matrix2];
                        [p_theta p_theta_n pn vec phase_mean kappa]=GaussianVonMisesDistribution2(matrix,pmatrix1,pmatrix2);
                        [I hy hyx hx]=Information(p_theta,pn,p_theta_n);
                        info(j)=I;
                    end
                    data{index}{i}.random.I_np=info;
                end
            end
        end
    end
    name=strcat(directory,'CA1','_firing_phase');
    save(name,'data','-mat')