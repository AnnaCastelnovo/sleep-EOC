clc
clear all
close all
%%
addpath('.../BrainNetViewer_20191031')
n_scouts=148;
surface='.../destrieux/BrainMesh_ICBM152.nv';
node='.../destrieux/destrieux_corrected.node';
conf='.../destrieux/cfg4.mat';
baseline_names = {'averageBaseline_delta','averageBaseline_theta','averageBaseline_alpha','averageBaseline_sigma','averageBaseline_beta','averageBaseline_gamma'};
pre_names = {'averagePre_delta','averagePre_theta','averagePre_alpha','averagePre_sigma','averagePre_beta','averagePre_gamma'};
%%
%conf
configuration=load(conf);
EC=configuration.EC;
EC.edg.draw_threshold=0.015;
%EC.edg.draw_sparsity=0.08;
EC.edg.directed=1;
cm=colormap(hot(64));
cm=flipud(cm);
EC.edg.CM=cm;
EC.edg.CMt=cm;
EC.edg.CM_custom=cm;
EC.edg.color_map_low=0.00;
EC.edg.color_map_high=0.50;
EC.edg.size_ratio=10; %scale factor
EC.edg.size_abs=1; %true
EC.edg.size_value=2;
save(conf,'EC')
%%
epochsDir = dir ('.../doa_templates_connectivity/*');
% Get a logical vector that tells which is a directory.
dirFlags = [epochsDir.isdir]
% Extract only those that are directories.
subFolders = epochsDir(dirFlags)
subFolders(ismember( {subFolders.name}, {'.', '..'})) = [];
epochs = struct2table(subFolders).name;
%

datapath='.../doa_templates_connectivity';
basename1= 'timefreq_connectn_pte_200609';
basename2= 'timefreq_connectn_pte_200610';
basename3='timefreq_connectn_pte_200611'; % no interp

%
ALLConnect = zeros(148,148);
ALLConnect_baseline = cell(1, 1);
Frequency_Band = 1;
for ep=1:length(epochs)
    epoch=epochs{ep};
    %
    list1=dir([datapath,'/',epoch,'/',basename1,'*.mat']);
    list2=dir([datapath,'/',epoch,'/',basename2,'*.mat']);
    list3=dir([datapath,'/',epoch,'/',basename3,'*.mat']);
    list=[list1;list2;list3];
    %
    M=load(fullfile(list(1).folder,list(1).name));
    Freqs=M.Freqs;
    if size(Freqs,2)==1
        freqs_name=cell(size(Freqs,1),1);
        for i=1:size(Freqs,1)
            freqs_name{i}=[num2str(Freqs(i),'%08.4f'),'Hz'];
        end
    elseif size(Freqs,2)==3
        freqs_name=Freqs(:,1);
    end
    %
    num_freq=length(freqs_name);
    %
    for w=1:length(list)
        %M.Time(1)
        %
        M=load(fullfile(list(w).folder,list(w).name));
        if M.Time(1)<60
%       Iterate over the first 10 epochs (10 x 6 seconds)
            TF=M.TF;
            %
            PTE=zeros(n_scouts,n_scouts);
            PTE_norm=zeros(n_scouts,n_scouts);
            PTE_fisher=zeros(n_scouts,n_scouts);
            Conn=zeros(n_scouts,n_scouts);
            %
            iBand = Frequency_Band; % select delta
            ik=1;
            for i=1:n_scouts
                for j=1:n_scouts
                    PTE(j,i)=TF(ik,1,iBand); 
                    ik=ik+1;
                end
            end
            %
            tmp = triu(PTE) + tril(PTE)';
            dPTE = triu(PTE./tmp,1) + tril(PTE./tmp',-1);
            PTE_norm1 = dPTE - 0.5; % Center result around 0

            diagM=zeros(n_scouts,n_scouts);
            PTE_norm2 = PTE_norm1 - diag(diag(PTE_norm1)) + diag(diagM);
            
            % make matrix asymmetric
            PTE_norm2(PTE_norm2<0)=0;
            
            %PTE_fisher = PTE_fisher - diag(diag(PTE_fisher)) + diag(diagM);
            %PTE_fisher = real(atanh(PTE_norm)); % Fisher transformation
            ALLConnect_baseline{ep, w} = PTE_norm; %rescale(PTE,0,1);
            
            %
            for i=1:n_scouts
                for j=1:i
                    val=PTE_norm(i,j);
                    if val>0
                        Conn(j,i)=0;
                    else
                        Conn(j,i)=abs(val);
                        Conn(i,j)=0;
                    end
                end
            end
            
        end
    end
end

%///////
% Baseline print edge File
B=cat(3,ALLConnect_baseline{:});
averageBaseline = mean(B,3);

datapath='.../surrogate_PTE2/';
s=fullfile(list(w).folder,list(w).name);
s=split(s,'.mat');
%
edge=[datapath,'PTE_edge','/','Dir_baseline','_',freqs_name{iBand},'.edge'];
fid=fopen(edge,'w');
for i=1:n_scouts
    for j=1:n_scouts
        fprintf(fid,'%f ',averageBaseline(i,j));
    end
    fprintf(fid,'\n');
end
fclose(fid);

% calculate threshold with bonferoni
threshold_bsl = prctile(averageBaseline_surr(:),95-(0.01/1700))

% ////////////////
% Average connectivity 6s pre-Event
datapath='.../doa_templates_connectivity';
basename1= 'timefreq_connectn_pte_200609';
basename2= 'timefreq_connectn_pte_200610';
basename3='timefreq_connectn_pte_200611'; % no interp

PreEventConnect = zeros(148,148);
PreEventConnect_fisher = cell(1, 1);

for ep=1:length(epochs)
    %
    epoch=epochs{ep};
    %
    list1=dir([datapath,'/',epoch,'/',basename1,'*.mat']);
    list2=dir([datapath,'/',epoch,'/',basename2,'*.mat']);
    list3=dir([datapath,'/',epoch,'/',basename3,'*.mat']);
    list=[list1;list2;list3];
    %

    %
    ep
    M=load(fullfile(list(1).folder,list(1).name));
    Freqs=M.Freqs;
    if size(Freqs,2)==1
        freqs_name=cell(size(Freqs,1),1);
        for i=1:size(Freqs,1)
            freqs_name{i}=[num2str(Freqs(i),'%08.4f'),'Hz'];
        end
    elseif size(Freqs,2)==3
        freqs_name=Freqs(:,1);
    end
    %
    num_freq=length(freqs_name);
    %
    for w=1:length(list)
        %
        M=load(fullfile(list(w).folder,list(w).name));
        if M.Time(1)==174 | M.Time(1)==180 | M.Time(1)==186 | M.Time(1)==192 | M.Time(1)==198
            %M.Time(1)
%       Iterate over the first 10 epochs (10 x 6 seconds)
            TF=M.TF;
            %
            PTE=zeros(n_scouts,n_scouts);
            PTE_norm=zeros(n_scouts,n_scouts);
            PTE_fisher=zeros(n_scouts,n_scouts);
            Conn=zeros(n_scouts,n_scouts);
            %
            iBand = Frequency_Band; % select delta
            ik=1;
            for i=1:n_scouts
                for j=1:n_scouts
                    PTE(j,i)=TF(ik,1,iBand); %(j,i)
                    ik=ik+1;
                end
            end
            %
            tmp = triu(PTE) + tril(PTE)';
            dPTE = triu(PTE./tmp,1) + tril(PTE./tmp',-1);
            PTE_norm = dPTE - 0.5; % Center result around 0
            
            diagM=zeros(n_scouts,n_scouts);
            PTE_norm = PTE_norm - diag(diag(PTE_norm)) + diag(diagM);
            
            % make matrix asymmetric
            PTE_norm(PTE_norm<0)=0;
            
            PreEventConnect_fisher{ep, w} = PTE_norm; %rescale(PTE,0,1);
            %
            for i=1:n_scouts
                for j=1:i
                    val=PTE_norm(i,j);
                    if val>0
                        Conn(j,i)=0;
                    else
                        Conn(j,i)=abs(val);
                        Conn(i,j)=0;
                    end
                end
            end
        end
    end
end


%/////////////////

newMat = cell(1,2);
for i=1:17
    tempRow = PreEventConnect_fisher(i,:);
    newMat{i} = tempRow(~cellfun(@isempty, tempRow));
end
PreEventConnect_fisher2 = vertcat(newMat{:});

% ////////////////
datapath='.../surrogate_PTE2/';
%%

for time=4:4 %(1:5)
    PreEventConnect_fisher3 = PreEventConnect_fisher2(:,time);
    B=cat(3,PreEventConnect_fisher3{:});
    %averageBaselinePre = mean(abs(B),3);
    averageBaselinePre = mean(B,3);

    s=fullfile(list(w).folder,list(w).name);
    s=split(s,'.mat');
    %
    edge=[datapath,'PTE_edge','/','Dir_preEvent_',num2str(time),'_',freqs_name{iBand},'.edge'];   
    fid=fopen(edge,'w');
    for i=1:n_scouts
        for j=1:n_scouts
            fprintf(fid,'%f ',averageBaselinePre(i,j));
        end
        fprintf(fid,'\n');
    end
    fclose(fid);
    threshold_pre = prctile(averageBaselinePre_surr(:),100-(0.01/170));

my_field = strcat('bls',freqs_name{iBand});
BSLmeanNetworks.(my_field) = averageBaseline;
my_field2 = strcat('bls',freqs_name{iBand});
PREmeanNetworks.(my_field2) = averageBaselinePre;
    
end
    
