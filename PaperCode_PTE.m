clear all; close all;
%% replace localpath with correct path
localpath= 'path/to/files'
addpath('localpath/BrainNetViewer_20191031')
n_scouts=148;

for fband=1:6
Frequency_Band = fband;

badblocks_bsl={};
badblocks_pre={};


epochsDir = dir ('localpath/DOA/data/P18_Average_reference/*');
% Get a logical vector that tells which is a directory.
dirFlags = [epochsDir.isdir]
% Extract only those that are directories.
subFolders = epochsDir(dirFlags)
subFolders(ismember( {subFolders.name}, {'.', '..'})) = [];
epochs = struct2table(subFolders).name(4:end);
%

datapath='localpath/DOA/data/P18_Average_reference';
basename1= 'timefreq_connectn_pte';

%
ALLConnect_baseline = cell(1, 1);
for ep=1:length(epochs)
    %
    epoch=epochs{ep};
    %
    list1=dir([datapath,'/',epoch,'/',basename1,'*.mat']);
    listt=[list1];

    M=load(fullfile(listt(1).folder,listt(1).name));
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
    for w=1:length(listt)
        %
        M=load(fullfile(listt(w).folder,listt(w).name));
        if M.Time(1)<60
        if not(contains(M.Comment, "PTE_norm"))
         
%       Iterate over the first 10 epochs (10 x 6 seconds)
            TF=M.TF;
            %
            PTE=zeros(n_scouts,n_scouts);
            PTE_norm=zeros(n_scouts,n_scouts);
            PTE_fisher=zeros(n_scouts,n_scouts);
            Conn=zeros(n_scouts,n_scouts);
            %
            iBand = Frequency_Band;
            ik=1;
            for i=1:n_scouts
                for j=1:n_scouts
                    PTE(j,i)=TF(ik,1,iBand);
                    ik=ik+1;
                end
            end
            
            % /// uncomment this code to average PTE to main regions (e.g. scouts)
%             M_reduced = [];
%             A=PTE;
%             for i=1:size(uniq)
%                 ind = find(strcmp(uniq.(1)(i),regions));
%                 M_reduced(i,:) = mean(A(ind, :));
%             end
%             M_reduced2 = [];
%             for i=1:size(uniq)
%                 ind = find(strcmp(uniq.(1)(i),regions));
%                 M_reduced2(:, i)= mean(M_reduced(:,ind),2);
%             end
%             PTE = M_reduced2;
            % /// avere PTE to main regions
            
            tmp = triu(PTE) + tril(PTE)';
            dPTE = triu(PTE./tmp,1) + tril(PTE./tmp',-1);
            PTE_norm = dPTE - 0.5;

            diagM=zeros(148,148);
            PTE_norm = PTE_norm - diag(diag(PTE_norm)) + diag(diagM);
            
            % make matrix asymmetric
            PTE_norm(PTE_norm<0)=0;
            
            ALLConnect_baseline{ep, w} = PTE_norm;
            
        end
        end
    end
end


%% ////////////////
% Average connectivity 6s pre-Event
datapath='localpath/DOA/data/P18_Average_reference';
basename1= 'timefreq_connectn_pte';

PreEventConnect = cell(1, 1);

for ep=1:length(epochs)
    %
    epoch=epochs{ep};
    %
    list1=dir([datapath,'/',epoch,'/',basename1,'*.mat']);
    listt=[list1];

    M=load(fullfile(listt(1).folder,listt(1).name));
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
    for w=1:length(listt)
        %
        M=load(fullfile(listt(w).folder,listt(w).name));
        mm = round(M.Time(1));
        if not(contains(M.Comment, "PTE_norm"))
        if mm==174 | mm==180 | mm==186 | mm==192
            TF=M.TF;
            %
            PTE=zeros(n_scouts,n_scouts);
            PTE_norm=zeros(n_scouts,n_scouts);
            PTE_fisher=zeros(n_scouts,n_scouts);
            Conn=zeros(n_scouts,n_scouts);
            %
            iBand = Frequency_Band;
            ik=1;
            for i=1:n_scouts
                for j=1:n_scouts
                    PTE(j,i)=TF(ik,1,iBand);
                    ik=ik+1;
                end
            end
            %
            % /// uncomment this code to average PTE to main regions (e.g. scouts)
%             M_reduced = [];
%             A=PTE;
%             for i=1:size(uniq)
%                 ind = find(strcmp(uniq.(1)(i),regions));
%                 M_reduced(i,:) = mean(A(ind, :));
%             end
%             M_reduced2 = [];
%             for i=1:size(uniq)
%                 ind = find(strcmp(uniq.(1)(i),regions));
%                 M_reduced2(:, i)= mean(M_reduced(:,ind),2);
%             end
%             PTE = M_reduced2;
            % /// avere PTE to main regions
            
            tmp = triu(PTE) + tril(PTE)';
            dPTE = triu(PTE./tmp,1) + tril(PTE./tmp',-1);
            PTE_norm = dPTE - 0.5; % Center result around 0
            
            diagM=zeros(148,148);
            PTE_norm = PTE_norm - diag(diag(PTE_norm)) + diag(diagM);
            
            % make matrix asymmetric
            PTE_norm(PTE_norm<0)=0;
            
            PreEventConnect{ep, w} = PTE_norm;
            %
        end
        end
    end
end


%/////////////////

newMat = cell(1,2);
for i=1:20
    tempRow = PreEventConnect(i,:);
    newMat{i} = tempRow(~cellfun(@isempty, tempRow));
end
PreEventConnect2 = vertcat(newMat{:});

newMat2 = cell(1,2);
for i=1:20
    tempRow2 = ALLConnect_baseline(i,:);
    newMat2{i} = tempRow2(~cellfun(@isempty, tempRow2));
end
ALLConnect_baseline2 = vertcat(newMat2{:});


% ////////////////
datapath='localpath/EEG_All_Conn/';

% store data
my_field = strcat('bsl_',freqs_name{iBand});
All_Conn_Real2.(my_field) = ALLConnect_baseline2;

my_field2 = strcat('pre_',freqs_name{iBand});
All_Conn_Real2.(my_field2) = PreEventConnect2(:,1);

my_field3 = strcat('event0_',freqs_name{iBand});
All_Conn_Real2.(my_field3) = PreEventConnect2(:,2);

my_field4 = strcat('event1_',freqs_name{iBand});
All_Conn_Real2.(my_field4) = PreEventConnect2(:,3);

my_field5 = strcat('event2_',freqs_name{iBand});
All_Conn_Real2.(my_field5) = PreEventConnect2(:,4);
     
end

%%
%save('All_Conn.mat', 'All_Conn_Real2');


