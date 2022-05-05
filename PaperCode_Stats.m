restoredefaultpath
%%
addpath .../fieldtrip-20210616
ft_defaults

%% PLot reduced brainnet viewer Diffs
addpath('.../BrainNetViewer_20191031')
surface='.../destrieux/BrainMesh_ICBM152.nv';
conf='.../destrieux/cfg_reduced1.mat';
node='.../destrieux/reducedNodesMean.node';

fn = fieldnames(All_Conn_Real2);
freqs_name = {'delta', 'theta', 'alpha', 'sigma', 'beta', 'gamma'};
mean_pre = [];

name1 = 'pre_';
name2 = 'event2_';

for iBand=1:6        
        field1 = strcat(name1,freqs_name{iBand});
        field2 = strcat(name2,freqs_name{iBand});
        bsl = All_Conn_Real2.(field1(:));
        bsl = bsl(:).';
        dev = All_Conn_Real2.(field2(:));
        dev = dev(:).';
        % calculate threshold
        dist_bsl_dev = vertcat(bsl{:}, dev{:});
        dist_bsl_dev = dist_bsl_dev(dist_bsl_dev~=0);
        threshold = prctile(dist_bsl_dev(:),85);
        p_mask_perm = zeros(148,148, "double");
  
        
        %calculate mean of indv
        bsl_t = cat(3,bsl{:});
        bsl_mean = mean(bsl_t,3);
        

        dev_t = cat(3,dev{:});
        dev_mean = mean(dev_t,3);
        
        %diff of means
        diff = dev_mean - bsl_mean;
        
        %!change file name: lookup significance
        signif_band = adj_p_values_pre_event2_t_f.(freqs_name{iBand});
        
        %apply significance threshold
        diff(signif_band==0) = 0;
        
        % keep max of directed weight pair
        j=1;
        for i=1:14
            
            j
            for j=1:j
                if abs(diff(i,j))>abs(diff(j,i))
                    diff_max(i,j) = diff(i,j);
                    diff_max(j,i) = 0;
                else
                    diff_max(i,j) = 0;
                    diff_max(j,i) = diff(j,i);
                end
            end
            j=j+1;
        end
        
       % store edge file
       edge=['.../EEG_All_Conn/reduced_edges','_',strcat(name1,name2),freqs_name{iBand},'.edge'];
       fid=fopen(edge,'w');
       %edge
       for i=1:14
           for j=1:14
               fprintf(fid,'%f ',diff_max(i,j));
           end
           fprintf(fid,'\n');
       end
       fclose(fid);
       
       %plot brainnet
       f_img=['.../EEG_All_Conn/imgs/Connectivity_reduced/',strcat(name1,name2),freqs_name{iBand},'.jpg'];
       BrainNet_MapCfg(surface,node,edge,conf,f_img)
        
end

%%

addpath('.../connectivity/BrainNetViewer_20191031')
surface='.../dsk/connectivity/destrieux/BrainMesh_ICBM152.nv';
conf='.../dsk/connectivity/destrieux/cfg_reduced1.mat';
node='.../dsk/connectivity/destrieux/reducedNodesMean.node';
edge='.../EEG_All_Conn/reduced_edges_bsl_event1_delta.edge';

BrainNet_MapCfg(surface,node, edge,conf)

edge_files_strength = {
                'reduced_edges_bsl_event1_delta.edge'
                };
            
for i=1:length(edge_files_strength)
    edge = edge_files_strength{i};
    f_img=['.../imgs/connectivity_strength2/','strength_', node_files_strength{i}(20:end),'.jpg'];
    BrainNet_MapCfg(surface,node,edge,conf,f_img);
end



%% destrieux mean reduced
A = [];
B = [];
C = [];
D = [];
for i=1:14
    ii = uniq.Region(i);
    
    ispresent = strcmp(regions.Region, ii{1,1});
    iidx=find(ispresent).';
    group = destrieuxcorrected(iidx, [1,2,3]);
    means = mean(table2array(group))
    %means
    A(i) = means(1,1);
    B(i) = means(1,2);
    C(i) = means(1,3);
    D(i) = 1;
    E(i) = 1;
end

T = table(A',B',C',D',D',uniq.Region);

%save
writetable(T, 'reducedNodesMean.node','Delimiter','\t', 'FileType', 'text', 'WriteVariableNames', 0);

%% transform brainstorm TF matrix strength!
all = [adj_p_values_pre_event1_t.delta(:); 
    adj_p_values_pre_event1_t.theta(:);
    adj_p_values_pre_event1_t.alpha(:);
    adj_p_values_pre_event1_t.sigma(:);
    adj_p_values_pre_event1_t.beta(:);
    adj_p_values_pre_event1_t.gamma(:)
    ];

pmap_n = reshape(all, [21904, 1, 6]);
pmap_n(pmap_n==0) = 1;  

pmap_n_pre_event1 = pmap_n;
save ("pmap_n_pre_event1.mat", "pmap_n_pre_event1");

% open brainstorm conn and replace pmap

pmap = pmap_n_pre_event1;

%% Reduce connectivity matrices to main regions
M_reduced = [];
A=All_Conn_Real2.bsl_delta{1,1};
for i=1:size(uniq)
    ind = find(strcmp(uniq.(1)(i),regin.Region));
    M_reduced(i,:) = mean(A(ind, :));
end
M_reduced2 = [];
for i=1:size(uniq)
    ind = find(strcmp(uniq.(1)(i),regin.Region));
    M_reduced2(:, i)= mean(M_reduced(:,ind),2);
end

%% Diff of means undirected 
fn = fieldnames(All_Conn_Real2);
freqs_name = {'delta', 'theta', 'alpha', 'sigma', 'beta', 'gamma'};
mean_pre = [];

name1 = 'bsl_';
name2 = 'event0_';

for iBand=1:6        
        field1 = strcat(name1,freqs_name{iBand});
        field2 = strcat(name2,freqs_name{iBand});
        bsl = All_Conn_Real2.(field1(:));
        bsl = bsl(:).';
        dev = All_Conn_Real2.(field2(:));
        dev = dev(:).';
        % calculate threshold
        dist_bsl_dev = vertcat(bsl{:}, dev{:});
        dist_bsl_dev = dist_bsl_dev(dist_bsl_dev~=0);
        threshold = prctile(dist_bsl_dev(:),99.9);
        p_mask_perm = zeros(148,148, "double");
        
        %copy triu to tril
        for i = 1:length(bsl)
            a = bsl{i};
            %apply threshold
            a(a<threshold)=0;
            
            al = triu(a)' + tril(a);
            bsl{i} = al;
        end
        for i = 1:length(dev)
            a = dev{i};
            %apply threshold
            a(a<threshold)=0;
            
            al = triu(a)' + tril(a);
            dev{i} = al;
        end
        
        %calculate mean of indv
        bsl_t = cat(3,bsl{:});
        bsl_mean = mean(bsl_t,3);
        
        dev_t = cat(3,dev{:});
        dev_mean = mean(dev_t,3);
        
        %diff of means
        diff = dev_mean - bsl_mean;
        
        %lookup significance
        %signif = ['adj_p_values_',name1, name2, 't_undir'];
        %signif = 'adj_p_values_bsl_pre_t_undir';
        signif_band = adj_p_values_bsl_event1_t_undir.(freqs_name{iBand});
        
        %apply significance threshold
        diff(signif_band==0) = 0;
        diff_sym = tril(diff)' + tril(diff);
        
        %half size
        mask = tril(true(size(diff_sym)));
        %// Use mask to select lower triangular elements from input array
        diff_sym = diff_sym(mask);
        mean_pre = [mean_pre, diff_sym(:)];

end

%mean_pre = reshape(mean_pre, [21904,1,6]);
% -> half size
mean_pre = reshape(mean_pre, [11026,1,6]);

%% Threshold between all segments I

fn = fieldnames(All_Conn_Real2);
freqs_name = {'delta', 'theta', 'alpha', 'sigma', 'beta', 'gamma'};
threshold_overAll = [];

name1 = 'pre_';
name2 = 'event1_';
name3 = 'event2_';
name4 = 'bsl_';
name5 = 'event0_';

for iBand=1:6        
        field1 = strcat(name1,freqs_name{iBand});
        field2 = strcat(name2,freqs_name{iBand});
        field3 = strcat(name3,freqs_name{iBand});
        field4 = strcat(name4,freqs_name{iBand});
        field5 = strcat(name5,freqs_name{iBand});
        bsl = All_Conn_Real2.(field1(:));
        bsl = bsl(:).';
        dev = All_Conn_Real2.(field2(:));
        dev = dev(:).';
        dev2 = All_Conn_Real2.(field3(:));
        dev2 = dev2(:).';
        dev3 = All_Conn_Real2.(field4(:));
        dev3 = dev3(:).';
        dev4 = All_Conn_Real2.(field5(:));
        dev4 = dev4(:).';
        % calculate threshold
        dist_bsl_dev1 = vertcat(bsl{:}, dev{:});
        dist_bsl_dev2 = vertcat(dev2{:}, dev3{:});
        dist_bsl_dev0 = vertcat(dist_bsl_dev1, dist_bsl_dev2);
        dist_bsl_dev = vertcat(dist_bsl_dev0, dev4{:});
        
        dist_bsl_dev = dist_bsl_dev(dist_bsl_dev~=0);
        threshold = prctile(dist_bsl_dev(:),85);
        %threshold
        threshold_overAll = [threshold_overAll, threshold];
        %end
end

%% mean of all networks with threshold over all II 
% change names of event!

fn = fieldnames(All_Conn_Real2);
mean_event0 = [];

for k=1:numel(fn)
        % manually change names of event !!!
        if contains(fn{k}, 'event0')
        
        m = All_Conn_Real2.(fn{k});
        m = m(:);
        for i=1:length(m)
            
            mi = m{i};
            
            %threshold = prctile(mi(:),85);
            
            freqs_name = {'delta', 'theta', 'alpha', 'sigma', 'beta', 'gamma'};
            t_i = find(contains(freqs_name,extractAfter(fn{k},"_")));
            threshold = threshold_overAll(t_i);
            threshold
            mi(mi<threshold)=0;
            mii = tril(mi) - triu(mi).';
            
            m{i} = mii;
        end
        B = cat(3,m{:});
        %apply threshold
        
        average = mean(B,3);
        fn{k}
        mean_event0 = [mean_event0, average(:)];
        
        end
end

%mean_pre = reshape(mean_pre, [21904,1,6]);

%% Threshold of means all segments III

threshold = [];

for i=1:6
    dist_bsl_dev = mean_event0(:,i);
    dist_bsl_dev = dist_bsl_dev(dist_bsl_dev~=0);
    thresholdi = prctile(dist_bsl_dev(:),85);
    threshold(i) = thresholdi;
    
    e0 = mean_event0(:,i);
    e0(abs(e0) < thresholdi) = 0;
    mean_event0(:,i) = e0;
    
%     e1 = mean_event1(:,i);
%     e1(abs(e1) < thresholdi) = 0;
%     mean_event1(:,i) = e1;
%     
%     e2 = mean_event2(:,i);
%     e2(abs(e2) < thresholdi) = 0;
%     mean_event2(:,i) = e2;
%     
%     p = mean_pre(:,i);
%     p(abs(p) < thresholdi) = 0;
%     mean_pre(:,i) = p;
end

mean_event0 = reshape(mean_event0, [21904,1,6]);
% mean_event1 = reshape(mean_event1, [21904,1,6]);
% mean_event2 = reshape(mean_event2, [21904,1,6]);
% mean_pre = reshape(mean_pre, [21904,1,6]);
% mean_bsl = reshape(mean_bsl, [21904,1,6]);

%mean_event0 = reshape(mean_event0, [196,1,6]);
% mean_event1 = reshape(mean_event1, [196,1,6]);
% mean_event2 = reshape(mean_event2, [196,1,6]);
% mean_pre = reshape(mean_pre, [196,1,6]);
% mean_bsl = reshape(mean_bsl, [196,1,6]);

%% add significance mask

fn = fieldnames(adj_p_values_pre_event2_t);
pre_event2 = [];
for k=1:numel(fn)
        if 	size(adj_p_values_pre_event2_t.(fn{k})) == 148, 148

        m = adj_p_values_pre_event2_t.(fn{k});

        pre_event2 = [pre_event2, m(:)];
        end
end

pre_event2 = reshape(pre_event2, [21904,1,6]);
TF(pre_event2==0)=0;

%% change brainstorm TF columns to conn Matrix
R = bst_memory("GetConnectMatrix", dif_pte_example);
size(R)
%% 
% Edit node file for BrainNet viewer strength!
freqs_name = {'delta', 'theta', 'alpha', 'sigma', 'beta', 'gamma'};

name2 = Strength_diff_pre_event0;

name1 = Strength_adj_Stats_pre_event0_t;

name_eps = 'Strength_adj_Stats_pre_event0_T_';

for iBand=1:6
     field1 = strcat('Nsignificant',freqs_name{iBand});
%     field2 = strcat(field1,name2);
%     stats_file = strcat(field2,freqs_name{iBand});
    if name1.(field1) > 0
    
        % edit node size -> diff in node strenght
        sig = name1.(freqs_name{iBand});

        % minus since mean(Bas)-mean(Dev)
        diff_file = name2.(freqs_name{iBand});
        diff_file_abs = abs(diff_file);
        diff_file_abs(sig==0) = 0;

        %color -> sign of diff in node strength
        diff_dir = sign( - diff_file) + 2;
        diff_dir(sig==0) = 0;

        destrieuxcorrected.('VarName4')(1:148) = diff_dir;

        % edit node size -> difference
        destrieuxcorrected.('VarName5')(1:148) = diff_file_abs;

        %save
        n = strcat(name_eps, freqs_name{iBand});
        save_name = strcat(n, '.node');
        writetable(destrieuxcorrected, save_name,'Delimiter','\t', 'FileType', 'text', 'WriteVariableNames', 0);

        save_name
    end
    
end


%% 
% Edit node file for BrainNet viewer strength DIR !!!

freqs_name = {'delta', 'theta', 'alpha', 'sigma', 'beta', 'gamma'};

name2 = Strength_diff_pre_event0;

name1 = Strength_adj_Stats_pre_event0_T_net;

name_file = 'Strength_adj_Stats_pre_event0_T_net_';

for iBand=1:6
     field1 = strcat('Nsignificant',freqs_name{iBand});
%     field2 = strcat(field1,name2);
%     stats_file = strcat(field2,freqs_name{iBand});
    if name1.(field1) > 0
    
        % edit node size -> diff in node strenght
        sig = name1.(freqs_name{iBand});

        % minus since [strength_diff, mean(Dev) - mean(Bas)];
        diff_file = name2.(freqs_name{iBand});
        diff_file_abs = abs(diff_file);
        diff_file_abs(sig==0) = 0;

        %color -> sign of diff in node strength
        diff_dir = sign( diff_file) + 2;
        diff_dir(sig==0) = 0;

        destrieuxcorrected.('VarName4')(1:148) = diff_dir;

        % edit node size -> difference
        destrieuxcorrected.('VarName5')(1:148) = diff_file_abs;

        %save
        n = strcat(name_file, freqs_name{iBand});
        save_name = strcat(n, '.node');
        writetable(destrieuxcorrected, save_name,'Delimiter','\t', 'FileType', 'text', 'WriteVariableNames', 0);

        save_name
    end
    
end

%% Print images BrainNet
addpath('.../connectivity/BrainNetViewer_20191031')
n_scouts=148;
surface='.../destrieux/BrainMesh_ICBM152.nv';
%%node='/Users/julian/Desktop/connectivity/destrieux/destrieux_corrected.node';
conf='.../destrieux/cfg5.mat';

node_files_strength = {
                'Strength_adj_Stats_bsl_event2_T_beta.node',
                ...
                };

node_files_meanPTE = {
                'Strength_adj_Stats_pre_event2_T_net_sigma.node',
                ...
                };

for i=1:length(node_files_strength)
    node = node_files_strength{i};
    f_img=['/.../connectivity_strength2_w0/','strength_', node_files_strength{i}(20:end),'.jpg'];
    BrainNet_MapCfg(surface,node,conf,f_img);
end

for i=1:length(node_files_meanPTE)
    node = node_files_meanPTE{i};
    f_img=['.../connectivity_strength2_w0/','meanPTE_', node_files_meanPTE{i}(20:end),'.jpg'];
    BrainNet_MapCfg(surface,node,conf,f_img);
end

%%
surface='...y/destrieux/BrainMesh_ICBM152.nv';
node='.../destrieux_corrected_bsl_event1_strenght.node';

%% combined thresholds
freqs_name = {'delta', 'theta', 'alpha', 'sigma', 'beta', 'gamma'};

for iBand=1:6
    field1 = strcat('bsl_',freqs_name{iBand});
    field2 = strcat('pre_',freqs_name{iBand});
    field3 = strcat('event1_',freqs_name{iBand});
    field4 = strcat('event2_',freqs_name{iBand});
    bsl = All_Conn_Real.(field1(:));
    pre = All_Conn_Real.(field2(:));
    event1 = All_Conn_Real.(field3(:));
    event2 = All_Conn_Real.(field4(:));

    dist_bsl_pre = vertcat(bsl{:}, pre{:});
    dist_bsl_event1 = vertcat(bsl{:}, event1{:});
    dist_bsl_event2 = vertcat(bsl{:}, event2{:});
    dist_pre_event1 = vertcat(pre{:}, event1{:});
    dist_pre_event2 = vertcat(pre{:}, event2{:});
    
    dist_bsl_pre = dist_bsl_pre(dist_bsl_pre~=0);
    dist_bsl_event1 = dist_bsl_event1(dist_bsl_event1~=0);
    dist_bsl_event2 = dist_bsl_event2(dist_bsl_event2~=0);
    dist_pre_event1 = dist_pre_event1(dist_pre_event1~=0);
    dist_pre_event2 = dist_pre_event2(dist_pre_event2~=0);

    threshold_bsl_pre.(freqs_name{iBand}) = prctile(dist_bsl_pre(:),85);
    threshold_bsl_event1.(freqs_name{iBand}) = prctile(dist_bsl_event1(:),85);
    threshold_bsl_event2.(freqs_name{iBand}) = prctile(dist_bsl_event2(:),85);
    threshold_pre_event1.(freqs_name{iBand}) = prctile(dist_pre_event1(:),85);
    threshold_pre_event2.(freqs_name{iBand}) = prctile(dist_pre_event2(:),85);
end

save 'threshold_bsl_pre.mat' 'threshold_bsl_pre';
save 'threshold_bsl_event1.mat' 'threshold_bsl_event1';
save 'threshold_bsl_event2.mat' 'threshold_bsl_event2';
save 'threshold_pre_event1.mat' 'threshold_pre_event1';
save 'threshold_pre_event2.mat' 'threshold_pre_event2';
%% /////////////////////
% Permutation test function
%  /////////////////////
% iterate over each scout, delete subthreshold & create stats mask

%trhesholded
[adj_p_values_bsl_pre_t, p_values_bsl_pre_t] = perm_stat_conn_dir_thresh(All_Conn_Real2, "bsl_", "pre_");
[adj_p_values_bsl_event1_t, p_values_bsl_event1_t] = perm_stat_conn_dir_thresh(All_Conn_Real, "bsl_", "event1_");
[adj_p_values_bsl_event2_t, p_values_bsl_event2_t] = perm_stat_conn_dir_thresh(All_Conn_Real, "bsl_", "event2_");
[adj_p_values_pre_event1_t, p_values_pre_event1_t] = perm_stat_conn_dir_thresh(All_Conn_Real, "pre_", "event1_");
[adj_p_values_pre_event2_t, p_values_pre_event2_t] = perm_stat_conn_dir_thresh(All_Conn_Real, "pre_", "event2_");


%% Dir trhesholded with FULL matrix
adj_p_values_bsl_pre_t_f = perm_stat_conn_dir_full_thresh(All_Conn_Real2, "bsl_", "pre_");
adj_p_values_bsl_event1_t_f = perm_stat_conn_dir_full_thresh(All_Conn_Real2, "bsl_", "event1_");
adj_p_values_bsl_event2_t_f = perm_stat_conn_dir_full_thresh(All_Conn_Real2, "bsl_", "event2_");
adj_p_values_pre_event1_t_f = perm_stat_conn_dir_full_thresh(All_Conn_Real2, "pre_", "event1_");
adj_p_values_pre_event2_t_f = perm_stat_conn_dir_full_thresh(All_Conn_Real2, "pre_", "event2_");

%% Dir trhesholded with FULL matrix | event 0
adj_p_values_bsl_event0_t_f = perm_stat_conn_dir_full_thresh(All_Conn_Real2, "bsl_", "event0_");
adj_p_values_pre_event0_t_f = perm_stat_conn_dir_full_thresh(All_Conn_Real2, "pre_", "event0_");

%% no threshold
[adj_p_values_bsl_pre, p_values_bsl_pre] = perm_stat_conn_dir(All_Conn_Real, "bsl_", "pre_");
[adj_p_values_bsl_event1, p_values_bsl_event1] = perm_stat_conn_dir(All_Conn_Real, "bsl_", "event1_");
[adj_p_values_bsl_event2,p_values_bsl_event2] = perm_stat_conn_dir(All_Conn_Real, "bsl_", "event2_");
[adj_p_values_pre_event1, p_values_pre_event1] = perm_stat_conn_dir(All_Conn_Real, "pre_", "event1_");
[adj_p_values_pre_event2, p_values_pre_event2] = perm_stat_conn_dir(All_Conn_Real, "pre_", "event2_");

%% undirected
% no threshold undir
[adj_p_values_bsl_pre_undir, p_values_bsl_pre_undir] = perm_stat_conn_undir(All_Conn_Real, "bsl_", "pre_");
[adj_p_values_bsl_event1_undir, p_values_bsl_event1_undir] = perm_stat_conn_undir(All_Conn_Real, "bsl_", "event1_");
[adj_p_values_bsl_event2_undir,p_values_bsl_event2_undir] = perm_stat_conn_undir(All_Conn_Real, "bsl_", "event2_");
[adj_p_values_pre_event1_undir, p_values_pre_event1_undir] = perm_stat_conn_undir(All_Conn_Real, "pre_", "event1_");
[adj_p_values_pre_event2_undir, p_values_pre_event2_undir] = perm_stat_conn_undir(All_Conn_Real, "pre_", "event2_");

%% undirected additional event0
[adj_p_values_pre_event0_undir, p_values_pre_event0_undir] = perm_stat_conn_undir(All_Conn_Real2, "pre_", "event0_");
[adj_p_values_bsl_event0_undir, p_values_bsl_event0_undir] = perm_stat_conn_undir(All_Conn_Real2, "bsl_", "event0_");

%%
%thresholded
[adj_p_values_bsl_pre_t_undir, p_values_bsl_pre_t_undir] = perm_stat_conn_undir_thresh(All_Conn_Real2, "bsl_", "pre_");
[adj_p_values_bsl_event1_t_undir, p_values_bsl_event1_t_undir] = perm_stat_conn_undir_thresh(All_Conn_Real2, "bsl_", "event1_");
[adj_p_values_bsl_event2_t_undir, p_values_bsl_event2_t_undir] = perm_stat_conn_undir_thresh(All_Conn_Real2, "bsl_", "event2_");
[adj_p_values_pre_event1_t_undir, p_values_pre_event1_t_undir] = perm_stat_conn_undir_thresh(All_Conn_Real2, "pre_", "event1_");
[adj_p_values_pre_event2_t_undir, p_values_pre_event2_t_undir] = perm_stat_conn_undir_thresh(All_Conn_Real2, "pre_", "event2_");

%% Strength (undirected)
[Strength_adj_Stats_bsl_event0_T, Strength_diff_bsl_event0] = perm_stat_strength_undir(All_Conn_Real2, "bsl_", "event0_");
[Strength_adj_Stats_pre_event0_t, Strength_diff_pre_event0] = perm_stat_strength_undir(All_Conn_Real2, "pre_", "event0_");

[Strength_adj_Stats_bsl_pre_T, Strength_diff_bsl_pre] = perm_stat_strength_undir(All_Conn_Real2, "bsl_", "pre_");
[Strength_adj_Stats_bsl_event1_T, Strength_diff_bsl_event1] = perm_stat_strength_undir(All_Conn_Real2, "bsl_", "event1_");
[Strength_adj_Stats_bsl_event2_T, Strength_diff_bsl_event2] = perm_stat_strength_undir(All_Conn_Real2, "bsl_", "event2_");
[Strength_adj_Stats_pre_event1_t, Strength_diff_pre_event1] = perm_stat_strength_undir(All_Conn_Real2, "pre_", "event1_");
[Strength_adj_Stats_pre_event2_T, Strength_diff_pre_event2] = perm_stat_strength_undir(All_Conn_Real2, "pre_", "event2_");

%% Strength net flow dir
[Strength_adj_Stats_bsl_event0_T_net, Strength_diff_bsl_event0] = perm_stat_strength(All_Conn_Real2, "bsl_", "event0_");
[Strength_adj_Stats_pre_event0_T_net, Strength_diff_pre_event0] = perm_stat_strength(All_Conn_Real2, "pre_", "event0_");

[Strength_adj_Stats_bsl_pre_T_net, Strength_diff_bsl_pre] = perm_stat_strength(All_Conn_Real2, "bsl_", "pre_");
[Strength_adj_Stats_bsl_event1_T_net, Strength_diff_bsl_event1] = perm_stat_strength(All_Conn_Real2, "bsl_", "event1_");
[Strength_adj_Stats_bsl_event2_T_net, Strength_diff_bsl_event2] = perm_stat_strength(All_Conn_Real2, "bsl_", "event2_");
[Strength_adj_Stats_pre_event1_T_net, Strength_diff_pre_event1] = perm_stat_strength(All_Conn_Real2, "pre_", "event1_");
[Strength_adj_Stats_pre_event2_T_net, Strength_diff_pre_event2] = perm_stat_strength(All_Conn_Real2, "pre_", "event2_");




%%
freqs_name = {'delta', 'theta', 'alpha', 'sigma', 'beta', 'gamma'};



for iBand=1:6
    F=Strength_adj_Stats_pre_event2_T.(freqs_name{iBand});
%     F(F>0.05)=0;
    field = strcat("Nsignificant",freqs_name{iBand});
    
    %adj_p_values_pre_event2.(freqs_name{iBand}) = F;
    Strength_adj_Stats_pre_event2_T.(field) = nnz(F)
end

%% mean TPE
for iBand=1:6
    iBand
    field1 = strcat('bsl_p_',freqs_name{iBand});
    field2 = strcat('pre_',freqs_name{iBand});
    bsl = All_Conn_Real_subthresh.(field1(:));
    pre = All_Conn_Real_subthresh.(field2(:));
    
    strength_bsl_p = [];
    for i=1:170
        try
            [is,os,str] = group_mean_TPE(bsl{i}, labels1);
        catch
        end
        strength_bsl_p = [strength_bsl_p; str];
    end

    strength_pre = [];
    for i=1:17
        try
            [is,os,str] = group_mean_TPE(pre{i}, labels1);
        catch
        end
        strength_pre = [strength_pre; str];
    end

    %permutation test and mean
    stats_perm_strength=[];
    strength_diff=[];
    for i=1:14
        Bas = strength_bsl_p(:,i);
        Pre = strength_pre(:,i);
        x=[Bas(:);Pre(:)];

        g1=ones(size(Bas,1),1);
        g2=2*ones(size(Pre,1),1);
        g=vertcat(g1,g2);

        stats_perm_strength = [stats_perm_strength, permutationTest(Pre,Bas,1000)];
        strength_diff = [strength_diff, mean(Bas)-mean(Pre)];
    end
    
    % ////event
    
    field3 = strcat('bsl_e_',freqs_name{iBand});
    field4 = strcat('event_',freqs_name{iBand});
    bsl = All_Conn_Real_subthresh.(field3(:));
    event = All_Conn_Real_subthresh.(field4(:));
    
    strength_bsl_e = [];
    for i=1:170
        try
            [is,os,str] = group_mean_TPE(bsl{i}, labels1);
        catch
        end
       
        strength_bsl_e = [strength_bsl_e; str];
    end

    strength_event = [];
    for i=1:17
        try
            [is,os,str] = group_mean_TPE(event{i}, labels1);
        catch
        end
        strength_event = [strength_event; str];
    end

    %permutation test and mean
    stats_perm_strength_event=[];
    strength_diff_event=[];
    for i=1:14
        Bas = strength_bsl_p(:,i);
        Pre = strength_event(:,i);
        x=[Bas(:);Pre(:)];

        g1=ones(size(Bas,1),1);
        g2=2*ones(size(Pre,1),1);
        g=vertcat(g1,g2);

        stats_perm_strength_event = [stats_perm_strength_event, permutationTest(Pre,Bas,1000)];
        strength_diff_event = [strength_diff_event, mean(Bas)-mean(Pre)];
    end
    
    Strength_Pre_Stats.(freqs_name{iBand}) = stats_perm_strength;
    Strength_Pre_Diff.(freqs_name{iBand}) = strength_diff;
    Strength_Event_Stats.(freqs_name{iBand}) = stats_perm_strength_event;
    Strength_Event_Diff.(freqs_name{iBand}) = strength_diff_event;
end

% FDR correction
%[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(Strength_Pre_Stats.delta,.05,'pdep','yes');
for iBand=1:6
    field = freqs_name{iBand};
    
    [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(Strength_Pre_Stats.(field),.05,'pdep','yes');
    Strength_Pre_adj_Stats.(freqs_name{iBand}) = adj_p
    
    [h, crit_p, adj_ci_cvrg, adj_p2]=fdr_bh(Strength_Event_Stats.(field),.05,'pdep','yes');
    Strength_Event_adj_Stats.(freqs_name{iBand}) = adj_p2
end


for iBand=1:6
    field = freqs_name{iBand};
    
    Strength_Pre_adj_Stats_masked.(freqs_name{iBand}) = Strength_Pre_adj_Stats.(freqs_name{iBand})(Strength_Pre_adj_Stats.(freqs_name{iBand})<0.05);
    Strength_Event_adj_Stats_masked.(freqs_name{iBand}) = Strength_Event_adj_Stats.(freqs_name{iBand})(Strength_Event_adj_Stats.(freqs_name{iBand})<0.05);
end
%% degree
for iBand=1:6
    iBand
    field1 = strcat('bsl_p_',freqs_name{iBand});
    field2 = strcat('pre_',freqs_name{iBand});
    bsl = All_Conn_Real_subthresh.(field1(:));
    pre = All_Conn_Real_subthresh.(field2(:));
    
    strength_bsl = [];
    for i=1:170
        [is,os,str] = degrees_dir(bsl{i});
        strength_bsl = [strength_bsl; str];
    end

    strength_pre = [];
    for i=1:17
        [is,os,str] = degrees_dir(pre{i});
        strength_pre = [strength_pre; str];
    end

    %permutation test and mean
    stats_perm_strength=[];
    strength_diff=[];
    for i=1:148
        Bas = strength_bsl(:,i);
        Pre = strength_pre(:,i);
        x=[Bas(:);Pre(:)];

        g1=ones(size(Bas,1),1);
        g2=2*ones(size(Pre,1),1);
        g=vertcat(g1,g2);

        stats_perm_strength = [stats_perm_strength, permutationTest(Pre,Bas,1000)];
        strength_diff = [strength_diff, mean(Bas)-mean(Pre)];
    end
    
    % ////event
    
    field3 = strcat('bsl_e_',freqs_name{iBand});
    field4 = strcat('event_',freqs_name{iBand});
    bsl = All_Conn_Real_subthresh.(field3(:));
    event = All_Conn_Real_subthresh.(field4(:));
    
    strength_bsl = [];
    for i=1:170
        [is,os,str] = degrees_dir(bsl{i});
        strength_bsl = [strength_bsl; str];
    end

    strength_event = [];
    for i=1:17
        [is,os,str] = degrees_dir(event{i});
        strength_event = [strength_event; str];
    end

    %permutation test and mean
    stats_perm_strength_event=[];
    strength_diff_event=[];
    for i=1:148
        Bas = strength_bsl(:,i);
        Pre = strength_event(:,i);
        x=[Bas(:);Pre(:)];

        g1=ones(size(Bas,1),1);
        g2=2*ones(size(Pre,1),1);
        g=vertcat(g1,g2);

        stats_perm_strength_event = [stats_perm_strength_event, permutationTest(Pre,Bas,1000)];
        strength_diff_event = [strength_diff_event, mean(Bas)-mean(Pre)];
    end
    
    Degree_Pre_Stats.(freqs_name{iBand}) = stats_perm_strength;
    Degree_Pre_Diff.(freqs_name{iBand}) = strength_diff;
    Degree_Event_Stats.(freqs_name{iBand}) = stats_perm_strength_event;
    Degree_Event_Diff.(freqs_name{iBand}) = strength_diff_event;
end

%FDR correction

for iBand=1:6
    field = freqs_name{iBand};
    
    [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(Degree_Pre_Stats.(field),.05,'pdep','yes');
    Degree_Pre_adj_Stats.(freqs_name{iBand}) = adj_p
    
    [h, crit_p, adj_ci_cvrg, adj_p2]=fdr_bh(Degree_Event_Stats.(field),.05,'pdep','yes');
    Degree_Event_adj_Stats.(freqs_name{iBand}) = adj_p2
end

for iBand=1:6
    field = freqs_name{iBand};
    
    Degree_Pre_adj_Stats_masked.(freqs_name{iBand}) = Degree_Pre_adj_Stats.(freqs_name{iBand})(Degree_Pre_adj_Stats.(freqs_name{iBand})<0.05);
    Degree_Event_adj_Stats_masked.(freqs_name{iBand}) = Degree_Event_adj_Stats.(freqs_name{iBand})(Degree_Event_adj_Stats.(freqs_name{iBand})<0.05);
end

%s_b1=strengths_dir(Conn_baseline_subThreshold{1});
%significance_mask1 = stats_perm_strength(stats_perm_strength<0.05);
%%
bls_mean = ALLConnect_baseline(:);
bls_mean_masked = {};
for i=1:170
    temp_matrix = cell2mat(bls_mean(i,1));
    if size(temp_matrix)>0
        v = temp_matrix .* significance_mask;
    end
    bls_mean_masked(i) = {v};%%%%
end

pre_mean = PreEventConnect_fisher3;
pre_mean_masked = {};
for i=1:17
    temp_matrix = cell2mat(pre_mean(i,1));
    v = temp_matrix .* significance_mask;
    pre_mean_masked(i) = {v};
end

%%% baseline delete below threshold
bls_mean_maskedB=cat(3,bls_mean_masked{:});
meanBsl_masked = mean(bls_mean_maskedB,3);
meanBsl_masked(abs(meanBsl_masked) < threshold_bsl) = 0;


%%% pre Event delete below threshold
pre_mean_maskedB=cat(3,pre_mean_masked{:});
meanPre_masked = mean(pre_mean_maskedB,3);
meanPre_masked(abs(meanPre_masked) < threshold_pre) = 0;

% % intersection of p surroage and p between
% meanBsl_masked(meanBsl_masked==0 | meanPre_masked==0)=0;
% meanPre_masked(meanBsl_masked==0 | meanPre_masked==0)=0;

% what is a good way to calculate diff???
mean_diff_connection =  abs(meanPre_masked) - abs(meanBsl_masked);

ts = mean_diff_connection(mean_diff_connection~=0);
tss = meanBsl_masked(meanBsl_masked~=0);
tsss = meanPre_masked(meanPre_masked~=0);



%% Function strength net flow (directed)
function [Strength_adj_Stats, Strength_diff] = perm_stat_strength(data, name1, name2)
    freqs_name = {'delta', 'theta', 'alpha', 'sigma', 'beta', 'gamma'};
    All_Conn_Real=data;
    for iBand=1:6
        iBand
        
        field1 = strcat(name1,freqs_name{iBand});
        field2 = strcat(name2,freqs_name{iBand});
        bsl = All_Conn_Real.(field1(:));
        bsl = bsl(:).';
        dev = All_Conn_Real.(field2(:));
        dev = dev(:).';
        % calculate threshold
        dist_bsl_dev = vertcat(bsl{:}, dev{:});
        dist_bsl_dev = dist_bsl_dev(dist_bsl_dev~=0);
        threshold = prctile(dist_bsl_dev(:),85);
        p_mask_perm = zeros(148,148, "double");

        strength_bsl = [];
        for i=1:size(bsl,2)
            bsl_tmp = bsl{i};
            bsl_tmp(bsl_tmp<threshold)=0;
            [is,os,str] = strengths_dir_mod(bsl_tmp);
            strength_bsl = [strength_bsl; str];
        end

        strength_dev = [];
        for i=1:size(dev,2)
            dev_tmp = dev{i};
            dev_tmp(dev_tmp<threshold)=0;
            [is,os,str] = strengths_dir_mod(dev_tmp);
            strength_dev = [strength_dev; str];
        end

        %permutation test and mean
        stats_perm_strength=[];
        strength_diff=[];
        for i=1:148
            Bas = strength_bsl(:,i);
            Dev = strength_dev(:,i);

            stats_perm_strength = [stats_perm_strength, permutationTest(Dev,Bas,1000)];
            %strength_diff = [strength_diff, mean(Bas)-mean(Dev)];
            strength_diff = [strength_diff, mean(Dev) - mean(Bas)];
        end

       
        [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(stats_perm_strength,.05,'pdep','yes');
        
        adj_p(adj_p>0.05) = 0;
        Strength_adj_Stats.(freqs_name{iBand}) = adj_p;
        Strength_diff.(freqs_name{iBand}) = strength_diff;
        
        field = strcat("Nsignificant",freqs_name{iBand});    
        Strength_adj_Stats.(field) = nnz(adj_p);
    end
end


%% Function strength (undirected)
function [Strength_adj_Stats, Strength_diff] = perm_stat_strength_undir(data, name1, name2)
    freqs_name = {'delta', 'theta', 'alpha', 'sigma', 'beta', 'gamma'};
    All_Conn_Real=data;
    for iBand=1:6
        iBand
        
        field1 = strcat(name1,freqs_name{iBand});
        field2 = strcat(name2,freqs_name{iBand});
        bsl = All_Conn_Real.(field1(:));
        bsl = bsl(:).';
        dev = All_Conn_Real.(field2(:));
        dev = dev(:).';
        % calculate threshold
        dist_bsl_dev = vertcat(bsl{:}, dev{:});
        dist_bsl_dev = dist_bsl_dev(dist_bsl_dev~=0);
        threshold = prctile(dist_bsl_dev(:),85);
        p_mask_perm = zeros(148,148, "double");

        strength_bsl = [];
        for i=1:size(bsl,2)
            bsl_tmp = bsl{i};
            bsl_tmp(bsl_tmp<threshold)=0;
            [is,os,str] = strengths_dir_mod(bsl_tmp);
            strength_bsl = [strength_bsl; str];
        end

        strength_dev = [];
        for i=1:size(dev,2)
            dev_tmp = dev{i};
            dev_tmp(dev_tmp<threshold)=0;
            [is,os,str] = strengths_dir(dev_tmp);
            strength_dev = [strength_dev; str];
        end

        %permutation test and mean
        stats_perm_strength=[];
        strength_diff=[];
        for i=1:148
            Bas = strength_bsl(:,i);
            Dev = strength_dev(:,i);

            stats_perm_strength = [stats_perm_strength, permutationTest(Dev,Bas,1000)];
            %strength_diff = [strength_diff, mean(Bas)-mean(Dev)];
            strength_diff = [strength_diff, mean(Dev) - mean(Bas)];
        end

       
        [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(stats_perm_strength,.05,'pdep','yes');
        
        adj_p(adj_p>0.05) = 0;
        Strength_adj_Stats.(freqs_name{iBand}) = adj_p;
        Strength_diff.(freqs_name{iBand}) = strength_diff;
        
        field = strcat("Nsignificant",freqs_name{iBand});    
        Strength_adj_Stats.(field) = nnz(adj_p);
    end
end


%% Threshold directed
function [adj_p_values, p_values] = perm_stat_conn_dir_thresh(data, name1, name2)
    freqs_name = {'delta', 'theta', 'alpha', 'sigma', 'beta', 'gamma'};
    All_Conn_Real = data;
    for iBand=1:1
        
        field1 = strcat(name1,freqs_name{iBand});
        field2 = strcat(name2,freqs_name{iBand});
        bsl = All_Conn_Real.(field1(:));
        bsl = bsl(:).';
        dev = All_Conn_Real.(field2(:));
        dev = dev(:).';
        % calculate threshold
        dist_bsl_dev = vertcat(bsl{:}, dev{:});
        dist_bsl_dev = dist_bsl_dev(dist_bsl_dev~=0);
        threshold = prctile(dist_bsl_dev(:),85);
        p_mask_perm = zeros(148,148, "double");
        kk=0;
        size(bsl)
        for j=size(bsl,2)
            for k=size(bsl,2)
                Bsl = [];
                for i=1:size(bsl,2)
                    if bsl{1,i}(j,k)>0
                        if bsl{1,i}(j,k) > threshold
                        Bsl(i) = bsl{1,i}(j,k);
                        else
                             Bsl(i) = 0;
                        end
                    else
                        if bsl{1,i}(k,j) < -threshold
                        Bsl(i) = - bsl{1,i}(k,j);
                        else
                             Bsl(i) = 0;
                        end
                    end
                end
                Dev = [];
                for ii=1:size(dev,2)
                    if dev{1,ii}(j,k)>0
                        if dev{1,ii}(j,k) > threshold
                        Dev(ii) = dev{1,ii}(j,k);
                        else
                             Dev(ii) = 0;
                        end
                    else
                        if dev{1,ii}(k,j) < -threshold
                        Dev(ii) = - dev{1,ii}(k,j);
                        else
                             Dev(ii) = 0;
                        end
                    end
                end
                % permutation test
                Bsl = Bsl.';
                Dev = Dev.';

                if sum(Dev)+sum(Bsl) == 0
                    p_mask_perm(j,k) = 1;
                else
                    p_mask_perm(j,k) = permutationTest(Dev,Bsl,1000);
                end
                
            end
            kk=kk+1;
            %kk
        end
    stats_band_dir = p_mask_perm;
    %FDR correction
    p_values.(freqs_name{iBand}) = p_mask_perm;
    
    m  = tril(true(size(stats_band_dir)), -1);
    v  = stats_band_dir(m).';
    [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(v,.05,'pdep','yes');

    A = tril(ones(148),-1);
    A(A > 0) = adj_p;
    A(A>0.05)=0;
    adj_p_values.(freqs_name{iBand}) = A;
    
    field = strcat("Nsignificant",freqs_name{iBand});    
    adj_p_values.(field) = nnz(A);
    end
end

%% Threshold directed individual nodes (Full matrix)
function adj_p_values = perm_stat_conn_dir_full_thresh(data, name1, name2)
    freqs_name = {'delta', 'theta', 'alpha', 'sigma', 'beta', 'gamma'};
    All_Conn_Real = data;
    for iBand=1:6
        
        field1 = strcat(name1,freqs_name{iBand});
        field2 = strcat(name2,freqs_name{iBand});
        bsl = All_Conn_Real.(field1(:));
        bsl = bsl(:).';
        dev = All_Conn_Real.(field2(:));
        dev = dev(:).';
        % calculate threshold
        dist_bsl_dev = vertcat(bsl{:}, dev{:});
        dist_bsl_dev = dist_bsl_dev(dist_bsl_dev~=0);
        threshold = prctile(dist_bsl_dev(:),85);
        p_mask_perm = zeros(14,14, "double");
        %kk=0;
        for j=1:148
            for k=1:148
                Bsl = [];
                for i=1:size(bsl,2)
                    if bsl{1,i}(j,k)>0
                        if bsl{1,i}(j,k) > threshold
                        Bsl(i) = bsl{1,i}(j,k);
                        else
                             Bsl(i) = 0;
                        end
                    end
                end
                Dev = [];
                for ii=1:size(dev,2)
                    if dev{1,ii}(j,k)>0
                        if dev{1,ii}(j,k) > threshold
                        Dev(ii) = dev{1,ii}(j,k);
                        else
                             Dev(ii) = 0;
                        end
                    end
                end
                % permutation test
                Bsl = Bsl.';
                Dev = Dev.';

                if sum(Dev)+sum(Bsl) == 0
                    p_mask_perm(j,k) = 1;
                else
                    p_mask_perm(j,k) = permutationTest(Dev,Bsl,1000);
                end
                
            end
            %kk=kk+1;
            %kk
        end
       
    
    [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p_mask_perm,.05,'pdep','yes');    
    adj_p(adj_p>0.05)=0;
    adj_p_values.(freqs_name{iBand}) = adj_p;

    A = adj_p;
    A(A>0.05)=0;
    adj_p_values.(freqs_name{iBand}) = A;
    
    field = strcat("Nsignificant",freqs_name{iBand});    
    adj_p_values.(field) = nnz(A);
    end
end

%% Threshold Undirected
%///////// permutation statistics function
function [adj_p_values, p_values] = perm_stat_conn_undir_thresh(data, name1, name2)
    freqs_name = {'delta', 'theta', 'alpha', 'sigma', 'beta', 'gamma'};
    All_Conn_Real = data;
    for iBand=1:6
        
        field1 = strcat(name1,freqs_name{iBand});
        field2 = strcat(name2,freqs_name{iBand});
        bsl = All_Conn_Real.(field1(:));
        bsl = bsl(:).';
        dev = All_Conn_Real.(field2(:));
        dev = dev(:).';
        % calculate threshold
        dist_bsl_dev = vertcat(bsl{:}, dev{:});
        dist_bsl_dev = dist_bsl_dev(dist_bsl_dev~=0);
        threshold = prctile(dist_bsl_dev(:),85);
        p_mask_perm = zeros(14,14, "double");
        kk=0;
        for j=1:14
            for k=1:kk
                Bsl = [];
                for i=1:size(bsl,2)
                    if bsl{1,i}(j,k)>0
                        if bsl{1,i}(j,k) > threshold
                        Bsl(i) = bsl{1,i}(j,k);
                        else
                             Bsl(i) = 0;
                        end
                    else
                        if bsl{1,i}(k,j) > threshold
                        Bsl(i) = bsl{1,i}(k,j);
                        else
                             Bsl(i) = 0;
                        end
                    end
                end
                Dev = [];
                for ii=1:size(dev,2)
                    if dev{1,ii}(j,k)>0
                        if dev{1,ii}(j,k) > threshold
                        Dev(ii) = dev{1,ii}(j,k);
                        else
                             Dev(ii) = 0;
                        end
                    else
                        if dev{1,ii}(k,j) > threshold
                        Dev(ii) = dev{1,ii}(k,j);
                        else
                             Dev(ii) = 0;
                        end
                    end
                end
                % permutation test
                Bsl = Bsl.';
                Dev = Dev.';

                if sum(Dev)+sum(Bsl) == 0
                    p_mask_perm(j,k) = 1;
                else
                    p_mask_perm(j,k) = permutationTest(Dev,Bsl,1000);
                end
                
            end
            kk=kk+1;
            %kk
        end
    stats_band_dir = p_mask_perm;
    %FDR correction
    p_values.(freqs_name{iBand}) = p_mask_perm;
    
    m  = tril(true(size(stats_band_dir)), -1);
    v  = stats_band_dir(m).';
    [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(v,.05,'pdep','yes');

    A = tril(ones(14),-1);
    A(A > 0) = adj_p;
    A(A>0.05)=0;
    adj_p_values.(freqs_name{iBand}) = A;
    
    field = strcat("Nsignificant",freqs_name{iBand});    
    adj_p_values.(field) = nnz(A);
    end
end

%% NO threshold
%///////// permutation statistics function
function [adj_p_values, p_values] = perm_stat_conn_dir(data, name1, name2)
    freqs_name = {'delta', 'theta', 'alpha', 'sigma', 'beta', 'gamma'};
    All_Conn_Real = data;
    for iBand=1:6
        
        field1 = strcat(name1,freqs_name{iBand});
        field2 = strcat(name2,freqs_name{iBand});
        bsl = All_Conn_Real.(field1(:));
        bsl = bsl(:).';
        dev = All_Conn_Real.(field2(:));
        dev = dev(:).';
        % calculate threshold
        dist_bsl_dev = vertcat(bsl{:}, dev{:});
        dist_bsl_dev = dist_bsl_dev(dist_bsl_dev~=0);
        threshold = prctile(dist_bsl_dev(:),85);
        p_mask_perm = zeros(148,148, "double");
        kk=0;
        for j=1:148
            for k=1:kk
                Bsl = [];
                for i=1:size(bsl,2)
                    if bsl{1,i}(j,k)>0
                        Bsl(i) = bsl{1,i}(j,k);

                    else
                        Bsl(i) = - bsl{1,i}(k,j);

                    end
                end
                Dev = [];
                for ii=1:size(dev,2)
                    if dev{1,ii}(j,k)>0
                        Dev(ii) = dev{1,ii}(j,k);

                    else
                        Dev(ii) = - dev{1,ii}(k,j);
                    end
                end
                % permutation test
                Bsl = Bsl.';
                Dev = Dev.';

                if sum(Dev)+sum(Bsl) == 0
                    p_mask_perm(j,k) = 1;
                else
                    p_mask_perm(j,k) = permutationTest(Dev,Bsl,1000);
                end
                
            end
            kk=kk+1;
            %kk
        end
    stats_band_dir = p_mask_perm;
    %FDR correction
    p_values.(freqs_name{iBand}) = p_mask_perm;
    m  = tril(true(size(stats_band_dir)), -1);
    v  = stats_band_dir(m).';
    [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(v,.05,'pdep','yes');

    A = tril(ones(148),-1);
    A(A > 0) = adj_p;
    A(A>0.05)=0;
    adj_p_values.(freqs_name{iBand}) = A;
    

    field = strcat("Nsignificant",freqs_name{iBand});    
    adj_p_values.(field) = nnz(A);
    
    end
end

%% NO threshold Undirected
function [adj_p_values, p_values] = perm_stat_conn_undir(data, name1, name2)
    freqs_name = {'delta', 'theta', 'alpha', 'sigma', 'beta', 'gamma'};
    All_Conn_Real = data;
    for iBand=1:6
        
        field1 = strcat(name1,freqs_name{iBand});
        field2 = strcat(name2,freqs_name{iBand});
        bsl = All_Conn_Real.(field1(:));
        bsl = bsl(:).';
        dev = All_Conn_Real.(field2(:));
        dev = dev(:).';
        % calculate threshold
        dist_bsl_dev = vertcat(bsl{:}, dev{:});
        dist_bsl_dev = dist_bsl_dev(dist_bsl_dev~=0);
        threshold = prctile(dist_bsl_dev(:),85);
        p_mask_perm = zeros(148,148, "double");
        kk=0;
        for j=1:148
            for k=1:kk
                Bsl = [];
                for i=1:size(bsl,2)
                    if bsl{1,i}(j,k)>0
%                         if bsl{1,i}(j,k) > threshold
                        Bsl(i) = bsl{1,i}(j,k);
%                         else
%                              Bsl(i) = 0;
%                         end
                    else
%                         if bsl{1,i}(k,j) < -threshold
                        Bsl(i) = bsl{1,i}(k,j);
%                         else
%                              Bsl(i) = 0;
%                         end
                    end
                end
                Dev = [];
                for ii=1:size(dev,2)
                    if dev{1,ii}(j,k)>0
%                         if dev{1,ii}(j,k) > threshold
                        Dev(ii) = dev{1,ii}(j,k);
%                         else
%                              Dev(ii) = 0;
%                         end
                    else
%                         if dev{1,ii}(k,j) < -threshold
                        Dev(ii) = dev{1,ii}(k,j);
%                         else
%                              Dev(ii) = 0;
%                         end
                    end
                end
                % permutation test
                Bsl = Bsl.';
                Dev = Dev.';

                if sum(Dev)+sum(Bsl) == 0
                    p_mask_perm(j,k) = 1;
                else
                    p_mask_perm(j,k) = permutationTest(Dev,Bsl,1000);
                end
                
            end
            kk=kk+1;
            %kk
        end
    stats_band_dir = p_mask_perm;
    %FDR correction
    p_values.(freqs_name{iBand}) = p_mask_perm;
    m  = tril(true(size(stats_band_dir)), -1);
    v  = stats_band_dir(m).';
    [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(v,.05,'pdep','yes');

    A = tril(ones(148),-1);
    A(A > 0) = adj_p;
    A(A>0.05)=0;
    adj_p_values.(freqs_name{iBand}) = A;
    end
end



