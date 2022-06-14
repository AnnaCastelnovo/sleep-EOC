restoredefaultpath
%%
addpath .../fieldtrip-20210616
ft_defaults

%%
addpath('/.../BrainNetViewer_20191031')
surface='/.../BrainMesh_ICBM152.nv';
node='/.../destrieux_corrected.node';
conf='/.../destrieux/cfg4.mat';

%edge='/.../reduced_edges_bsl_event1_delta.edge';
%BrainNet_MapCfg(surface,node, edge,conf)

%% /////////////////////
% Permutation test for Node strength and net flow information 
%  /////////////////////
%% Node Strength
[Strength_adj_Stats_bsl_event0_T, Strength_diff_bsl_event0] = perm_stat_strength_undir(All_Conn, "bsl_", "event0_");
[Strength_adj_Stats_pre_event0_t, Strength_diff_pre_event0] = perm_stat_strength_undir(All_Conn, "pre_", "event0_");
[Strength_adj_Stats_bsl_pre_T, Strength_diff_bsl_pre] = perm_stat_strength_undir(All_Conn, "bsl_", "pre_");
[Strength_adj_Stats_bsl_event1_T, Strength_diff_bsl_event1] = perm_stat_strength_undir(All_Conn, "bsl_", "event1_");
[Strength_adj_Stats_bsl_event2_T, Strength_diff_bsl_event2] = perm_stat_strength_undir(All_Conn, "bsl_", "event2_");
[Strength_adj_Stats_pre_event1_t, Strength_diff_pre_event1] = perm_stat_strength_undir(All_Conn, "pre_", "event1_");
[Strength_adj_Stats_pre_event2_T, Strength_diff_pre_event2] = perm_stat_strength_undir(All_Conn, "pre_", "event2_");

%% Net information flow
[Strength_adj_Stats_bsl_event0_T_net, Strength_diff_bsl_net_event0] = perm_stat_strength(All_Conn, "bsl_", "event0_");
[Strength_adj_Stats_pre_event0_T_net, Strength_diff_pre_net_event0] = perm_stat_strength(All_Conn, "pre_", "event0_");
[Strength_adj_Stats_bsl_pre_T_net, Strength_diff_bsl_net_pre] = perm_stat_strength(All_Conn, "bsl_", "pre_");
[Strength_adj_Stats_bsl_event1_T_net, Strength_diff_bsl_net_event1] = perm_stat_strength(All_Conn, "bsl_", "event1_");
[Strength_adj_Stats_bsl_event2_T_net, Strength_diff_bsl_net_event2] = perm_stat_strength(All_Conn, "bsl_", "event2_");
[Strength_adj_Stats_pre_event1_T_net, Strength_diff_pre_net_event1] = perm_stat_strength(All_Conn, "pre_", "event1_");
[Strength_adj_Stats_pre_event2_T_net, Strength_diff_pre_net_event2] = perm_stat_strength(All_Conn, "pre_", "event2_");

%% This section edits the node file for BrainNetViewer
% 1. Node Strength
freqs_name = {'delta', 'theta', 'alpha', 'sigma', 'beta', 'gamma'};

%edit names manually !!
name2 = Strength_diff_pre_event0;
name1 = Strength_adj_Stats_pre_event0_t;

name_eps = 'Strength_adj_Stats_pre_event0_T_';

for iBand=1:6
     field1 = strcat('Nsignificant',freqs_name{iBand});

    if name1.(field1) > 0
    
        sig = name1.(freqs_name{iBand});

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

        %save_name
    end
    
end


%% This section edits the node file for BrainNetViewer
% 1. Node Strength
freqs_name = {'delta', 'theta', 'alpha', 'sigma', 'beta', 'gamma'};

%load names manually!!
name2 = Strength_diff_pre_event0;
name1 = Strength_adj_Stats_pre_event0_T_net;

name_file = 'Strength_adj_Stats_pre_event0_T_net_';

for iBand=1:6
     field1 = strcat('Nsignificant',freqs_name{iBand});
    if name1.(field1) > 0
    
        % edit node size -> diff in node strenght
        sig = name1.(freqs_name{iBand});

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

        %save_name
    end
    
end


%% Plot images with BrainNetViewer
%add correct path names
addpath('.../connectivity/BrainNetViewer_20191031')
surface='.../destrieux/BrainMesh_ICBM152.nv';
conf='.../destrieux/cfg5.mat';
n_scouts=148;

%load all node files
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

            stats_perm_strength = [stats_perm_strength, permutationTest(Dev,Bas,5000)];
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


function [is,os,str] = strengths_dir_mod(CIJ)
% compute net flow
is = sum(CIJ,1);    % instrength = column sum of CIJ
os = sum(CIJ,2)';   % outstrength = row sum of CIJ
str = os - is;
end

function [is,os,str] = strengths_dir(CIJ)
% compute strengths
is = sum(CIJ,1);    % instrength = column sum of CIJ
os = sum(CIJ,2)';   % outstrength = row sum of CIJ

str = is+os;        % strength = instrength+outstrength

