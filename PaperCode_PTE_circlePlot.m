freqs_name = {'delta', 'theta', 'alpha', 'sigma', 'beta', 'gamma'};
epi = mean_pre;
for f=1:6
    name = 'pre_';
    TF1 = epi.TF(:,:,f);

    %threshold = prctile(TF(:),99.5);
    threshold_k = min(maxk(abs(TF1(:)),100));

    TF2 = TF1;
    TF2(abs(TF2)<threshold_k)=0;
    %
    TFM = reshape(TF2, [148, 148]);

    to_del = [];
    for i = 1:148
        nnz(TFM(i,:)) + nnz(TFM(:,i))
        if nnz(TFM(i,:)) + nnz(TFM(:,i)) == 0
            to_del = [to_del i];
        end
    end

    Batlas = epi.Atlas.Scouts;
    Batlas(to_del) = [];

    for i=1:148
        for j=1:148
            if TFM(i,j) < 0
                TFM(j,i) = -TFM(i,j);
                TFM(i,j) = 0;
            end
        end
    end

    TFM(TFM>0)=1;
    regions=struct2cell(Batlas);
    labels1=regions(6,:).';
    scoutl=regions(4,:).';

    uniques={'LL';'LO';'LP';'LT';'LC';'LF';'LPF';'RPF';'RF';'RC';'RT';'RP';'RO';'RL' };
    rng(1);
    ii=1;
    reord_label = cell(148,1);
    reord_scout = cell(148,1);
    reord_idx = [];
    for i=1:length(uniques);
        for j=1:length(scoutl);
            if ismember(labels1(j),uniques(i));
                reord_label(ii) = uniques(i);
                reord_idx(ii) = j;
                reord_scout(ii) = scoutl(j);
                ii=ii+1;
            end
        end     
    end

    B = TFM(:,reord_idx);
    reord_M = B(reord_idx,:);
end

save ("reord_labelHigh.mat", "reord_label");
save ("reord_scouts.mat", "reord_scout");
save ("TFM_reord.mat", "reord_M");