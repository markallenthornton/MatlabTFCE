% fwe_control.m
% demo of package's control over familywise error rate
nsim = 100;
%% one sample ttest, one-sided
fp = NaN(nsim,1);
for i = 1:nsim
    imgs = randn(4,4,4,20);
    pcorr = matlab_tfce('onesample',1,imgs);
    fp(i) = sum(pcorr(:)<.05);
end
sum(fp>0)/nsim % false positive rate

%% one sample ttest, two-sided
fp = NaN(nsim,1);
for i = 1:nsim
    imgs = randn(4,4,4,20);
    [pcorr_pos,pcorr_neg] = matlab_tfce('onesample',2,imgs);
    fp(i) = sum(pcorr_pos(:)<.05)+sum(pcorr_neg(:)<.05);
end
sum(fp>0)/nsim % false positive rate

%% independent sample ttest, one-sided
fp = NaN(nsim,1);
for i = 1:nsim
    imgs = randn(4,4,4,20);
    imgs2 = randn(4,4,4,20);
    [pcorr_pos] = matlab_tfce('independent',1,imgs,imgs2);
    fp(i) = sum(pcorr_pos(:)<.05)+sum(pcorr_neg(:)<.05);
end
sum(fp>0)/nsim % false positive rate

%% independent sample ttest, two-sided
fp = NaN(nsim,1);
for i = 1:nsim
    imgs = randn(4,4,4,20);
    imgs2 = randn(4,4,4,20);
    [pcorr_pos,pcorr_neg] = matlab_tfce('independent',2,imgs,imgs2);
    fp(i) = sum(pcorr_pos(:)<.05)+sum(pcorr_neg(:)<.05);
end
sum(fp>0)/nsim % false positive rate
%% correlation, one-sided
fp = NaN(nsim,1);
for i = 1:nsim
    imgs = randn(4,4,4,20);
    covariate = randn(20,1);
    pcorr = matlab_tfce('correlation',1,imgs,[],covariate,500);
    fp(i) = sum(pcorr(:)<.05);
end
sum(fp>0)/nsim % false positive rate

%% correlation, two-sided
fp = NaN(nsim,1);
for i = 1:nsim
    imgs = randn(4,4,4,20);
    covariate = randn(20,1);
    [pcorr_pos,pcorr_neg] = matlab_tfce('correlation',2,imgs,[],covariate);
    fp(i) = sum(pcorr_pos(:)<.05)+sum(pcorr_neg(:)<.05);
end
sum(fp>0)/nsim % false positive rate
%% one-way ANOVA
fp = NaN(nsim,1);
for i = 1:nsim
    imgs = {randn(4,4,4,20);randn(4,4,4,20);randn(4,4,4,20)};
    pcorr = matlab_tfce('rm_anova1',1,imgs);
    fp(i) = sum(pcorr(:)<.05);
end
sum(fp>0)/nsim % false positive rate
%% two-way ANOVA
fp1 = NaN(nsim,1);
fp2 = NaN(nsim,1);
fpi = NaN(nsim,1);
for i = 1:nsim
    imgs = {randn(4,4,4,20) randn(4,4,4,20);randn(4,4,4,20) randn(4,4,4,20);randn(4,4,4,20) randn(4,4,4,20)};
    [pcorr_fac1,pcorr_fac2,pcorr_int] = matlab_tfce('rm_anova2',1,imgs);
    fp1(i) = sum(pcorr_fac1(:)<.05);
    fp2(i) = sum(pcorr_fac2(:)<.05);
    fpi(i) = sum(pcorr_int(:)<.05);
end
sum(fp1>0)/nsim % false positive rate
sum(fp2>0)/nsim % false positive rate
sum(fpi>0)/nsim % false positive rate

%% multiple regression
fp = NaN(nsim,5);
for i = 1:nsim
    imgs = randn(4,4,4,20);
    covariate = [ones(20,1) randn(20,4)];
    [pcorr_pos,pcorr_neg] = matlab_tfce('regression',2,imgs,[],covariate,[],[],[],[],[],[],[1 0 0 1 1]);
    for j = 1:5
        fp(i,j) = sum(pcorr_pos{j}(:)<.05)+sum(pcorr_neg{j}(:)<.05);
    end
end
sum(fp>0)/nsim