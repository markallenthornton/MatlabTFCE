% fwe_control.m
% demo of package's control over familywise error rate
nsim = 1000;
%% one sample ttest, one-sided
fp = NaN(nsim,1);
for i = 1:nsim
    imgs = randn(4,4,4,20);
    pcorr = matlab_tfce('onesample',1,imgs);
    fp(i) = sum(pcorr(:)<.05);
end
sum(fp)/nsim % false positive rate

%% one sample ttest, two-sided
fp = NaN(nsim,1);
for i = 1:nsim
    imgs = randn(4,4,4,20);
    [pcorr_pos,pcorr_neg] = matlab_tfce('onesample',2,imgs);
    fp(i) = sum(pcorr_pos(:)<.05)+sum(pcorr_neg(:)<.05);
end
sum(fp)/nsim % false positive rate

%% independent sample ttest, one-sided
fp = NaN(nsim,1);
for i = 1:nsim
    imgs = randn(4,4,4,20);
    imgs2 = randn(4,4,4,20);
    [pcorr_pos] = matlab_tfce('independent',1,imgs,imgs2);
    fp(i) = sum(pcorr_pos(:)<.05)+sum(pcorr_neg(:)<.05);
end
sum(fp)/nsim % false positive rate

%% independent sample ttest, two-sided
fp = NaN(nsim,1);
for i = 1:nsim
    imgs = randn(4,4,4,20);
    imgs2 = randn(4,4,4,20);
    [pcorr_pos,pcorr_neg] = matlab_tfce('independent',2,imgs,imgs2);
    fp(i) = sum(pcorr_pos(:)<.05)+sum(pcorr_neg(:)<.05);
end
sum(fp)/nsim % false positive rate
%% correlation, one-sided
fp = NaN(nsim,1);
for i = 1:nsim
    imgs = randn(4,4,4,20);
    covariate = randn(20,1);
    pcorr = matlab_tfce('correlation',1,imgs,[],covariate);
    fp(i) = sum(pcorr(:)<.05)+sum(pcorr(:)<.05);
end
sum(fp)/nsim % false positive rate

%% correlation, two-sided
fp = NaN(nsim,1);
for i = 1:nsim
    imgs = randn(4,4,4,20);
    covariate = randn(20,1);
    [pcorr_pos,pcorr_neg] = matlab_tfce('correlation',2,imgs,[],covariate);
    fp(i) = sum(pcorr_pos(:)<.05)+sum(pcorr_neg(:)<.05);
end
sum(fp)/nsim % false positive rate
