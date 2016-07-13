function [varargout] = matlab_tfce_rm_anova2(imgs,nperm,H,E,C,dh)
%MATLAB_TFCE_RM_ANOVA2 performs maximal statistic permutation testing
%   [varargout] = matlab_tfce_rm_anova2(imgs,nperm,H,E,C,dh) corrects
%   a two-way factorial ANOVA (mean ~=) for multiple comparisons via a 
%   permutation procedure (random condition switching). Maximal f-stats of
%   permuted data are compared with f-stats in the unshuffled data.
%
%   Arguments:
%   imgs -- an m x 1 cell array of 4D (x,y,z,subject) matrices of images
%   nperm -- number of permutations to perform. More permutations yield
%   more precise correct p-values.
%   -- img the 3D image to be transformed
%   -- H height exponent
%   -- E extent exponent
%   -- C connectivity
%   -- ndh step number for cluster formation
%
%   Output:
%	pcorr_fac1 -- corrected p-values for factor one
%	pcorr_fac2 -- corrected p-values for factor two
%	pcorr_int  -- corrected p-values for interaction term

% calculate matrix size
bsize = size(imgs{1});
nsub = bsize(4);
bsize = bsize(1:3);
levels = size(imgs);

% set tranform function
transform = @matlab_tfce_transform;

% calculate implicit mask
cellmeans = cellfun(@(x) sum(x,4),imgs,'UniformOutput',false);
implicitmask = sum(cat(4,cellmeans{:}),4);
implicitmask = ~(implicitmask == 0 | isnan(implicitmask));
nvox = sum(implicitmask(:));

% extract occupied voxels for permutation test
occimgs = NaN(nsub,levels(1),levels(2),nvox);
for i = 1:levels(1)
    for j = 1:levels(2)
        for s = 1:nsub
            curimg = imgs{i,j}(:,:,:,s);
            occimgs(s,i,j,:) = curimg(implicitmask);
        end
    end
end
clear imgs curimg

% calculate means
smeans = mean(mean(occimgs,2),3);
c1means = mean(mean(occimgs,1),3);
c2means = mean(mean(occimgs,1),2);
c1smeans = mean(occimgs,3);
c2smeans = mean(occimgs,2);
ccmeans = mean(occimgs,1);
gmeans = mean(mean(mean(occimgs)));
gmeans_c1 = repmat(gmeans,[1,levels(1),1,1]);
gmeans_c2 = repmat(gmeans,[1,1,levels(2),1]);
gmeans_c1s = repmat(gmeans,[nsub,levels(1),1,1]);
gmeans_c2s = repmat(gmeans,[nsub,1,levels(2),1]);
gmeans_cc = repmat(gmeans,[1,levels(1),levels(2),1]);
gmeans_s = repmat(gmeans,[nsub,1,1,1]);

% calculate SS
SSfac1 = sum((c1means-gmeans_c1).^2,2)*levels(2)*nsub;
SSfac2 = sum((c2means-gmeans_c2).^2,3)*levels(1)*nsub;
SScell = sum(sum((ccmeans-gmeans_cc).^2,2),3)*nsub;
SSsub = sum((smeans-gmeans_s).^2,1)*levels(1)*levels(2);
SSc1s = sum(sum((c1smeans-gmeans_c1s).^2,2),1)*levels(2);
SSc2s = sum(sum((c2smeans-gmeans_c2s).^2,3),1)*levels(1);
SStot = NaN([1,1,1,nvox]);
for v = 1:nvox
    diffs = occimgs(:,:,:,v)-gmeans(v);
    SStot(v) = sum(diffs(:).^2);
end
SSerr1 = SSc1s - SSfac1 - SSsub; 
SSerr2 = SSc2s - SSfac2 - SSsub; 
SSint = SScell - SSfac1 - SSfac2;
SSerr = SStot-SSfac1-SSfac2-SSint-SSerr1-SSerr2-SSsub;

% calculate MS
MSfac1 = SSfac1./(levels(1)-1);
MSfac2 = SSfac2./(levels(2)-1);
MSint = SSint./((levels(1)-1)*(levels(2)-1));
MSerr1 = SSerr1./((levels(1)-1)*(nsub-1));
MSerr2 = SSerr2./((levels(2)-1)*(nsub-1));
MSerr = SSerr./((levels(1)-1)*(levels(2)-1)*(nsub-1));

% calculate F
fvals1 = MSfac1./MSerr1;
fvals2 = MSfac2./MSerr2;
fvalsi = MSint./MSerr;

% perform TFCE
truestat1 = zeros(bsize);
truestat2 = zeros(bsize);
truestati = zeros(bsize);
truestat1(implicitmask) = fvals1;
truestat2(implicitmask) = fvals2;
truestati(implicitmask) = fvalsi;
tfcestat1 = transform(truestat1,H,E,C,dh);
tfcestat2 = transform(truestat2,H,E,C,dh);
tfcestati = transform(truestati,H,E,C,dh);
tfcestat1 = tfcestat1(implicitmask);
tfcestat2 = tfcestat2(implicitmask);
tfcestati = tfcestati(implicitmask);

% initialize progress indicator
parfor_progress(nperm);
global parworkers

% cycle through permutations
exceedances1 = zeros(nvox,1);
exceedances2 = zeros(nvox,1);
exceedancesi = zeros(nvox,1);
parfor(p = 1:nperm,parworkers)
    
    % permute labels
    roccimgs = NaN(size(occimgs));
    for s = 1:nsub
        relabeling1 = randperm(levels(1))';
        relabeling2 = randperm(levels(2))';
        roccimgs(s,relabeling1,relabeling2,:) = occimgs(s,:,:,:);
    end
    
    % calculate permuted means
    rc1means = mean(mean(roccimgs,1),3);
    rc2means = mean(mean(roccimgs,1),2);
    rc1smeans = mean(roccimgs,3);
    rc2smeans = mean(roccimgs,2);
    rccmeans = mean(roccimgs,1);
    
    % calculate SS
    rSSfac1 = sum((rc1means-gmeans_c1).^2,2)*levels(2)*nsub;
    rSSfac2 = sum((rc2means-gmeans_c2).^2,3)*levels(1)*nsub;
    rSScell = sum(sum((rccmeans-gmeans_cc).^2,2),3)*nsub;
    rSSc1s = sum(sum((rc1smeans-gmeans_c1s).^2,2),1)*levels(2);
    rSSc2s = sum(sum((rc2smeans-gmeans_c2s).^2,3),1)*levels(1);
    rSSerr1 = rSSc1s - rSSfac1 - SSsub; 
    rSSerr2 = rSSc2s - rSSfac2 - SSsub; 
    rSSint = rSScell - rSSfac1 - rSSfac2;
    rSSerr = SStot-rSSfac1-rSSfac2-rSSint-rSSerr1-rSSerr2-SSsub;

    % calculate MS
    rMSfac1 = rSSfac1./(levels(1)-1);
    rMSfac2 = rSSfac2./(levels(2)-1);
    rMSint = rSSint./((levels(1)-1)*(levels(2)-1));
    rMSerr1 = rSSerr1./((levels(1)-1)*(nsub-1));
    rMSerr2 = rSSerr2./((levels(2)-1)*(nsub-1));
    rMSerr = rSSerr./((levels(1)-1)*(levels(2)-1)*(nsub-1));

    % calculate F
    rfvals1 = rMSfac1./rMSerr1;
    rfvals2 = rMSfac2./rMSerr2;
    rfvalsi = rMSint./rMSerr;
    
    % perform TFCE
    rbrain1 = zeros(bsize);
    rbrain2 = zeros(bsize);
    rbraini = zeros(bsize);
    rbrain1(implicitmask) = rfvals1;
    rbrain2(implicitmask) = rfvals2;
    rbraini(implicitmask) = rfvalsi;
    rbrain1 = transform(rbrain1,H,E,C,dh);
    rbrain2 = transform(rbrain2,H,E,C,dh);
    rbraini = transform(rbraini,H,E,C,dh);
    rstats1 = rbrain1(implicitmask);
    rstats2 = rbrain2(implicitmask);
    rstatsi = rbraini(implicitmask);
    
    % compare maxima to true statistic and increment as appropriate
    curexceeds1 = max(rstats1) >= tfcestat1;
    curexceeds2 = max(rstats2) >= tfcestat2;
    curexceedsi = max(rstatsi) >= tfcestati;
    exceedances1 = exceedances1 + curexceeds1;
    exceedances2 = exceedances2 + curexceeds2;
    exceedancesi = exceedancesi + curexceedsi;
    
    % update progress indicator (only does so 1 in 5 to minimize overhead)
    if ~randi([0 4]);
        parfor_progress;
    end
    
end

% create corrected p-value image
corrected1 = exceedances1./nperm;
corrected2 = exceedances2./nperm;
correctedi = exceedancesi./nperm;
pcorr_fac1 = ones(bsize);
pcorr_fac2 = ones(bsize);
pcorr_int = ones(bsize);
pcorr_fac1(implicitmask) = corrected1;
pcorr_fac2(implicitmask) = corrected2;
pcorr_int(implicitmask) = correctedi;

% assign output to varargout
varargout{1} = pcorr_fac1;
varargout{2} = pcorr_fac2;
varargout{3} = pcorr_int;

end

