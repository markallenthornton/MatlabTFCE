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

% calculate true f-stat image
S = repmat((1:nsub)',levels(1)*levels(2),1);
F1 = repmat(sort(repmat((1:levels(1))',nsub,1)),levels(2),1);
F2 = sort(repmat((1:levels(2))',nsub*levels(1),1));
fvals1 = NaN(nvox,1);
fvals2 = NaN(nvox,1);
fvalsi = NaN(nvox,1);
for v = 1:nvox
    curimg = occimgs(:,:,:,v);
    Y = curimg(:);
    stats = rm_anova2(Y,S,F1,F2,{'F1','F2'});
    fvals1(v) = stats{2,5};
    fvals2(v) = stats{3,5};
    fvalsi(v) = stats{4,5};
end
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

% cycle through permutations
exceedances1 = zeros(nvox,1);
exceedances2 = zeros(nvox,1);
exceedancesi = zeros(nvox,1);
for p = 1:nperm
    
    % permute labels
    roccimgs = occimgs;
    for s = 1:nsub
        relabeling1 = randsample(levels(1),levels(1));
        relabeling2 = randsample(levels(2),levels(2));
        roccimgs(s,relabeling1,relabeling2,:) = roccimgs(s,:,:,:);
    end
    
    % calculate permutation statistic
    rfvals1 = NaN(nvox,1);
    rfvals2 = NaN(nvox,1);
    rfvalsi = NaN(nvox,1);
    for v = 1:nvox
        curimg = roccimgs(:,:,:,v);
        Y = curimg(:);
        stats = rm_anova2(Y,S,F1,F2,{'F1','F2'});
        rfvals1(v) = stats{2,5};
        rfvals2(v) = stats{3,5};
        rfvalsi(v) = stats{4,5};
    end
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

