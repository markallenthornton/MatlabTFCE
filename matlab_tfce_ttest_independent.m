function [varargout] = matlab_tfce_ttest_independent(imgs1,imgs2,tails,nperm,H,E,C,dh)
%MATLAB_TFCE_TTEST_INDEPENDENT tests means difference, independent samples
%   [varargout] = matlab_tfce_ttest_independent(imgs1,imgs2,tails,nperm)
%   Independent (two-sample) t-tests. The alternative hypothesis is 
%   imgs1>imgs2. Maximal t-statistics from tests on permuted data are 
%   compared with t-sstatistics in the original data.
%
%   Arguments:
%   imgs1 -- images from group 1 with dimensions x,y,z,nsubject1
%   imgs2 -- images from group 2 with dimensions x,y,z,nsubject2
%	tails -- 1 or 2 tailed test
%   nperm -- number of permutations
%   -- img the 3D image to be transformed
%   -- H height exponent
%   -- E extent exponent
%   -- C connectivity
%   -- ndh step number for cluster formation
%
%   Output:
%	If tails == 1:
%   pcorr -- wholebrain map of corrected p-values
%	If tails == 2:
%	pcorr_pos -- corrected p-values for positive effects
%	pcorr_neg -- corrected p-values for negative effects

% calculate matrix size
bsize = size(imgs1);
nsub1 = bsize(4);
bsize = size(imgs2);
nsub2 = bsize(4);
bsize = bsize(1:3);
nsub = nsub1+nsub2;

% set tranform function
if tails == 1
    transform = @matlab_tfce_transform;
else
    transform = @matlab_tfce_transform_twotailed;
end

% calculate true mean image
truestat = (mean(imgs1,4)-mean(imgs2,4))./sqrt(var(imgs1,0,4)/nsub1+var(imgs2,0,4)/nsub2);
implicitmask = ~isnan(truestat);
tfcestat = transform(truestat,H,E,C,dh);

% p-values for comparison
tvals = tfcestat(implicitmask);
if tails == 2
	tvals=abs(tvals);
end
nvox = length(tvals);

% extract occupied voxels for permutation test
imgs = cat(4,imgs1,imgs2);
glabs = [ones(nsub1,1);ones(nsub2,1)*2];
occimgs = NaN(nvox,nsub);
for s = 1:nsub
    curimg = imgs(:,:,:,s);
    occimgs(:,s) = curimg(implicitmask);
end

% initialize progress indicator
parfor_progress(nperm);
global parworkers

% cycle through permutations
exceedances = zeros(nvox,1);
parfor(p = 1:nperm,parworkers)
    
    % permute signs
    rsel = randperm(nsub)';
    relabeling = glabs(rsel);
    rimgs1 = occimgs(:,relabeling==1);
    rimgs2 = occimgs(:,relabeling==2);
    
    % calculate permutation means
    rstats = (mean(rimgs1,2)-mean(rimgs2,2))./sqrt(var(rimgs1,0,2)/nsub1+var(rimgs2,0,2)/nsub2);
    rbrain = zeros(bsize);
    rbrain(implicitmask) = rstats;
    rbrain = transform(rbrain,H,E,C,dh);
    rstats = rbrain(implicitmask);
    if tails == 2
        rstats = abs(rstats);
    end
    
    % compare maxima to t-values and increment as appropriate
    curexceeds = max(rstats) >= tvals;
    exceedances = exceedances + curexceeds;
    
    % update progress indicator (only does so 1 in 5 to minimize overhead)
    if ~randi([0 4]);
        parfor_progress;
    end
    
end

% create corrected p-value image
corrected = exceedances./nperm;
pcorr = ones(bsize);
pcorr(implicitmask) = corrected;

% split into positive and negative effects (if needed)
if tails == 2
	pos = truestat>0;
	pcorr_pos = pcorr;
	pcorr_pos(~pos) = 1;
	pcorr_neg = pcorr;
	pcorr_neg(pos) = 1;
end

% assigne output to varargout
if tails == 1
	varargout{1} = pcorr;
else
	varargout{1} = pcorr_pos;
	varargout{2} = pcorr_neg;
end

end

