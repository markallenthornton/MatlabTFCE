function [varargout] = matlab_tfce_correlation(imgs,covariate,tails,nperm,H,E,C,dh)
% MATLAB_TFCE_CORRELATION computes TFCE corrected p-values for
% individual difference correlation between a covariate and brain activity.
% Note that actual inference is performed on Fisher r-to-z transforms of
% the Pearson correlation coefficents for linearity.
%
%   Arguments:
%   imgs -- a 4D (x,y,z,subject) matrix of images.
%   covariate -- a vector of length = number of subjects containing values
%   to be correlated with brain activity
%	tails -- 1 or 2 tailed test
%   nperm -- number of permutations to perform. More permutations yield
%   more precise correct p-values.
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
bsize = size(imgs);
nsub = bsize(4);
covariate = covariate(:);
bsize = bsize(1:3);

% set tranform function
if tails == 1
    transform = @matlab_tfce_transform;
else
    transform = @matlab_tfce_transform_twotailed;
end

% calculate implicit mask
sumimg = sum(imgs,4);
implicitmask = ~isnan(sumimg) & sumimg~=0;
nvox = sum(implicitmask(:));

% extract occupied voxels for permutation test
occimgs = NaN(nsub,nvox);
for s = 1:nsub
    curimg = imgs(:,:,:,s);
    occimgs(s,:) = curimg(implicitmask);
end

% calculate true correlation image
truestat = corr(occimgs,covariate);
truestat =.5*log((1+truestat)./(1-truestat));
trueimg=NaN(bsize);
trueimg(implicitmask) = truestat;
trueimg = transform(trueimg,H,E,C,dh);
tfcestat = trueimg(implicitmask);
cvals = tfcestat;
if tails == 2
	cvals = abs(tfcestat);
end

% initialize progress indicator
parfor_progress(nperm);
global parworkers

% cycle through permutations
exceedances = zeros(nvox,1);
parfor(p = 1:nperm,parworkers)
    
    % permute covariates
    rsel = randperm(nsub);
    rcov = covariate(rsel);
    
    % calculate permutation correlations
    rstats = corr(occimgs,rcov);
    rstats =.5*log((1+rstats)./(1-rstats));
    rbrain = zeros(bsize);
    rbrain(implicitmask) = rstats;
    rbrain = transform(rbrain,H,E,C,dh);
    rstats = rbrain(implicitmask);
    if tails == 2
        rstats = abs(rstats);
    end
    
    % compare maxima to t-values and increment as appropriate
    curexceeds = max(rstats) >= cvals;
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
    btruestat = NaN(bsize);
    btruestat(implicitmask) = truestat;
	pos = btruestat>0;
	pcorr_pos = pcorr;
	pcorr_pos(~pos) = 1;
	pcorr_neg = pcorr;
	pcorr_neg(pos) = 1;
end

% assign output to varargout
if tails == 1
	varargout{1} = pcorr;
else
	varargout{1} = pcorr_pos;
	varargout{2} = pcorr_neg;
end

end