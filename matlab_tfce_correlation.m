function [varargout] = matlab_tfce_correlation(imgs,covariate,varargin)
% TFCE_CORRELATION computes threshold free cluster enhancement for a an
% individual difference correlation between a covariate and brain activity.
% The maximal statistic technique is augmented with sequential stepdown to
% maximize power.
%
%   Arguments:
%   imgs -- a 4D (x,y,z,subject) matrix of images. Like tfce_permutation,
%   these should have had tfce_transform applied already.
%   covariate -- a vector of length = number of subjects containing values
%   to be correlated with brain activity
%   nperm -- number of permutations to perform. More permutations yield
%   more precise correct p-values.
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
if tails == 1
	cvals = truestat;
else
	cvals = abs(truestat);
end

% cycle through permutations
exceedances = zeros(nvox,1);
for p = 1:nperm
    
    % permute covariates
    rcov = randsample(covariate,nsub);
    
    % calculate permutation correlations
    rstats = corr(occimgs,rcov);
	if tails == 1
		rcvals = rstats;
	else
		rcvals = abs(rstats);
	end
    
    % compare maxima to t-values and increment as appropriate
    curexceeds = max(rstats) >= cvals;
    exceedances = exceedances + curexceeds;
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

% assigne output to varargout
if tails == 1
	varargout{1} = pcorr;
else
	varargout{1} = pcorr_pos;
	varargout{2} = pcorr_neg;
end

end