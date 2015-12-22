function [varargout] = matlab_tfce_ttest_independent(imgs1,imgs2,tails,varargin)
%TFCE_PERMUTATION_INDEPENDENT tests means difference, independent samples
%   [varargout] = tfce_permutation_independent(imgs1,imgs2,nperm) performs a
%   comparison of the voxelwise image means of two independent groups (as
%   in an independent t-test). This version performs 1-tailed tests. The 
%   alternative hypothesis is imgs1>imgs2. Maximal statistics from tests on
%   permuted data are compared with statistics in the original
%   data using a stepdown procedure described by Holmes, Blair, Watson, and
%   Ford (1996) with an implementation due to Westfall and Young (1993).
%   This procedure is analogous to the more commonly known sequentially
%   rejective test due to Holm (1979). Note that independent sample tests
%   require homogeneity of variance.
%
%   Arguments:
%   imgs1 -- images from group 1 with dimensions x,y,z,nsubject1
%   imgs2 -- images from group 2 with dimensions x,y,z,nsubject2
%	tails -- 1 or 2 tailed test
%   nperm -- number of permutations (1000 default)
%
%   Output:
%	If tails == 1:
%   pcorr -- wholebrain map of corrected p-values
%	If tails == 2:
%	pcorr_pos -- corrected p-values for positive effects
%	pcorr_neg -- corrected p-values for negative effects


% set defaults
nperm = 1000;
if nargin > 1
    nperm = varargin{1};
end

% calculate matrix size
bsize = size(imgs1);
nsub1 = bsize(4);
bsize = size(imgs2);
nsub2 = bsize(4);
bsize = bsize(1:3);
nsub = nsub1+nsub2;

% calculate true mean image
truestat = (mean(imgs1,4)-mean(imgs2,4))/sqrt(var(imgs1,0,4)/nsub1+var(imgs2,0,4)/nsub1);
implicitmask = ~isnan(truestat);

% p-values for comparison
tvals = truestat(implicitmask);
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

% cycle through permutations
exceedances = zeros(nvox,1);
for p = 1:nperm
    
    % permute signs
    relabeling = randsample(glabs,nsub);
    rimgs1 = occimgs(:,relabeling==1);
    rimgs2 = occimgs(:,relabeling==2);
    
    % calculate permutation means
    rstats = (mean(rimgs1,2)-mean(rimgs2,2))/sqrt(var(rimgs1,0,2)/nsub1+var(rimgs2,0,2)/nsub2);
    
    % compare maxima to t-values and increment as appropriate
    if tails == 1
		curexceeds = max(rstats) >= tvals;
	else
		curexceeds = max(abs(rstats)) >= tvals;
    end
    exceedances = exceedances + curexceeds;
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

