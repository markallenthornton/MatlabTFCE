function [varargout] = stepdown_tfce_ttest_onesample(imgs,tails,nperm,H,E,C,dh)
%STEPDOWN_TFCE_TTEST_ONESAMPLE performs maximal statistic permutation testing
%   [varargout] = matlab_tfce_ttest_onesample(imgs,tails,nperm) corrects
%   a one-sample t-test (mean > 0) for multiple comparisons via a 
%   permutation procedure (random sign flipping). Maximal t-stats of
%   permuted data are compared with t-stats in the unshuffled data.
%
%   Arguments:
%   imgs -- a 4D (x,y,z,subject) matrix of images
%	tails -- 1 or 2 tailed test
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
bsize = bsize(1:3);

% set tranform function
if tails == 1
    transform = @stepdown_tfce_transform;
else
    transform = @stepdown_tfce_transform_twotailed;
end

% calculate true mean image
truestat = mean(imgs,4)./(std(imgs,0,4)./sqrt(nsub));
implicitmask = ~isnan(truestat);
truestat = transform(truestat,H,E,C,dh);


% p-values for comparison
tvals = truestat(implicitmask);
if tails == 2
	tvals=abs(tvals);
end
nvox = length(tvals);

% extract occupied voxels for permutation test
occimgs = NaN(nvox,nsub);
for s = 1:nsub
    curimg = imgs(:,:,:,s);
    occimgs(:,s) = curimg(implicitmask);
end

% cycle through permutations
exceedances = zeros(nvox,1);
for p = 1:nperm
    
    % permute signs
    relabeling = randsample([-1 1],nsub,'true');
    roccimgs = occimgs;
    for s = 1:nsub
        if relabeling(s) == -1;
            roccimgs(:,s) = -occimgs(:,s);
        end
    end
    
    % calculate permutation statistic
    rstats =  mean(roccimgs,2)./(std(roccimgs,0,2)./sqrt(nsub));
    rbrain = zeros(bsize);
    rbrain(implicitmask) = rstats;
    rbrain = transform(rbrain,H,E,C,dh);
    rstats = rbrain(implicitmask);

    % compare maxima to true statistic and increment as appropriate
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

