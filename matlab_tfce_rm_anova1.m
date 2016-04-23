function pcorr = matlab_tfce_rm_anova1(imgs,nperm,H,E,C,dh)
%MATLAB_TFCE_RM_ANOVA1 performs maximal statistic permutation testing
%   pcorr = matlab_tfce_rm_anova1(imgs,nperm,H,E,C,dh) corrects
%   a one-way omnibus ANOVA (means ~=) for multiple comparisons via a 
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
%   pcorr -- wholebrain map of corrected p-values

% calculate matrix size
bsize = size(imgs);
nsub = bsize(4);
bsize = bsize(1:3);

% set tranform function
if tails == 1
    transform = @matlab_tfce_transform;
else
    transform = @matlab_tfce_transform_twotailed;
end

% calculate true mean image
truestat = mean(imgs,4)./(std(imgs,0,4)./sqrt(nsub));
implicitmask = ~isnan(truestat);
tfcestat = transform(truestat,H,E,C,dh);


% p-values for comparison
tvals = tfcestat(implicitmask);
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
    if tails == 2
        rstats = abs(rstats);
    end
    
    % compare maxima to true statistic and increment as appropriate
    curexceeds = max(rstats) >= tvals;
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

