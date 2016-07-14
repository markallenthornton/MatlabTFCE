function [varargout] = matlab_tfce_regression(imgs,preds,tails,nperm,H,E,C,dh,nuisance)
% MATLAB_TFCE_REGRESSION computes TFCE corrected p-values for each
% parameter in a multiple regression
%
%   Arguments:
%   imgs -- a 4D (x,y,z,subject) matrix of images.
%   preds -- an m x n matrix, where m = no. of subjects and n = the number
%   of predictors in the model (the first should be the interecept)
%	tails -- 1 or 2 tailed test
%   nperm -- number of permutations to perform. More permutations yield
%   more precise correct p-values.
%   -- img the 3D image to be transformed
%   -- H height exponent
%   -- E extent exponent
%   -- C connectivity
%   -- ndh step number for cluster formation
%   -- nuisance specifies which variables are not of interest
%   -- omnibusf requests overall model significance test
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
npred = size(preds,2);

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

% calculate true regression images
xsol = (preds'*preds)\preds';
bs = xsol*occimgs;  % undstandardized regression coefficients

yty = NaN(1,nvox);
for i = 1:nvox
    yty(i) = occimgs(:,i)'*occimgs(:,i);
end

bpy = NaN(1,nvox);
for i = 1:nvox
    bpy(i) = bs(:,i)'*preds'*occimgs(:,i);
end

sxy = sqrt((yty-bpy)/(nsub-npred));
sqp = sqrt(inv(preds'*preds));
serr = NaN(npred,nvox);
for i = 1:nvox
    serr(:,i) = diag(sxy(i)*sqp);
end

truestat = bs./serr;
tvals = NaN(npred,nvox);
for i = 1:npred
    trueimg=NaN(bsize);
    trueimg(implicitmask) = truestat(i,:);
    trueimg = transform(trueimg,H,E,C,dh);
    tfcestat = trueimg(implicitmask);
    tvals(i,:) = tfcestat;
    if tails == 2
        tvals(i,:) = abs(tfcestat);
    end  
end

% initialize progress indicator
parfor_progress(nperm);
global parworkers

% cycle through permutations
exceedances = zeros(nvox,npred);
parfor(p = 1:nperm,parworkers)
    
    % permute predictors
    rsel = randperm(nsub);
    rpreds = preds(rsel,:);
    
    % calculate permuted t-values
    xsol = (rpreds'*rpreds)\rpreds';
    bs = xsol*occimgs;

    bpy = NaN(1,nvox);
    for i = 1:nvox
        bpy(i) = bs(:,i)'*rpreds'*occimgs(:,i);
    end
    
    sxy = sqrt((yty-bpy)/(nsub-npred));
    sqp = sqrt(inv(rpreds'*rpreds));
    serr = NaN(npred,nvox);
    for i = 1:nvox
        serr(:,i) = diag(sxy(i)*sqp);
    end

    rtstats = bs./serr;
    curexceeds = zeros(nvox,npred);
    for i = 1:npred
        if ~nuisance(i)
            rbrain=NaN(bsize);
            rbrain(implicitmask) = rtstats(i,:);
            rbrain = transform(rbrain,H,E,C,dh);
            rstats = rbrain(implicitmask);
            if tails == 2
                rstats = abs(rstats);
            end  
            % compare maxima to t-values and increment as appropriate
            curexceeds(:,i) = max(rstats) >= tvals(i,:);
        end
    end
    exceedances = exceedances + curexceeds;
    
    % update progress indicator (only does so 1 in 5 to minimize overhead)
    if ~randi([0 4]);
        parfor_progress;
    end
    
end

if tails == 1
    cormaps = {};
else
    cormaps_pos = {};
    cormaps_neg = {};
end

for i = 1:npred
    if ~nuisance(i)
        % create corrected p-value image
        corrected = exceedances(:,i)./nperm;
        pcorr = ones(bsize);
        pcorr(implicitmask) = corrected;

        % split into positive and negative effects (if needed)
        if tails == 2
            btruestat = NaN(bsize);
            btruestat(implicitmask) = truestat(i,:);
            pos = btruestat>0;
            pcorr_pos = pcorr;
            pcorr_pos(~pos) = 1;
            pcorr_neg = pcorr;
            pcorr_neg(pos) = 1;
        end

        % assign to cormaps
        if tails == 1
            cormaps{i} = pcorr;
        else
            cormaps_pos{i} = pcorr_pos;
            cormaps_neg{i} = pcorr_neg;
        end
    else
        % assign to cormaps
        if tails == 1
            cormaps{i} = NaN;
        else
            cormaps_pos{i} = NaN;
            cormaps_neg{i} = NaN;
        end
    end
end

% assign output to varargout
if tails == 1
	varargout{1} = cormaps;
else
	varargout{1} = cormaps_pos;
	varargout{2} = cormaps_neg;
end

end