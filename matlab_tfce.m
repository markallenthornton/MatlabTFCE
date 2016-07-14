function [varargout] = matlab_tfce(analysis,tails,imgs,varargin)
%MATLAB_TFCE general wrapper for specifying analyses. Receives arguments
% from matlab_tfce_gui.m if specificed manually. All intended
% functionality in the package can be accessed via this function or the gui
% input version. This matlab_tfce and the functions it calls are all
% standalone - i.e. they do not rely on functions from other packages. In
% contrast, matlab_tfce_gui uses (included) functions from the 
% 'NIfTI and ANALYZE tools' package to facilitate file io and the gui.
% Using matlab_tfce directly allows for headless sessions and more
% customization of the file io.
%
% This package offers a standalone implemetation of multiple comparison
% correction for fMRI data. It achieves this through a permutation testing
% approach which controls familywise error rate by comparing voxelwise
% statistics to the maximal statistics obtained from repeating the analysis
% with randomized data. See Nichols & Holmes (2002) for a detailed
% treatment of this approach. 
%
% The maximal statistic technique is combined
% with the threshold free cluster enhancement (TFCE) transformation due to
% Smith & Nichols (2009), which obviates the need for arbitrary voxelwise
% cluster-forming thresholds and instead produces continuous correct
% p-values for all voxels. Although some spatial specifity is lost
% relative to purely voxelwise approach, this approach, like cluster
% corrections, is substantially less conservative due to the fact that
% it capitalizes on spatial dependency in the data. 
%
% [varargout] = matlab_tfce(analysis,tails,imgs,imgs2,covariate,nperm,H,E,C,dh,parworkers,nuisance)
% [pcorr] = matlab_tfce(analysis,1,imgs,imgs2,covariate,nperm,H,E,C,dh,parworkers,nuisance)
% [pcorr_pos,pcorr_neg] = matlab_tfce(analysis,2,imgs,imgs2,covariate,nperm,H,E,C,dh,parworkers,nuisance)
% [pcorr_fac1,pcorr_fac2,pcorr_int] = matlab_tfce(analysis,1,imgs,imgs2,covariate,nperm,H,E,C,dh,parworkers,nuisance)
%
% Arguments:
%
% analysis -- type of analysis to perform. Options include:
%   -- 'onesample' -- tests one sample hypothesis mean > 0
%   -- 'paired' -- paired (dependent samples) test mean(imgs) > mean(imgs2)
%   -- 'independent' -- independent (two sample) test mean(imgs) > mean(imgs2)
%   -- 'correlation' -- correlation across subjects of imgs with covariate
%   -- 'rm_anova1' -- one-factor repeated measures ANOVA
%   -- 'rm_anova2' -- two-factor repeated measures ANOVA
%   -- 'regression' -- multiple linear regression with covariate matrix
%
% tails -- specify a 1 or 2 tailed test (unidirectional or bidirectional)
% that can be combined with t-tests, correlations, and regressions. Ignored
% for repeated measures ANOVAs.
%
% imgs -- a 4D matrix of imaging data for analysis. Dimensions are expected
% to be x,y,z,subject. Alternatively, for repeated measures ANOVAs, a cell
% array in which each cell contains the 4D (x,y,z,subject) matrix for one
% experimental cell in the ANOVA. Should be m x 1 for one-factor ANOVAs, or
% m x n for two-factor ANOVAs, where m and n are the number of factor
% levels of the first and second factor respectively.
%
% Optional arguments (can supply [] to skip):
%
% imgs2 -- a second 4D matrix as above, required for paired and independent
% analysis options. Must have same xyz dimensions as imgs, and if the
% analysis is paired, subject number must match as well. Note that this 
% argument is not used for repeated measures anovas.
%
% covariate -- If analysis type is 'correlation', a subject x 1 matrix 
% containing an individual difference covariate for correlation across 
% subjects with voxelwise activity. If analysis is 'regression', a
% subject x predictor matrix consisting of multiple covariates. First
% column must be constant 1s for intercept. Compatible with categorical
% predictors (e.g. between-subjects ANOVA), but coding/contrasts are not
% handled internally and must be built into predictor matrix a priori.
% 
% nperm -- number of permutations to perform. default = 5000
%
% H -- height exponent, default = 2
%
% E -- extent exponent, default = 0.5
%
% C -- connectivity, default = 26 (6 = surface, 18 = edge, 26 = corner)
%
% dh -- step size for cluster formation, default = .1
%
%   More steps will be more precise but will require more time and memory.
%   The H & E default parameter settings match FSL's randomise/fslmaths.
%   The C default setting matches FSL's ramdomise default setting. To
%   match SPM's default cluster forming, use 18 instead.
%
% parworkers -- number of parallel workers to use, default = 0 (no
% parallelization). If a parallel pool is already open, its size will
% override this argument. If a new pool is opened for MatlabTFCE, it will
% be closed at the end of script. Opening and closing parallel pools is
% time consuming so it may be more efficient to do outside the package if
% one is running several analyses.
%
% nuisance -- commandline-only argument for specifying which variables in
% a multiple regression are nuisance regressors. Accepts a (1 x predictor)
% binary (0/1 true/false) vector. 1s/trues indicate variables that should
% be treated as nuisance predictors (vs. variables of interest). These
% variables will be included in the model, but will not be subjected to
% TFCE, accelerating calculations. Default is all 0s (i.e. all variables
% treated as of non-nuisance). Ignorned if analysis not regression.
%
% Output: 
% If tails == 1, a single output image with the same xyz dimensions as imgs
% consisting of corrected p-values with be returned.
%
% If tails == 2, two such output images will be returned, one for the
% 'positive' tail and one for the 'negative' tail of the test,
% respectively.
%
% If analysis == rm_anova2, three output images will be produced, one for
% each main effect and the interaction term. The first main effect
% corresponds to the rows of the input cell array, and the second main
% effect corresponds to the columns of the input cell array.
%
% If analysis == regression, output will be 1 or 2 cell arrays (depending
% on 'tails'). If 1 tail is requested, the array will consist of images
% equal in number to the predictors in the 'covariate' matrix (and in the
% same order). Each corrected p-value image will represent the hypothesis
% that the given coefficient > 0. If 2 tails are requested, the first will
% reflect the same as the above, and the second will reflect the opposite 
% direction of the test (coefficients < 0).
%
% Note that for convenience, if using matlab_tfce_gui.m, result images will
% be written out as 1-pcorr instead.


%% input checks and default setting

% setting defaults
imgs2 = [];
covariate = [];
nperm = 5000;
H = 2;
E = .5;
C = 26;
dh = .1;
global parworkers;
parworkers = 0;

% adjusting optional arguments based on input
fixedargn = 3;
if nargin > (fixedargn + 0)
    imgs2 = varargin{1};
end
if nargin > (fixedargn + 1)
    covariate = varargin{2};
end
nuisance = zeros([1 size(covariate,2)]);
if nargin > (fixedargn + 2)
    if ~isempty(varargin{3})
        nperm = varargin{3};
    end
end
if nargin > (fixedargn + 3)
    if ~isempty(varargin{4})
        H = varargin{4};
    end
end
if nargin > (fixedargn + 4)
    if ~isempty(varargin{5})
        E = varargin{5};
    end
end
if nargin > (fixedargn + 5)
    if ~isempty(varargin{6})
        C = varargin{6};
    end
end
if nargin > (fixedargn + 6)
    if ~isempty(varargin{7})
        dh = varargin{7};
    end
end
if nargin > (fixedargn + 7)
    if ~isempty(varargin{8})
        parworkers = varargin{8};
    end
end
if nargin > (fixedargn + 8)
    if ~isempty(varargin{9})
        nuisance = varargin{9};
    end
end

% check that tails are appropriate
if ~(sum(tails==[1 2]))
    error('Inappropriate number of tails (must be 1 or 2)');
end

% check class and size of input images
bsize = size(imgs);
if iscell(imgs)
    if ~(strcmp(analysis,'rm_anova1') || strcmp(analysis,'rm_anova2'))
        error('Only ANOVAs expect cell array imgs input')
    end
    if strcmp(analysis,'rm_anova1')
        if bsize(2)~=1
            if bsize(1) == 1
                warning('m x 1 cell array expected for one-way anova; correcting orientation')
                imgs = imgs';
            else
                error('m x 1 cell array expected for one-way anova')
            end 
        end
    else
        if sum(bsize==1)>0
            error('m x n (m & n >2) cell array expected for two-way anova')
        end
    end
    if sum(sum(cellfun(@(x) length(size(x)),imgs)~=4))>0
        error('Image data not 4D - each cell must be x-y-z-subject')
    end
elseif (strcmp(analysis,'rm_anova1') || strcmp(analysis,'rm_anova2'))
    error('imgs must be cell array for ANOVAs')
elseif length(bsize) ~= 4
    error('Image data not 4D - must be x-y-z-subject')
elseif bsize(4) < 13
    warning('Low N limits number of unique permutations. Approximate (vs. exact) permutation may not be appropriate.');
    nsub = bsize(4);
else
    nsub = bsize(4);
end

% check properties of imgs2
if ~isempty(imgs2)
    imgs1 = imgs;
    bsize2 = size(imgs2);
    if ~(strcmp(analysis,'independent') || strcmp(analysis,'paired'))
        warning(['The analysis ' analysis ' ignores the imgs2 argument']);
    end
    if sum(bsize(1:3) == bsize(1:3)) ~= 3
        error('XYZ dimensions of imgs and imgs2 do not match');
    end
    nsub = bsize(4) + bsize2(4);
else
    if (strcmp(analysis,'independent') || strcmp(analysis,'paired'))
        error('The imgs2 argument must be specified for this analysis type.');
    end
end

% check the subject number is the same for paired tests
if strcmp(analysis,'paired')
    if bsize(4) ~= bsize2(4)
        error('The 4th (subject number) dimension of imgs1 and imgs2 must match for paired tests');
    else
        nsub = bsize(4);
    end
end

% check the subject number is the same for repeated measures anovas
if strcmp(analysis,'rm_anova1') || strcmp(analysis,'rm_anova2')
    nsub = cellfun(@(x) size(x,4),imgs);
    if length(unique(nsub)) > 1
        error('Not an equal number of images in each cell')
    end
    nsub = nsub(1);
end

% check covariate
if strcmp(analysis,'correlation')
    covariate = covariate(:);
    if length(covariate) ~= bsize(4)
        error('Covariate length does not equal 4th dimension of images');
    else
        nsub = bsize(4);
    end
elseif strcmp(analysis,'regression')
    if size(covariate,1) ~= bsize(4)
        error('Covariate length does not equal 4th dimension of images');
    else
        nsub = bsize(4);
    end  
end

% check if parallel pool exists; use it if it does, create it if not (but
% is requested).
pp = gcp('nocreate');
try 
    parworkers = pp.NumWorkers;
    parexist = 1;
    warning('Existing parallel pool discovered with %d workers. Overriding parworkers argument to use this pool.',parworkers);
catch
    parexist = 0;
    if parworkers > 0
        parpool(parworkers);
    end
end

% check nuisance
if strcmp(analysis,'regression')
    ns = size(nuisance);
    if ~isequal(ns,[1 size(covariate,2)])
        error('nuisance size does not equal predictor number');
    end
end

%% analysis calls

% indicate analysis start
fprintf('Starting %i-tailed %s analysis...\n',[tails,analysis]);
fprintf('MatlabTFCE started at %s \n',datestr(now));
fprintf('Permuting data from %i subjects %i times\n',[nsub,nperm])


% select appropriate analysis
switch analysis
    
    % one sample test (mean > 0)
    case 'onesample'
        if tails == 1
            pcorr = matlab_tfce_ttest_onesample(imgs,tails,nperm,H,E,C,dh);
        else
            [pcorr_pos,pcorr_neg]= matlab_tfce_ttest_onesample(imgs,tails,nperm,H,E,C,dh); 
        end
    
    % paired (repeated measures) test (imgs1>imgs2)
    case 'paired'
        imgs = imgs1-imgs2;
        if tails == 1
            pcorr = matlab_tfce_ttest_onesample(imgs,tails,nperm,H,E,C,dh);
        else
            [pcorr_pos,pcorr_neg]= matlab_tfce_ttest_onesample(imgs,tails,nperm,H,E,C,dh); 
        end
    
    % independent (two sample) samples test (imgs1>imgs2)
    case 'independent'
        if tails == 1
            pcorr = matlab_tfce_ttest_independent(imgs1,imgs2,tails,nperm,H,E,C,dh);
        else
            [pcorr_pos,pcorr_neg] = matlab_tfce_ttest_independent(imgs1,imgs2,tails,nperm,H,E,C,dh);
        end
    
    % covariate-img correlation (R>0)
    case 'correlation'
        if tails == 1
            pcorr = matlab_tfce_correlation(imgs,covariate,tails,nperm,H,E,C,dh);
        else
            [pcorr_pos,pcorr_neg] = matlab_tfce_correlation(imgs,covariate,tails,nperm,H,E,C,dh);
        end
        
    % repeated measure (within subject) omnibus one-way ANOVA
    case 'rm_anova1'
        pcorr = matlab_tfce_rm_anova1(imgs,nperm,H,E,C,dh);
    
    % repeated measure (within subject) two-way factorial ANOVA
    case 'rm_anova2'
        [pcorr_fac1,pcorr_fac2,pcorr_int] = matlab_tfce_rm_anova2(imgs,nperm,H,E,C,dh);
        
    % multiple regression 
    case 'regression'
        if tails == 1
            pcorr = matlab_tfce_regression(imgs,covariate,tails,nperm,H,E,C,dh,nuisance);
        else
            [pcorr_pos,pcorr_neg] = matlab_tfce_regression(imgs,covariate,tails,nperm,H,E,C,dh,nuisance);
        end
        
    % unrecognized analysis input
    otherwise
        error('Analysis type not recognized');
end

%% assign output
if strcmp(analysis,'rm_anova2')
    varargout{1} = pcorr_fac1;
    varargout{2} = pcorr_fac2;
    varargout{3} = pcorr_int;
elseif strcmp(analysis,'rm_anova1')
    varargout{1} = pcorr;
elseif tails == 1
    varargout{1} = pcorr;
else
    varargout{1} = pcorr_pos;
    varargout{2} = pcorr_neg;
end
fprintf('Permutation complete!\n');
fprintf('MatlabTFCE finished at %s \n',datestr(now));

%% clean up parallel pool if created in-script
if ~parexist
    delete(gcp('nocreate'));
end
end

