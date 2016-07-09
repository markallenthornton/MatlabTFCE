# MatlabTFCE
Standalone MATLAB implementation of permutation TFCE correction

This package offers a standalone implemetation of multiple comparison
correction for fMRI data. It achieves this through a permutation testing
approach which controls familywise error rate by comparing voxelwise
statistics to the maximal statistics obtained from repeating the analysis
with randomized data. See Nichols & Holmes (2002) for a detailed
treatment of this approach. 

This maximal permuted statistic correction technique is combined
with the threshold free cluster enhancement (TFCE) transformation due to
Smith & Nichols (2009), which obviates the need for arbitrary voxelwise
cluster-forming thresholds and instead produces continuous correct
p-values for all voxels. Although some spatial specifity is lost
relative to purely voxelwise approach, this approach, like cluster
corrections, is substantially less conservative due to the fact that
it capitalizes on spatial dependency in the data.

All functionality in the package can be accessed via matlab_tfce.m and
that file also provides a description of the relevant parameters. The
file matlab_tfce_gui.m provides a convenient user interface.

The following statistical tests can be computed:

'onesample' -- tests one sample hypothesis mean > 0

'paired' -- paired (dependent samples) test mean(imgs) > mean(imgs2)

'independent' -- independent (two sample) test mean(imgs) > mean(imgs2)

'correlation' -- correlation across subjects of imgs with covariate

'rm_anova1' -- one-factor repeated measures ANOVA

'rm_anova2' -- two-factor (with interaction) repeated measures ANOVA

'regression' -- multiple linear regression

Other than the ANOVAs, all tests permit both one- and two-tailed modes. See the demo script to observe proof of familywise error rate control under trivial conditions. The TFCE transformation itself and the paired/one-sample t-test routines have been validated directly against FSL's randomise using real fMRI data. Default parameter settings were also adopted from FSL.

See here for additional information: http://markallenthornton.com/blog/matlab-tfce/
