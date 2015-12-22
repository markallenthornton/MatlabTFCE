%% output to nii
% write out corrected images
if tails == 1
    niiout.img = pcorr;
    save_nii(niiout,[ananame '_' analysis '_pcorr.nii']);
else
    niiout.img = pcorr_pos;
    save_nii(niiout,[ananame '_' analysis '_pcorr_positive.nii']);
    niiout.img = pcorr_neg;
    save_nii(niiout,[ananame '_' analysis '_pcorr_negative.nii']);
end
