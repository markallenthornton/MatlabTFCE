function matlab_tfce_gui()
% MATLAB_TFCE_GUI convenient GUI interface for matlab_tfce
% See matlab_tfce.m for details.

anafig = figure('Name','MatlabTFCE',...
                      'MenuBar','none', ...
                      'Toolbar','none', ...
                      'NumberTitle','off', ...
                      'Resize', 'off', ...
                      'Position',[450,400,400,200]); %X,Y then Width, Height

set(anafig, 'Color', ([128,0,0] ./255)); 

uicontrol('Style','PushButton','HorizontalAlignment','left','String','One-sample t-test','Position',[33,120,150,40],'Callback',{@analysis, 'onesample'});
uicontrol('Style','PushButton','HorizontalAlignment','left','String','Paired-sample t-test','Position',[216,120,150,40],'Callback',{@analysis, 'paired'});
uicontrol('Style','PushButton','HorizontalAlignment','left','String','Indepedent-sample t-test','Position',[33,40,150,40],'Callback',{@analysis, 'independent'});
uicontrol('Style','PushButton','HorizontalAlignment','left','String','Correlation test','Position',[216,40,150,40],'Callback',{@analysis, 'correlation'});

% %% output to nii
% % write out corrected images
% if tails == 1
%     niiout.img = pcorr;
%     save_nii(niiout,[ananame '_' analysis '_pcorr.nii']);
% else
%     niiout.img = pcorr_pos;
%     save_nii(niiout,[ananame '_' analysis '_pcorr_positive.nii']);
%     niiout.img = pcorr_neg;
%     save_nii(niiout,[ananame '_' analysis '_pcorr_negative.nii']);
% end
end