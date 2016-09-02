function matlab_tfce_gui()
% MATLAB_TFCE_GUI convenient GUI interface for matlab_tfce
% See matlab_tfce.m for details.

% generate analysis selection figure
anafig = figure('Name','MatlabTFCE',...
                      'MenuBar','none', ...
                      'Toolbar','none', ...
                      'NumberTitle','off', ...
                      'Resize', 'off', ...
                      'Position',[450,400,360,270]); %X,Y then Width, Height
set(anafig, 'Color', ([128,0,0] ./255)); 
uipanel('Title','Select analysis',...
             'BackgroundColor','white',...
             'Units','Pixels',...
             'Position',[10 10 340 260]);
uicontrol('Style','PushButton','HorizontalAlignment','left','String','One-sample t-test','Position',[20,200,150,40],'Callback',{@anadiag, 'onesample'});
uicontrol('Style','PushButton','HorizontalAlignment','left','String','Paired-sample t-test','Position',[190,200,150,40],'Callback',{@anadiag, 'paired'});
uicontrol('Style','PushButton','HorizontalAlignment','left','String','Indepedent-sample t-test','Position',[20,140,150,40],'Callback',{@anadiag, 'independent'});
uicontrol('Style','PushButton','HorizontalAlignment','left','String','Correlation test','Position',[190,140,150,40],'Callback',{@anadiag, 'correlation'});
uicontrol('Style','PushButton','HorizontalAlignment','left','String','One-way RM ANOVA','Position',[20,80,150,40],'Callback',{@anadiag, 'rm_anova1'});
uicontrol('Style','PushButton','HorizontalAlignment','left','String','Two-way RM ANOVA','Position',[190,80,150,40],'Callback',{@anadiag, 'rm_anova2'});
uicontrol('Style','PushButton','HorizontalAlignment','left','String','Multiple linear regression','Position',[105,20,150,40],'Callback',{@anadiag, 'regression'});
end

function anadiag(~,~,analysis)
    close all;
    % select input images
    if strcmp(analysis,'independent') || strcmp(analysis,'paired')
        [filenames1, pathnames1] = uigetfile({'*.nii';'*.img'},'Select the 1st set of images (in order)','MultiSelect', 'on');
        [filenames2, pathnames2] = uigetfile({'*.nii';'*.img'},'Select the 2nd set of images (in order)','MultiSelect', 'on');
        if ischar(filenames1) || ischar(filenames2)
            error('Must specify multiple files - one per participant');
        end
        [imgs1,niiout] = imgreader(filenames1,pathnames1);
        [imgs2,~] = imgreader(filenames2,pathnames2);
        imgs = {imgs1,imgs2};
    else
        if strcmp(analysis,'correlation') || strcmp(analysis,'regression')
            [filenames, pathnames] = uigetfile({'*.nii';'*.img'},'Select images (in order of covariate)','MultiSelect', 'on');
            [covname, covpath] = uigetfile('*.csv','Select csv with covariate(s)','MultiSelect', 'off');
            if ~(ischar(covname) && ischar(covpath))
                error('Only one covariate file permitted');
            end
            covariate = csvread([covpath covname]);
        elseif strcmp(analysis,'rm_anova1')
            anovaui = figure('Name','MatlabTFCE',...
                      'MenuBar','none', ...
                      'Toolbar','none', ...
                      'NumberTitle','off', ...
                      'Resize', 'off', ...
                      'Position',[450,400,220,110]); %X,Y then Width, Height
            set(anovaui, 'Color', ([128,0,0] ./255)); 
            yp = 20;
            uipanel('Title','Number of levels in factor',...
                         'BackgroundColor','white',...
                         'Units','Pixels',...
                         'Position',[30 40+yp 160 40]);
            yp = 15;
            u1 = uicontrol('Style','edit','HorizontalAlignment','left','String','3','Position',[50,50+yp,120,20]);
            uicontrol('Style','PushButton','HorizontalAlignment','center','String','Next','Position',[30,10,160,40],'Callback',{@anovacallback,[u1]});
            uiwait(gcf); 
            close(anovaui);
            global levels
            filenames = cell(levels,1);
            pathnames = cell(levels,1);
            for i = 1:levels
                [filenames{i}, pathnames{i}] = uigetfile({'*.nii';'*.img'},['Select images for ANOVA cell (' num2str(i) ')'],'MultiSelect', 'on');
            end
        elseif strcmp(analysis,'rm_anova2')
            anovaui = figure('Name','MatlabTFCE',...
                      'MenuBar','none', ...
                      'Toolbar','none', ...
                      'NumberTitle','off', ...
                      'Resize', 'off', ...
                      'Position',[450,400,220,160]); %X,Y then Width, Height
            set(anovaui, 'Color', ([128,0,0] ./255)); 
            yp = 20;
            uipanel('Title','Number of levels in factor 1',...
                         'BackgroundColor','white',...
                         'Units','Pixels',...
                         'Position',[30 90+yp 160 40]);
            uipanel('Title','Number of levels in factor 2',...
                         'BackgroundColor','white',...
                         'Units','Pixels',...
                         'Position',[30 40+yp 160 40]);
            yp = 15;
            u1 = uicontrol('Style','edit','HorizontalAlignment','left','String','2','Position',[50,100+yp,120,20]);
            u2 = uicontrol('Style','edit','HorizontalAlignment','left','String','2','Position',[50,50+yp,120,20]);
            uicontrol('Style','PushButton','HorizontalAlignment','center','String','Next','Position',[30,10,160,40],'Callback',{@anovacallback,[u1 u2]});
            uiwait(gcf); 
            close(anovaui);
            global levels
            for i = 1:levels(2)
                for j = 1:levels(1)
                    [filenames{j,i}, pathnames{j,i}] = uigetfile({'*.nii';'*.img'},['Select images for ANOVA cell (' num2str(j) ',' num2str(i) ')'],'MultiSelect', 'on');
                end
            end
        else
            [filenames, pathnames] = uigetfile({'*.nii';'*.img'},'Select images','MultiSelect', 'on');
        end
        if ischar(filenames)
            error('Must specify multiple files - one per participant');
        end
        [imgs,niiout] = imgreader(filenames,pathnames);
    end
if ~exist('covariate', 'var')
    covariate = [];
end
 
% ui appearance
inctails = ~(strcmp(analysis,'rm_anova1') || strcmp(analysis,'rm_anova2'));
anafig = figure('Name','MatlabTFCE',...
                      'MenuBar','none', ...
                      'Toolbar','none', ...
                      'NumberTitle','off', ...
                      'Resize', 'off', ...
                      'Position',[450,450,220,470]); %X,Y then Width, Height
set(anafig, 'Color', ([128,0,0] ./255)); 
yp = 20;

uipanel('Title','Contrast name',...
             'BackgroundColor','white',...
             'Units','Pixels',...
             'Position',[30 390+yp 160 40]);
uipanel('Title','Tails of test',...
             'BackgroundColor','white',...
             'Units','Pixels',...
             'Position',[30 340+yp 160 40]);
uipanel('Title','Permutation number',...
             'BackgroundColor','white',...
             'Units','Pixels',...
             'Position',[30 290+yp 160 40]);
uipanel('Title','TFCE H parameter',...
             'BackgroundColor','white',...
             'Units','Pixels',...
             'Position',[30 240+yp 160 40]);
uipanel('Title','TFCE E parameter',...
             'BackgroundColor','white',...
             'Units','Pixels',...
             'Position',[30 190+yp 160 40]);
uipanel('Title','Connectivity',...
             'BackgroundColor','white',...
             'Units','Pixels',...
             'Position',[30 140+yp 160 40]);
uipanel('Title','TFCE dh parameter',...
             'BackgroundColor','white',...
             'Units','Pixels',...
             'Position',[30 90+yp 160 40]);
uipanel('Title','Size of parallel pool',...
             'BackgroundColor','white',...
             'Units','Pixels',...
             'Position',[30 40+yp 160 40]);
yp = 15;

% ui controls
u1 = uicontrol('Style','edit','HorizontalAlignment','left','String','new_contrast','Position',[50,400+yp,120,20]);
if inctails
    u2 = uicontrol('Style','popupmenu','HorizontalAlignment','left','String',{'1','2'},'Position',[50,350+yp,120,20]);
else
    u2 = uicontrol('Style','popupmenu','HorizontalAlignment','left','String',{'1'},'Position',[50,350+yp,120,20]);
end
u3 = uicontrol('Style','edit','HorizontalAlignment','left','String','5000','Position',[50,300+yp,120,20]);
u4 = uicontrol('Style','edit','HorizontalAlignment','left','String','2','Position',[50,250+yp,120,20]);
u5 = uicontrol('Style','edit','HorizontalAlignment','left','String','0.5','Position',[50,200+yp,120,20]);
u6 = uicontrol('Style','popupmenu','HorizontalAlignment','left','String',{'26','18','6'},'Position',[50,150+yp,120,20]);
u7 = uicontrol('Style','edit','HorizontalAlignment','left','String','0.1','Position',[50,100+yp,120,20]);
u8 = uicontrol('Style','edit','HorizontalAlignment','left','String','0','Position',[50,50+yp,120,20]);
uicontrol('Style','PushButton','HorizontalAlignment','center','String','Begin','Position',[30,10,160,40],'Callback',{@run_tfce,[u1 u2 u3 u4 u5 u6 u7,u8],analysis,niiout,imgs,covariate});
    
end

function [imgs,niiout] = imgreader(filenames,pathnames)
% reads in nifti images
if iscell(filenames{1})
    levsize = size(filenames);
    imgs = cell(levsize);
    nsub = length(filenames{1});
    for i = 1:levsize(1)
        for j = 1:levsize(2)
            curfiles = filenames{i,j};
            curpath = pathnames{i,j};
            for s = 1:nsub
                curnii = load_nii([curpath curfiles{s}]);
                if s == 1
                    curimgs = NaN([size(curnii.img) nsub]);
                    niiout = curnii;
                end
                curimgs(:,:,:,s) = curnii.img;
            end
            imgs{i,j} = curimgs;
        end
    end  
else
    nsub = length(filenames);
    if ischar(pathnames)
        for s = 1:nsub
            curnii = load_nii([pathnames filenames{s}]);
            if s == 1
                imgs = NaN([size(curnii.img) nsub]);
                niiout = curnii;
            end
            imgs(:,:,:,s) = curnii.img;
        end
    elseif length(pathnames) == length(filenames)
        for s = 1:nsub
            curnii = load_nii([pathnames{s} filenames{s}]);
            if s == 1
                imgs = NaN([size(curnii.img) nsub]);
                niiout = curnii;
            end
            imgs(:,:,:,s) = curnii.img;
        end
    else
        error('Files must all originate from same directory or a different directory for each file');
    end
end
end

function anovacallback(~,~,uis)
strings = arrayfun(@(x) get(x,'String'),uis,'UniformOutput',false);
global levels
levels = arrayfun(@(x) str2double(x),strings);
uiresume(gcf);
end

function run_tfce(~,~,uis,analysis,niiout,imgs,covariate)

% arrange gui outputs
strings = arrayfun(@(x) get(x,'String'),uis,'UniformOutput',false);
values = arrayfun(@(x) get(x,'Value'),uis);
ananame = strings{1};
tails = str2double(strings{2}{values(2)});
nperm = str2double(strings{3});
H = str2double(strings{4});
E = str2double(strings{5});
C = str2double(strings{6}{values(6)});
dh = str2double(strings{7});
parworkers = str2double(strings{8});
if ~(strcmp(analysis,'rm_anova1') || strcmp(analysis,'rm_anova2')) && length(imgs) == 2 
    imgs2 = imgs{2};
    imgs = imgs{1};
else
    imgs2 = [];
end
close all

% run analysis
disp(['Analysis ' ananame ' ' analysis ' started at ' datestr(now)]);
if strcmp(analysis,'rm_anova2')
    [pcorr_fac1,pcorr_fac2,pcorr_int] = matlab_tfce(analysis,tails,imgs,imgs2,covariate,nperm,H,E,C,dh,parworkers);
elseif tails == 1
    pcorr = matlab_tfce(analysis,tails,imgs,imgs2,covariate,nperm,H,E,C,dh,parworkers);
else
    [pcorr_pos,pcorr_neg] = matlab_tfce(analysis,tails,imgs,imgs2,covariate,nperm,H,E,C,dh,parworkers);
end
disp(['Analysis ' ananame ' ' analysis ' completed at ' datestr(now)]);

% write output
disp('Writing output to file...')
if strcmp(analysis,'rm_anova2')
    fname = [ananame '_' analysis '_pcorr_factor1.nii'];
    write_res(fname,niiout,pcorr_fac1);
    fname = [ananame '_' analysis '_pcorr_factor2.nii'];
    write_res(fname,niiout,pcorr_fac2);
    fname = [ananame '_' analysis '_pcorr_interaction.nii'];
    write_res(fname,niiout,pcorr_int);
elseif tails == 1
    fname = [ananame '_' analysis '_pcorr.nii'];
    write_res(fname,niiout,pcorr);
elseif strcmp(analysis,'regression')
    for i = 1:size(covariate,2)
        if tails == 1
            fname = [ananame '_' analysis '_b' num2str(i-1) '_pcorr.nii'];
            write_res(fname,niiout,pcorr{i})
        else
            fname = [ananame '_' analysis '_b' num2str(i-1) '_pcorr_positive.nii'];
            write_res(fname,niiout,pcorr_pos{i})
            fname = [ananame '_' analysis '_b' num2str(i-1) '_pcorr_negative.nii'];
            write_res(fname,niiout,pcorr_neg{i})
        end
    end
else
    fname = [ananame '_' analysis '_pcorr_positive.nii'];
    write_res(fname,niiout,pcorr_pos);
    fname = [ananame '_' analysis '_pcorr_negative.nii'];
    write_res(fname,niiout,pcorr_neg);
end
disp('Writing complete!')

end

function write_res(fname,niiout,pcorr)
% output writing function (1-corrected p-value)
niiout.img = 1-pcorr;
save_nii(niiout,fname);
end