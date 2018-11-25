clear;
clc;

% ------------ INPUTS -------------------

global mcRoot
mcRoot = '/home/slab/users/mangstad/repos/MethodsCorePsych';
addpath(fullfile(mcRoot,'matlabScripts'));

fid = fopen('/net/parasite/saige/ABIDE/Scripts/CPM/sublist','r'); %read in Nx1 subject list
s = textscan(fid,'%s');
SubjDir = s{1};

fid = fopen('/net/parasite/saige/ABIDE/Scripts/CPM/ADOS_total','r'); %read in Nx1 ADOS total score
s2 = textscan(fid,'%f');
ADOS_total = s2{1};

fid = fopen('/net/parasite/saige/ABIDE/Scripts/CPM/meanFD','r'); %read in Nx1 meanFD 
s3 = textscan(fid,'%f');
meanFD = s3{1};

[R_test, P_test] = corr(meanFD, ADOS_total) %check that meanFD & ADOS score are not correlated

CorrTemplate = '/net/parasite/saige/ABIDE/Subjects/[Subject]/Power/Power_corr.mat';

p = 264; %number of ROIs in Power atlas
n = size(SubjDir,1); %number of subjects

%load connectivity matrix for all subjects in sublist varaible (SubjDir)
all_mats = zeros(p,p,n);
for i = 1:size(SubjDir,1)
    Subject = SubjDir{i};
    path = mc_GenPath(CorrTemplate);
    data = load(path);
    z = mc_FisherZ(data.rMatrix);
    z(logical(eye(p))) = 0; %sets diagonal of connectivity matrix = 0, not NaN
    sr = sum(isnan(z)); %lines 37-40 set any column/row of all_mats = 0 if all values in the column/row are NaN. This would mean the subject did not have any data at that node (found 1 subject with this problem due to cerebellum cutoff).
    br = sr>=263;
    z(br,:) = 0;
    z(:,br) = 0;
    
    all_mats(:,:,i) = z;
end

% threshold for feature selection
thresh = 0.01;

% ---------------------------------------
all_behav = ADOS_total;

no_sub = size(all_mats,3);
no_node = size(all_mats,1);

behav_pred_pos = zeros(no_sub,1);
behav_pred_neg = zeros(no_sub,1);


for leftout = 1:no_sub;
    fprintf('\n Leaving out subj # %6.3f',leftout);
    
    % leave out subject from matrices and behavior
    
    train_mats = all_mats;
    train_mats(:,:,leftout) = [];
    train_vcts = reshape(train_mats,[],size(train_mats,3));
    
    train_behav = all_behav;
    train_behav(leftout) = [];
    
    % correlate all edges with behavior

    [r_mat,p_mat] = corr(train_vcts',train_behav);
    vals.p = 1;
    
    FD = meanFD;
    FD(leftout) = [];
    
    stats = mc_CovariateCorrectionFast(train_vcts', [train_behav FD],2,vals); %Controls for motion by using subject's meanFD value from SPM12 realignment preprocessing step as a covariate. 
    r_mat = stats.b(2,:);
    p_mat = stats.p(2,:);
    
    r_mat = reshape(r_mat,no_node,no_node);
    p_mat = reshape(p_mat,no_node,no_node);
    
    % set threshold and define masks
    
    pos_mask = zeros(no_node,no_node);
    neg_mask = zeros(no_node,no_node);
    
    pos_edges = find(r_mat > 0 & p_mat < thresh);
    neg_edges = find(r_mat < 0 & p_mat < thresh);
    
    pos_mask(pos_edges) = 1;
    neg_mask(neg_edges) = 1;
    
    % get sum of all edges in TRAIN subs (divide by 2 to control for the
    % fact that matrices are symmetric)
    
    train_sumpos = zeros(no_sub-1,1);
    train_sumneg = zeros(no_sub-1,1);
    
    for ss = 1:size(train_sumpos);
        train_sumpos(ss) = sum(sum(train_mats(:,:,ss).*pos_mask))/2;
        train_sumneg(ss) = sum(sum(train_mats(:,:,ss).*neg_mask))/2;
    end
    
    % build model on TRAIN subs
    
    fit_pos = polyfit(train_sumpos, train_behav,1);
    fit_neg = polyfit(train_sumneg, train_behav,1);
    
    % run model on TEST sub
    
    test_mat = all_mats(:,:,leftout);
    test_sumpos = sum(sum(test_mat.*pos_mask))/2;
    test_sumneg = sum(sum(test_mat.*neg_mask))/2;
    
    behav_pred_pos(leftout) = fit_pos(1)*test_sumpos + fit_pos(2);
    behav_pred_neg(leftout) = fit_neg(1)*test_sumneg + fit_neg(2);
    
end

% compare predicted and observed scores

[R_pos, P_pos] = corr(behav_pred_pos,all_behav)
[R_neg, P_neg] = corr(behav_pred_neg,all_behav)

figure(1); plot(behav_pred_pos,all_behav,'r.'); lsline
figure(2); plot(behav_pred_neg,all_behav,'b.'); lsline