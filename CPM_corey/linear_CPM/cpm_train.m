function [r,p,pmask,mdl]=cpm_train(x,y,pthresh,varargin)
% Train a Connectome-based Predictive Model
% x            Predictor variable
% y            Outcome variable
% pthresh      p-value threshold for feature selection
% r            Correlations between all x and y
% p            p-value of correlations between x and y
% pmask        Mask for significant features
% mdl          Coefficient fits for linear model relating summary features to y

% Parse inputs
p=inputParser;

addRequired(p,'x',@isnumeric);
addRequired(p,'y',@isnumeric); % must be n x nsubs
addRequired(p,'pthresh',@isnumeric);
addParameter(p,'corrtype',@ischar);
addParameter(p,'z',@isnumeric);

parse(p,x,y,pthresh,varargin{:});

corrtype=p.Results.corrtype;

if strcmp(corrtype, 'partial')
    z=p.Results.z;
end


clearvars p


switch corrtype
    
    case 'pearsons'
        [r,p]=corr(x',y);
    case 'partial'
        [r,p]=partialcorr(x',y,z);
end
        

% Select significant features

pmask=(+(r>0))-(+(r<0)); 
pmask=pmask.*(+(p<pthresh));

% For each subject, summarize selected features 
for i=1:size(x,2)
    summary_feature(i)=nanmean(x(pmask>0,i))-nanmean(x(pmask<0,i));
end

% Fit y to summary features
mdl=robustfit(summary_feature,y');    