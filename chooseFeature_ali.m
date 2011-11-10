% chooseFeature() - choose and format the features for classification
%
% Usage:  [Xtrain,Xtest] = chooseFeature(EEGf,EXP,crossval,varargin)
%
% Requierd input:     meaning (size or options) {default}
%
%   EEGf            = structure containing the features
%   EXP             = structure resluting from classification over several
%                     cross validation
%   crossval        = indice of the crossvalidation
%   
%
% Optional input:     meaning (size or options) {default}
%
%   normalization   = normalize each features ('minMax','off') {'Zscore'}
%   features        = feature to use (see EEGf.featClassif(indBin)) {'mean1'}
%   timepoint       = indice of a particula time point at which you want to
%                     run the classification, by default work on all time
%                     point {'all'}
%   area            = selection electrode from a particular area ('V2','PIT',...) {'all'}
%   indBin          = indice of binsize to take, usefull when EEGf contain features comuted on different bin size {1}
%
% Outputs:            meaning (size or options)
%
%   Xtrain          = matrix for training the classifieur (trials*features)
%   Xtest           = matrix for test the classifieur (trials*features)
%
%
% Autors: Maxime, Serre Lab, 2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Xtrain,Xtest] = chooseFeature_ali(EEGf,EXP,crossval,varargin)

%-----------------------------------------

if nargin < 2
    error('chooseFeature : TooFewInputs Requires an EEGf, EXP, crossval and timepoint');
end

defaults.normalization             = 'Zscore';
defaults.features                  = 'mean';
defaults.area                      = 'all';
defaults.indBin                    = 1;
defaults.timepoint                 = 1:EEGf.pnts;

args                               = propval(varargin,defaults);

normalization                      = args.normalization;
features                           = args.features;
area                               = args.area;
indBin                             = args.indBin;
timepoint                          = args.timepoint;
%-----------------------------------------

%-------------------------%
% Select channels         %
%-------------------------%
if  ~strcmp(area,'all')&&isfield(EXP.chanlocs,'area')&&~isnumeric(area)
    channels = elecByArea (EXP,'area',area);
elseif isnumeric(area)
    channels = area;
else
    channels = 1:size(EEGf.data,1);
end

%-------------------------%
% Select features         %
%-------------------------%
switch features
    case 'all'
        %--------------%
        % all features %
        %--------------%
        Xtrain = [];
        Xtest  = [];
        allFeat=fieldnames(EEGf.featClassif);%{'data','range','contrast','std','der1','der2','der3','raw'};
        for f=1:length(allFeat)
            [XtrainF,XtestF]= chooseFeature_ali(EEGf,EXP,crossval,'timepoint',timepoint,'indBin',indBin,'features',allFeat{f},'area',area,'normalization',normalization);
            % Normalize each feature
            %------------------------
            XtrainF = reshape(zscore(XtrainF(:)),size(XtrainF,1),size(XtrainF,2));
            Xtrain = [Xtrain, XtrainF];
            XtestF = reshape(zscore(XtestF(:)),size(XtestF,1),size(XtestF,2));
            Xtest  = [Xtest, XtestF];
        end
    case {'mean','mean1','mean2','range','contrast','std','der1','der2','der3'...
            'powDelta', 'powTheta', 'powAlpha', 'powBeta', 'powGamma', 'powHigamma', 'phsDelta', 'phsTheta', 'phsAlpha', 'phsBeta', 'phsGamma', 'phsHigamma',...
            'kurto','variance','slope','skew','integral','pow'}
        %------------------------%
        % one value by electrode %
        %------------------------%
        Xtrain = squeeze(getfield(EEGf.featClassif(indBin),features,{channels,timepoint,find(EXP.selector(crossval,:)==1)}));
        Xtrain = reshape(Xtrain,length(channels)*length(timepoint),length(find(EXP.selector(crossval,:)==1)))';
        Xtest  = squeeze(getfield(EEGf.featClassif(indBin),features,{channels,timepoint,find(EXP.selector(crossval,:)==2)}));
        Xtest  = reshape(Xtest,length(channels)*length(timepoint),length(find(EXP.selector(crossval,:)==2)))';
        
        % Normalization
        %--------------
        if strcmp(normalization,'Zscore')
            % Zscore
            %--------
            [Xtrain,Xtest] = normalizeTrainTest(Xtrain,Xtest);
            
        elseif strcmp(normalization,'minMax')
            % MinMax
            %--------
            maxx = max(Xtrain); minxx = min(Xtrain);
            Xtrain  = (Xtrain - repmat(minxx,size(Xtrain,1),1)) ./ repmat(maxx-minxx,size(Xtrain,1),1);
            maxx = max(Xtest); minxx = min(Xtest);
            Xtest  = (Xtest - repmat(minxx,size(Xtest,1),1)) ./ repmat(maxx-minxx,size(Xtest,1),1);
        end
    case {'raw','poly1','poly2','poly3','poly4','poly5','poly6'}
        %-----------------%
        % bin of raw data %
        %-----------------%
        Xtrain = squeeze(getfield(EEGf.featClassif(indBin),features,{channels,1:size(EEGf.featClassif(indBin).(features),2),timepoint,find(EXP.selector(crossval,:)==1)}));
        Xtrain = reshape(Xtrain,length(channels)*size(EEGf.featClassif(indBin).(features),2),length(timepoint),length(find(EXP.selector(crossval,:)==1)));
        Xtrain = reshape(Xtrain,length(channels)*size(EEGf.featClassif(indBin).(features),2)*length(timepoint),length(find(EXP.selector(crossval,:)==1)))';
        Xtest  = squeeze(getfield(EEGf.featClassif(indBin),features,{channels,1:size(EEGf.featClassif(indBin).(features),2),timepoint,find(EXP.selector(crossval,:)==2)}));
        Xtest = reshape(Xtest,length(channels)*size(EEGf.featClassif(indBin).(features),2),length(timepoint),length(find(EXP.selector(crossval,:)==2)));
        Xtest = reshape(Xtest,length(channels)*size(EEGf.featClassif(indBin).(features),2)*length(timepoint),length(find(EXP.selector(crossval,:)==2)))';
        % Normalization
        %--------------
        if strcmp(normalization,'Zscore')
            % Zscore
            %--------
            [Xtrain,Xtest] = normalizeTrainTest(Xtrain,Xtest);
        elseif strcmp(normalization,'minMax')
            % MinMax
            %--------
            maxx = max(Xtrain); minxx = min(Xtrain);
            Xtrain  = (Xtrain - repmat(minxx,size(Xtrain,1),1)) ./ repmat(maxx-minxx,size(Xtrain,1),1);
            maxx = max(Xtest); minxx = min(Xtest);
            Xtest  = (Xtest - repmat(minxx,size(Xtest,1),1)) ./ repmat(maxx-minxx,size(Xtest,1),1);
        end
    case {'PW','Phase'}
        %-------------------------%
        % spectrale analysis      %
        %-------------------------%
    otherwise
        Xtrain = squeeze(getfield(EEGf.featClassif(indBin),features,{channels,timepoint,find(EXP.selector(crossval,:)==1)}));
        Xtrain = reshape(Xtrain,length(channels)*length(timepoint),length(find(EXP.selector(crossval,:)==1)))';
        Xtest  = squeeze(getfield(EEGf.featClassif(indBin),features,{channels,timepoint,find(EXP.selector(crossval,:)==2)}));
        Xtest  = reshape(Xtest,length(channels)*length(timepoint),length(find(EXP.selector(crossval,:)==2)))'; 
                % Normalization
        %--------------
        if strcmp(normalization,'Zscore')
            % Zscore
            %--------
            [Xtrain,Xtest] = normalizeTrainTest(Xtrain,Xtest);
            
        elseif strcmp(normalization,'minMax')
            % MinMax
            %--------
            maxx = max(Xtrain); minxx = min(Xtrain);
            Xtrain  = (Xtrain - repmat(minxx,size(Xtrain,1),1)) ./ repmat(maxx-minxx,size(Xtrain,1),1);
            maxx = max(Xtest); minxx = min(Xtest);
            Xtest  = (Xtest - repmat(minxx,size(Xtest,1),1)) ./ repmat(maxx-minxx,size(Xtest,1),1);
        end
end
        
        %display(['-------- Training: ' num2str(size(Xtrain,1))  ' trials for ' num2str(size(Xtrain,2)) ' features --------'])
        %display(['-------- Test: ' num2str(size(Xtest,1)) ' trials for ' num2str(size(Xtest,2)) ' features --------'])
end