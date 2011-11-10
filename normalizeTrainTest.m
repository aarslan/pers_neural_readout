% normalizeTrainTest() - normalize train set and testing set by features using the mean and sd
%                        calculated on the train set
%
% Usage:  [XtrainNorm,XtestNorm] = normalizeTrainTest(Xtrain,Xtest)
%
% Requierd input:     meaning (size or options) {default}
%
%   Xtrain          = matrix for training the classifier (trials*features)
%   Xtest           = matrix for testing the classifier (trials*features)
%
% Outputs:            meaning (size or options)
%
%   XtrainNorm      = normalized matrix for training the classifier (trials*features)
%   XtestNorm       = normalized matrix for testing the classifier (trials*features)
%
%
% Modified by Sebastien, Serre lab, 2011
% Authors: Maxime, cerco, 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [XtrainNorm,XtestNorm] =  normalizeTrainTest(Xtrain,Xtest)

m = mean(Xtrain,1);
s = std(Xtrain,[],1);

XtrainNorm = (Xtrain - repmat(m,size(Xtrain,1),1)) ./ repmat(s,size(Xtrain,1),1);
XtestNorm  = (Xtest - repmat(m,size(Xtest,1),1)) ./ repmat(s,size(Xtest,1),1);

% if some features are always 0, then s=0 and then there is inf values
XtrainNorm(isinf(XtrainNorm)) = 0;
XtestNorm(isinf(XtestNorm))   = 0;

% also remove nans, not sure why there is some sometimes 
XtrainNorm(isnan(XtrainNorm)) = 0;
XtestNorm(isnan(XtestNorm))   = 0;

end