function [generalParam chooseFeatParam genExpParam classifParam permLabParam] = paramMaker

%% assign parameters for generalParam
generalParam.timeMode           = 'pointByPoint';
generalParam.timeInter          = [0,550];

%% assign parameters for generateExperiment

genExpParam.eventField          = {'correct'};
genExpParam.eventCode           = {1};

%maxime put a conditional here to distinguish 'correct' and
%'correct','type'. That also exists for the regressor. integrate that
%conditional in some way in the future.

genExpParam.nbTrials                = 1; %Number of trials by images you want to keep, {-1}
genExpParam.selectM                 = 'first';
    
genExpParam.regressorsEventField    = 'type';
genExpParam.regressorsEventCode     = {1,0};
    
genExpParam.preClassif              = 'on';
genExpParam.classifType             = 'generalization';
genExpParam.trainTestProp           = 1;
genExpParam.crossValid              = 10;
genExpParam.maxTrainTrials          = 170;
    
%% assign parameters for chooseFeature
chooseFeatParam.area             = 'all';
chooseFeatParam.normalization    = 'Zscore';

%% Assign parameters for the main classifier
classifParam.classifieur         = 'logReg'; % 'L1'
classifParam.optimization        = 'off';

%% Assign parameters for permutLabelClassif
permLabParam.BSiter              = 300;
permLabParam.threshold           = 0.975;
   

end