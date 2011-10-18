% runClassif() - main function to run a classification on neural data
%
% Usage: [EXP] = mainClassification(EEG, paramStruct, feats, rand, timeMode)
%
% Optional input:     meaning (size or options) {default}
%
%   EEG             = to work directly on an EEG structure {''};
%   paramStruct     = Structure containing the parameters. It can be either
%                     a structure, a filepath pointing to a matlab variable, or [] to invoke
%                     the user interface(not implemented yet)
%   feats           = which feats to be used for classification (cell
%                       array). Currently just sequential
%   rand            = allow to run the classification with random label to
%                     compute statistical threshlod, save a .mat with the same name
%                     but finished by rand in ./rand {'off'};
%
%   timeMode        = different way of doing the temporal classification can be:
%                      - ('pointByPoint'): run a classifier for each point of time
%                      - ('allEpoch')    : run one classifier using the allepoch data points
%                      - {'both'}        : run both to plot the accuracy of the first and the weights of the second
% Outputs:            meaning (size or options)
%
%   EXP             = structure all parameters about classification
%                     experiment, saved as [ dir_out name '.mat']
%
% Author: Ali, Maxime, Serre Lab, 2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [EXP] = mainClassification(EEG, paramStruct, feats, rand, timeMode, outPath )

if ischar(paramStruct)
    load(paramStruct)
else
    % create the parameters
    [generalParam chooseFeatParam genExpParam classifParam permLabParam] = paramMaker;
end

% parse the parameters from the structure and call generateExperiment in
% eval()
genExpString = [];
genExpFields = fieldnames(genExpParam);
for ff=1:length(genExpFields)
    genExpString = [genExpString ', ''' genExpFields{ff} ''', genExpParam.' genExpFields{ff}];
end
genExpString = ['generateExperiment(EEG, ' genExpString(2:end) ');'];
[EXP] = eval(genExpString);

EXP.timeMode        = generalParam.timeMode;
EXP.timeInter       = generalParam.timeInter;
EXP.chooseFeatParam = chooseFeatParam;
EXP.genExpParam     = genExpParam;
EXP.classifParam    = classifParam;
EXP.permLabParam    = permLabParam;

% Time working interval %

indt1     = max(find(EEG.times<=EXP.timeInter(1)));
indt2     = min(find(EEG.times>=EXP.timeInter(2)));
t1        = EEG.times(indt1);
t2        = EEG.times(indt2);
timeInter = [t1,t2];

% Store variables %
EXP.timesOfClassif              = EEG.times(indt1:indt2);
EXP.EEGname                     = EEG.filename;
EXP.EEGfname                    = EEG.filepath;
EXP.machine                     = computer;
EXP.subject                     = EEG.filename(1:3);
EXP.rand                        = rand;

% put Everything in a big allParam structure for convenience (preserving the backwards compatibility)
EXP.allParam.generalParam    = generalParam;
EXP.allParam.chooseFeatParam = chooseFeatParam;
EXP.allParam.genExpParam     = genExpParam;
EXP.allParam.classifParam    = classifParam;
EXP.allParam.permLabParam    = permLabParam;

%----------------%
% Classification %
%----------------%

% Loop over cross validation(train test combination)
%---------------------------------------------------
for i=1:EXP.genExpParam.crossValid
    
    % Labels
    %-------
    Ytrain    = EXP.regressors(EXP.selector(i,:)==1)';
    Ytest     = EXP.regressors(EXP.selector(i,:)==2)';
    
    % Labels for binary classification
    %---------------------------------
    Ytrain(Ytrain==0) = -1;
    Ytest(Ytest==0)   = -1;
    
    % Randomized label for statistcal threshold
    %------------------------------------------
    if strcmp(rand,'on')
        Ytrain = Ytrain(randperm(length(Ytrain)));
        Ytest  = Ytest(randperm(length(Ytest)));
    end
    
    
    % loop over features
    %-------------------
    for fe=1:length(feats)

        % Deal with bin shift
        % ALI: this is wrong.
        %---------------------------
                  times     = EEG.times;%+(EEG.binSizes(1)/2)*(1000/EEG.srate);
                  indt1     = max(find(times<=timeInter(1)));
                  indt2     = min(find(times>=timeInter(2)));

        % PointByPoint
        %-------------
        if strcmp(timeMode,'pointByPoint') || strcmp(timeMode,'both')
            time = indt1:indt2;
            % loop over point of time
            %-------------------------
            for j=1:length(time)
                display(['PointByPoint Classification: ' num2str(j) ' over ' num2str(length(time)) '/ crossvalidation ' num2str(i) ' over ' num2str(genExpParam.crossValid)]);
                % chooseFeature_ali
                %--------------
                [Xtrain,Xtest] = chooseFeature_ali(EEG,EXP,i,'timepoint',time(j),'area',chooseFeatParam.area,'features', feats{fe},'normalization', chooseFeatParam.normalization);
                %classification
                %---------------

                %[YTrainPredict,YTestPredict,B,bestGamma, mnrStats] = classification(Xtrain,Xtest,Ytrain,Ytest, classifParam.classifieur, classifParam.optimization);
                [YTrainPredict,YTestPredict,B,bestGamma] = classification(Xtrain,Xtest,Ytrain,Ytest, classifParam.classifieur, classifParam.optimization);

                %mesure performances
                %-------------------
                perfTrain           = measure_perf(YTrainPredict,Ytrain);
                perfTest            = measure_perf(YTestPredict,Ytest);
                perfTest.B          = B;
                perfTest.bestGamma  = bestGamma;
                %permutation for statistical threshold estimation
                %-------------------------------------------------
                BSthreshold          = permutLabelClassif(Ytest,YTestPredict, permLabParam.BSiter, permLabParam.threshold);
                perfTest.BSthreshold = BSthreshold;
                %store classification performances
                %---------------------------------
                EXP.classifResults.(feats{fe}).pointByPoint.perfTrain(i,j)  = perfTrain;
                EXP.classifResults.(feats{fe}).pointByPoint.perfTest(i,j)   = perfTest;
                %EXP.classifResults.(feats{fe}).pointByPoint.mnrStats(i,j)   = mnrStats;
            end
        end
        
        % allEpoch
        %---------
        if strcmp(timeMode,'allEpoch') || strcmp(timeMode,'both')
            display(['Classif allEpoch:  crossvalidation ' num2str(i) ' over ' num2str(genExpParam.crossValid)]);
            % chooseFeature_ali
            %--------------
            [Xtrain,Xtest] = chooseFeature_ali(EEG,EXP,i,'timepoint', [indt1:indt2],'area', chooseFeatParam.area,'features', feats{fe},'normalization',chooseFeatParam.normalization);
            %classification
            %---------------
            [YTrainPredict,YTestPredict,B,bestGamma] = classification(Xtrain,Xtest,Ytrain,Ytest, classifParam.classifieur, classifParam.optimization);
            %mesure performances
            %-------------------
            perfTrain           = measure_perf(YTrainPredict,Ytrain);
            perfTest            = measure_perf(YTestPredict,Ytest);
            perfTest.B          = B;
            perfTest.bestGamma  = bestGamma;
            %permutation for statistical threshold estimation
            %-------------------------------------------------
            BSthreshold          = permutLabelClassif(Ytest,YTestPredict, permLabParam.BSiter, permLabParam.threshold);
            perfTest.BSthreshold = BSthreshold;
            %store classification performances
            %---------------------------------
            EXP.classifResults.(feats{fe}).allEpoch.perfTrain(i)  = perfTrain;
            EXP.classifResults.(feats{fe}).allEpoch.perfTest(i)   = perfTest;
            %EXP.classifResults.(feats{fe}).allEpoch.mnrStats(i)   = mnrStats;
        end
        
    end
end


fileName = [EEG.filename(1:3) 'rand' upper(rand) '_' sprintf('%dmsTo%dms_', generalParam.timeInter) genExpParam.regressorsEventField ...
    sprintf('_%s_%s_%s_%s_%s_%s_%s', feats{:}) chooseFeatParam.normalization '_' classifParam.classifieur];

%Include creation timestamps of the parameter and the EXP file inside the
%file in case it's needed in the future.

EXP.fileTimeStamp.results = datestr(clock,'HHMMmmddyy');
EXP.fileTimeStamp.params = datestr(clock,'HHMMmmddyy');

if outPath ~= 0
    try
        save([outPath '/' fileName '.mat'],'EXP')
    catch
        warning('couldnt find the outout path, saving in current directory')
        save([fileName '.mat'],'EXP')
    end
end
end
