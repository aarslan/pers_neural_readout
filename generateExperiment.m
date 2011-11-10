% generateExperiment() -  general function that deal with the conditions to analyse data.
%
%                 * From an EEG structure of eeglab create a EXP structure that gives the
%                   regressors to compare conditions and the selector to select trials indices
%                   based either on event field or stimuli name list
%
%                 * This fonction can be used for any type of analysis (ERP, timeFrequency, Classification ...)
%
%                 * In the case o classification it also return a field classifParma with
%                   information relative to classification disign
%
% Usage: [EXP] = generateExperiment(EEG,'option1','valueOp1',...)
%
% Requierd input:     meaning (size or options) {default}
%
%   EEG                  = an EEGlab structure
%
% Optional input:     meaning (size or options or exemple) {default}
%   experiement          = for spetial case, experiment difine in runClassif
%
%
%  Select by event field
%   eventField           = cell of string containing EEG.event you want to use for selection ({'correct','mask', ... });
%   eventCode            = cell of double/string you want to select ({1,'no',...});
%
%  Select by image list
%   imgList              = a cell of string containing the names of images you want to select
%   selectM              = mode of selection for repetition by images:
%                           ('first') : select nbTrials first repetions
%                           ('rand')  : select nbTrials random trials
%                           {'all'}   : select all images
%   nbTrials             = Number of trials by images you want to keep, {-1};
%   nbRepetTrain         = Number of repetitions for training for each image for stability decoding
%   testSameRepet        = to evaluate the effect of number of repetition by images test allways on
%                          same repetion previously saved in the directory (in runClassif) {'off'}
%
%
%  Contrast by eventfield
%   regressorsEventField = EEG.event you want to contrast ('type');
%   regressorsEventCode  = cell of double/string you want to contrast ({1,0});
%
%  Parameters for classification
%   preClassif           = either you run preclassif stuf or not {'on'};
%   classifType          = type of classification
%                           ('stability')       : train and test on same images differents trials
%                           ('diffTrainTest')   : train and test on differents conditions
%                           {'generalization'}  : train and test on differents images
%
%   trainTestProp        = proportion of trials for train/ trials test {1};
%   crossValid           = number of crossvalidation {30};
%   maxTrainTrials       = maximum trials for training {10000};
%   selectorQuantile     = vector of proportion you want to seperate (Ex: [0,33 0,33] will separeate
%                          the data in 3 categories containing each 1/3 of the value)
%
%  For classifType = 'diffTrainTest'
%
%  trainEventField      = cell of string containing EEG.event you want to use for train selection ({'correct','mask', ... });
%  trainEventCode       = cell of double/string you want to select for training({1,'no',...});
%  testEventField       = cell of string containing EEG.event you want to use for test selection ({'correct','mask', ... });
%  testEventCode        = cell or cellarray of double/string you want to select for test({1,'no',...});
%  testEventFieldRelevantTri = to work on differant relevant trials for the test
%  testEventCodeRelevantTri  = to work on differant relevant trials for the test
%  nbTestDiffcateg           = nb of stim you want to use in the test category that are different
%                               from the train one
%
%
% Outputs:            meaning (size or options)
%
%   EXP                  = an EEGlab structure
%
%   EXP.relevantTrials   = vector of integer, indices of trials selected either by eventField
%                          combinaison of eventField,stimuli names, combinaison of eventField and
%                          stimuli name
%   EXP.regressors       = vectors of integer indicating conditions to contrast (size=EEG.nbtrials)
%                          each condition is represented by a number condition1 =0, condition 2 = 1 ...
%   EXP.selector         = if (preClassif='off') this is a vector (size=EEG.nbtrials)
%                          of 1 for relevant trials and 0 for irrelevant
%                          if (preClassif='on') this a matrix (size=nbCrossValidation*EEG.nbtrials)
%                          of 1 for train, 2 for test and 0 for irrelevant trials [0,1,2 ... 1;0,..,2]
%   EXP.listImageName    = cell of string contaning the unique of stiumliname {'stim1.jpg',...,'stimN.jpg'}
%   EXP.classifParametr  = classification parameters
%
%
% Autors: Ali, Maxime, Serre Lab, 2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [EXP] = generateExperiment(EEG, varargin)
%-----------------------------------------
if nargin < 0
    error('generateExperiment:TooFewInputs','Requires an EEG structure');
end
defaults.eventField             = '';
defaults.eventCode              = '';


defaults.selectM                = 'all';
defaults.nbTrials               = -1;
defaults.imgList                = '';
defaults.nbRepetTrain           = '';
defaults.nbRepetTest            = '';
defaults.testSameRepet          ='off';


defaults.regressorsEventField   = '';
defaults.regressorsEventCode    = '';

defaults.preClassif             = 'off';
defaults.classifType            = 'generalization';
defaults.trainTestProp          = 1;
defaults.crossValid             = 30;
defaults.trainEventField        = '';
defaults.trainEventCode         = '';
defaults.testEventField         = '';
defaults.testEventCode          = '';
defaults.maxTrainTrials         = 10000;
defaults.selectorQuantile       = '';
defaults.experiment             = '';
defaults.nbTestDiffcateg        = '';
defaults.testOnSameCateg        = 'on';

args = propval(varargin,defaults);

eventField              = args.eventField;
eventCode               = args.eventCode;


nbTrials                = args.nbTrials;
selectM                 = args.selectM;
imgList                 = args.imgList;
nbRepetTrain            = args.nbRepetTrain;
nbRepetTest             = args.nbRepetTest;

testSameRepet           = args.testSameRepet;

regressorsEventField    = args.regressorsEventField;
regressorsEventCode     = args.regressorsEventCode;

preClassif              = args.preClassif;
classifType             = args.classifType;
trainTestProp           = args.trainTestProp;
crossValid              = args.crossValid ;
trainEventField         = args.trainEventField ;
trainEventCode          = args.trainEventCode ;
testEventField          = args.testEventField ;
testEventCode           = args.testEventCode ;
maxTrainTrials          = args.maxTrainTrials ;
selectorQuantile        = args.selectorQuantile ;
experiment              = args.experiment;
nbTestDiffcateg         = args.nbTestDiffcateg;
testOnSameCateg         = args.testOnSameCateg;

%-------------------------------------------------------------------------
% EXP field common with EEG field
%--------------------------------
EEG.featClassif   = [];
EXP               = EEG;
EXP.data          = [];

%---------------------------------------------------------%
% Select trials based on event: create EEG.relevantTrials %
%---------------------------------------------------------%

if ~isempty(eventField)&(~isempty(eventCode))
    EXP.relevantTrials = selectRelevantTrials(EEG,eventField,eventCode);
else
    EXP.relevantTrials = 1:EEG.trials;
    
end

%---------------------------------------------------------------------%
% Select trials based on images repetitons: update EEG.relevantTrials %
%---------------------------------------------------------------------%
if nbTrials > 0 & isfield(EEG.event,'imageName')
    
    [listImageName, ind] = selectRepetIm(EEG,EXP.relevantTrials,imgList,selectM,nbTrials);
    
    % EXP field
    %----------
    EXP.relevantTrials = ind;
    EXP.listImageName = listImageName;
end



%---------------------------------------------------%
% Code conditons to contrast: Create EXP.regressors %
%---------------------------------------------------%

if ~isempty(regressorsEventField) & ~isempty(regressorsEventCode)
    
    % EXP field
    %----------
    EXP.regressors           = ones(1,EEG.trials)*-1;
    
    for i=1:length(regressorsEventCode)
        ind = findStruct(EEG.event,regressorsEventField,regressorsEventCode{i});
        display (['------Code regressors based on code : ' num2str(regressorsEventCode{i}) 'of field : ' regressorsEventField  ' : corresponding to ' num2str(length(ind))  ' trials-----'])
        EXP.regressors(ind)=i-1;
    end
    
    clear ind
end

%------------%
% PreClassif %
%------------%

if strcmp(preClassif,'on')
    
    display (['------Compute classification of ' classifType ' type-----'])
    if ~isfield(EXP,'regressors')
        error (' PreClassif requierd and EXP.regressors field')
    end
    
    % General classif parameters
    %---------------------------
    EXP.classifParam                    = {};
    EXP.classifParam.trainTestProp      = trainTestProp;
    EXP.classifParam.crossValid         = crossValid;
    EXP.selector                        = [];
    EXP.descriptor                      = [];
    
    % Indices positive and negatives
    %-------------------------------
    positiveInd = find(EXP.regressors==0);
    negativeInd = find(EXP.regressors==1);
    if isfield(EXP,'relevantTrials')
        positiveInd = intersect(EXP.relevantTrials,positiveInd);
        negativeInd = intersect(EXP.relevantTrials,negativeInd);
    end
    
    %----------%
    % Selector %
    %----------%
    % Tranform quantitative variable in discret one
    %----------------------------------------------
    if ~isempty(selectorQuantile)
        val         = [EEG.event.(trainEventField{1})];
        valNoNan    = val(intersect(find(~isnan(val)),EXP.relevantTrials));
        lim         = quantile(valNoNan,selectorQuantile);
        % Attribute a discrete code
        %---------------------------
        val2 = val;
        val2(~isnan(val)&val<lim(1))=1;
        
        for li = 1:length(lim)
            val2(~isnan(val)&val>=lim(li))=li+1;
        end
        val2(isnan(val)) = 0;
    end
    
    % Stability
    %----------
    if strcmp(classifType,'stability')&&strcmp(testSameRepet,'on')
        
        load ([EEG.filepath  EXP.filename(1:3) 'StabTestInd.mat']);
        
        positiveIndRandTest = intersect(StabTestInd,positiveInd);
        negativeIndRandTest = intersect(StabTestInd,negativeInd);
        
        %----------------------%
        % Different train test
        %----------------------%
        
    elseif strcmp(classifType,'diffTrainTest')
        
        if isempty(selectorQuantile)
            % Indeces for training
            %---------------------
            indice={};
            for i=1:length(trainEventField) %loop on events
                display (['------Select ' trainEventField{i} ' ' num2str(trainEventCode{i}) ' for training-----'])
                ind = findStruct(EEG.event,trainEventField{i},trainEventCode{i});
                indice{i}    = ind;
            end
            if length(trainEventField)==1
                indiceTrain = indice{1};
            end
            % If you want to select by a combination of two events
            %-----------------------------------------------------
            if length(trainEventField)==2
                indiceTrain = intersect(indice{1},indice{2});
                % If you want to select by a combination of  more than two events
                %----------------------------------------------------------------
            elseif length(trainEventField)>2
                error('Selection by more than 2 event Field has not been implemented yet')
            end
            
            positiveIndTrain    = intersect(positiveInd,indiceTrain);
            negativeIndTrain    = intersect(negativeInd,indiceTrain);
            
            % Indeces for testing
            %---------------------
            % Special test relevantTrials
            %----------------------------
            if strcmp(experiment,'distanceEffectFirstIM')
                tri    = [EEG.event.trial];
                li     = [EEG.event.list];
                typ    = str2num(char(EEG.event.type));
                type   = typ';
                indFirstPres = find((type==120|type==220)&tri<321&li==1);
                
                positiveIndTest     = intersect(positiveInd,indFirstPres);
                negativeIndTest     = intersect(negativeInd,indFirstPres);
            else
                positiveIndTest     = positiveInd;
                negativeIndTest     = negativeInd;
            end
        end
        %------------------------------------------------------
        
        
        if length(testEventCode)>0 % to allow difftrain test on same subCat with several repet on train
            indice={};
            for i=1:length(testEventField) %loop on events
                if length(testEventCode{i})==1
                    display (['------Select ' testEventField{i} ' ' num2str(testEventCode{i}) ' for testing-----'])
                    %[EEG, ind]   = pop_selectevent( EEG,testEventField{i},testEventCode{i},'deleteepochs','off');
                    ind = findStruct(EEG.event,testEventField{i},testEventCode{i});
                    indice{i}    = ind;
                else
                    % Loop on different test eventCode
                    %---------------------------------
                    nbCodePo    = [];
                    nbCodeNe    = [];
                    indiceTest  = [];
                    for co=1:length(testEventCode{i})
                        code       = testEventCode{i};
                        display (['------Select ' testEventField{i} ' ' num2str(code(co)) ' for testing-----'])
                        %[EEG, ind]         = pop_selectevent( EEG,testEventField{i},code(co),'deleteepochs','off');
                        ind = findStruct(EEG.event,testEventField{i},code(co));
                        
                        positiveIndice{co} = intersect(positiveIndTest,ind);
                        negativeIndice{co} = intersect(negativeIndTest,ind);
                        nbCodePo           = [nbCodePo, length(positiveIndice{co})];
                        nbCodeNe           = [nbCodeNe, length(negativeIndice{co})];
                        indiceTest         = [indiceTest,ind];
                    end
                end
            end
            
            if length(testEventField)==1&length(testEventCode{1})==1
                indiceTest = indice{1};
                % If you want to select by a combination of two events
                %-----------------------------------------------------
            elseif length(testEventField)==2&length(testEventCode{1})==1&length(testEventCode{2})==1
                indiceTest = intersect(indice{1},indice{2});
                % If you want to select by a combination of  more than two events
                %----------------------------------------------------------------
            elseif length(testEventField)>2
                error('Selection by more than 2 event Field has not been implemented yet')
            elseif length(testEventField)==2&(length(testEventCode{1})>1|length(testEventCode{2})>1)
                error('Selection by 2 event Field and several code for at least one of them has not been implemented yet')
            end 
        end
    elseif ~isempty(selectorQuantile)
        % Selection only by positive
        %----------------------------
        positiveIndTrain    = intersect(positiveInd,find(val2==trainEventCode{1}));
        negativeIndTrain    = negativeInd;
        positiveIndTest     = intersect(positiveInd,find(val2==testEventCode{1}));
        negativeIndTest     = negativeInd;
        
        
    end
    if nbRepetTrain>1
        % nb images with nbTrain repetition
        %----------------------------------
        if strcmp(testOnSameCateg,'on')
            [listImTrainPos, StabTrainInd] = selectRepetIm(EEG,positiveIndTrain,imgList,selectM,nbRepetTrain+nbRepetTest);
            display(['Difftrain test with several image rep for training: select in total ' num2str(nbRepetTrain+nbRepetTest) ' repetition on ' num2str(length(unique(listImTrainPos))) ' positives images '])
            
            [listImTrainNeg, StabTrainInd] = selectRepetIm(EEG,negativeIndTrain,imgList,selectM,nbRepetTrain+nbRepetTest);
            display(['Difftrain test with several image rep for training: select in total ' num2str(nbRepetTrain+nbRepetTest) ' repetition on ' num2str(length(unique(listImTrainNeg))) 'negatives images '])
            
            nReal= min(length(listImTrainPos),length(listImTrainNeg))/2;
            if nReal<maxTrainTrials
                maxTrainTrials = nReal;
            end
        else
            [listImTrainPos, StabTrainInd] = selectRepetIm(EEG,positiveIndTrain,imgList,selectM,nbRepetTrain);
            display(['Difftrain test with several image rep for training: select in total ' num2str(nbRepetTrain) ' repetition on ' num2str(length(unique(listImTrainPos))) ' positives images '])
            
            [listImTrainNeg, StabTrainInd] = selectRepetIm(EEG,negativeIndTrain,imgList,selectM,nbRepetTrain);
            display(['Difftrain test with several image rep for training: select in total ' num2str(nbRepetTrain) ' repetition on ' num2str(length(unique(listImTrainNeg))) 'negatives images '])
        end
            
    end
    
    %-----------------%
    % Crossvalidation %
    %-----------------%
    for i=1:crossValid
        %rand('state', sum(100*clock)); %initialisation du randLabel

        % 0 for non selected trials
        %--------------------------
        crossValidVector        = zeros(1,EEG.trials);
        descriptorVector        = ones(1,EEG.trials)*-1;%for complex paradigm
        % Generalisation
        %---------------
        if strcmp(classifType,'generalization')
            % Number initiale of trials
            %--------------------------
            EXP.classifParam.nbInitialPostive   = length(positiveInd);
            EXP.classifParam.nbInitialNegative  = length(negativeInd);
            % number of train and test trials to select (same for positive and negative)
            %---------------------------------------------------------------------------
            nbTrain                       = floor(min([EXP.classifParam.nbInitialPostive,EXP.classifParam.nbInitialNegative])/2);
            if nbTrain < maxTrainTrials/2
                nbTrain                   = nbTrain;
            else
                nbTrain                   = maxTrainTrials/2;
            end
            EXP.classifParam.nbTrain            = nbTrain;
            EXP.classifParam.nbTest             = floor(nbTrain*trainTestProp);
            % Randomize indices of two type of regressors
            %--------------------------------------------
            %rand('state', sum(100*clock)); %initialisation du randLabel
            positiveIndRand         = positiveInd(randperm(length(positiveInd)));
            negativeIndRand         = negativeInd(randperm(length(negativeInd)));
            % Indices for training
            %---------------------
            positiveIndRandTrain    = positiveIndRand(1:EXP.classifParam.nbTrain);
            negativeIndRandTrain    = negativeIndRand(1:EXP.classifParam.nbTrain);
            
            % Indices for testing
            %--------------------
            positiveIndRandTest     = setdiff(positiveInd,positiveIndRandTrain);
            negativeIndRandTest     = setdiff(negativeInd,negativeIndRandTrain);
            positiveIndRandTest     = positiveIndRandTest(1:EXP.classifParam.nbTest);
            negativeIndRandTest     = negativeIndRandTest(1:EXP.classifParam.nbTest);
            
        elseif strcmp(classifType,'diffTrainTest')&nbRepetTrain<2
            % Number initiale of trials
            %--------------------------
            EXP.classifParam.nbInitialpositiveIndTrain   = length(positiveIndTrain);
            EXP.classifParam.nbInitialnegativeIndTrain   = length(negativeIndTrain);
            EXP.classifParam.nbInitialpositiveIndTest    = length(positiveIndTest);
            EXP.classifParam.nbInitialnegativeIndTest    = length(negativeIndTest);
            % number of train and test trials to select (same number for positive and negative)
            %---------------------------------------------------------------------------
            if length(testEventCode{1})==1
                nbTrain = min([length(positiveIndTrain),length(negativeIndTrain),length(positiveIndTest),length(negativeIndTest)]);
            else
                % One eventField several eventCode
                %---------------------------------
                nbTrain = min([length(positiveIndTrain)/2,length(negativeIndTrain)/2,nbCodePo,nbCodeNe]);
            end
            % Two select an arbitrary number of trials = maxTrainTrials
            %----------------------------------------------------------
            if nbTrain < maxTrainTrials/2
                nbTrain                   = nbTrain;
            else
                nbTrain                   = maxTrainTrials/2;
            end
            
            EXP.classifParam.nbTrain  = nbTrain;
            EXP.classifParam.nbTest   = floor(nbTrain*trainTestProp);
            
            display (['------Number of trials for training by regressor: ' num2str(nbTrain) '-----'])
            display (['------Number of trials for testing by regressor: ' num2str(EXP.classifParam.nbTest) '-----'])
            
            %--------------------------------%
            % Train indeces random selection
            %--------------------------------%
            % Randomize indices
            %-------------------
            positiveIndRandTrain        = positiveIndTrain(randperm(length(positiveIndTrain)));
            negativeIndRandTrain        = negativeIndTrain(randperm(length(negativeIndTrain)));
            % Selection
            %----------
            positiveIndRandTrain        = positiveIndRandTrain(1:EXP.classifParam.nbTrain);
            negativeIndRandTrain        = negativeIndRandTrain(1:EXP.classifParam.nbTrain);
            
            %--------------------------------%
            % Test indeces random selection
            %--------------------------------%
            positiveIndRandTest         = [];
            negativeIndRandTest         = [];
            for co = 1:length(testEventCode{1})
                % One eventField several eventCode
                %---------------------------------
                if length(testEventCode{1})>1
                    positiveIndTest = setdiff(positiveIndice{co},positiveIndRandTrain);
                    negativeIndTest = setdiff(negativeIndice{co},negativeIndRandTrain);
                end
                % Randomize indices
                %-------------------
                positiveIndRandomTest       = positiveIndTest(randperm(length(positiveIndTest)));
                negativeIndRandomTest       = negativeIndTest(randperm(length(negativeIndTest)));
                % Selection
                %----------
                positiveIndRandTest         = [positiveIndRandTest positiveIndRandomTest(1:EXP.classifParam.nbTest)];
                negativeIndRandTest         = [negativeIndRandTest negativeIndRandomTest(1:EXP.classifParam.nbTest)];
                
            end
            
        elseif strcmp(classifType,'stability')
            
            if strcmp(testSameRepet,'off')
                [listImTot, StabTestInd] = selectRepetIm(EEG,EXP.relevantTrials,imgList,selectM,nbRepetTrain+nbRepetTest);
                display(['Satbility preClassification: select in total ' num2str(nbRepetTrain+nbRepetTest) ' repetition on ' num2str(length(unique(listImTot))) ' images for CV ' num2str(i)])
                
                [listImageNameTest, StabTestInd] = selectRepetIm(EEG,EXP.relevantTrials,unique(listImTot),selectM,nbRepetTest);
                
                positiveIndRandTest = intersect(StabTestInd,positiveInd);
                negativeIndRandTest = intersect(StabTestInd,negativeInd);
            end
            
            % To select repetition different from test
            %------------------------------------------
            relevantTrials = setdiff(EXP.relevantTrials,StabTestInd);
            
            [listImageNameTrain, StabTrainInd] = selectRepetIm(EEG,relevantTrials,listImageNameTest,selectM,nbRepetTrain);
            if ~isequal(listImageNameTrain,listImageNameTest)
                error('You try to do a stability classification with different images name, this is not possible!')
            end
            
            % Selection
            %-----------
            positiveIndRandTrain = intersect(StabTrainInd,positiveInd);
            negativeIndRandTrain = intersect(StabTrainInd,negativeInd);
            
        elseif strcmp(classifType,'diffTrainTest')&nbRepetTrain>1
            display('Use different train and test mode with several repetion by train images')
            % Train
            %------
            if strcmp(testOnSameCateg,'off')
                [listImTrainPosTr, positiveIndRandTrain] = selectRepetIm(EEG,positiveIndTrain,listImTrainPos,selectM,nbRepetTrain);
                display(['Difftrain test with several image rep for training: select in total ' num2str(nbRepetTrain) ' repetition on ' num2str(length(unique(listImTrainPos))) ' positives images '])
                
                [listImTrainNegTr, negativeIndRandTrain] = selectRepetIm(EEG,negativeIndTrain,listImTrainNeg,selectM,nbRepetTrain);
                display(['Difftrain test with several image rep for training: select in total ' num2str(nbRepetTrain) ' repetition on ' num2str(length(unique(listImTrainNeg))) ' negatives images '])
                positiveIndRandTest = [];
                negativeIndRandTest = [];
                
            else
                
                listImTrainPosrand      = listImTrainPos(randperm(length(listImTrainPos)));
                listImTrainNegrand      = listImTrainNeg(randperm(length(listImTrainNeg)));
                listImTrainPosTr    = listImTrainPosrand(1:maxTrainTrials);
                listImTrainNegTr    = listImTrainNegrand(1:maxTrainTrials);
                
                [listImTrainPosTr, positiveIndRandTrain] = selectRepetIm(EEG,positiveIndTrain,listImTrainPosTr,selectM,nbRepetTrain);
                display(['Difftrain test with several image rep for training: select in total ' num2str(nbRepetTrain) ' repetition on ' num2str(length(unique(listImTrainPosTr))) ' positives images '])
                
                [listImTrainNegTr, negativeIndRandTrain] = selectRepetIm(EEG,negativeIndTrain,listImTrainNegTr,selectM,nbRepetTrain);
                display(['Difftrain test with several image rep for training: select in total ' num2str(nbRepetTrain) ' repetition on ' num2str(length(unique(listImTrainNegTr))) ' negatives images '])
                
                % Test
                %-----
                
                positiveIndRandTest = [];
                negativeIndRandTest = [];
                
                % on same images from same type than training
                positiveIndTrainTe = setdiff(positiveIndTrain,positiveIndRandTrain);
                negativeIndTrainTe = setdiff(negativeIndTrain,negativeIndRandTrain);
                
                [listImTrainPosTr, po] = selectRepetIm(EEG,positiveIndTrainTe,listImTrainPosTr,selectM,nbRepetTest);
                display(['Difftrain test with several image rep for testing: select in total ' num2str(nbRepetTest) ' repetition on ' num2str(length(unique(listImTrainPosTr))) ' same positives images '])
                
                [listImTrainNegTr, ne] = selectRepetIm(EEG,negativeIndTrainTe,listImTrainNegTr,selectM,nbRepetTest);
                display(['Difftrain test with several image rep for testing: select in total ' num2str(nbRepetTest) ' repetition on ' num2str(length(unique(listImTrainNegTr))) ' same negatives images '])
                
                positiveIndRandTest = [positiveIndRandTest po];
                negativeIndRandTest = [negativeIndRandTest ne];
                descriptorVector([po ne])= 1;
                
                
                % on different images from same type than training
                listImTestPos = setdiff(listImTrainPos,listImTrainPosTr);
                listImTestPos = listImTestPos(1:maxTrainTrials);
                listImTestNeg = setdiff(listImTrainNeg,listImTrainNegTr);
                listImTestNeg = listImTestNeg(1:maxTrainTrials);
                
                [listImTestPos, po] = selectRepetIm(EEG,positiveIndTrain,listImTestPos,selectM,nbRepetTest);
                display(['Difftrain test with several image rep for testing: select in total ' num2str(nbRepetTest) ' repetition on ' num2str(length(unique(listImTestPos))) ' different positives images '])
                
                [listImTestNeg, ne] = selectRepetIm(EEG,negativeIndTrain,listImTestNeg,selectM,nbRepetTest);
                display(['Difftrain test with several image rep for testing: select in total ' num2str(nbRepetTest) ' repetition on ' num2str(length(unique(listImTrainNegTr))) 'different negatives images '])
                
                positiveIndRandTest = [positiveIndRandTest po];
                negativeIndRandTest = [negativeIndRandTest ne];
                descriptorVector([po ne])= 2;
            end
            % on different type
            if length(testEventCode)>0
                for co = 1:length(testEventCode{1})
                    positiveIndTest = positiveIndice{co};
                    negativeIndTest = negativeIndice{co};
                    % Randomize indices
                    %-------------------
                    positiveIndRandomTest       = positiveIndTest(randperm(length(positiveIndTest)));
                    negativeIndRandomTest       = negativeIndTest(randperm(length(negativeIndTest)));
                    % Selection
                    %----------
                    for u=1:length(positiveIndice);taillepo(u)=length(positiveIndice{u});end
                    for u=1:length(negativeIndice);taillene(u)=length(negativeIndice{u});end
                    mini= min([taillepo taillene]);
                    
                    if length(nbTestDiffcateg)<1
                        
                        if mini>=maxTrainTrials
                            
                            positiveIndRandTest         = [positiveIndRandTest positiveIndRandomTest(1:maxTrainTrials)];
                            negativeIndRandTest         = [negativeIndRandTest negativeIndRandomTest(1:maxTrainTrials)];
                            descri=testEventCode{1}(co);
                            descriptorVector([positiveIndRandomTest(1:maxTrainTrials) negativeIndRandomTest(1:maxTrainTrials)]) = descri;
                            
                        elseif mini<maxTrainTrials
                            maxTestTrials = mini;
                            
                            positiveIndRandTest         = [positiveIndRandTest positiveIndRandomTest(1:maxTestTrials)];
                            negativeIndRandTest         = [negativeIndRandTest negativeIndRandomTest(1:maxTestTrials)];
                            descri=testEventCode{1}(co);
                            descriptorVector([positiveIndRandomTest(1:maxTestTrials) negativeIndRandomTest(1:maxTestTrials)]) = descri;
                        end
                    else
                        positiveIndRandTest         = [positiveIndRandTest positiveIndRandomTest(1:nbTestDiffcateg)];
                        negativeIndRandTest         = [negativeIndRandTest negativeIndRandomTest(1:nbTestDiffcateg)];
                        descri=testEventCode{1}(co);
                        descriptorVector([positiveIndRandomTest(1:nbTestDiffcateg) negativeIndRandomTest(1:nbTestDiffcateg)]) = descri;
                        
                    end
                end
            end
            
        end
        
        % 1 for training trials
        %----------------------
        crossValidVector([positiveIndRandTrain,negativeIndRandTrain])  = 1;
        % 2 for testing trials
        %---------------------
        crossValidVector([positiveIndRandTest,negativeIndRandTest])    = 2;
        % Selector matrix: crossValid*nbtrials
        %-------------------------------------
        EXP.selector   = [EXP.selector ; crossValidVector];
        EXP.descriptor = [EXP.descriptor ; descriptorVector];
        
        
    end
    
end

end




function   [indiceRelevant] = selectRelevantTrials(EEG,eventField,eventCode)
% Select relevant trials
%-----------------------
indice={};
for i=1:length(eventField) %loop on eventField
    if ~iscell(eventCode{i})
        ind = findStruct(EEG.event,eventField{i},eventCode{i});
        display (['------Select ' eventField{i} ' ' num2str(eventCode{i}) ': keep ' num2str(length(ind))  ' trials-----'])
        indice{i}=ind;
    else % To select one EventField with 2 eventCode
        eventCo  = eventCode{i};
        indice2 = [];
        for j=1:length(eventCo) %loop on eventCode
            ind = findStruct(EEG.event,eventField{i},eventCo{j});
            display (['------Select ' eventField{i} ' ' num2str(eventCo{j}) ': keep ' num2str(length(ind))  ' trials-----'])
            if size(ind,1)>size(ind,2)
                indice2 =[indice2; ind];
            else
                indice2 =[indice2, ind];
            end
        end
        indice{i}=indice2;
    end
end

if length(eventField)==1
    indiceRelevant = indice{1};
end
% If you want to select by a combination of two events
%-----------------------------------------------------
if length(eventField)==2
    indiceRelevant = intersect(indice{1},indice{2});
    
    % If you want to select by a combination of  more than two events
    %----------------------------------------------------------------
elseif length(eventField)==3
    indiceRelevant = intersect(indice{1},indice{2});
    indiceRelevant = intersect(indiceRelevant,indice{3});
elseif length(eventField)==4
    indiceRelevant = intersect(indice{1},indice{2});
    indiceRelevant = intersect(indiceRelevant,indice{3});
    indiceRelevant = intersect(indiceRelevant,indice{4});
end

display(['Keeping ' num2str(length(indiceRelevant)) ' relevant trials'])

end




function     [listImageName, ind] = selectRepetIm(EEG,relevantTrials,imgList,selectM,nbTrials)

% All images names over all repetitions
%--------------------------------------
for j=1:EEG.trials
    nameall(j) = getfield(EEG.event,{j},'imageName');
end

% If no image list given in input use images names present in the EEG
%--------------------------------------------------------------------
if isempty(imgList)
    name = unique(nameall);
else
    name = imgList;
end
% Selection process
%------------------


ind                 = [];% indices to select
listImageName       = [];%list of images names that are keeped

display (['------Select ' num2str(nbTrials) ' repetion by images with the ' selectM ' mode-----']);

for j = 1:length(name) %loop by image name
    nbTrialsReal(j) = sum(strcmp(name(j),nameall));% number of repetion by images
    
    if nbTrialsReal<nbTrials
        display(['------You can not select more than ' num2str(nbTrialsReal(j)) ' repetions for image ' name{j} ' you are going to skip images------'])
    end
    
    indice = find(strcmp(name(j),nameall));% indices of trials coresponding to image j
    
    indice = intersect(relevantTrials,indice);% combination with EXP.relevantTrials
    
    % Select nbTrials images
    %-----------------------
    if nbTrials > 0
        if length(indice)>=nbTrials
            
            % Use only images that are reapeted enough
            %-----------------------------------------
            if strcmp(selectM,'first')
                
                % Select first nbTrial images
                %----------------------------
                ind     = [ind indice(1:nbTrials)];
            else
                % Select nbTrial random images
                %-----------------------------
                indice  = indice(randperm(length(indice)));
                ind     = [ind indice(1:nbTrials)];
            end
            listImageName = [listImageName name(j)];
        end
    end
end

end



