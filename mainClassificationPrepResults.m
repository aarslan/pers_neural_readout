function [performance decInfo] = mainClassificationPrepResults(EXP)

% EXP:      classification outputfile

%Initialize
%----------
perfValue = 'accuracy'; %hardcoded for convenience
wheightsOrder = 'area'; %'contribution'; %'area';
errorB = sqrt(EXP.genExpParam.crossValid); %this can also be 1 depending on your error bar calc
weights   = [];
trial     = 1:size(EXP.selector,2);
featList = fields(EXP.classifResults);


%loop over each feature classified
for ff=1:length(featList)
    
    %loop over cross validation
    %--------------------------
    for cv=1:EXP.genExpParam.crossValid
        if cv==1
            EXP.trialUsed  = [];
        end
        % For accuracy
        %-------------
        
        mode = fields(EXP.classifResults.(featList{ff})); mode = mode{:};
        if strcmp(perfValue,'accuracy')
            EXP.classifResults.(featList{ff}).acc(cv,:)     = [EXP.classifResults.(featList{ff}).(mode).perfTest(cv,:).accuracy]*100;
            EXP.classifResults.(featList{ff}).accRand(cv,:) = [EXP.classifResults.(featList{ff}).(mode).perfTest(cv,:).BSthreshold]*100;
        end
        
        % Find ind of trials used for test
        %---------------------------------
        EXP.trialUsed = [EXP.trialUsed, trial(EXP.selector(cv,:)==2)];
        
        % For weights
        %------------
        wei       = [];
        if strcmp(EXP.timeMode,'pointByPoint')
            wei = [EXP.classifResults.(featList{ff}).pointByPoint.perfTest(cv,:).B];
            wei = wei(1:(size(wei,1)-1),:);
        elseif strcmp(EXP.timeMode,'allEpoch')
            wei = [EXP.classifResults.(featList{ff}).allEpoch.perfTest(cv).B];
            wei = wei(1:(length(wei)-1)); % remove biais
            wei = reshape(wei,length(EXP.chanlocs.(validLocLabel)),length(EXP.timesOfClassif));
        end
        weights = cat(3,weights,wei);
    end
    
    % Average performances on all the CV
    %-----------------------------------
    EXP.classifResults.(featList{ff}).accMean = mean(EXP.classifResults.(featList{ff}).acc,1);
    EXP.classifResults.(featList{ff}).accSD   = std(EXP.classifResults.(featList{ff}).acc,[],1)/errorB;
    EXP.classifResults.(featList{ff}).accPop = EXP.classifResults.(featList{ff}).acc;
    % Random thershold computed on test only
    %---------------------------------------
    EXP.classifResults.(featList{ff}).accRandMean = mean(EXP.classifResults.(featList{ff}).accRand,1);
    EXP.classifResults.(featList{ff}).accRandSD   = std(EXP.classifResults.(featList{ff}).accRand,[],1)/errorB;
    EXP.classifResults.(featList{ff}).accRandMean = EXP.classifResults.(featList{ff}).accRand;
    % Weights
    %---------
    EXP.classifResults.(featList{ff}).weightsMean = abs(squeeze(mean(weights,3))); % we take the abs of the weights
    EXP.classifResults.(featList{ff}).weightsStd = abs(squeeze(std(weights,[],3)));
    % Tials used for Test
    %--------------------
    EXP.trialUsed = unique([EXP.trialUsed]);
    
    if ~isempty(EXP.chanlocs)
        

        
        

            weights_mean{ff} = EXP.classifResults.(featList{ff}).weightsMean;
            weights_std{ff} = EXP.classifResults.(featList{ff}).weightsStd;

    else
        warning('YOUR DATASET DOES NOT CONTAIN ANY LOCATIONS TO ORGANIZE THE GRAPH')
        weights_mean{ff} = EXP.classifResults.(featList{ff}).weightsMean;
    end
    

end

% Delimitation
%-------------
try
    if strcmp(wheightsOrder,'area')
        for i=1:length(areas)
            limit = max(find(sort(or)==i));
            limit = limit+.5;
            med   = median(find(sort(or)==i));
            line(EXP.timesOfClassif([1 end]), [limit limit], 'Color','w', 'LineWidth',1, 'LineStyle', '--' )
            text(xlimit(2)+10,med,areas{i},'Color',[0,0,0],'FontSize',10);
        end
        % Electrode label
        %----------------
        for el = 1 : length(order)
            text(xlimit(2)+5,el, EXP.chanlocs.(validLocLabel){order(el)},'Color',[0,0,0],'FontSize',5);
        end
    end
catch
end

% Title
%-------
if strcmp(EXP.timeMode,'pointByPoint')
    name   = ['Subj' EXP.subject(1:3) ' pointByPoint classification weights, ' (featList{ff}) ' feature'];
else
    name   = ['Subj' EXP.subject(1:3) ' whole epoch classif weights, ' (featList{ff}) ' feature'];
end

for ff=1:length(featList)
    performanceCurves.classifPopulationSD{ff} = EXP.classifResults.(featList{ff}).accSD;
    performanceCurves.timesOfClassif{ff} = EXP.timesOfClassif;
    performanceCurves.classifCurve{ff} = EXP.classifResults.(featList{ff}).accMean;
    performanceCurves.classifPop{ff} = EXP.classifResults.(featList{ff}).accPop;
    
    performanceWeights.weights_mean{ff} = weights_mean{ff};
    performanceWeights.weights_std{ff} = weights_std{ff};
end

performance.performanceCurves = performanceCurves;
performance.performanceWeights = performanceWeights;

decInfo.Subj = EXP.subject(1:3);
decInfo.feats  = featList;
%decInfo.binLen = num2str(EXP.binSizes*1000/EXP.srate);
decInfo.chanlocs = EXP.chanlocs;
decInfo.nbChan = size(performanceWeights.weights_mean{1},1);
decInfo.name = name;
decInfo.xlimit = EXP.timeInter;

end
