function [weights_mean performanceCurves decInfo] = mainClassificationPlotPerf(EXP, EXPrand, graphMode)
% EXP:      actual classification
% EXPrand : classification with randomized labels
% graphMode : 0 = both weigths and perf
%           1 = only perf
%           3 = only weights

ALLEXP = [EXP EXPrand];
figure

%Initialize
%----------
perfValue = 'accuracy'; %hardcoded for convenience
wheightsOrder = 'none'; %'contribution'; %'area';
errorB = sqrt(ALLEXP(1).genExpParam.crossValid); %this can also be 1 depending on your error bar calc
weights   = [];
trial     = 1:size(ALLEXP(1).selector,2);
Nbin = 7;
featList = fields(EXP.classifResults);

%subplot parameters
if (graphMode == 0)
    nRows = 2; % bottom: weights, top: performance
else
    nRows = 1;
end
nCol = length(fields(EXP.classifResults));
weightHandles = [];
accuHandles = [];

validLocLabels = {'area', 'gyrus', 'labels'};
for ii=1:length(validLocLabels)
    if isfield(ALLEXP(1).chanlocs, validLocLabels{ii})
        validLocLabel = validLocLabels{ii};
        plotLabels = 0;
    end
end

%loop over each feature classified
for ff=1:length(featList)
    
    %loop over cross validation
    %--------------------------
    for i=1:ALLEXP(1).genExpParam.crossValid
        % Loop over EXP (for Rand and Cond2)
        %-----------------------------------
        for iexp = 1:length(ALLEXP)
            if i==1
                % Initialisation
                %---------------
                ALLEXP(iexp).trialUsed  = [];
            end
            % For accuracy
            %-------------
            
            mode = fields(EXP.classifResults.(featList{ff})); mode = mode{:};
            if strcmp(perfValue,'accuracy')
                ALLEXP(iexp).classifResults.(featList{ff}).acc(i,:)     = [ALLEXP(iexp).classifResults.(featList{ff}).(mode).perfTest(i,:).accuracy]*100;
                ALLEXP(iexp).classifResults.(featList{ff}).accRand(i,:) = [ALLEXP(iexp).classifResults.(featList{ff}).(mode).perfTest(i,:).BSthreshold]*100;
                chance = 50;
                % For Auc
                %--------
            elseif strcmp(perfValue,'auc')
                ALLEXP(iexp).classifResults.(featList{ff}).acc(i,:)     = [ac;ALLEXP(iexp).classifResults.(featList{ff}).pointByPoint.perfTest(i,:).auc];
                ALLEXP(iexp).classifResults.(featList{ff}).accRand(i,:) = [accRan;ALLEXP(iexp).classifResults.(featList{ff}).pointByPoint.perfTest(i,:).BSthreshold];
                chance = 0.5;
            end
            % Find ind of trials used for test
            %---------------------------------
            ALLEXP(iexp).trialUsed = [ALLEXP(iexp).trialUsed, trial(ALLEXP(iexp).selector(i,:)==2)];
        end
        
        % For weights
        %------------
        wei       = [];
        if strcmp(ALLEXP(1).timeMode,'pointByPoint')
            wei = [ALLEXP(1).classifResults.(featList{ff}).pointByPoint.perfTest(i,:).B];
            wei = wei(1:(size(wei,1)-1),:);
        elseif strcmp(ALLEXP(1).timeMode,'allEpoch')
            wei = [ALLEXP(1).classifResults.(featList{ff}).allEpoch.perfTest(i).B];
            wei = wei(1:(length(wei)-1)); % remove biais
            wei = reshape(wei,length(ALLEXP(1).chanlocs.(validLocLabel)),length(ALLEXP(1).timesOfClassif));
        end
        weights = cat(3,weights,wei);
    end
    
    for iexp = 1:length(ALLEXP)
        % Average performances on all the CV
        %-----------------------------------
        ALLEXP(iexp).classifResults.(featList{ff}).accMean = mean(ALLEXP(iexp).classifResults.(featList{ff}).acc,1);
        ALLEXP(iexp).classifResults.(featList{ff}).accSD   = std(ALLEXP(iexp).classifResults.(featList{ff}).acc,[],1)/errorB;
        % Random thershold computed on test only
        %---------------------------------------
        ALLEXP(iexp).classifResults.(featList{ff}).accRandMean = mean(ALLEXP(iexp).classifResults.(featList{ff}).accRand,1);
        ALLEXP(iexp).classifResults.(featList{ff}).accRandSD   = std(ALLEXP(iexp).classifResults.(featList{ff}).accRand,[],1)/errorB;
        % Weights
        %---------
        ALLEXP(iexp).classifResults.(featList{ff}).weightsMean = abs(squeeze(mean(weights,3))); % we take the abs of the weights
        % Tials used for Test
        %--------------------
        ALLEXP(iexp).trialUsed = unique([ALLEXP(iexp).trialUsed]);
        
    end
    
    plotLabels = 1;
    validLocLabels = {'area', 'gyrus', 'labels'};
    xlimit = EXP.timeInter;
    
    if ~isempty(ALLEXP(1).chanlocs)
        for ii=1:length(validLocLabels)
            if isfield(ALLEXP(1).chanlocs, validLocLabels{ii})
                validLocLabel = validLocLabels{ii};
                plotLabels = 0;
            end
        end
        
        for i=1:size(ALLEXP(1).chanlocs,2)
            if ~iscell(ALLEXP(1).chanlocs(i).(validLocLabel))
                ALLEXP(1).chanlocs(i).(validLocLabel) = {ALLEXP(1).chanlocs(i).(validLocLabel)};
            end
        end
        
        elec   = [ALLEXP(1).chanlocs.(validLocLabel)];
        elec   = elec(~strcmp(elec,'NULL'));
        areas  = unique(elec);
        
        if strcmp(wheightsOrder,'contribution')
            weights_meanm = mean (ALLEXP(1).classifResults.(featList{ff}).weightsMean,2);% work on workInt
            [~,order]  = sort(weights_meanm,'descend');
        elseif strcmp(wheightsOrder,'area')
            or    = [];
            order = [];
            el    = 1 :length(elec);
            for ar=1:length(areas)
                or(strcmp(elec,areas{ar})) = ar;
            end
            for ar=1:length(areas)
                order = [order, el(or==ar)];
            end
        else
            order = 1:length(elec);
        end
        
        if strcmp(EXP.chooseFeatParam.area,'all')
            weights_mean = ALLEXP(1).classifResults.(featList{ff}).weightsMean(order,:);
        else
            weights_mean = ALLEXP(1).classifResults.(featList{ff}).weightsMean;
        end
    else
        warning('YOUR DATASET DOES NOT CONTAIN ANY LOCATIONS TO ORGANIZE THE GRAPH')
        weights_mean = ALLEXP(1).classifResults.(featList{ff}).weightsMean;
    end
    
    % Zscore by colon for plot
    %-------------------------
    if strcmp(ALLEXP(1).timeMode,'pointByPoint')
        weights_mean = zscore(weights_mean);
    end
    
    % Plot only weights in xlimit
    %----------------------------
    if (graphMode == 0) || (graphMode == 3)
        subplot(nRows,nCol, ff)
        imagesc(ALLEXP(1).timesOfClassif,1:size(weights_mean,1),weights_mean)
        if strcmp(wheightsOrder,'area')
            for i=1:length(areas)
                limit = max(find(sort(or)==i));
                limit = limit+.5;
                med   = median(find(sort(or)==i));
                hold on
                line(ALLEXP(1).timesOfClassif([1 end]), [limit limit], 'Color','w', 'LineWidth',1, 'LineStyle', '--' )
                
            end
        end
        hold off
    end
end

if (graphMode == 0) || (graphMode == 3)
    % Delimitation
    %-------------
    try
        if strcmp(wheightsOrder,'area')
            for i=1:length(areas)
                limit = max(find(sort(or)==i));
                limit = limit+.5;
                med   = median(find(sort(or)==i));
                hold on
                line(ALLEXP(1).timesOfClassif([1 end]), [limit limit], 'Color','w', 'LineWidth',1, 'LineStyle', '--' )
                text(xlimit(2)+10,med,areas{i},'Color',[0,0,0],'FontSize',10);
            end
            % Electrode label
            %----------------
            if plotLabels
                for el = 1 : length(order)
                    text(xlimit(2)+5,el, EXP.chanlocs.(validLocLabel){order(el)},'Color',[0,0,0],'FontSize',5);
                end
            end
        end
    catch
    end
    
    % Plot labels
    %------------
    xlabel('Time (ms)')
    % Title
    %-------
    if strcmp(ALLEXP(1).timeMode,'pointByPoint')
        name   = ['Subj' ALLEXP(1).subject(1:3) ' pointByPoint classification weights, ' (featList{ff}) ' feature, bin ' num2str(ALLEXP(1).binSizes*1000/ALLEXP(1).srate) ' ms'];
    else
        name   = ['Subj' ALLEXP(1).subject(1:3) ' whole epoch classif weights, ' (featList{ff}) ' feature, bin ' num2str(ALLEXP(1).binSizes*1000/ALLEXP(1).srate) ' ms'];
    end
    
    title(name)
    
end
%-----------------------------------------
col1    = [1,0,0];
colRand = [0,1,0];
col2    = [0,0,1];




for ff=1:length(featList)
    
    % Random classification accuracy curve
    %--------------------------------------
    performanceCurves.randPopulation{ff, 1} = [ALLEXP(1).timesOfClassif,ALLEXP(1).timesOfClassif(end:-1:1)];
    performanceCurves.randPopulation{ff, 2} = [ALLEXP(2).classifResults.(featList{ff}).accMean-ALLEXP(2).classifResults.(featList{ff}).accSD,...
        ALLEXP(2).classifResults.(featList{ff}).accMean(end:-1:1)+ALLEXP(2).classifResults.(featList{ff}).accSD(end:-1:1)];
    
    performanceCurves.randCurve{ff, 1} = ALLEXP(1).timesOfClassif;
    performanceCurves.randCurve{ff, 2} = ALLEXP(2).classifResults.(featList{ff}).accMean;
    
    if (graphMode == 0) || (graphMode == 1)
        accuHandles(ff) = subplot(nRows,nCol, ff);%% spaghetti coding change length(featList)+ ff);
        % Classif Value accuracy curve
        %-----------------------------
        performanceCurves.classifPopulation{ff, 1} = [ALLEXP(1).timesOfClassif,ALLEXP(1).timesOfClassif(end:-1:1)];
        performanceCurves.classifPopulation{ff, 2} = [ALLEXP(1).classifResults.(featList{ff}).accMean-ALLEXP(1).classifResults.(featList{ff}).accSD,...
            ALLEXP(1).classifResults.(featList{ff}).accMean(end:-1:1)+ALLEXP(1).classifResults.(featList{ff}).accSD(end:-1:1)];
        
        performanceCurves.classifCurve{ff, 1} = ALLEXP(1).timesOfClassif;
        performanceCurves.classifCurve{ff, 2} = ALLEXP(1).classifResults.(featList{ff}).accMean;
        
        
        %Plot both curves
        fill(performanceCurves.randPopulation{ff, 1}, performanceCurves.randPopulation{ff, 2},colRand), hold on;
        %alpha(0.2);
        plot(performanceCurves.randCurve{ff, 1}, performanceCurves.randCurve{ff, 2},'Color',colRand,'LineWidth',2);
        
        fill(performanceCurves.classifPopulation{ff, 1}, performanceCurves.classifPopulation{ff, 2},col1), hold on;
        %alpha(0.2);
        plot(performanceCurves.classifCurve{ff, 1}, performanceCurves.classifCurve{ff, 2},'Color',col1,'LineWidth',2);

        % Chance level
        %--------------
        plot(ALLEXP(1).timesOfClassif,repmat(chance,length(ALLEXP(1).timesOfClassif),1),'--k')
        
        %-----------%
        % workInt   %
        %-----------%
        indt1 = max(find(ALLEXP(1).timesOfClassif<0));
        indt2 = min(find(ALLEXP(1).timesOfClassif>600));
        t1    = ALLEXP(1).timesOfClassif(indt1);
        t2    = ALLEXP(1).timesOfClassif(indt2);
        
        %---------%
        % PlotMax %
        %---------%
        % Loop over EXP (for Rand and Cond2)
        %-----------------------------------
        
        for iexp = 1:length(ALLEXP)
            [ALLEXP(iexp).classifResults.(featList{ff}).Max,inLat] = max(ALLEXP(iexp).classifResults.(featList{ff}).accMean(indt1:indt2));
            ALLEXP(iexp).classifResults.(featList{ff}).indMax      = indt1+inLat-1;
            ALLEXP(iexp).classifResults.(featList{ff}).latMax      = ALLEXP(1).timesOfClassif(ALLEXP(iexp).classifResults.(featList{ff}).indMax);
        end
        % Plot a cross : condition1
        %--------------------------
        plot(repmat(ALLEXP(1).classifResults.(featList{ff}).latMax,1,2),...
            [ALLEXP(1).classifResults.(featList{ff}).accMean(ALLEXP(1).classifResults.(featList{ff}).indMax)-2....
            ,ALLEXP(1).classifResults.(featList{ff}).accMean(ALLEXP(1).classifResults.(featList{ff}).indMax)+2],'k','LineWidth',1.2)
        plot(ALLEXP(1).timesOfClassif(ALLEXP(1).classifResults.(featList{ff}).indMax)-4:ALLEXP(1).timesOfClassif(ALLEXP(1).classifResults.(featList{ff}).indMax)+4,...
            repmat(ALLEXP(1).classifResults.(featList{ff}).Max,1,9),'k','LineWidth',1.2)
        text(ALLEXP(1).classifResults.(featList{ff}).latMax + 6,ALLEXP(1).classifResults.(featList{ff}).Max + 2,...
            [num2str(round(ALLEXP(1).classifResults.(featList{ff}).Max)) ' % at ' num2str(ALLEXP(1).classifResults.(featList{ff}).latMax) ' ms' ],'Color',col1 )
        
        %------------%
        % plotLatSig %
        %------------%
        % Loop over EXP (for Rand and Cond2)
        %-----------------------------------
        for iexp = 1:length(ALLEXP)
            
            % Test to decide when classif is different from rand bedore train
            %----------------------------------------------------------------
            for t=1:length(EXP.timesOfClassif)
                [h,pvals] = ttest2(ALLEXP(iexp).classifResults.(featList{ff}).acc(:,t),ALLEXP(2).classifResults.(featList{ff}).acc(:,t));
                p(t) = pvals;
            end
            ALLEXP(iexp).classifResults.(featList{ff}).indLat = min(find(conv(double(p<0.05),ones(Nbin,1),'same')==Nbin))-floor(Nbin/2);
            
            % Indices
            %--------
            ALLEXP(iexp).classifResults.(featList{ff}).lat    = ALLEXP(1).timesOfClassif(ALLEXP(iexp).classifResults.(featList{ff}).indLat);
            ALLEXP(iexp).classifResults.(featList{ff}).latVal = ALLEXP(1).classifResults.(featList{ff}).accMean(ALLEXP(iexp).classifResults.(featList{ff}).indLat);
            
        end
        % Plot a line at the significance time : Condition1
        %--------------------------------------------------
        plot(repmat(ALLEXP(1).classifResults.(featList{ff}).lat-1,1,2),...
            [ALLEXP(1).classifResults.(featList{ff}).accMean(ALLEXP(1).classifResults.(featList{ff}).indLat)-2,...
            ALLEXP(1).classifResults.(featList{ff}).accMean(ALLEXP(1).classifResults.(featList{ff}).indLat)+2],'K','LineWidth',1.2)
        text(ALLEXP(1).classifResults.(featList{ff}).lat + 6,...
            ALLEXP(1).classifResults.(featList{ff}).latVal - 1,...
            [ num2str(ALLEXP(1).classifResults.(featList{ff}).lat) ' ms' ], 'Color', col1 )
        
        set(gca, 'ytick', [])
        if max(performanceCurves.randCurve{ff, 2}) >75
        ylim([45 90])
        else
            ylim ([45 80])
        end
        xlim(xlimit)
        % Title
        %-------
        name   = ['Subj ' ALLEXP(1).subject(1:3) ' ' (featList{ff}) ' feature, bin ' num2str(ALLEXP(1).binSizes*1000/ALLEXP(1).srate) ' ms'];
        title(name)
        xlabel('Time (ms)')
    end
end

if (graphMode == 0) || (graphMode == 1)
    set(accuHandles(1), 'ytickMode', 'auto')
    if strcmp(perfValue,'accuracy')
        ylabel(accuHandles(1) , 'Accuracy (%)')
        
    elseif strcmp(perfValue,'auc')
        ylabel(ccuHandles(1) , 'AUC')
    end
end

decInfo.Subj = ALLEXP(1).subject(1:3);
decInfo.feats  = (featList{ff});
decInfo.binLen = num2str(ALLEXP(1).binSizes*1000/ALLEXP(1).srate);
decInfo.chanlocs = {EXP.chanlocs.(validLocLabel){order}};
decInfo.nbChan = numel(elec); 

end