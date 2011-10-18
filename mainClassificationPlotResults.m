function mainClassificationPlotResults(performance, decInfo)

zsc = 1;
xlimit = decInfo.xlimit;

validLocLabels = {'area', 'gyrus', 'labels'};
for ii=1:length(validLocLabels)
    if isfield(decInfo.chanlocs, validLocLabels{ii})
        validLocLabel = validLocLabels{ii};
    end
end

wheightsOrder = 'area';
elec   = [decInfo.chanlocs.(validLocLabel)];
elec   = elec(~strcmp(elec,'NULL'));
areas  = unique(elec);


for ff=1:numel(decInfo.feats);
    
    timesOfClassif = performance.performanceCurves.timesOfClassif{ff};
    weights_mean = performance.performanceWeights.weights_mean{ff};
    
    
    % Zscore by colon for plot
    %-------------------------
    if zsc
        weights_mean = zscore(weights_mean);
    end
    
    
    if strcmp(wheightsOrder,'contribution')
        weights_meanm = mean (EXP.classifResults.(featList{ff}).weightsMean,2);% work on workInt
        [~,order]  = sort(weights_meanm{ff},'descend');
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
    
    
    for i=1:size(decInfo.chanlocs,2)
        if ~iscell(decInfo.chanlocs(i).(validLocLabel))
            decInfo.chanlocs(i).(validLocLabel) = {decInfo.chanlocs(i).(validLocLabel)};
        end
    end
    
    
    %----------------------------
    imagesc(timesOfClassif,1:size(weights_mean,1),weights_mean)
    if strcmp(wheightsOrder,'area')
        for i=1:length(areas)
            limit = max(find(sort(or)==i));
            limit = limit+.5;
            med   = median(find(sort(or)==i));
            line(timesOfClassif([1 end]), [limit limit], 'Color','w', 'LineWidth',1, 'LineStyle', '--' )
        end
    end
    
    % Delimitation
    %-------------
    
    if strcmp(wheightsOrder,'area')
        for i=1:length(areas)
            limit = max(find(sort(or)==i));
            limit = limit+.5;
            med   = median(find(sort(or)==i));
            line(timesOfClassif([1 end]), [limit limit], 'Color','w', 'LineWidth',1, 'LineStyle', '--' )
            text(xlimit(2)+10,med,areas{i},'Color',[0,0,0],'FontSize',10);
        end
        % Electrode label
        %----------------
        for el = 1 : length(order)
            text(xlimit(2)+5,el, decInfo.chanlocs.(validLocLabel){order(el)},'Color',[0,0,0],'FontSize',5);
        end
    end
    
end

end