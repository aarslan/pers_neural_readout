function [EXP varargout] = electrodeSelector(EEG, paramStruct, feats, timeMode, outPath, electrodes, manualInd )
%this is a wrapper function for mainClassification. What it does is to pick
%out electrodes in certain ways and runs the classification on these subset
%of electrodes (or individual electrodes) separately.

%ARGUMENTS:
%For the first few arguments, please see the argument descriptions in the mainClassification function
%For the rest of the arguments:
%   electrodes: is a string with a few options:
%               'individual': run each electrode separately
%               'group' : look at channel locs, ask user which grouping
%               structure to use, and run classification on each group of
%               electrodes
%               'subset': manually provide the electrode indices in the
%               next argument (as an array), all of them to be run in one
%               go

%determines if the EXP files will be saved to outPath as well (figures are saved in their individual directories)
saveSwitch = 0;

%get names of the electrode labels
try
    allGroups = cellfun(@(x) [x ' '], fieldnames(EEG.chanlocs)', 'UniformOutput', false);
catch
    warning('set does not have channel locations, run this function with "individual" option')
end
if strcmp(electrodes, 'group')
    taxonomy = 'regions'
else
    
taxonomy = input(['found these fields in chanlocs, pick one for grouping: ' [allGroups{:}] '\n'], 's');
end
if ~iscell([EEG.chanlocs.(taxonomy)]) %replace the selected labels with cells if they're not already
    temp = {EEG.chanlocs.(taxonomy)};
    EEG = rmfield(EEG, 'chanlocs');
    EEG.chanlocs.(taxonomy) = temp;
end
groups = unique(EEG.chanlocs.(taxonomy));

for ff=1:numel(feats)
    
    switch electrodes
        case 'individual'
            pathName = [outPath '/indiv_elec_deco/' feats{ff} '/'];
            mkdir(pathName)
            for ee=1:EEG.nbchan
                tempEEG = sliceSet(EEG, ee);
                if saveSwitch; decodeOutPath = pathName; else decodeOutPath = 0; end
                EXP  =  mainClassification(tempEEG, paramStruct, feats(ff), 'off', timeMode, decodeOutPath);
                rand =  mainClassification(tempEEG, paramStruct, feats(ff), 'on', timeMode, decodeOutPath);
                h = figure;
                mainClassificationPlotPerf(EXP, rand);
                if ~isempty(EEG.chanlocs); loc = [EEG.chanlocs.(taxonomy){ee}]; else loc = []; end
                saveas(h, [pathName 'chn_' num2str(ee) loc '_' feats{ff}], 'pdf')
                saveas(h, [pathName 'chn_' num2str(ee) loc '_' feats{ff}], 'epsc')
                close(h);
            end
        case 'group'
            %put a conditional here to direct user to someother case
            
            for gg=1:numel(groups)
                grpInd = strcmpi(EEG.chanlocs.(taxonomy), groups{gg});
                tempEEG = sliceSet(EEG, grpInd);
                if saveSwitch;
                    decodeOutPath = pathName;
                    pathName = [outPath '/group_chn_deco/' feats{ff} '/'];
                    mkdir(pathName);
                else
                    decodeOutPath = 0;
                end
                EXP{gg}  =  mainClassification(tempEEG, paramStruct, feats(ff), 'off', decodeOutPath);
                %rand =  mainClassification(tempEEG, paramStruct, feats(ff), 'on', timeMode, decodeOutPath);
                
                %[rien performanceCurves{gg} decInfo{gg}] = mainClassificationPlotPerf(EXP, rand, 1);
                %h = max(get(0, 'child'));
                %saveas(h, [pathName 'grp_' groups{gg} '_' feats{ff} 'numEl' num2str(decInfo{gg}.nbChan)], 'pdf')
                %saveas(h, [pathName 'grp_' groups{gg} '_' feats{ff} 'numEl' num2str(decInfo{gg}.nbChan)], 'epsc')
                %close(h);
                
            end
            varargout{1} = groups;
        case 'subset'
            pathName = [outPath '/subset_chn_deco/' feats{ff} '/'];
            mkdir(pathName);
            tempEEG = sliceSet(EEG, manualInd);
            if saveSwitch; decodeOutPath = pathName; else decodeOutPath = 0; end
            EXP  =  mainClassification(tempEEG, paramStruct, feats(ff), 'off', timeMode, decodeOutPath);
            rand =  mainClassification(tempEEG, paramStruct, feats(ff), 'on', timeMode, decodeOutPath);
            %h = figure;
            %mainClassificationPlotPerf(EXP, rand);
            %saveas(h, [pathName 'chn_' num2str(ee) '_' feats{ff}], 'pdf')
            %saveas(h, [pathName 'chn_' num2str(ee) '_' feats{ff}], 'epsc')
            %close(h);
    end
end
end

function newEEG = sliceSet(EEG, elecList)
newEEG = EEG;
newEEG.nbchan = numel(elecList);
newEEG.data = newEEG.data(elecList, :,:);
fieldz = fieldnames(newEEG.chanlocs);
for aa=1:numel(fieldz) %go over each field and pick out the chanlocs of interest
    tomp = newEEG.chanlocs.(fieldz{aa});
    newEEG.chanlocs.(fieldz{aa}) = tomp(elecList);
end

allFields = fieldnames(newEEG.featClassif);
for aa=1:numel(allFields) %go over each field and pick out the chanlocs of interest
    tomp = newEEG.featClassif.(allFields{aa});
    newEEG.featClassif.(allFields{aa}) = tomp(elecList,:,:);
end


end