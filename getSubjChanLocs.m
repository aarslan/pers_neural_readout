% getSubjChanLocs() - add chanlocs field to EEG
%
% Usage: [EEGnew] = getSubjChanLocs (EEG,indSubject)
%
% Requierd input:     meaning (size or options) {default}
%
%   EEG             = EEG structures
%   indSubject      = subject's indece (in subjectList lign 54 runClassif)
%
% Outputs:            meaning (size or options)
%
%   EEG             = EEG structures with electrode localisation in
%                     EEG.chanlocs
% 
% Autors: Ali, Serre Lab, 2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [EEGnew] = getSubjChanLocs (EEG,indSubject)

% Directory containning .loc depending on the machine
%----------------------------------------------------
if strcmp(computer,'MACI64')            % check if it's maxime's local machine
    if exist('/Users/maximecauchoix/Documents/clusterCopy/data/')
        directory = '/Users/maximecauchoix/Documents/clusterCopy/data/monkey_human/';
    else
        directory = '/Users/aarslan/Brown/Readout/ecog/Data-Analysis/';
    end
else                                           % cluster
    directory = '/gpfs/data/tserre/data/neural_readout/human';
end

subjectList = {'Dicky','Roucky','004','005','006'};

%-----------%
% For Human %
%-----------%
if indSubject>2
    fid = fopen([directory '/' subjectList{indSubject} '_locs.txt']);
    locations = textscan(fid, '%s %s %s %s %s %s', 'delimiter','	');
    fclose(fid);
    temp = cell2struct(locations(1), {'hemis'},1 );    locs.hemis    = temp.hemis;
    temp = cell2struct(locations(2), {'lobe'},1 );     locs.lobe     = temp.lobe;
    temp = cell2struct(locations(3), {'gyrus'},1 );    locs.gyrus    = temp.gyrus;
    temp = cell2struct(locations(4), {'matter'},1 );   locs.matter   = temp.matter;
    temp = cell2struct(locations(5), {'ba'},1 );       locs.ba       = temp.ba;
    temp = cell2struct(locations(6), {'distance'},1 ); locs.distance = temp.distance;
    clear temp fid;
    chanlocs = locs;
    
    %------------%
    % For Monkey %
    %------------%
    
else
    %Dicky
    %-----
    if indSubject==1
        % Manual localisation
        %--------------------
        % Anatomic criterium
        %-------------------
        dicky.v2    = [4 7 8];
        dicky.v4    = [2 3 5 6];
        dicky.pit   = [1 10 11];
        dicky.fef   = [13];
        dicky.lip   = [ ];
        dicky.unc   = [12 9];
        % Latency criterium
        %-------------------
%         dicky.low    = [6 7 8];
%         dicky.med    = [9 4 3];
%         dicky.hig    = [1 10 11];
%         dicky.unc   = setdiff(1:13,[dicky.low dicky.med dicky.hig]);
       
        
        for i=1:size(EEG.data,1) % iterate over all electrodes
            fields = fieldnames(dicky);
            for r=1:length(fields)
                if sum(ismember(dicky.(fields{r}), i)) == 1
                    chanlocs.area{i} = fields{r};
                end
            end
            chanlocs.labels{i} = EEG.chanlocs(i).labels;
        end
        chanlocs.area   = chanlocs.area';
        chanlocs.labels = chanlocs.labels';
        
        % Hemispher field
        %----------------
        dicky = [];
        dicky.left      = [1:6];
        dicky.right     = [7:12];
        dicky.frontal   = 13;

        for i=1:size(EEG.data,1) % iterate over all electrodes
            fields = fieldnames(dicky);
            for r=1:length(fields)
                if sum(ismember(dicky.(fields{r}), i)) == 1
                    chanlocs.hemi{i} = fields{r};
                end
            end
        end
        chanlocs.hemi   = chanlocs.hemi';
        
        
    end
    
    %Roucky
    %------
    if indSubject==2
        % Anatomic criterium
        %-------------------
        roucky.v2   = [8 14 15 12];
        roucky.v4   = [16 4 6 7];
        roucky.pit  = [2 11];
        roucky.fef  = [21 22];
        roucky.lip  = [17 18 19 20];
        roucky.unc  = [9 10 13 5 3 1];
        % Latency criterium
        %-------------------
%         roucky.low   = [8 14 15 12];
%         roucky.med   = [3 5 9 4];
%         roucky.hig   = [2 11 6 7];
%         roucky.unc   = setdiff(1:22,[roucky.low roucky.med roucky.hig]);


        
        for i=1:size(EEG.data,1) % iterate over all electrodes
            fields = fieldnames(roucky);
            for r=1:length(fields)
                if sum(ismember(roucky.(fields{r}), i)) == 1
                    chanlocs.area{i} = fields{r};
                end
            end
            chanlocs.labels{i} = EEG.chanlocs(i).labels;
        end
        chanlocs.area   = chanlocs.area';
        chanlocs.labels = chanlocs.labels';
        
        % Hemispher field
        %----------------
        roucky = [];
        roucky.left     = [1:9];
        roucky.right    = [10:16];
        roucky.parietal = [17:20];
        roucky.frontal  = [21:22];

        
        
        for i=1:size(EEG.data,1) % iterate over all electrodes
            fields = fieldnames(roucky);
            for r=1:length(fields)
                if sum(ismember(roucky.(fields{r}), i)) == 1
                    chanlocs.hemi{i} = fields{r};
                end
            end
        end
        chanlocs.hemi   = chanlocs.hemi';
        
    end
end
EEG.chanlocs = chanlocs;
EEGnew = EEG;
end