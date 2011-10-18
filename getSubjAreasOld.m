% getSubjAreas() - %% takes a directory where the talairach locations for the patients are,
% and returns the structure containing electrodes that fall into the regions you specify as seperate
% strings. The field named ".all" has all the electrode numbers of interest without
% duplicates.
%
% functional areas (v1, v2 etc) are hard coded, so it's limited to the
% assignments I made with BA's (therefore not pretty accurate). Abbrevations for anatomical locations (gyri) are
% read from the talairach labels.

%
% Usage: [areasStruct] = getSubjChanLocs (EEG,indSubject)
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
% Autors: Ali, Maxime, Serre Lab, 2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [areasStruct] = getSubjAreas (directory, varargin)
if isempty(dir([directory '/*_locs.txt']))
    directory = '/gpfs/data/tserre/data/neural_readout/human';
    display(['the path you gave doesnt have any loc files in it, switching to: ' directory])
end

funcLabels = {'v1', 'v2', 'v3', 'v4', 'v8', 'mt'};
funcLabelsDict.v1 = {'Brodmann area 17'};
funcLabelsDict.v2 = {'Brodmann area 18'};
funcLabelsDict.v3 = {'Brodmann area 18'};
funcLabelsDict.v4 = {'Brodmann area 19'};
funcLabelsDict.v8 = {'Brodmann area 19'};
funcLabelsDict.mt = {'Brodmann area 19'};

locFiles = dir([directory '/*_locs.txt']);

for theFile=1:size(locFiles,1)
    fid = fopen([directory  '/' locFiles(theFile).name]);
    locations = textscan(fid, '%s %s %s %s %s %s', 'delimiter','	');
    fclose(fid);
    locs(theFile).name = locFiles(theFile).name;
    temp =  cell2struct(locations(1), {'hemis'},1 ); locs(theFile).hemis = temp.hemis;
    temp = cell2struct(locations(2), {'lobe'},1 );locs(theFile).lobe = temp.lobe;
    temp = cell2struct(locations(3), {'gyrus'},1 );locs(theFile).gyrus = temp.gyrus;
    temp = cell2struct(locations(4), {'matter'},1 );locs(theFile).matter = temp.matter;
    temp = cell2struct(locations(5), {'ba'},1 );locs(theFile).ba = temp.ba;
    temp = cell2struct(locations(6), {'distance'},1 );locs(theFile).distance = temp.distance;
    clear temp fid;
end

currFuncLabels = varargin(ismember(varargin, funcLabels));
currAnatLabels = varargin(~ismember(varargin, funcLabels));
currBALabels = [];

for y=1:length(currAnatLabels)
    bik = currAnatLabels(y);
    if strcmpi(bik{1}(1:2), 'ba')
        currBALabels = [currBALabels currAnatLabels(y)];
    end
end

elecsArea = [];

%% this loop handles anatomical landmark (gyri) assignments
for lab=1:length(currAnatLabels) 
    for pat=1:size(locs,2)
        areasStruct(pat).name = locs(pat).name(1:3); %put name in the patient's structure
        for elec=1:size(locs(pat).gyrus, 1)
            s = locs(pat).gyrus{elec};
            s = s(regexp(s, '[A-Z]', 'start'));
            if strcmpi(s, currAnatLabels(lab))
                elecsArea = [elecsArea elec];
                areasStruct(pat).(currAnatLabels{lab}) = elecsArea;
            end
        end
        elecsArea = [];
    end
end

%% this loop handles Brodmann Area assignments
for lab=1:length(currBALabels)
    for pat=1:size(locs,2)
        areasStruct(pat).name = locs(pat).name(1:3); %put name in the patient's structure
        for elec=1:size(locs(pat).gyrus, 1)
            elec
            pat
            if length(locs(pat).ba{elec})>4 %strcmpi(locs(pat).ba{elec}(1), 'B') && isempty(locs(pat).ba{elec}(1))
                    s = locs(pat).ba{elec};
                    s
                    s = s([1 10 15:end]);
                    if strcmpi(s, currBALabels(lab))
                        elecsArea = [elecsArea elec];
                        areasStruct(pat).(currBALabels{lab}) = elecsArea;
                    end
            end
        end
        elecsArea = [];
    end
end

%% this loop handles functional landmark assignments using the dictionary
%%(funcLabelsDict) above
for lab=1:length(currFuncLabels) 
    for pat=1:size(locs,2)
        areasStruct(pat).name = locs(pat).name(1:3); %put name in the patient's structure
        for elec=1:size(locs(pat).gyrus, 1)
            s = locs(pat).ba{elec};
            if sum(strcmpi(s, funcLabelsDict.(currFuncLabels{lab}))) %find if there's a match in any where of the cell array
                elecsArea = [elecsArea elec];
                areasStruct(pat).(currFuncLabels{lab}) = elecsArea;
            end
        end
        elecsArea = [];
    end
end



%% monkey assignments
% these assignments are made by hand, desired locations are read from the
% varargin and the corresponding regions are written by hand

areasStruct(size(locs,2)+1).name = 'dicky';
areasStruct(size(locs,2)+2).name = 'roucky';

if ismember('v2', varargin)
    areasStruct(size(locs,2)+1).v2 = [4 7 8];
    areasStruct(size(locs,2)+2).v2 = [14 15];
end

if ismember('v4', varargin)
    areasStruct(size(locs,2)+1).v4 = [2 3 5 6 9 11];
    areasStruct(size(locs,2)+2).v4 = [11 12 16 15 6 7 4];
end

if ismember('pit', varargin)
    areasStruct(size(locs,2)+1).pit = [1 10 12];
    areasStruct(size(locs,2)+2).pit = [10 13 1 2 3 5];
end

if ismember('fef', varargin)
    areasStruct(size(locs,2)+1).fef = [13];
    areasStruct(size(locs,2)+2).fef = [21 22];
end

if ismember('lip', varargin)
    areasStruct(size(locs,2)+1).lip = [ ];
    areasStruct(size(locs,2)+2).lip = [17 18 19 20];
end

%% this loop puts every electrode in one array in the field '.all' for
% convenience (without duplicate electrode numbers)
for pat=1:length(areasStruct)
    all = [];
    fields = fieldnames(areasStruct(pat));
    for field=2:length(fieldnames(areasStruct(pat)))
        all = [all areasStruct(pat).(fields{field})];
        all = unique(all);
        areasStruct(pat).all = all; 
    end
end

end