% getSubjElec() - takes individual subj/monkey and returns all electrodes
%
% Usage: [elec] = getSubjElec(subj)
% Requierd input: patient number ('00#') or monkey name (currently 'dicky' or 'roucky')
%
% Outputs:
%    for patients : a cell array containing all the electrode locations in this format: 'hemisphere-lobe-gyrus-BA'
%    for monkeys  : a structure array containing locations as seperate
%       field names. All of these are hand written.
function [elecs] = getSubjElec (directory, subj)
if isempty(dir([directory '/*_locs.txt']))
    

        directory = '/Users/maximecauchoix/Documents/clusterCopy/data/monkey_human';

        %directory = '/gpfs/data/tserre/data/neural_readout/human';

    
    display(['the path you gave doesnt have any loc files in it, switching to: ' directory])
end
monkeys = {'dicky' 'roucky'};

if isempty(dir([directory '/' subj '_locs.txt']));
    if sum(ismember(monkeys, subj)) == 0
        error('no such location data or patient (or monkey)')
    end
end

if sum(ismember(monkeys, subj)) == 0 %%%%%%%%%%%%%%%%%%%%%%%%human parts
    fid = fopen([directory '/' subj '_locs.txt']);
    locations = textscan(fid, '%s %s %s %s %s %s', 'delimiter','	');
    fclose(fid);
    locs.name = subj;
    temp =  cell2struct(locations(1), {'hemis'},1 ); locs.hemis = temp.hemis;
    temp = cell2struct(locations(2), {'lobe'},1 );locs.lobe = temp.lobe;
    temp = cell2struct(locations(3), {'gyrus'},1 );locs.gyrus = temp.gyrus;
    temp = cell2struct(locations(4), {'matter'},1 );locs.matter = temp.matter;
    temp = cell2struct(locations(5), {'ba'},1 );locs.ba = temp.ba;
    temp = cell2struct(locations(6), {'distance'},1 );locs.distance = temp.distance;
    clear temp fid;
    
    elecs = {};
    
    for y=1:length(locs.hemis)
        if length(locs.hemis{y}) <=5 || length(locs.hemis{y}) == 20
            elecs{y} = 'inval';
        else
            g= locs.gyrus{y};
            g = g(regexp(g, '[A-Z]', 'start'));
            if length(locs.ba{y})>4
                elecs{y} = sprintf('%s-%s-%s-%s', locs.hemis{y}(1), locs.lobe{y}(1), g, locs.ba{y}(15:end));
            else
                elecs{y} = sprintf('%s-%s-%s-%s', locs.hemis{y}(1), locs.lobe{y}(1), g, '00');
            end
        end
        elecs = elecs';
    end
    
else %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%monkey parts
    if strcmpi('dicky', subj)
        dicky.v2 = [4 7 8];
        dicky.v4 = [2 3 5 6 9 11];
        dicky.pit = [1 10 12];
        dicky.fef = [13];
        dicky.lip = [ ];
        for i=1:EEG.nbchan % iterate over all electrodes
            fields = fieldnames(dicky);
            for r=1:length(fields)
                if sum(ismember(dicky.(fields{r}), i)) == 1
                    elecs{i} = fields{r};
                end
            end
        end
        elecs = elecs';
    end
    if strcmpi('roucky', subj)
        roucky.v2 = [8 14 15];
        roucky.v4 = [11 12 16 15 6 7 9 4];
        roucky.pit = [10 13 1 2 3 5];
        roucky.fef = [21 22];
        roucky.lip = [17 18 19 20];
        for i=1:EEG.nbchan % iterate over all electrodes
            fields = fieldnames(roucky);
            for r=1:length(fields)
                if sum(ismember(roucky.(fields{r}), i)) == 1
                    elecs{i} = fields{r};
                end
            end
        end
        elecs = elecs';
    end
end
end