function allSpectPlotter
subjList = {'006'};


myEXPDir = '/gpfs/data/tserre/aarslan/NIPS_alter/EXPfiles/';
myResDir = '/gpfs/home/aarslan/scratch/NIPS_alter/results/';

for subjNo = 1:numel(subjList)
    
    subjName = subjList{subjNo};
    mkdir([myResDir subjName]);
    myStemDir = [subjName '_grouped/'];
    areaFolders = dir([myEXPDir myStemDir]);
    areaFolders = {areaFolders.name};
    areaFolders = areaFolders(3:end);
    [r c] = calcSubplot(numel(areaFolders));
    %hands = zeros(1,numel(areaFolders));
    
    for areaNo=1:numel(areaFolders)
        area = areaFolders{areaNo};
        
        hands(areaNo) = spectDecPlotter(subjName, area);
        title(hands(areaNo).main, regexprep(area, '_', '\\_'))
                print(hands(areaNo).fig, [myResDir subjName '/' area] , '-djpeg')
                print(hands(areaNo).fig, [myResDir subjName '/' area] , '-dpdf')
        close(hands(areaNo).fig)
    end
%     
%     for areaNo=1:numel(areaFolders)
%         area = areaFolders{areaNo};
%         
%     end
    
end


end

function [r c] = calcSubplot(numFold)

r = ceil(sqrt(numFold));
c = ceil(sqrt(numFold));
end