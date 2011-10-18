function allSpectPlotter
subjList = {'006', '007', '008' };


myEXPDir = '/gpfs/home/aarslan/work/NIPS_alter/EXPfiles/';
myResDir = '/gpfs/home/aarslan/work/NIPS_alter/results/';

for subjNo = 1:numel(subjList)
    
    subjName = subjList{subjNo};
    mkdir([myResDir subjName]);
    myStemDir = [subjName '_grouped/'];
    areaFolders = dir([myEXPDir myStemDir]);
    areaFolders = {areaFolders.name};
    areaFolders = areaFolders(3:end);
    [r c] = calcSubplot(numel(areaFolders));
    hands = zeros(1,numel(areaFolders));
    
    for areaNo=1:numel(areaFolders)
        area = areaFolders{areaNo};
        
        hands(areaNo) = spectDecPlotter(subjName, area);
        title(area)
    end
    
    for areaNo=1:numel(areaFolders)
        area = areaFolders{areaNo};
        print(hands(areaNo), [myResDir subjName '/' area] , '-djpeg')
        close(hands(areaNo))
    end
    
end


end

function [r c] = calcSubplot(numFold)

r = ceil(sqrt(numFold));
c = ceil(sqrt(numFold));
end