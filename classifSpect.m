function classifSpect(monNom)

myDir = ['./sets/tfSets/' monNom '/'];
myResDir = './EXPfiles/';
myFiles = dir([myDir '*.mat']);

randMode = 'off';
timeMode = 'pointByPoint';

parfor ff=1:numel(myFiles)
    ff
    name = [myDir myFiles(ff).name];
    fName = [myFiles(ff).name(1:end-4)];
    classifThis(name, randMode, timeMode, myResDir,fName)
end

end



function classifThis(name, randMode, timeMode, myResDir, fName)
    load(name);
    myFields = fields(EEGf.featClassif);
    %EXP_pow = mainClassification(EEGf, [], myFields(1), randMode, timeMode, []);
    [EXPs grps] = electrodeSelector(EEGf, [], myFields(1), timeMode, [], 'group');
    %EXP_phs = mainClassification(EEGf, [], myFields(2), randMode, timeMode, []);
    if numel(EXPs) == 1
        save( [ myResDir fName '_pow.mat'], 'EXP_pow')
    else
        for ar=1:numel(EXPs)
            arDir = [myResDir fName(1:3) '_grouped/' grps{ar} '/'];
            if ~exist(arDir)
            mkdir(arDir)
            end
            EXP = EXPs{ar};
        save( [ arDir fName '_pow.mat'], 'EXP')
        end
    end
    
    %save( [ resName '_phs.mat'], 'EXP_phs')
end