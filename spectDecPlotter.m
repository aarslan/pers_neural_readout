function hand = spectDecPlotter(subjName, area)
myDir = './sets/tfSets/';
myEXPDir = '/gpfs/home/aarslan/work/NIPS_alter/EXPfiles/';
myStemDir = [subjName '_grouped/'];
myRandStemDir = '007_grouped_rand/';
myAreaDir = area; %'G_oc-temp_med-Parahip/';
myResDir = [myEXPDir myStemDir myAreaDir '/'];
myRandResDir = myResDir; %[myEXPDir myRandStemDir myAreaDir];

myFiles = dir([myResDir '*pow.mat']);
mrRandFiles = dir([myRandResDir '*pow.mat']);

enableRandTest = 0;
regularT = 0;

myPerfGraph_pow = [];
myPerfPop = [];
myPerfOrd_pow = [];

for fi=1:numel(myFiles)
    load([myResDir myFiles(fi).name])
    myFiles(fi).name
    name = regexp( myFiles(fi).name, '_', 'split');
    myFreqs(fi) = str2num(name{3});
    [performance decInfo] = mainClassificationPrepResults(EXP);
    myPerfGraph_pow(fi,:) = [performance.performanceCurves.classifCurve{:}];
    %myWeiGraph(fi,:,:) = [performance.performanceWeights.weights_mean];
    
    %myPerfPop(fi,:,:) = performance.performanceCurves.classifPop{:};
    
    if enableRandTest && (numel(mrRandFiles) ~= 0) && regularT
        load([myRandResDir myFiles(fi).name])
        [performance decInfo] = mainClassificationPrepResults(EXP);
        myPerfPop_rand(fi,:,:) = performance.performanceCurves.classifPop{:};
    end
    
end

[new ord] = sort(myFreqs);
myPerfOrd_pow = myPerfGraph_pow(ord,:);

clims = [45 75];
hand = figure; 
imagesc(EXP.timesOfClassif, new, myPerfOrd_pow, clims)
xlabel('time in ms')
ylabel('frequency in Hz')
set(gca, 'YDir', 'normal')

% figure
% myFiles = dir([myResDir '*phs.mat']);
% 
% for fi=1:numel(myFiles)
%     load([myResDir myFiles(fi).name])
%     myFiles(fi).name
%     name = regexp( myFiles(fi).name, '_', 'split');
%     myFreqs(fi) = str2num(name{3});
%     [performance decInfo] = mainClassificationPrepResults(EXP);
%     myPerfGraph_phs(fi,:) = [performance.performanceCurves.classifCurve{:}];
%    %myWeiGraph(fi,:,:) = performance.performanceWeights.performance.performanceWeights
% end
% [new ord] = sort(myFreqs);
% myPerfOrd_phs = myPerfGraph_phs(ord,:);
% 
% imagesc(EXP_pow.timesOfClassif, new, myPerfOrd_phs)
% 
% end
