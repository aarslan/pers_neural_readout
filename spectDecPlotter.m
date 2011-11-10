function hand = spectDecPlotter(subjName, area)
myDir = './sets/tfSets/';
myEXPDir = '/gpfs/data/tserre/aarslan/NIPS_alter/EXPfiles/';
myStemDir = [subjName '_grouped/'];
myRandStemDir = '007_grouped_rand/';
myAreaDir = area; %'G_oc-temp_med-Parahip/';
myResDir = [myEXPDir myStemDir myAreaDir '/'];
myRandResDir = myResDir; %[myEXPDir myRandStemDir myAreaDir];

myFiles = dir([myResDir '*pow.mat']);
mrRandFiles = dir([myRandResDir '*pow.mat']);

enableRandTest = 0;
regularT = 0;

%check out what's the size
load([myResDir myFiles(1).name])
[performance ~] = mainClassificationPrepResults(EXP);
thisLong = numel([performance.performanceCurves.classifCurve{:}]);
%%%

myPerfGraph_pow = zeros(numel(myFiles), thisLong);
myPerfPop = [];
myPerfOrd_pow = [];

myFreqs = zeros(numel(myFiles),1);

parfor fi=1:numel(myFiles)
    thisOne = load([myResDir myFiles(fi).name]);
    myFiles(fi).name
    name = regexp( myFiles(fi).name, '_', 'split');
    myFreqs(fi) = str2double(name{3});
    [performance decInfo] = mainClassificationPrepResults(thisOne.EXP);
    myPerfGraph_pow(fi,:) = [performance.performanceCurves.classifCurve{:}];
    %myWeiGraph(fi,:,:) = [performance.performanceWeights.weights_mean];
    
    %myPerfPop(fi,:,:) = performance.performanceCurves.classifPop{:};
    
%     if enableRandTest && (numel(mrRandFiles) ~= 0) && regularT
%         load([myRandResDir myFiles(fi).name])
%         [performance decInfo] = mainClassificationPrepResults(EXP);
%         myPerfPop_rand(fi,:,:) = performance.performanceCurves.classifPop{:};
%     end
    
end
locFlds = fields(EXP.chanlocs);
elecCnt = numel(EXP.chanlocs.(locFlds{1}));

[new ord] = sort(myFreqs);
myPerfOrd_pow = myPerfGraph_pow(ord,:);

clims = [45 75];
hand.fig = figure;

or1 = [0.2 0.2];
msize = [0.75 0.70];
h_main = axes('position', [or1 msize]);
imagesc(EXP.timesOfClassif, new, myPerfOrd_pow, clims)
text(EXP.timesOfClassif(4), new(4), [num2str(elecCnt) ' elect'], 'color', [1 1 1 ], 'fontsize', 20)
set(h_main, 'YDir', 'normal')
set(h_main, 'xtick', [])
set(h_main, 'ytick', [])
%xlabel('time in ms')
%ylabel('frequency in Hz')
colorbar
sizes = get(h_main, 'position');
strtx = sizes(1); strty = sizes(2); h = sizes(4); w = sizes(3);

powSum = mean(myPerfOrd_pow,2);
timeSum = mean(myPerfOrd_pow);

sideGraphSize = 0.1;
han_pow = axes;
plot(new, powSum)
xlim([new(1), new(end)])
ylim([50 70])
wpow = sideGraphSize;
set(han_pow, 'View', [90 270])
set(han_pow, 'position', [or1(1)-wpow or1(2) wpow msize(2)])
set(han_pow, 'YDir', 'reverse')
xlabel('frequency in Hz')

han_time = axes;
plot(EXP.timesOfClassif,timeSum)
htime = sideGraphSize;
set(han_time, 'position', [or1(1) or1(2)-htime sizes(3) htime])
xlim([EXP.timesOfClassif(1) EXP.timesOfClassif(end)])
ylim([50 70])
xlabel('time in ms')
hand.main = h_main;
hand.vert = han_pow;
hand.hori = han_time;

end
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
