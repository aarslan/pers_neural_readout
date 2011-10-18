function [EEGf] =  featureMaker4(EEG, varargin)
% receives an EEG set and adds another field containted .featClassif
%%which contains several sub fields to describe the data. All of these
%%statistics are produced within a moving window which moves by every
%%consecutive sample. Created for every single trial in
%%every channel.

% Inputs:
% EEG = the source EEG structure from eeglab
% binSize = the size of the bin in samples, this can be an array as ylong as you want, but use only even numbers
% domain = 'time'|'spectral'|'mean'|'bands' you need to specify which features will
% be calculated with this switch. Structures will be created accordingly.

%example use:
% EEG = featureMaker4(EEG, [2 4 6 8 10], 'time', 'spectral' )

% Features:
% mean
% range = |min-max|
% contrast = range/|max+min|
% std = sexually transmitted disease
% der1, der2, der3 = 1st, 2nd, 3rd derivatives
% raw = unprocessed data bins

display(sprintf('Processing patient %s', EEG.filename(1:3)))

% initiating the structure and the arrays
nfreqBins=20; %this was the default number given by eeglab for our 250hz, [-200 650] epochs
FFTwinSize = 32;
nTimePoints = EEG.pnts-FFTwinSize;
freqs = linspace(5, 120, nfreqBins);
cycles = [linspace(1, 8, 8) ones(1, 12)*8];
spectralDone = 0; %this tag is used to prevent repetitive spectral analysis when the time analysis is done for more than one bin size
bandsDone = 0; %this tag is used to prevent repetitive bands analysis when the time analysis is done for more than one bin size

%%this is just one trial of analysis to find out EEGLAB's default array size for the
%%spectral stuff
[ersp,itc,powbase,times,freqs,erspboot,itcboot] = ...
    newtimef(EEG.data(1,:,1), EEG.pnts, [EEG.times(1) EEG.times(end)], EEG.srate, cycles,...
    'freqs', freqs, 'timesout', nTimePoints,  ...
    'plotersp', 'off', 'plotitc', 'off', 'plotphasesign', 'off', 'plotphaseonly', 'off', 'newfig', 'off' );
nTimePoints = length(times);

zeroArr = zeros(1, EEG.pnts);
zeroVect = zeros(EEG.nbchan,EEG.pnts, EEG.trials);
zeroVect2 = zeros(EEG.nbchan, nfreqBins, nTimePoints, EEG.trials);




for k=1:length(varargin{1})
    binSizes = varargin{1};
    binSize = binSizes(k);
    if ~isempty(cell2mat(strfind(varargin(2:end), 'time')))
        meanArr = zeroArr;
        rangeArr = zeroArr;
        mmean = zeroArr;
        contrastArr = zeroArr;
        stdArr = zeroArr;
        der1 = zeroArr;
        der2 = zeroArr;
        der3 = zeroArr;
        raw = zeros(binSize, EEG.pnts); %%this field contains an array for each sample, so it's a 2d array
        
        mean1s = zeroVect;
        mean2s = zeroVect;
        ranges = zeroVect;
        contrasts = zeroVect;
        stds = zeroVect;
        der1s = zeroVect;
        der2s = zeroVect;
        der3s = zeroVect;
        raws = zeros(EEG.nbchan, binSize, EEG.pnts, EEG.trials);
    end
    
    if ~isempty(cell2mat(strfind(varargin(2:end), 'spectral')))
        ersps = zeroVect2;
        itcs = zeroVect2;
    end
    
    if ~isempty(cell2mat(strfind(varargin(2:end), 'mean')))
        mean1s = zeroVect;
        meanArr = zeroArr;
    end
    
    if ~isempty(cell2mat(strfind(varargin(2:end), 'bands'))) && ~bandsDone
        %%filtering the data for 'bands' features and hilbert transform
        filtOrder = round(EEG.pnts/3-1);
        delta = reshape(eegfilt(EEG.data, EEG.srate, 1, 4, 0, filtOrder),EEG.nbchan,EEG.pnts, EEG.trials);
        theta = reshape(eegfilt(EEG.data, EEG.srate, 4, 8, 0, filtOrder), EEG.nbchan,EEG.pnts, EEG.trials);
        alpha = reshape(eegfilt(EEG.data, EEG.srate, 8, 15, 0, filtOrder), EEG.nbchan,EEG.pnts, EEG.trials);
        beta = reshape(eegfilt(EEG.data, EEG.srate, 15, 30, 0, filtOrder), EEG.nbchan,EEG.pnts, EEG.trials);
        gamma = reshape(eegfilt(EEG.data, EEG.srate, 30, 60, 0, filtOrder), EEG.nbchan,EEG.pnts, EEG.trials);
        higamma = reshape(eegfilt(EEG.data, EEG.srate, 60, 0, 0, filtOrder), EEG.nbchan,EEG.pnts, EEG.trials);
        display('eeg bandpass filtering is completed');
        
        deltaPow = zeroVect;
        thetaPow = zeroVect;
        alphaPow = zeroVect;
        betaPow = zeroVect;
        gammaPow = zeroVect;
        higammaPow = zeroVect;
        
        deltaPhs = zeroVect;
        thetaPhs = zeroVect;
        alphaPhs = zeroVect;
        betaPhs = zeroVect;
        gammaPhs = zeroVect;
        higammaPhs = zeroVect;
    end
    
    % Filling the .featClassif structure
   
    for chn=1: EEG.nbchan
        tic
        for  tri=1: EEG.trials
            if ~isempty(cell2mat(strfind(varargin(2:end), 'time')))
                %padding the extremes of the array with repetition
                paddedArr = [fliplr(EEG.data(chn,1:binSize/2,tri)) EEG.data(chn,:,tri) fliplr(EEG.data(chn,end-binSize/2:end,tri))];
                %mean
                meanArr = wkeep(convn(paddedArr, ones(1,binSize)/binSize,'same'), EEG.pnts,'c');
                
                i=0;
                for samp=binSize/2+1:EEG.pnts+binSize/2
                    i=i+1;
                    currBin = paddedArr(samp-binSize/2:samp+binSize/2-1); %get the relevant interval
                    rangeArr(i) = abs(min(currBin)-max(currBin));
                    mmean(i) = mean(currBin);
                    contrastArr(i) = rangeArr(i)/abs(min(currBin)+max(currBin));
                    stdArr(i) = std(currBin);
                    raw(:,i) = currBin';
                end
                %derivatives
                der1 = diff(paddedArr);
                der1 = wkeep(der1(binSize/2:end),EEG.pnts,'l') ;
                der2 = diff(paddedArr,2);
                der2 = wkeep(der2(binSize/2:end),EEG.pnts,'l') ;
                der3 = diff(paddedArr, 3);
                der3 = wkeep(der3(binSize/2:end),EEG.pnts,'l') ;
                
                mean1s(chn,:,tri) = meanArr;
                mean2s(chn,:,tri) = mmean;
                ranges(chn,:,tri) = rangeArr;
                contrasts(chn,:,tri) = contrastArr;
                stds(chn,:, tri) = stdArr;
                der1s(chn,:,tri) = der1;
                der2s(chn,:,tri) = der2;
                der3s(chn,:,tri) = der3;
                raws(chn,:,:,tri) = raw;
            end
            %tic
            if ~isempty(cell2mat(strfind(varargin(2:end), 'spectral'))) && ~spectralDone
                [ersp,itc,powbase,times,freqs,erspboot,itcboot] = ...
                    newtimef(EEG.data(chn,:,tri), EEG.pnts, [EEG.times(1) EEG.times(end)], EEG.srate, cycles,...
                    'freqs', freqs, 'timesout', nTimePoints,  ...
                    'plotersp', 'off', 'plotitc', 'off', 'plotphasesign', 'off', 'plotphaseonly', 'off', 'newfig', 'off' );
                ersps(chn,:,:,tri) = ersp;
                itcs(chn,:,:,tri) = itc;
                
            end
            
            if ~isempty(cell2mat(strfind(varargin(2:end), 'mean')))
                paddedArr = [fliplr(EEG.data(chn,1:binSize/2,tri)) EEG.data(chn,:,tri) fliplr(EEG.data(chn,end-binSize/2:end,tri))];
                meanArr = wkeep(convn(paddedArr, ones(1,binSize)/binSize,'same'), EEG.pnts,'c');
                mean1s(chn,:,tri) = meanArr;
            end
            
            %toc
            %display(sprintf('Trial........................................%d', tri))
        end
        
        if ~isempty(cell2mat(strfind(varargin(2:end), 'bands')))&& ~bandsDone
            deltaChn = hilbert(delta(chn,:,:));
            thetaChn = hilbert(theta(chn,:,:));
            alphaChn = hilbert(alpha(chn,:,:));
            betaChn = hilbert(beta(chn,:,:));
            gammaChn = hilbert(gamma(chn,:,:));
            higammaChn = hilbert(higamma(chn,:,:));
            
            deltaPow(chn,:,:) = abs(delta(chn,:,:)).^2;
            thetaPow(chn,:,:) = abs(theta(chn,:,:)).^2;
            alphaPow(chn,:,:) = abs(alpha(chn,:,:)).^2;
            betaPow(chn,:,:) = abs(beta(chn,:,:)).^2;
            gammaPow(chn,:,:) = abs(gamma(chn,:,:)).^2;
            higammaPow(chn,:,:) = abs(higamma(chn,:,:)).^2;
            
            deltaPhs(chn,:,:) = angle(delta(chn,:,:)).^2;
            thetaPhs(chn,:,:) = angle(theta(chn,:,:)).^2;
            alphaPhs(chn,:,:) = angle(alpha(chn,:,:)).^2;
            betaPhs(chn,:,:) = angle(beta(chn,:,:)).^2;
            gammaPhs(chn,:,:) = angle(gamma(chn,:,:)).^2;
            higammaPhs(chn,:,:) = angle(higamma(chn,:,:)).^2;
            bandsDone = 1;
        end
        display(sprintf('Channel %d', chn))
    end
    toc
    if ~isempty(cell2mat(strfind(varargin(2:end), 'time')))
        EEG.featClassif(k).mean1 = mean1s;
        EEG.featClassif(k).mean2 = mean2s;
        EEG.featClassif(k).range = ranges;
        EEG.featClassif(k).contrast = contrasts;
        EEG.featClassif(k).std = stds;
        EEG.featClassif(k).der1 = der1s;
        EEG.featClassif(k).der2 = der2s;
        EEG.featClassif(k).der3 = der3s;
        EEG.featClassif(k).raw = raws;
    end
    
    if ~isempty(cell2mat(strfind(varargin(2:end), 'spectral')))
        EEG.featClassifSpect.ersp = ersps;
        EEG.featClassifSpect.itc = itcs;
        EEG.featClassifSpect.freqs = freqs;
        EEG.featClassifSpect.times = times;
        spectralDone = 1;
    end
    
    if ~isempty(cell2mat(strfind(varargin(2:end), 'mean')))
        EEG.featClassif(k).mean1 = mean1s;
    end
    
    if ~isempty(cell2mat(strfind(varargin(2:end), 'bands')))
        EEG.featClassif.pow.delta = deltaPow;
        EEG.featClassif.pow.theta = thetaPow;
        EEG.featClassif.pow.alpha = alphaPow;
        EEG.featClassif.pow.beta = betaPow;
        EEG.featClassif.pow.gamma = gammaPow;
        EEG.featClassif.pow.higamma = higammaPow;
        
        EEG.featClassif.phs.delta = deltaPhs;
        EEG.featClassif.phs.theta = thetaPhs;
        EEG.featClassif.phs.alpha = alphaPhs;
        EEG.featClassif.phs.beta = betaPhs;
        EEG.featClassif.phs.gamma = gammaPhs;
        EEG.featClassif.phs.higamma = higammaPhs;
        bandsDone = 1;
    end
end
tic
EEGf=EEG;
EEGf.binSizes = varargin{1};
%naming and saving
tag=[];
if length(varargin{1})<5
    for i=1:length(varargin{1})
        tag = [tag int2str(varargin{1}(i)) '_'];
    end
else
    tag = ['binSize_' int2str(varargin{1}(1)) 'to' int2str(varargin{1}(end)) '_'];
end

if length(varargin)~=2
    dom = [upper(varargin{2}(1)), upper(varargin{3}(1))];
else
    dom = upper(varargin{2}(1));
end
subj   = EEG.filename(1:findstr(EEG.filename,'merged.set')-1);
save([subj sprintf('_feats_%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s', tag, dom) '.mat' ], 'EEGf', '-v7.3')
toc
end
