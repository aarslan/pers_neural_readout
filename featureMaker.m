
% featureMaker () - create features for classification
%
% Usage: EEG = featureMaker(EEG, binSizes, feats, save)
%
% Requierd input:     meaning (size or options) {default}
%
%   EEG             = EEG structures
%
% Optional input:     meaning (size or options) {default}
%   save             = save the EEG in a .mat with the features {'off'}
%   binSizes         = vector of sizes of bin in nb of point {[2]}
%   feats            = cell array of string features to extract ('time','timef','mean','bands')
%                         time features (computed by bin of BinSizes)
%                               - mean
%                               - range = |min-max|
%                               - contrast = range/|max+min|
%                               - std = sexually transmitted disease
%                               - der1, der2, der3 = 1st, 2nd, 3rd
%                               derivatives
%                               - raw = unprocessed data bins
%
% Outputs:            meaning (size or options)
%
%   EEG             = EEG structures with EEG.featClassif
%
% Autors: Ali, Serre Lab, 2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [EEGf] =  featureMaker(EEG, binSizes, feats, saveFlag)

display(sprintf('Processing patient %s', EEG.filename(1:3)))

%--------------%
% Initializing %
%--------------%
    bandsDone = 0; %this tag is used to prevent repetitive bands analysis when the time analysis is done for more than one bin size

% Newtimef parameter
%-------------------
if ismember('timef',feats)
    nfreqBins   = 20; %this was the default number given by eeglab for our 250hz, [-200 650] epochs
    FFTwinSize  = 32;
    nTimePoints = EEG.pnts-FFTwinSize;
    freqs       = linspace(5, 120, nfreqBins);
    cycles      = [linspace(1, 8, 8) ones(1, 12)*8];
    timefDone = 0; %this tag is used to prevent repetitive timef analysis when the time analysis is done for more than one bin size
    
    %this is just one trial of analysis to find out EEGLAB's default array size for the timef stuff
    [ersp,itc,powbase,times,freqs,erspboot,itcboot] = ...
        newtimef(EEG.data(1,:,1), EEG.pnts, [EEG.times(1) EEG.times(end)], EEG.srate, cycles,...
        'freqs', freqs, 'timesout', nTimePoints,  ...
        'plotersp', 'off', 'plotitc', 'off', 'plotphasesign', 'off', 'plotphaseonly', 'off', 'newfig', 'off' );
    nTimePoints = length(times);
    zeroVect2   = zeros(EEG.nbchan, nfreqBins, nTimePoints, EEG.trials);
end
zeroVect    = zeros(EEG.nbchan,EEG.pnts, EEG.trials);
zeroArr     = zeros(1, EEG.pnts);



for k=1:length(binSizes)
    binSize      = binSizes(k);
    % Times features
    %---------------
    if ismember('time',feats)
        meanArr     = zeroArr;
        rangeArr    = zeroArr;
        mmean       = zeroArr;
        contrastArr = zeroArr;
        stdArr      = zeroArr;
        der1        = zeroArr;
        der2        = zeroArr;
        der3        = zeroArr;
        raw         = zeros(binSize, EEG.pnts); %%this field contains an array for each sample, so it's a 2d array
        
        mean1s      = zeroVect;
        mean2s      = zeroVect;
        ranges      = zeroVect;
        contrasts   = zeroVect;
        stds        = zeroVect;
        der1s       = zeroVect;
        der2s       = zeroVect;
        der3s       = zeroVect;
        raws        = zeros(EEG.nbchan, binSize, EEG.pnts, EEG.trials);
    end
    
    % Timef features
    %---------------
    if ismember('timef',feats)
        ersps = zeroVect2;
        itcs  = zeroVect2;
    end
    
    % mean features
    %--------------
    
    if ismember('mean',feats)
        mean1s  = zeroVect;
        meanArr = zeroArr;
    end
    
    % bands features
    %---------------
    if ismember('bands',feats) & ~bandsDone
        %%filtering the data for 'bands' features and hilbert transform
        filtOrder = round(EEG.pnts/3-1);
        delta     = reshape(eegfilt(EEG.data, EEG.srate, 3, 5, 0, filtOrder),EEG.nbchan,EEG.pnts, EEG.trials);
        theta     = reshape(eegfilt(EEG.data, EEG.srate, 5, 8, 0, filtOrder), EEG.nbchan,EEG.pnts, EEG.trials);
        alpha     = reshape(eegfilt(EEG.data, EEG.srate, 8, 15, 0, filtOrder), EEG.nbchan,EEG.pnts, EEG.trials);
        beta      = reshape(eegfilt(EEG.data, EEG.srate, 15, 30, 0, filtOrder), EEG.nbchan,EEG.pnts, EEG.trials);
        gamma      = reshape(eegfilt(EEG.data, EEG.srate, 30, 60, 0, filtOrder), EEG.nbchan,EEG.pnts, EEG.trials);
        higamma   = reshape(eegfilt(EEG.data, EEG.srate, 60, 0, 0, filtOrder), EEG.nbchan,EEG.pnts, EEG.trials);
        display('eeg bandpass filtering is completed');
        
        deltaPow    = zeroVect;
        thetaPow    = zeroVect;
        alphaPow    = zeroVect;
        betaPow     = zeroVect;
        gammaPow    = zeroVect;
        higammaPow  = zeroVect;
        
        deltaPhs    = zeroVect;
        thetaPhs    = zeroVect;
        alphaPhs    = zeroVect;
        betaPhs     = zeroVect;
        gammaPhs    = zeroVect;
        higammaPhs  = zeroVect;
    end
    
    %------------------%
    % Extract Features %
    %------------------%
    
    for chn=1: EEG.nbchan
        tic
        for  tri=1: EEG.trials
            % Times features
            %---------------
            if ismember('time',feats)
                % Padding the extremes of the array with repetition
                paddedArr = [fliplr(EEG.data(chn,1:binSize/2,tri)) EEG.data(chn,:,tri) fliplr(EEG.data(chn,end-binSize/2:end,tri))];
                %mean
                meanArr   = wkeep(convn(paddedArr, ones(1,binSize)/binSize,'same'), EEG.pnts,'c');
                
                i=0;
                for samp            = binSize/2+1:EEG.pnts+binSize/2
                    i=i+1;
                    currBin         = paddedArr(samp-binSize/2:samp+binSize/2-1); %get the relevant interval
                    rangeArr(i)     = abs(min(currBin)-max(currBin));
                    mmean(i)        = mean(currBin);
                    contrastArr(i)  = rangeArr(i)/abs(min(currBin)+max(currBin));
                    stdArr(i)       = std(currBin);
                    raw(:,i)        = currBin';
                end
                %derivatives
                der1 = diff(paddedArr);
                der1 = wkeep(der1(binSize/2:end),EEG.pnts,'l') ;
                der2 = diff(paddedArr,2);
                der2 = wkeep(der2(binSize/2:end),EEG.pnts,'l') ;
                der3 = diff(paddedArr, 3);
                der3 = wkeep(der3(binSize/2:end),EEG.pnts,'l') ;
                
                mean1s(chn,:,tri)       = meanArr;
                mean2s(chn,:,tri)       = mmean;
                ranges(chn,:,tri)       = rangeArr;
                contrasts(chn,:,tri)    = contrastArr;
                stds(chn,:, tri)        = stdArr;
                der1s(chn,:,tri)        = der1;
                der2s(chn,:,tri)        = der2;
                der3s(chn,:,tri)        = der3;
                raws(chn,:,:,tri)       = raw;
            end
            
            % Timef features
            %---------------
            if ismember('timef',feats) & ~timefDone
                [ersp,itc,powbase,times,freqs,erspboot,itcboot] = ...
                    newtimef(EEG.data(chn,:,tri), EEG.pnts, [EEG.times(1) EEG.times(end)], EEG.srate, cycles,...
                    'freqs', freqs, 'timesout', nTimePoints,  ...
                    'plotersp', 'off', 'plotitc', 'off', 'plotphasesign', 'off', 'plotphaseonly', 'off', 'newfig', 'off' );
                ersps(chn,:,:,tri) = ersp;
                itcs(chn,:,:,tri)  = itc;
                
            end
            
            % mean features
            %---------------
            if ismember('mean',feats)
                paddedArr           = [fliplr(EEG.data(chn,1:binSize/2,tri)) EEG.data(chn,:,tri) fliplr(EEG.data(chn,end-binSize/2:end,tri))];
                meanArr             = wkeep(convn(paddedArr, ones(1,binSize)/binSize,'same'), EEG.pnts,'c');
                mean1s(chn,:,tri)   = meanArr;
            end
            
        end
        
        % Bands features
        %---------------
        if ismember('bands',feats)& ~bandsDone
            % Hilbert
            deltaChn        = squeeze(hilbert(delta(chn,:,:)));
            thetaChn        = squeeze(hilbert(theta(chn,:,:)));
            alphaChn        = squeeze(hilbert(alpha(chn,:,:)));
            betaChn         = squeeze(hilbert(beta(chn,:,:)));
            gammaChn        = squeeze(hilbert(gamma(chn,:,:)));
            higammaChn      = squeeze(hilbert(higamma(chn,:,:)));
            % Power
            deltaPow(chn,:,:)   = abs(deltaChn(:,:)).^2;
            thetaPow(chn,:,:)   = abs(thetaChn(:,:)).^2;
            alphaPow(chn,:,:)   = abs(alphaChn(:,:)).^2;
            betaPow(chn,:,:)    = abs(betaChn(:,:)).^2;
            gammaPow(chn,:,:)   = abs(gammaChn(:,:)).^2;
            higammaPow(chn,:,:) = abs(higammaChn(:,:)).^2;
            % Phase
            deltaPhs(chn,:,:)   = angle(deltaChn(:,:)).^2;
            thetaPhs(chn,:,:)   = angle(thetaChn(:,:)).^2;
            alphaPhs(chn,:,:)   = angle(alphaChn(:,:)).^2;
            betaPhs(chn,:,:)    = angle(betaChn(:,:)).^2;
            gammaPhs(chn,:,:)   = angle(gammaChn(:,:)).^2;
            higammaPhs(chn,:,:) = angle(higammaChn(:,:)).^2;
            
        end
        display(sprintf('Channel %d', chn))
    end
    toc
    
    %----------------%
    % Store Features %
    %----------------%
    
    % Times features
    %---------------
    if ismember('time',feats)
        EEG.featClassif(k).mean1    = mean1s;
        EEG.featClassif(k).mean2    = mean2s;
        EEG.featClassif(k).range    = ranges;
        EEG.featClassif(k).contrast = contrasts;
        EEG.featClassif(k).std      = stds;
        EEG.featClassif(k).der1     = der1s;
        EEG.featClassif(k).der2     = der2s;
        EEG.featClassif(k).der3     = der3s;
        EEG.featClassif(k).raw      = raws;
    end
    
    % Timef features
    %---------------
    if ismember('timef',feats)
        EEG.featClassifSpect.ersp   = ersps;
        EEG.featClassifSpect.itc    = itcs;
        EEG.featClassifSpect.freqs  = freqs;
        EEG.featClassifSpect.times  = times;
        timefDone = 1;
    end
    
    % mean features
    %---------------
    if ismember('mean',feats)
        EEG.featClassif(k).mean1    = mean1s;
    end
    
    % Bands features
    %---------------
    if ismember('bands',feats)
        % Power
        EEG.featClassif.powDelta   = deltaPow;
        EEG.featClassif.powTheta   = thetaPow;
        EEG.featClassif.powAlpha   = alphaPow;
        EEG.featClassif.powBeta    = betaPow;
        EEG.featClassif.powGamma   = gammaPow;
        EEG.featClassif.powHigamma = higammaPow;
        % Phase
        EEG.featClassif.phsDelta   = deltaPhs;
        EEG.featClassif.phsTheta   = thetaPhs;
        EEG.featClassif.phsAlpha   = alphaPhs;
        EEG.featClassif.phsBeta    = betaPhs;
        EEG.featClassif.phsGamma   = gammaPhs;
        EEG.featClassif.phsHigamma = higammaPhs;
        bandsDone = 1;
    end
end
toc
EEGf=EEG;
EEGf.binSizes = binSizes;

% Naming and saving
%------------------
if saveFlag
    tag=[];
    if length(binSizes)<5
        for i=1:length(binSizes)
            tag = [tag int2str(binSizes(i))];
        end
    else
        tag = ['binSize_' int2str(binSizes(1)) 'to' int2str(binSizes(end))];
    end
    
    subj   = EEG.filename(1:findstr(EEG.filename,'merged.set')-1);
    save([subj sprintf('_feats_%s_%s_%s_%s_%s_%s_%s', tag, feats{:}) '.mat' ], 'EEGf', '-v7.3')
    toc
end
end
