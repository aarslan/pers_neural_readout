setName = EEG.filename(1:3);
myDir = ['./sets/tfSets/' setName '/'];
if ~exist(myDir)
    mkdir(myDir)
end

%% get the tf
totaltf = [];
for sec=1:50:size(EEG.data, 3)
    sec
    if sec+50 < size(EEG.data, 3)
    stopPnt = sec+49;
    else
      stopPnt = size(EEG.data, 3);
    end
[temptf, freqs, times] = timefreq(EEG.data(:,:,sec:stopPnt), EEG.srate, ... 
    'wletmethod', 'dftfilt3', 'padratio', 4, 'cycles', [0 2], 'freqs', [4 120], ...
    'nfreqs', 60 ,'tlimits', [EEG.times(1) EEG.times(end)], 'timesout', linspace(-130,580, 180)); %'ntimesout', 180,
totaltf = cat(4, totaltf, temptf);
end

%% disguise as a new EEG
EEGf = EEG;
EEGf.pnts = size(totaltf,3);
EEGf.times = times;
EEGf.xmin = times(1);
EEGf.xmax = times(end);
EEGf.data = false(size(totaltf,1), size(totaltf,3), size(totaltf,4));

%% prep  phs n pow
% figure 
% 
 regularPhs = double(angle(totaltf));
% subplot(2,2,1); imagesc(regularPhs); colorbar
% 
 unwrpPhs = unwrap(regularPhs,[], 3);
% subplot(2,2,2); imagesc(unwrpPhs); colorbar
% 
 diffPhs = diff(unwrpPhs,[],3);
% subplot(2,2,3); imagesc(diffPhs); colorbar  
powz = abs(totaltf);

%% create features
for fr=1:size(totaltf,2)
featClassif.(['hz_pow_' mat2str(round(freqs(fr))) ]) = squeeze(powz(:,fr,:,:));
featClassif.(['hz_phs_' mat2str(round(freqs(fr))) ]) = squeeze(diffPhs(:,fr,:,:));
EEGf.featClassif = featClassif;
clear featClassif
fileName = [EEG.filename(1:3) '_hz_' mat2str(round(freqs(fr))) '.mat'];
save( [myDir fileName], 'EEGf','-v7.3');
end

