
%%This script takes the Subject structures  for different sets (loads from the given path) and
% 1) imports them to EEGLAB
% 2) extract epochs
% 3) concatanates the 2 files

%Subject parameters--------------------------------------------------------
subjNumber = '005';
subjInd = 4; %Maxime's stupid subject indices
chnSkip = 1; %Number of channels to skip at the bottom of the array
filePath = './'; %'/gpfs/data/tserre/data/neural_readout/human/';
setCount = 2;
setNames = ['a' 'b']; %this is used only when LOADING the data
eeglab
newSampRate = 250;
epochLenght = [-0.20 0.650];

for i=1:setCount;
    
    load ([filePath 'Subj_alter_' subjNumber setNames(i) '.mat'])
    
    %parameters----------------------------------------------------------------
    nbchan=size(Subj.ieeg.dat, 1);
    srate=Subj.ieeg.sr;
    eegArray = Subj.ieeg.dat(1:end-chnSkip,:);
    name = Subj.beha.PARAMS.userID;
    name = name(1:3);
    session = Subj.beha.PARAMS.set;
    
    %mat2set-------------------------------------------------------------------
    
    CreateEventTableforEEGLAB([name(1:3) session '_events.txt'], Subj);
    display('event table generated')
    
    EEG = pop_importdata( 'dataformat', 'array', 'data', eegArray, 'setname', [name(1:3) session], 'srate', srate);
    
    EEG = pop_importevent(EEG, 'event', [name session '_events.txt'], 'fields',{'type','latency','trial','distance', 'imageName', 'rt', 'isi', 'correct', 'mask', 'session'}, 'skipline',1, 'timeunit',1);
   
    %Maxime's change to match with monkey imageName as to be cell of string;
    im={EEG.event.imageName};
    EEG = pop_editeventfield(EEG, 'typeinfo', [], 'latencyinfo', [], 'imageName', im);
    %----------
    EEG = getSubjChanLocs(EEG, subjInd);
    EEG = pop_resample( EEG, newSampRate);
    EEG = pop_epoch( EEG, {  }, epochLenght , 'epochinfo', 'yes');
    
    if length(EEG.event) ~= length(EEG.epoch)
        error('the number of epochs and events do not match, check the epoch length')
    end
    
    [ALLEEG EEG i] = eeg_store(ALLEEG, EEG, i);
    
end
EEGmerge = pop_mergeset( ALLEEG, [1 2], 0);
 
EEGmerge = pop_eegfilt( EEGmerge, 55, 65, [], [1]);
[ALLEEG EEGmerge] = eeg_store(ALLEEG, EEGmerge);
EEGmerge = getSubjChanLocs(EEG, subjInd);
EEG = EEGmerge;
EEG = pop_saveset( EEGmerge,  'filename', [name 'merged.set'], 'filepath', filePath);