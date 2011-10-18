if 1

%this script plots ERP responses for each of the electrodes 
%launches eeglab and takes a set file.

%parametres
dir      = '/gpfs/data/tserre/data/neural_readout/monkey_human';%path
dataset     = '006merged.set';%name of the dataset
humanTag = 1; % put 1 for human, 0 for monkey
comp ={'anim', 'nanim'}; 


nfreqBins=20; 
freqs = linspace(5, 120, nfreqBins);
cycles = [linspace(1, 8, 8) ones(1, 12)*8];

[path subj] = fileparts(dataset);
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;


%Init 
EEG    = pop_loadset( 'filename', dataset , 'filepath', dir);
ERPs = zeros(length(EEG.times), EEG.nbchan);

%load, select correct, baseline extract
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 1 );
EEG   = pop_selectevent( EEG, 'correct', 1);
EEG         = pop_rmbase(EEG, [EEG.xmin find(EEG.times==0)]);

if humanTag
    fid = fopen( [dataset(1:3) '_locs.txt']);%name of the dataset
    locations = textscan(fid, '%s %s %s %s %s %s', 'delimiter','	');
    fclose(fid);
    temp = cell2struct(locations(1), {'hemis'},1 ); locs.hemis = temp.hemis;
    temp = cell2struct(locations(2), {'lobe'},1 );locs.lobe = temp.lobe;
    temp = cell2struct(locations(3), {'gyrus'},1 );locs.gyrus = temp.gyrus;
    temp = cell2struct(locations(5), {'ba'},1 );locs.ba = temp.ba;
    clear temp fid;

end

EEG1   = pop_selectevent( EEG, 'type' , 1);
EEG1         = pop_rmbase(EEG1, [EEG1.xmin find(EEG1.times==0)]);
[ALLEEG, EEG1, CURRENTSET] = eeg_store( ALLEEG, EEG1, 2 );

EEG2   = pop_selectevent( EEG, 'type' , 0);
EEG2         = pop_rmbase(EEG2, [EEG2.xmin find(EEG2.times==0)]);
[ALLEEG, EEG2, CURRENTSET] = eeg_store( ALLEEG, EEG2, 3 );

end


for i=1:EEG.nbchan
    figure
    
    if humanTag
        subplot(4, 1, 1)
        title(sprintf('elec %d', i))
        text(00, 0.3, sprintf('%s \n%s, %s \n%s ', locs.hemis{i}, locs.lobe{i}, locs.gyrus{i}, locs.ba{i}), 'FontSize', 15 )
        axis off;
    else
        subplot(4, 1, 1)
        title(sprintf('elec %d', i))
        text(00, 0.3, sprintf('electrode %s', EEG.chanlocs(i).labels), 'FontSize', 15 )
        axis off;
    end;
    
    subplot(4, 1, 2)
    ERP = ERPdrawer(EEG, i, 'r');
    
    
    subplot(4, 1, 3)
    ERP = ERPdrawer(EEG1, i, 'b');
    
    hold on
    legend(comp{1}, 'location','NorthWest')
    
    ERP = ERPdrawer(EEG2, i, 'g');
    hold off
    
    subplot(4, 1, 4)
    [ersp,itc,powbase,times,freqs,erspboot,itcboot] = ...
        newtimef(EEG.data(i,:,:), EEG.pnts, [EEG.times(1) EEG.times(end)], EEG.srate, cycles,...
        'freqs', freqs,  ...
        'plotitc', 'off', 'plotphasesign', 'off', 'plotphaseonly', 'off');
    
    set(gcf,'PaperUnits', 'centimeters');
    set(gcf,'PaperOrientation', 'portrait');
    set(gcf,'PaperPosition',[0 0 22 30]);
    set(gcf,'PaperType','usletter');
    
    %save
    figName = [EEG.filename(1:3) '_electrode_' int2str(i)];

    print(gcf, '-dpdf', ['ERPfigures/' figName '.pdf']);
    delete(gcf);
end
