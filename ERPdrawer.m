function ERP = ERPdrawer(EEG, elec, color)
%color is the color code appropriate for plot func

data = zeros( EEG.pnts, EEG.trials);

for i=1:EEG.trials
    data(:,i) = EEG.data(elec,:,i);
end
data = data';
ERP=mean(data);
ESM=std(data)/sqrt(EEG.trials);
plot(EEG.times, ERP,'Color', color,'LineWidth',2)%plot ERP
hold on;
fill([EEG.times,EEG.times(end:-1:1)],[ERP-ESM,ERP(end:-1:1)+ESM(end:-1:1)],color)

alpha(0.2);

y=ylim;
ymi=y(1);
yma=y(2);

xlim=([EEG.times(1) EEG.times(end)]);

%axes
plot(EEG.times,repmat(0,EEG.pnts,1),'k','LineWidth',1)
plot([0 0],[ymi yma],'k','LineWidth',1)

hold off
end
