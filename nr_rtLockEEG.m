function newEEG = nr_rtLockEEG(EEG, x, l)

RT = [EEG.event(:).rt];
%figure;
%hist(RT, 100)

tlimit = 0.650;
RT1 = RT < tlimit+l;
a = RT;
a(~RT1) = 0;
RT2 = 0 < a - (x+l);
newRT = a(RT2);

%figure;
%hist(newRT,50); % hist of selected trials
triInx = RT2;

featFields = fields(EEG.featClassif);

newEEG = EEG;
for fi=1:numel(featFields)
    newEEG.featClassif.(featFields{fi}) = EEG.featClassif.(featFields{fi})(:,:,triInx);
end

newEEG.trials = numel(find(triInx));
newEEG.event = newEEG.event(triInx);
newEEG.epoch = newEEG.epoch(triInx);
display(['selected ' num2str(newEEG.trials) ' trials'])

lenInMS = (x+l)*1000; % length in ms
[lel lol] = hist(diff(newEEG.times),100);
[ ~, most ]= max(lel);
samInt = lol(most); % find sampling interval
lenInSmp = floor(lenInMS/samInt); %length in samples

nudgeAmnt = (newRT*1000)-lenInMS;
nudgeAmntInSmp = round(nudgeAmnt/samInt);
%figure;
%hist(nudgeAmntInSmp,50);

newTimes = (-x)*1000:samInt:-l*1000;
for fi=1:numel(featFields)
    datSize = size(newEEG.featClassif.(featFields{fi}));
    newData = zeros(datSize(1), numel(newTimes), datSize(3));
    for tri=1:newEEG.trials
        [~, endInx] = min(abs(EEG.times-(newRT(tri)-l)*1000)); %end point 
        newData(:,:,tri) = newEEG.featClassif.(featFields{fi})(:, endInx-numel(newTimes):endInx-1, tri);
    end
    newEEG.featClassif.(featFields{fi}) = newData;
end

% for fi=1:numel(featFields)
%     datSize = size(newEEG.featClassif.(featFields{fi}));
%     newData = zeros(datSize + [0 max(nudgeAmntInSmp) 0]);
%     for tri=1:numel(nudgeAmntInSmp)
%         newData(:, nudgeAmntInSmp(tri)+1:nudgeAmntInSmp(tri)+datSize(2), tri) = newEEG.featClassif.(featFields{fi})(:, :, tri);
%     end
%     newTimes = (-x)*1000:samInt:-l*1000;
%     newEEG.featClassif.(featFields{fi}) = newData(:,end+1-numel(newTimes):end,:);
% end
%%%%%%%%%
newEEG.newTimes = newTimes;

end
