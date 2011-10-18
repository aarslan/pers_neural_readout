function CreateEventTableforEEGLAB(filename, Subj);
%% a function which takes output file name and Subject structure as inputs and outputs a space-delimited txt file. The text file contains all the events in the experiments with the following column names


fid = fopen(filename, 'w');
fprintf(fid, 'type	latency	trial	distance	imageName	rt	isi	correct	mask	session\n');

% type: (1=animal/ 0=nonanimal)
% trialType: trial number
% latency: in sec
% distance: (1=head/ 2=close body/3= medium/4= far body)
% image: name of the image used
% rt: in seconds
% isi: in seconds
% mask (1=masked/ 0=nonmasked)
% session: monkey = day of recording; human=session(A/B)

ind = find(Subj.ieeg.trig.stim);
nStim = length(ind);

if nStim ~= length(Subj.beha.resp)
    error('mismatch in num of stim presentations');
end
    
sr = Subj.ieeg.sr;

trialType   = Subj.beha.PARAMS.label;
latency     = ind/sr;
trial       = 1:1:nStim; 
distance	= Subj.beha.PARAMS.dist;
imageName   = Subj.beha.PARAMS.stim;
rt          = Subj.beha.DATA(:,1);
isi         = Subj.beha.PARAMS.isi;
correct     = Subj.beha.resp;
mask        = Subj.beha.mask;
session     = Subj.beha.PARAMS.set;

for i=1:nStim
    theImageName = char(imageName(i));
fprintf(fid, '%d\t%.3f\t%d\t%d\t%s\t%.3f\t%d\t%d\t%d\t%c\n', ...
    trialType(i), latency(i), trial(i),  distance(i), theImageName(6:end), rt(i), isi(1), correct(i), mask(i), session);
end

fclose(fid);