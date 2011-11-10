% findStruct() - find for a structure in the type EEG.event(1:n).correct
% Usage:  [indices] = findStruct(S,field,condition)
%
% Requierd input:     meaning (size or options) {default}
%
%   S               = structure 
%   field           = field on which you want to perform the finf
%   condition       = condition of selection
%
% Optional input:     meaning (size or options) {default}

% Outputs:            meaning (size or options)
%
%   indices         = indeces corresponding to the condition
%
% Autors: Maxime, Serre Lab, 2010
% Sebastien, Nov 9 2010 : rewrite with dynamic fields 
% Maxime, Dec 23 2010: adapt to char event
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [indices] = findStruct(S,field,condition)

if ~ischar(field); field = char(field); end

if ischar([S.(field)]);array=cellstr(char(S.(field)));else array=[S.(field)];end
    
if ischar(condition)
    indices = find(strcmp(array,condition));
else indices = find(array == condition);
end
    
    


% % Get field value
% %----------------
% if isstr(getfield(S,{1},field))
%     for i= 1: length(S)
%         fieldValue{i} = getfield(S,{i},field);
%     end
% else
%     for i= 1: length(S)
%         fieldValue{i} = num2str(getfield(S,{i},field));
%     end
% end
% % find
% %-----
%     if ~isstr(condition)
%         condition = num2str(condition);
%     end
% indices = find(strcmp(fieldValue,condition));
% end