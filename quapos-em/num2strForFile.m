function filenamefriendlyString = num2strForFile(aNumber, formatSpec)
% Adapt the beviour of MATLAB's internal num2str( )
% to produce strings with the dot '.' replaced by underscore '_' to avoid
% complications in file names
%
% Karl Hoffmann, Max Planck Institute of Molecular Cell Biology and
% Genetics, Dresden, Germany
% generated: 2021-10-20
% last edit: 2024-03-06

if  exist('formatSpec', 'var')
    tempStr = num2str(aNumber,formatSpec); 
else
    tempStr = num2str(aNumber);
end
filenamefriendlyString = strrep(tempStr, '.', '_');