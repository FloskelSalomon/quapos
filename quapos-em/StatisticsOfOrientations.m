% StatisticsOfOrientations
% Compute statistics for orientation angles respecting their periodicity of pi radians (or 180 degrees)
% by the Q tensor, which is equivalent to circular statistics for a doubled-angle representation.
%
% Dr. Karl Hoffmann, Max Planck Institute of Molecular Cell Biology and Genetics, Dresden, Germany,
% now at FI Freiberg Institut für Energie- und Klimaoekonomie GmbH, Freiberg, Germany
% file generated 2024-01-26
% last edit 2024-03-06

% read orientation angles, p.e. from .xlsx or .csv
FileWithOrientationAnglesInRadians = readmatrix('<your/filename/here.csv') % details to be filled

disp('the matrix has the following number of rows')
rows = size(FileWithOrientationAnglesInRadians,1)

%flatten to one-dimensional array by selecting the correct column of the matrix (here: column 2)
arrayOfOrientationAnglesInRadians = FileWithOrientationAnglesInRadians(1:rows,2);

% check data (disable for large data sets)
disp('arrayOfOrientationAnglesInRadians')
arrayOfOrientationAnglesInRadians

%optionally give weights
% weights = xlsxread( ) % details to be filled
% otherwise all orientation angles get the same weight
weights = ones(size(arrayOfOrientationAnglesInRadians));

[ dominantOrientation , coherency ] = ...
            ComputeQTensor( arrayOfOrientationAnglesInRadians, ...
                            weights, ...
                            1 );
assert(numel(dominantOrientation) == 2 , 'result should be a single two-component vector!')

% give dominant orientation also as angle (in radians)
dominantOrientationAngleInRadians = atan2(dominantOrientation(2), dominantOrientation(1));
dominantOrientationAngleInRadians

disp('average orientation angle (in degree)')
dominantOrientationAngleInDegree = mod(rad2deg(dominantOrientationAngleInRadians), 180)

disp('coherency')
coherency
