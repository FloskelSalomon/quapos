function [dominantOrientationField, coherencyField] = ComputeLocalAlignment(theta, magn, boxRadius_LocalAlignment, powerOfMagn_LocalAlignment)
% Compute local degree of alignment (also known as coherency) and locally 
% dominant orientations for orientation field given by magnitude magn and
% angle theta. Calculation is within square subpatches of side length
% 1+2*boxRadius, with result assigned to the center
%
% Iteration over the subpatches is with parallelisation (parfor)
%
% Karl Hoffmann, Max Planck Institute of Molecular Cell Biology and
% Genetics, Dresden, Germany
% last edit: 2024-03-06

%% preparations
% sanity checks
assert( all(size(theta) == size(magn)) , 'Input array of ComputeLocalAlignment must be of equal size' )
assert(boxRadius_LocalAlignment > 0 && floor(boxRadius_LocalAlignment)==boxRadius_LocalAlignment, 'Input boxRadius to ComputeLocalAlignment must be positive integer')

% set optional arguments if not provided
if (~exist('PowerOfMagn', 'var'))
	powerOfMagn_LocalAlignment = 1;
end

% prepare output arrays 
i_max = size(theta,1);
j_max = size(theta,2);
dominantOrientationField   = zeros( i_max, j_max, 2 );
dominantOrientationField_1 = zeros( i_max, j_max );
dominantOrientationField_2 = zeros( i_max, j_max );
coherencyField             = zeros( i_max, j_max );

%% calculation
parfor ii=boxRadius_LocalAlignment+1:i_max-boxRadius_LocalAlignment
    dominantOrientationField_1_Row = zeros(1,j_max);
    dominantOrientationField_2_Row = zeros(1,j_max);
    coherencyField_Row             = zeros(1,j_max);
    for jj=boxRadius_LocalAlignment+1:j_max-boxRadius_LocalAlignment
        %[ dominantOrientationField(ii,jj,:) , coherencyField(ii,jj) ] = ...
        [ dominantOrientation , coherency ] = ...
            ComputeQTensor( theta(ii-boxRadius_LocalAlignment:ii+boxRadius_LocalAlignment, jj-boxRadius_LocalAlignment:jj+boxRadius_LocalAlignment), ...
                            magn (ii-boxRadius_LocalAlignment:ii+boxRadius_LocalAlignment, jj-boxRadius_LocalAlignment:jj+boxRadius_LocalAlignment), ...
                            powerOfMagn_LocalAlignment ) ;
        dominantOrientationField_1_Row(jj) = dominantOrientation(1);
        dominantOrientationField_2_Row(jj) = dominantOrientation(2);
        coherencyField_Row(jj)             = coherency;
    end
    dominantOrientationField_1(ii,:) = dominantOrientationField_1_Row;
    dominantOrientationField_2(ii,:) = dominantOrientationField_2_Row;
    coherencyField(ii,:)             = coherencyField_Row;
end

dominantOrientationField(:,:,1) = dominantOrientationField_1;
dominantOrientationField(:,:,2) = dominantOrientationField_2;
