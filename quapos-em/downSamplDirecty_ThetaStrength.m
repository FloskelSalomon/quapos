function [theta_downSampl, strength_downSampl] = downSamplDirecty_ThetaStrength(theta, strength, downSamplFactor, skipRows, skipCols, strength_exponent)
% Downsampling of 2D or 3D orientation field (also known as directionality
% field) in the FIRST TWO dimensions,
% based on averaging in the domain of doubled angle, which is appropriate 
% for orientations with their periodicity of pi (180 degree)
% Input and output are of the form (angle, strength).
% [For input in the form (x_component, y_component) a simple convolution
% conv2( <x or y>, 1/(downSamplFactor^2)*ones(downSamplFactor) ) suffices.]
%
% inputs:
% 1.equal-size scalar fields theta and strength give the orientation angle, 
%   and a strength of the orientation. Strength must be nonnegative. 
%   strength 0 can be used to mask out some data.
% 2.downSamplFactor - blocks of square size downSamplFactor are used
% 3.optional scalar values skipRows and skipCols (default: 0). The first
%   block starts after this number of rows (skipRows) and columns
%   (skipCols). 
% 4.optional scalar value strength_exponent controls the weight of strength 
%   for the averaging within each block. Typically either 1 or 2 [from what
%   is found in the literature].  Default is 1.
%
% Parts of the input array, that do not fill a complete block at the last
% rows and/or columns, are dropped. Hence the size of the output array is
%   floor(  (size(input arrays)-[skipRows skipCols] )/downSamplFactor )
%    
%
% Karl Hoffmann, Max Planck Institute of Molecular Cell Biology and
% Genetics, Dresden, Germany
% generated: 2021-03-26
% last edit: 2022-08-25

% CHANGELOG
% 2021-07-08 include conv2in2or3D to allow 3d arrays to be processed


%% preparations
% sanity checks
assert(all(size(theta)==size(strength)), 'Input fields "theta" and "strength" must have equal shape')
assert(int8(downSamplFactor) == downSamplFactor, 'Input "downSamplFactor" must be integer')

% set optional arguments if not provided
if (~exist('skipRows', 'var'))
	skipRows = 0;
else
    assert(int8(skipRows) == skipRows  &&  skipRows > 0, 'Input "skipRows" must be non-negative integer')
end

if (~exist('skipCols', 'var'))
	skipCols = 0;
else
    assert(int8(skipCols) == skipCols  &&  skipCols > 0, 'Input "skipCols" must be non-negative integer')
end

if (~exist('strength_exponent', 'var'))
	strength_exponent = 1;
end

%% calculations
% convert to u and v (Cartesian coordinates)
prefactor = strength.^strength_exponent;
u_doubleAngle = prefactor .* cos( 2*theta );
v_doubleAngle = prefactor .* sin( 2*theta );

u_doubleAngle_blocks = 1/(downSamplFactor^2) ...
    * conv2in2or3D(u_doubleAngle, ones(downSamplFactor), 'valid') ;
v_doubleAngle_blocks = 1/(downSamplFactor^2) ...
    * conv2in2or3D(v_doubleAngle, ones(downSamplFactor), 'valid') ;

u_doubleAngle_downSampl = u_doubleAngle_blocks(1+skipRows:downSamplFactor:end, 1+skipCols:downSamplFactor:end, :);
v_doubleAngle_downSampl = v_doubleAngle_blocks(1+skipRows:downSamplFactor:end, 1+skipCols:downSamplFactor:end, :);

% transform back into strength and angle theta
strength_downSampl = (u_doubleAngle_downSampl .^2 + v_doubleAngle_downSampl .^2) .^0.5 ;
theta_downSampl    = 0.5*atan2( v_doubleAngle_downSampl, u_doubleAngle_downSampl );
