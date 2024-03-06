function convolved = conv2in2or3D(array, kernel, optargs)
% apply 2D convolution to each of the slices of a 3D array
%
%
% Karl Hoffmann, Max Planck Institute of Molecular Cell Biology and
% Genetics, Dresden, Germany
% file generated: 2021-07-08
% last edit: 2021-07-08


% for preallocation, run the convolution once for the first slice
convolved_slice1 = conv2( array(:,:,1), kernel, optargs) ;

NumSlices = size(array,3);
% if there are more than one slice, preallocate result array and perform
% convolution in each of the remaining slices
if NumSlices > 1
    convolved = zeros( [size(convolved_slice1) NumSlices] );
    convolved(:,:,1) = convolved_slice1;
    for kk = 2:NumSlices
         convolved(:,:,kk) = conv2( array(:,:,kk), kernel, optargs) ;
    end
else
    convolved = convolved_slice1;
end
