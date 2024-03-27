function [Directy] = ComputeDirectyGradient(anArray, kernel_size)
% Compute an orientation estimate (also known as directionality) for input
% matrix of grayscale image intensities, 
% using Scharr gradient estimation kernels of sizes 5 [recommended], 3 or 2 
% that are optimized for low orientational error, see Tabelle B.11 in
%
%   Hanno Scharr (2010): Optimale Operatoren in der digitalen
%   Bildverarbeitung, PhD Thesis at University Heidelberg, Germany,
%   DOI=10.11588/heidok.00000962, http://www.ub.uni-heidelberg.de/archiv/962 
%
% The gradient gives the direction of steepest increase (here: of image
% intensity), therefore the direction perpendicular to it is returned as
% orientation of objeccts in the image.
%
% Karl Hoffmann, Max Planck Institute of Molecular Cell Biology and
% Genetics, Dresden, Germany
% last edit: 2022_08_25
%
%
% CHANGELOG
% 2022_08_25 clean-up
% 2021_03_15 computation of u and v simplified


if min(anArray(:)) < 0
    warning('Negative entries in input array of ComputeDirecty detected, but gradient estimation not impeded.')
end
% do not enforce a normalization here, leave it to higher code levels

image_size = size(anArray);
   
if (kernel_size == 5)
    %{
    % Compare Directionality plugin in Fiji 
    filtering_kernel_x = [-2, -1, 0, 1, 2; -3, -2, 0, 2, 3; -4, -3, 0, 3, 4]/46;
    filtering_kernel_x = [filtering_kernel_x ; flipud(filtering_kernel_x(1:2,:))];
    filtering_kernel_y = rot90(filtering_kernel_x);
    %}
    % Use "5x5-opt" from [Scharr], see Tabelle B.11, and Abbildung 6.2. 
    % Note error of whopping 6.7e-05 (!)
    derivative_kernel = [21.27, 85.46,      0, -85.46, -21.27]/256;
    smoothing_kernel  = [ 5.91, 61.77, 120.64,  61.77,   5.91]/256;
    filtering_kernel_x= derivative_kernel .* smoothing_kernel' ;
    filtering_kernel_y = rot90(filtering_kernel_x);
elseif (kernel_size == 2)
    % The simplest kernel possible. 
    % See [Scharr], Tabelle B.11 and Abbildung 6.2 for error of 1.8e-02
    % We here directly combine derivative ("Ableitung") and smoothing ("Glättung")
    % NOTE that in the original work, the normalisation factor /2 is missing in the smoothing ("Glättung")
    filtering_kernel_x = [ 1, -1; 1, -1] /2;
    filtering_kernel_y = rot90(filtering_kernel_x);
elseif (kernel_size == 3)
    % See [Scharr], Tabelle B.11 and Abbildung 6.2 for error of 2.2e-03
    % both for "3x3-int" (which is optimised for low angle error in fixed point arithmetics)
    % and for "3x3-opt"  (which is optimised for low angle error in floating point arithmetics)
    % We here directly combine derivative ("Ableitung") and smoothing ("Glättung")
    % 3x3-int:  
    %filtering_kernel_x = [47, 0, -47; 162, 0, -162; 47, 0, -47]/512;
    % 3x3-opt
    filtering_kernel_x = [46.84, 0, -46.84; 162.32, 0, -162.32; 46.84, 0, -46.84]/512;
    filtering_kernel_y = rot90(filtering_kernel_x);
else
   error(['Please choose a kernel size 2, 3 or 5. Chosen kernel size was ' mat2str(kernel_size) '.']); 
end
    


%{
% the commented version yields derivative approximations even at the very
% boundary pixels where the filter is not fully within the input image
convolved_image_x = conv2(anArray, filtering_kernel_x);
convolved_image_y = conv2(anArray, filtering_kernel_y);

radius_kernel = floor((size(filtering_kernel_x,1) - 1)/2);

%prune
convolved_image_x = convolved_image_x(radius_kernel + 1 : image_size(1) + radius_kernel, ...
    radius_kernel + 1 : image_size(2) + radius_kernel);
convolved_image_y = convolved_image_y(radius_kernel + 1 : image_size(1) + radius_kernel, ...
    radius_kernel + 1 : image_size(2) + radius_kernel);
%}

convolved_image_x_oversized = conv2(anArray, filtering_kernel_x);
convolved_image_y_oversized = conv2(anArray, filtering_kernel_y);

convolved_image_x = NaN(size(anArray));
convolved_image_y = NaN(size(anArray));

%leave values to NaN where the kernel was not fully within the input image
kernel_oversize_tl = ceil ((size(filtering_kernel_x,1) - 1)/2); %how many pixels the kernel has to the top or left of its center ("center" = pixel that the filtering result is assigned to)
kernel_oversize_br = floor((size(filtering_kernel_x,1) - 1)/2);

convolved_image_x(kernel_oversize_tl+1:end-kernel_oversize_br, kernel_oversize_tl+1:end-kernel_oversize_br) = ...
    convolved_image_x_oversized(kernel_oversize_tl+kernel_oversize_br+1 : end-kernel_oversize_tl-kernel_oversize_br, kernel_oversize_tl+kernel_oversize_br+1 : end-kernel_oversize_tl-kernel_oversize_br);
convolved_image_y(kernel_oversize_tl+1:end-kernel_oversize_br, kernel_oversize_tl+1:end-kernel_oversize_br) = ...
    convolved_image_y_oversized(kernel_oversize_tl+kernel_oversize_br+1 : end-kernel_oversize_tl-kernel_oversize_br, kernel_oversize_tl+kernel_oversize_br+1 : end-kernel_oversize_tl-kernel_oversize_br);


% L2 norm
magnitudesL2 = sqrt(convolved_image_x.^2 + convolved_image_y.^2);

% normalized components of the identified direction
%{
angles = atan2(convolved_image_y , convolved_image_x);
orthogonal_angles = angles + pi/2;
u = cos(orthogonal_angles);
v = sin(orthogonal_angles);
%}
% more direct computation of normalized components u and v
u = -convolved_image_y ./magnitudesL2 ; %normalization and rotation by pi/2 in one step
v =  convolved_image_x ./magnitudesL2 ; % therefore "-y" assigned to u, and "+x" assigned to v

Directy(:,:,1) = u;
Directy(:,:,2) = v;
Directy(:,:,3) = magnitudesL2;

end

