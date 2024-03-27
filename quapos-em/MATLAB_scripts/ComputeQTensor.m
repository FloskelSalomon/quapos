function [dominantOrientation, coherency] = ComputeQTensor(theta, magn, PowerOfMagn)
% Compute the Q-tensor of a set of orientations, given as magnitude magn
% and angle theta. The function is similar to ComputeDirecty_FourierAndQ.m
% but considers the complete input arrays rather than sliding subpatches of
% defined size. 
% The shape of the input does *not* matter, and NaN values are ignored (the
% output is NaN only when all input is NaN)
%
% Karl Hoffmann, Max Planck Institute of Molecular Cell Biology and
% Genetics, Dresden, Germany
% last edit: 2024-03-06


%% use magn to form a distribution (called  power_to_consider)
power_to_consider = magn .^ PowerOfMagn ;
% turn it into a distribution from which we compute the Q-tensor (= the second moment tensor of the distribution, except for -0.5*unitmatrix)
power_to_consider = power_to_consider / nansum(power_to_consider(:)) ;  %normalise

Q_xx =  nansum(nansum( ( cos(theta).^2  - 0.5 )      .*  power_to_consider  )) ;
Q_xy =  nansum(nansum(  cos(theta)  .* sin(theta)  .*  power_to_consider  )) ;


%% eigen decomposition of the Q-tensor [Q_xx Q_xy ; Q_xy Q_yy]  where Q_yy = - Q_xx
% eigenvalues are  \pm \sqrt{ Q_{xx}^2 + Q_{xy}^2 }
%eigVal_large =   sqrt(Q_xx ^2 + Q_xy ^2) ;
%eigVal_small = - sqrt(Q_xx ^2 + Q_xy ^2) ;
eigVec_small = [ -Q_xy ; Q_xx+sqrt(Q_xx^2+Q_xy^2) ];
if  (Q_xy==0 && Q_xx<0)  %arithmetically, this case is equal to eigVec_small == [ 0; Q_xx+sqrt(Q_xx^2)] == [0 ; 0] for Q_xx<0
    % but there might be roundoff errors
    eigVec_small = [ 1 ; 0 ];
else
    eigVec_small = eigVec_small / norm(eigVec_small) ;
end
%eigVec_large = [ eigVec_small(2) ; -eigVec_small(1) ] ;

strength = 2 * sqrt( Q_xx^2 + Q_xy^2 ); %

dominantOrientation = [ eigVec_small(2) ; -eigVec_small(1) ] ;
%dominantOrientation = eigVec_large;

coherency = strength;