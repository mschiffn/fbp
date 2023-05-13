% exemplary application of
% the filtered backpropagation (FBP) algorithm
%
% -------------------------------------------------------------------------
% REFERENCES:
% -------------------------------------------------------------------------
%   [1] M. F. Schiffner and G. Schmitz, “Plane wave pulse-echo ultrasound diffraction tomography with a fixed linear transducer array,” in
%       Acoust. Imaging, ser. Acoust. Imaging, A. Nowicki, J. Litniewski, and T. Kujawska, Eds., vol. 31, Springer Netherlands, 2012, pp. 19–30.
%       DOI : 10.1007/978-94-007-2619-2_3
%   [2] A. J. Devaney, “A filtered backpropagation algorithm for diffraction tomography,”
%       Ultrasonic Imaging, vol. 4, no. 4, pp. 336–350, Oct. 1982.
%       DOI : 10.1016/0161-7346(82)90017-7
%
% -------------------------------------------------------------------------
% ABOUT:
% -------------------------------------------------------------------------
%   author: Martin F. Schiffner
%   date: 2023-05-03
%   modified: 2023-05-13

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% clear workspace
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 0.) parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load RF data acquired from tissue phantom
load( 'data_RF.mat' );

% bandwidth
f_bounds = [ 2.25, 6.75 ] * 1e6;

% dependent parameters
positions_z = ( 64 + (0:511) ) * element_pitch / 4;

% dynamic range for illustration
dynamic_range_dB = 60;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1.) beamform RF data and show results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% call filtered backpropagation (FBP)
[ image_compound, images_single, positions_x ] = fbp_pw( element_pitch, positions_z, data_RF, f_s, theta_incident, c_avg, f_bounds );

%--------------------------------------------------------------------------
% 1.) show results
%--------------------------------------------------------------------------
c_limits = [ -dynamic_range_dB, 0 ];

figure( 1 );
for index_theta = 1:numel( theta_incident )
    subplot( 1, numel( theta_incident ), index_theta );
    imagesc( positions_x * 1e3, positions_z * 1e3, 20 * log10( abs( images_single( :, :, index_theta ) ) / max( max( abs( images_single( :, :, index_theta ) ) ) ) ), c_limits );
    title( sprintf( 'Steering angle: %.1f°', rad2deg( theta_incident( index_theta ) ) ) );
    xlabel( 'Lateral position (mm)' );
    ylabel( 'Axial position (mm)' );
end

figure( 2 );
imagesc( positions_x * 1e3, positions_z * 1e3, 20 * log10( abs( image_compound ) / max( abs( image_compound( : ) ) ) ), c_limits );
title( 'Coherent compounding' );
xlabel( 'Lateral position (mm)' );
ylabel( 'Axial position (mm)' );
colormap gray;
