function [ gamma_kappa_recon, gamma_kappa_recon_theta, image_pos_x ] = fbp_pw( element_pitch, positions_z, data_RF, f_s, theta_incident, c_0, f_bounds, A_in_td, A_in_analy_threshold_dB, factor_interp, factor_zero_pad )
%
% Filtered backpropagation (FBP) algorithm for
% the RF data obtained from
% steered plane waves (PWs) in
% the two-dimensional space.
%
% The function applies
% the FBP algorithm [1, 2] to
% RF data.
% The RF data is acquired in
% multiple sequential pulse-echo measurements.
% The FBP algorithm recovers
% the relative spatial fluctuations in
% compressibility.
%
% minimal usage:
% gamma_kappa_recon = fbp_pw( element_pitch, positions_z, data_RF, f_s, theta_incident, c_0 );
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
%   01.) element_pitch:    element pitch (m)
%   02.) positions_z:      axial positions in image (m)
%   03.) data_RF:          matrix of RF data
%   04.) f_s:              sampling frequency (Hz)
%   05.) theta_incident:   steering angles (rad)    -> -pi/2, <- pi/2
%   06.) c_0:              speed of sound (m/s)
%
% OPTIONAL
%   07.) f_bounds:         frequency bounds (Hz)
%   08.) A_in_td:          incident pulse (waveform of incident PW)
%   09.) A_in_analy_threshold_dB: threshold for pseudo-inverse filter in dB (used for deconvolution with A_in_td)
%   10.) factor_interp:
%   11.) factor_zero_pad:
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
%   01.) gamma_kappa_recon = recovered compressibility fluctuations (1)
%   02.) gamma_kappa_recon_theta = recovered compressibility fluctuations for each propagation direction (1)
%   03.) image_pos_x = lateral positions (m)
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
% REMARKS:
% -------------------------------------------------------------------------
% - The two-dimensional space underestimates the diffraction-induced decay of the signals ( 1 / sqrt(R) vs 1 / R ).
% - The finite aperture of the linear array acts as an undesired window function.
%   The acoustic pressure outside the aperture is erroneously set to zero.
%	- The compensation eliminates the effects of the finite aperture on the modulus of the image pixels.
%	- It is unclear how to deal with the overlapping spectra in the compound image.
%
% -------------------------------------------------------------------------
% ABOUT:
% -------------------------------------------------------------------------
%   author: Martin F. Schiffner
%   date: 2013-01-26
%   modified: 2023-05-12

% print status
time_start = tic;
str_date_time = sprintf( '%04d-%02d-%02d: %02d:%02d:%02d', fix( clock ) );
fprintf( '\t %s: Filtered Backpropagation (FBP) [ steered plane wave ]...\n', str_date_time );

%--------------------------------------------------------------------------
% 1.) check arguments
%--------------------------------------------------------------------------
% ensure at least 6 and at most 11 arguments
narginchk( 6, 11 );

% ensure positive element_pitch
mustBePositive( element_pitch );

% ensure positive positions_z
mustBePositive( positions_z );

% ensure two- or three-dimensional array for data_RF
if ndims( data_RF ) > 3
	errorStruct.message = 'data_RF must be a two- or three-dimensional array!';
    errorStruct.identifier = 'fbp_pw:NoArray';
    error( errorStruct );
end

% ensure at least two spatial positions and two samples
if ~( size( data_RF, 1 ) > 1 && size( data_RF, 2 ) > 1  )
    errorStruct.message = 'data_RF must provide at least two temporal and two spatial samples!';
    errorStruct.identifier = 'fbp_pw:NoSamples';
    error( errorStruct );
end

% ensure positive f_s
mustBePositive( f_s );

% ensure valid theta_incident
mustBeInRange( theta_incident, - pi / 2, pi / 2 );

% ensure positive c_0
mustBePositive( c_0 );

% ensure existence of nonempty f_bounds
if nargin < 7 || isempty( f_bounds )
    f_bounds = [ 0, f_s / 2 ];
end

% ensure existence of nonempty A_in_td
if nargin < 8 || isempty( A_in_td )
    A_in_td = 1;
end

% ensure existence of nonempty A_in_analy_threshold_dB
if nargin < 9 || isempty( A_in_analy_threshold_dB )
    A_in_analy_threshold_dB = 6;
end

% ensure existence of nonempty factor_interp
if nargin < 10 || isempty( factor_interp )
    factor_interp = 4;
end

% ensure existence of nonempty factor_zero_pad
if nargin < 11 || isempty( factor_zero_pad )
    factor_zero_pad = 2;
end

%--------------------------------------------------------------------------
% 2.) time-frequency decomposition
%--------------------------------------------------------------------------
% number of samples in the time domain
N_samples_t = size( data_RF, 1 );

% number of points used in temporal DFT
% TODO: zero padding?
N_points_dft_t = N_samples_t;
if mod( N_points_dft_t, 2 ) == 0
    N_points_dft_t = N_points_dft_t + 1;
end
% assertion: N_points_dft_t is odd

% boundary frequency indices
index_Omega_lb = ceil( f_bounds( 1 ) * N_points_dft_t / f_s ) + 1;
index_Omega_ub = floor( f_bounds( 2 ) * N_points_dft_t / f_s ) + 1;
indices_Omega = (index_Omega_lb:index_Omega_ub).';

% frequency axis
axis_f_bp = f_s * ( indices_Omega - 1 ) / N_points_dft_t;
axis_k_bp = 2 * pi * axis_f_bp / c_0;

% number of relevant discrete frequencies
N_samples_f = numel( indices_Omega );

% compute frequency-continuous phasors
data_RF_phasors = fft( data_RF, N_points_dft_t,  1 ) / ( sqrt( N_points_dft_t ) * f_s );
data_RF_phasors_cropped = data_RF_phasors( indices_Omega, :, : );
clear data_RF_phasors;

% expected frequency-continuous spectrum of pressure wave
A_in = fft( A_in_td, N_points_dft_t ) / ( sqrt( N_points_dft_t ) * f_s );
A_in_analy_cropped = 2 * A_in( indices_Omega );

%--------------------------------------------------------------------------
% 3.) space-frequency decomposition
%--------------------------------------------------------------------------
% number of array elements
N_elements = size( data_RF, 2 );
M_elements = ( N_elements - 1 ) / 2;

% number of points used in lateral DFT
N_points_dft_x = factor_zero_pad * N_elements;
if mod( N_points_dft_x, 2 ) == 0
    N_points_dft_x = N_points_dft_x + 1;
end
% assertion: N_points_dft_x is odd

N_zeros_padded = N_points_dft_x - N_elements;

% centroids of element faces
%(for even N_elements: no element in pos_rx_x == 0)
% positions_ctr_x = (-M_elements:M_elements) * element_pitch;
pos_rx_x = [ (-M_elements:M_elements), M_elements + (1:N_zeros_padded) ] * element_pitch;

% axis of angular spatial frequencies
M_points_dft_x = ( N_points_dft_x - 1 ) / 2;
axis_k_x = 2 * pi * (-M_points_dft_x:M_points_dft_x) / ( N_points_dft_x * element_pitch );

% normalized spatial frequencies (direction vectors)
mat_k_x_norm = axis_k_x ./ axis_k_bp;
mat_k_z_norm = sqrt( 1 - mat_k_x_norm.^2 );

% propagable plane waves (PWs)
indicator_propagable = abs( mat_k_x_norm ) < 1;

% compute frequency-continuous angular spectrum (correct phase shift caused by pos_rx_x(1) ~= 0)
data_RF_phasors_cropped_ang_spec = fft( data_RF_phasors_cropped, N_points_dft_x, 2 ) * element_pitch / sqrt( N_points_dft_x );
data_RF_phasors_cropped_ang_spec = fftshift( data_RF_phasors_cropped_ang_spec, 2 ) .* exp( -1j * axis_k_x * pos_rx_x( 1 ) );
clear data_RF_phasors_cropped;

% angles of incident plane waves
N_theta = numel( theta_incident );
e_theta = [ sin( theta_incident( : )' ); cos( theta_incident( : )' ) ];

% modify angular spectrum to account for steering angle (element that fires first)
% check steering direction to find the element that fires first
pos_rx_x_ref = repmat( pos_rx_x( 1 ), [ 1, N_theta ] ); % steering direction: first quadrant
indicator = e_theta( 1, : ) < 0;
pos_rx_x_ref( indicator ) = pos_rx_x( N_elements ); % steering direction: second quadrant

A_in_analy_cropped_abs_max = max( abs( A_in_analy_cropped ) );
A_in_analy_threshold = 10^(-A_in_analy_threshold_dB / 20) * A_in_analy_cropped_abs_max;

filter_A_in_inv_pseudo = zeros( N_samples_f, 1 );
indicator = abs( A_in_analy_cropped ) > A_in_analy_threshold;
filter_A_in_inv_pseudo( indicator ) = 1 ./ A_in_analy_cropped(indicator);
filter_A_in_inv_pseudo( ~indicator ) = abs(A_in_analy_cropped(~indicator)) ./ (A_in_analy_cropped(~indicator) * A_in_analy_threshold);

% figure(1);
% subplot(2,2,1);
% plot((1:N_samples_f), abs(filter_A_in_inv_pseudo));
% subplot(2,2,2);
% plot((1:N_samples_f), unwrap(angle(filter_A_in_inv_pseudo)));
% subplot(2,2,3);
% plot((1:N_samples_f), abs(filter_A_in_inv_pseudo .* A_in_analy_cropped));
% subplot(2,2,4);
% plot((1:N_samples_f), unwrap(angle(filter_A_in_inv_pseudo .* A_in_analy_cropped)));

%--------------------------------------------------------------------------
% 4.) compute relative spatial fluctuations in compressibility
%--------------------------------------------------------------------------
% setup image geometry
N_image_axis = [ factor_interp * N_points_dft_x, numel( positions_z ) ]; % number of lattice points on each axis
M_image_axis = ( N_image_axis - 1 ) / 2;

image_delta_x = element_pitch / factor_interp; % spatial sampling interval in lateral direction
image_pos_x = (-M_image_axis( 1 ):M_image_axis( 1 )) * image_delta_x; % lateral positions

% allocate memory for results
gamma_kappa_recon_theta = zeros( N_image_axis( 2 ), N_image_axis( 1 ), N_theta );

% iterate steering angles
for index_theta = 1:N_theta

    % print status
    time_start_theta = tic;
    fprintf( '\t\tsteering angle: %6.2f° (%d of %d)... ', rad2deg( theta_incident( index_theta ) ), index_theta, N_theta );

    % current propagation direction
    e_theta_x_act = e_theta( 1, index_theta );
    e_theta_z_act = e_theta( 2, index_theta );

    % current reference element
	pos_rx_x_ref_act = pos_rx_x_ref( index_theta );

    % phasors
    data_RF_phasors_cropped_ang_spec_act = data_RF_phasors_cropped_ang_spec( :, :, index_theta );

    %----------------------------------------------------------------------
    % a) filter fields
    %----------------------------------------------------------------------
    % compute filter kernel
    filter = zeros( N_samples_f, N_points_dft_x );
    filter( indicator_propagable ) = e_theta_x_act * mat_k_x_norm( indicator_propagable ) + e_theta_z_act * mat_k_z_norm( indicator_propagable ) + 1;

    % apply filter
    data_RF_phasors_cropped_ang_spec_filtered = zeros( N_samples_f, N_points_dft_x );
    data_RF_phasors_cropped_ang_spec_filtered( indicator_propagable ) = data_RF_phasors_cropped_ang_spec_act( indicator_propagable ) .* filter( indicator_propagable );

    % apply pseudo inverse filter
    data_RF_phasors_cropped_ang_spec_filtered = data_RF_phasors_cropped_ang_spec_filtered .* filter_A_in_inv_pseudo ./ axis_k_bp;

    % iterate axial positions
    for index_pos_z = 1:N_image_axis( 2 )

        % print progress in percent
        if mod( index_pos_z, 50 ) == 0
            fprintf( '%5.1f %%', ( index_pos_z - 1 ) / N_image_axis( 2 ) * 1e2 );
        end

        %------------------------------------------------------------------
        % b) backpropagate filtered fields
        %------------------------------------------------------------------
        % compute backpropagator
        propagator = exp( 1j * axis_k_bp .* mat_k_z_norm * positions_z( index_pos_z ) );

        % backpropagate
        gamma_kappa_recon_dft_shift = zeros( N_samples_f, N_points_dft_x );
        gamma_kappa_recon_dft_shift( indicator_propagable ) = data_RF_phasors_cropped_ang_spec_filtered( indicator_propagable ) .* propagator( indicator_propagable );
        gamma_kappa_recon_dft = ifftshift( gamma_kappa_recon_dft_shift, 2 );

        %------------------------------------------------------------------
        % c) inverse lateral Fourier transform
        %------------------------------------------------------------------
        % lateral interpolation by zero insertion
        gamma_kappa_recon_interp_dft = zeros( N_samples_f, N_image_axis( 1 ) );
        gamma_kappa_recon_interp_dft( :, 1:(M_points_dft_x + 1) ) = gamma_kappa_recon_dft( :, 1:(M_points_dft_x + 1) );
        gamma_kappa_recon_interp_dft( :, ( N_image_axis( 1 ) - M_points_dft_x + 1 ):end ) = gamma_kappa_recon_dft( :, (M_points_dft_x + 2):end );

        % phase of incident plane wave
        phase_pw = exp( 1j * axis_k_bp .* ( e_theta_x_act * ( image_pos_x - pos_rx_x_ref_act ) + e_theta_z_act * positions_z( index_pos_z ) ) );

        temp = fftshift( ifft( gamma_kappa_recon_interp_dft, [], 2 ), 2 ) .* phase_pw;
        gamma_kappa_recon_theta( index_pos_z, :, index_theta ) = sum( temp, 1 );

        % erase progress in percent
        if mod( index_pos_z, 50 ) == 0
            fprintf( '\b\b\b\b\b\b\b' );
        end

    end % for index_pos_z = 1:N_image_axis( 2 )

    % infer and print elapsed time
    time_elapsed_theta = toc( time_start_theta );
    fprintf( 'done! (%.3f s)\n', time_elapsed_theta );

end % for index_theta = 1:N_theta

% compound image
gamma_kappa_recon = sum( gamma_kappa_recon_theta, 3 ) / N_theta;

% infer and print elapsed time
time_elapsed = toc( time_start );
fprintf( '\t\tdone! (%.3f s)\n', time_elapsed );

end % function [ gamma_kappa_recon, gamma_kappa_recon_theta, image_pos_x ] = fbp_pw( element_pitch, positions_z, data_RF, f_s, theta_incident, c_0, f_bounds, A_in_td, A_in_analy_threshold_dB, factor_interp, factor_zero_pad )
