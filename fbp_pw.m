function [ gamma_kappa_recon, gamma_kappa_recon_theta, image_pos_x ] = fbp_pw( element_pitch, positions_z, data_RF, f_s, theta_incident, c_0, f_bounds, F_number, A_in_td, A_in_analy_threshold_dB, factor_interp, N_zeros_pad )
% fbp_pw  Filtered backpropagation algorithm for steered plane waves.
%
% The filtered backpropagation (FBP) algorithm [1, 2] uses
% the radio frequency (RF) data acquired in
% multiple sequential pulse-echo measurements to recover
% the relative spatial fluctuations in
% compressibility.
%
% The FBP algorithm assumes
% the two-dimensional space,
% weak scattering (Born approximation),
% the transmission of steered plane waves (PWs), and
% a linear array.
%
% minimal usage:
% gamma_kappa_recon = fbp_pw( element_pitch, positions_z, data_RF, f_s, theta_incident, c_0 );
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
%   01.) element_pitch:    element pitch of the linear array (m)
%   02.) positions_z:      axial positions in image (m)
%   03.) data_RF:          array of RF data (a.u.)
%   04.) f_s:              sampling frequency (Hz)
%   05.) theta_incident:   steering angles (rad)    -> -pi/2, <- pi/2
%   06.) c_0:              speed of sound (m/s)
%
% OPTIONAL
%   07.) f_bounds:         frequency bounds (Hz)
%   08.) F_number:         receive F-number (1)
%   09.) A_in_td:          waveform of transmitted PW (a.u.)
%   10.) A_in_analy_threshold_dB: threshold for pseudo-inverse filter in dB (used for deconvolution with A_in_td)
%   11.) factor_interp:    lateral interpolation factor (1)
%   12.) N_zeros_pad:      number of zeros to pad in lateral discrete Fourier transform (1)
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
%   01.) gamma_kappa_recon:       recovered compressibility fluctuations (1)
%   02.) gamma_kappa_recon_theta: recovered compressibility fluctuations for each steering angle (1)
%   03.) image_pos_x:             lateral positions in image (m)
%
% -------------------------------------------------------------------------
% REFERENCES:
% -------------------------------------------------------------------------
%   [1] M. F. Schiffner and G. Schmitz, “Plane wave pulse-echo ultrasound diffraction tomography with a fixed linear transducer array,” in
%       Acoust. Imaging, ser. Acoust. Imaging, A. Nowicki, J. Litniewski, and T. Kujawska, Eds., vol. 31, Springer Netherlands, 2012, pp. 19–30.
%       DOI : <a href="matlab:web('https://dx.doi.org/10.1007/978-94-007-2619-2_3')">10.1007/978-94-007-2619-2_3</a>
%   [2] A. J. Devaney, “A filtered backpropagation algorithm for diffraction tomography,”
%       Ultrasonic Imaging, vol. 4, no. 4, pp. 336–350, Oct. 1982.
%       DOI : <a href="matlab:web('https://dx.doi.org/10.1016/0161-7346(82)90017-7')">10.1016/0161-7346(82)90017-7</a>
% 
% -------------------------------------------------------------------------
% REMARKS:
% -------------------------------------------------------------------------
% - The two-dimensional space underestimates the diffraction-induced decay of the signals ( 1 / sqrt(R) vs 1 / R ).
% - The finite-sized aperture of the linear array acts as an undesired window function.
%   The acoustic pressure field outside this aperture is erroneously assumed to be zero.
% - It is unclear how to deal with the overlapping spectra in the compound image.
%
% -------------------------------------------------------------------------
% ABOUT:
% -------------------------------------------------------------------------
%   author: Martin F. Schiffner
%   date: 2013-01-26
%   modified: 2023-05-14

% print status
time_start = tic;
str_date_time = sprintf( '%04d-%02d-%02d: %02d:%02d:%02d', fix( clock ) );
fprintf( '%s: Filtered Backpropagation (FBP) [ 2-space, Born approximation, steered plane waves ]...\n', str_date_time );

%--------------------------------------------------------------------------
% 1.) check arguments
%--------------------------------------------------------------------------
% ensure at least 6 and at most 12 arguments
narginchk( 6, 12 );

% ensure positive element_pitch
mustBeScalarOrEmpty( element_pitch );
mustBePositive( element_pitch );

% ensure positive positions_z
mustBePositive( positions_z );

% ensure strictly monotonic increase
if ~issorted( positions_z, 'strictascend' )
    errorStruct.message = 'positions_z must be strictly monotonically increasing!';
    errorStruct.identifier = 'fbp_pw:NoStrictIncreasePosZ';
    error( errorStruct );
end

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
    f_bounds = [ eps, f_s / 2 ];
end

% ensure existence of nonempty F_number
if nargin < 8 || isempty( F_number )
    F_number = 0;
end

% ensure existence of nonempty A_in_td
if nargin < 9 || isempty( A_in_td )
    A_in_td = 1;
end

% ensure existence of nonempty A_in_analy_threshold_dB
if nargin < 10 || isempty( A_in_analy_threshold_dB )
    A_in_analy_threshold_dB = 6;
end

% ensure existence of nonempty factor_interp
if nargin < 11 || isempty( factor_interp )
    factor_interp = 4;
end

% ensure existence of nonempty N_zeros_pad
if nargin < 12 || isempty( N_zeros_pad )
    N_zeros_pad = size( data_RF, 2 ) + 1;
end

%--------------------------------------------------------------------------
% 2.) time-frequency decomposition
%--------------------------------------------------------------------------
% number of samples in the time domain
N_samples_t = size( data_RF, 1 );

% number of points used in the temporal DFT
N_points_dft_t = N_samples_t;
if mod( N_points_dft_t, 2 ) == 0
    N_points_dft_t = N_points_dft_t + 1;
end
% assertion: N_points_dft_t is odd

% boundary frequency indices
index_Omega_lb = ceil( f_bounds( 1 ) * N_points_dft_t / f_s ) + 1;
index_Omega_ub = floor( f_bounds( 2 ) * N_points_dft_t / f_s ) + 1;
indices_Omega = (index_Omega_lb:index_Omega_ub).';

% number of relevant discrete frequencies
N_samples_f = numel( indices_Omega );

% frequency axis
axis_f_bp = f_s * ( indices_Omega - 1 ) / N_points_dft_t;
axis_k_bp = 2 * pi * axis_f_bp / c_0;

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

% number of points used in lateral DFT
N_points_dft_x = N_elements + N_zeros_pad;

% center of leftmost array element
pos_rx_x_left = - ( N_elements - 1 ) * element_pitch / 2;

% axis of angular spatial frequencies
index_dft_x_shift = ceil( N_points_dft_x / 2 );
axis_k_x = 2 * pi * (( index_dft_x_shift - N_points_dft_x ):(index_dft_x_shift - 1)) / ( N_points_dft_x * element_pitch );

% normalized spatial frequencies (direction vectors)
mat_k_x_norm = axis_k_x ./ axis_k_bp;
mat_k_z_norm = sqrt( 1 - mat_k_x_norm.^2 );

% axial frequencies
mat_k_z = axis_k_bp .* mat_k_z_norm;

% angular aperture
indicator_aperture = abs( mat_k_x_norm ) < 1 / sqrt( 1 + ( 2 * F_number )^2 );

% compute angular spectrum
data_RF_phasors_cropped_ang_spec = fft( data_RF_phasors_cropped, N_points_dft_x, 2 ) * element_pitch / sqrt( N_points_dft_x );

% correct phase shift required by pos_rx_x_left
phase_shift = exp( -1j * axis_k_x * pos_rx_x_left );
if mod( N_points_dft_x, 2 ) == 0
    phase_shift( 1 ) = real( phase_shift( 1 ) );
end
data_RF_phasors_cropped_ang_spec = fftshift( data_RF_phasors_cropped_ang_spec, 2 ) .* phase_shift;
clear data_RF_phasors_cropped;

%--------------------------------------------------------------------------
% 4.) regularization
%--------------------------------------------------------------------------
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
% 5.) compute relative spatial fluctuations in compressibility
%--------------------------------------------------------------------------
% angles of incident plane waves
N_theta = numel( theta_incident );
e_theta = [ sin( theta_incident( : )' ); cos( theta_incident( : )' ) ];

% check steering direction to find the element that fires first
pos_rx_x_ref = repmat( pos_rx_x_left, [ 1, N_theta ] );
indicator = e_theta( 1, : ) < 0;
pos_rx_x_ref( indicator ) = -pos_rx_x_left;

% setup image geometry
N_image_axis = [ factor_interp * N_points_dft_x, numel( positions_z ) ]; % number of lattice points on each axis
M_image_axis = ( N_image_axis - 1 ) / 2;

% spatial sampling interval in lateral direction
image_delta_x = element_pitch / factor_interp;

% lateral and axial positions
image_pos_x = (-M_image_axis( 1 ):M_image_axis( 1 )) * image_delta_x;
positions_z_diff = diff( [ 0, positions_z ] );
positions_z_diff_diff_abs = abs( diff( positions_z_diff ) );

% allocate memory for results
gamma_kappa_recon_theta = zeros( N_image_axis( 2 ), N_image_axis( 1 ), N_theta );

% iterate steering angles
for index_theta = 1:N_theta

    % print status
    time_start_theta = tic;
    fprintf( '\tsteering angle: %6.2f° (%d of %d)... ', rad2deg( theta_incident( index_theta ) ), index_theta, N_theta );

    % current propagation direction
    e_theta_x_act = e_theta( 1, index_theta );
    e_theta_z_act = e_theta( 2, index_theta );

    % current reference element
	pos_rx_x_ref_act = pos_rx_x_ref( index_theta );

    % lateral phase shift
    lateral = e_theta_x_act * ( image_pos_x - pos_rx_x_ref_act );

    % phasors
    data_RF_phasors_cropped_ang_spec_act = data_RF_phasors_cropped_ang_spec( :, :, index_theta );

    %----------------------------------------------------------------------
    % a) filter fields
    %----------------------------------------------------------------------
    % compute filter kernel
    filter = zeros( N_samples_f, N_points_dft_x );
    filter( indicator_aperture ) = e_theta_x_act * mat_k_x_norm( indicator_aperture ) + e_theta_z_act * mat_k_z_norm( indicator_aperture ) + 1;

    % apply filter
    data_RF_phasors_cropped_ang_spec_filtered = zeros( N_samples_f, N_points_dft_x );
    data_RF_phasors_cropped_ang_spec_filtered( indicator_aperture ) = data_RF_phasors_cropped_ang_spec_act( indicator_aperture ) .* filter( indicator_aperture );

    % apply pseudo inverse filter
    data_RF_phasors_cropped_ang_spec_filtered = data_RF_phasors_cropped_ang_spec_filtered .* filter_A_in_inv_pseudo ./ axis_k_bp;

    %----------------------------------------------------------------------
    % b) backpropagate filtered fields
    %----------------------------------------------------------------------
    % compute backpropagator
    propagator = exp( 1j * mat_k_z * positions_z_diff( 1 ) );

    % initialize
    gamma_kappa_recon_dft_shift = data_RF_phasors_cropped_ang_spec_filtered;

    % iterate axial positions
    for index_pos_z = 1:N_image_axis( 2 )

        % print progress in percent
        if mod( index_pos_z, 50 ) == 0
            fprintf( '%5.1f %%', ( index_pos_z - 1 ) / N_image_axis( 2 ) * 1e2 );
        end

        %------------------------------------------------------------------
        % i.) backpropagate by difference in positions_z
        %------------------------------------------------------------------
        % recompute backpropagator if necessary
        if ( index_pos_z > 1 ) && ( positions_z_diff_diff_abs( index_pos_z - 1 ) > eps )
            propagator = exp( 1j * mat_k_z * positions_z_diff( index_pos_z ) );
        end

        % backpropagate
        gamma_kappa_recon_dft_shift( indicator_aperture ) = gamma_kappa_recon_dft_shift( indicator_aperture ) .* propagator( indicator_aperture );
        gamma_kappa_recon_dft = ifftshift( gamma_kappa_recon_dft_shift, 2 );

        %------------------------------------------------------------------
        % ii.) lateral interpolation by zero insertion
        %------------------------------------------------------------------
        gamma_kappa_recon_interp_dft = zeros( N_samples_f, N_image_axis( 1 ) );
        gamma_kappa_recon_interp_dft( :, 1:index_dft_x_shift ) = gamma_kappa_recon_dft( :, 1:index_dft_x_shift );
        gamma_kappa_recon_interp_dft( :, ( N_image_axis( 1 ) - ( N_points_dft_x - index_dft_x_shift ) + 1 ):end ) = gamma_kappa_recon_dft( :, (index_dft_x_shift + 1):end );

        %------------------------------------------------------------------
        % iii.) inverse lateral Fourier transform
        %------------------------------------------------------------------
        % phase of incident plane wave
        phase_pw = exp( 1j * axis_k_bp .* ( lateral + e_theta_z_act * positions_z( index_pos_z ) ) );

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
fprintf( 'done! (%.3f s)\n', time_elapsed );

end % function [ gamma_kappa_recon, gamma_kappa_recon_theta, image_pos_x ] = fbp_pw( element_pitch, positions_z, data_RF, f_s, theta_incident, c_0, f_bounds, F_number, A_in_td, A_in_analy_threshold_dB, factor_interp, N_zeros_pad )
