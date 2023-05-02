function [ gamma_kappa_recon, gamma_kappa_recon_theta, gamma_kappa_recon_compensated, image_pos_x ] = fbp_pw( data_RF, A_in_td, f_lb, f_ub, A_in_analy_threshold_dB, samples_delta_x, image_pos_z, theta_incident, c_0, f_s, factor_interp, factor_zero_pad )
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
% INPUTS:
%   data_RF 	% matrix of RF data
%   A_in_td 	% incident pulse (waveform of incident PW)
%   f_lb        % lower frequency bound
%   f_ub        % upper frequency bound
%   A_in_analy_threshold_dB % threshold for pseudo-inverse filter in dB (used for deconvolution with A_in_td)
%   samples_delta_x % element pitch
%   image_pos_z     % axial positions in image (m)
%   theta_incident  % steering angles (rad)
%   c_0
%   f_s
%   factor_interp
%   factor_zero_pad
%
% OUTPUTS:
%   gamma_kappa_recon = recovered compressibility fluctuations (1)
%   gamma_kappa_recon_theta = recovered compressibility fluctuations for each propagation direction (1)
%   gamma_kappa_recon_compensated = compensated recovered compressibility fluctuations (1)
%   image_pos_x = lateral positions (m)
%
% REFERENCES:
%	[1] M. F. Schiffner and G. Schmitz, “Plane wave pulse-echo ultrasound diffraction tomography with a fixed linear transducer array,” in
%       Acoust. Imaging, ser. Acoust. Imaging, A. Nowicki, J. Litniewski, and T. Kujawska, Eds., vol. 31, Springer Netherlands, 2012, pp. 19–30.
%       DOI : 10.1007/978-94-007-2619-2_3
%   [2] A. J. Devaney, “A filtered backpropagation algorithm for diffraction tomography,”
%       Ultrasonic Imaging, vol. 4, no. 4, pp. 336–350, Oct. 1982.
%       DOI : 10.1016/0161-7346(82)90017-7
%
% REMARKS:
%	- It is unclear how to deal with the overlapping spectra in the compound image.
%	- The compensation eliminates the effects of the finite aperture.
%   - The two-dimensional space underestimates the diffraction-induced decay of the signals.
%
% author: Martin F. Schiffner
% date: 2013-01-26
% modified: 2020-09-09

%% compute angular spectra for all relevant frequencies

% number of samples in time domain (use odd number)
N_samples_t = size(data_RF, 1);
if mod(N_samples_t, 2) == 0
    N_samples_t = N_samples_t + 1;
end
% assertion: N_samples_t is odd

% number of elements in spatial domain
N_elements = size(data_RF, 2);
M_elements = (N_elements - 1) / 2;

% position of transmit receive elements
%(for even N_elements: no element in pos_rx_x == 0)
pos_rx_x = (-M_elements:M_elements) * samples_delta_x;

% logical number of samples used in DFT
N_samples = factor_zero_pad * N_elements;
M_samples = (N_samples - 1) / 2;
N_zeros_padded = N_samples - N_elements;

if mod(N_samples, 2) == 0
    
    N_samples = N_samples + 1;
    M_samples = ceil(M_samples);
    N_zeros_padded = N_zeros_padded + 1;
end
% assertion: N_samples is odd
pos_rx_x = [pos_rx_x, pos_rx_x(end) + (1:N_zeros_padded) * samples_delta_x];

% axis of angular spatial frequencies
axis_k_x = 2 * pi * (-M_samples:M_samples) / (N_samples * samples_delta_x);

% angles of incident plane waves
N_theta = numel( theta_incident );
e_theta = [ cos( theta_incident( : )' ); sin( theta_incident( : )' ) ];

%time frequency decomposition, frequency-continuous phasors
data_RF_phasors = fft(data_RF, N_samples_t, 1) / (sqrt(N_samples_t) * f_s);

%compute corresponding frequency axis
axis_omega =  2 * pi * (0:(N_samples_t - 1)) * f_s / N_samples_t;

% bandpass filter (only evaluate relevant frequencies)
omega_lb = 2 * pi * f_lb;
omega_ub = 2 * pi * f_ub;

indicator_k = (axis_omega >= omega_lb) & (axis_omega <= omega_ub);
axis_omega = axis_omega(indicator_k);
axis_k_norm = axis_omega / c_0;

N_samples_k = numel(axis_k_norm);

% expected frequency-continuous spectrum of pressure wave
A_in = fft(A_in_td, N_samples_t) / (sqrt(N_samples_t) * f_s);
A_in_analy_cropped = 2 * A_in(indicator_k);

A_in_analy_cropped_abs_max = max( abs( A_in_analy_cropped ) );
A_in_analy_threshold = 10^(-A_in_analy_threshold_dB / 20) * A_in_analy_cropped_abs_max;

filter_A_in_inv_pseudo = zeros(1, N_samples_k);
indicator = abs(A_in_analy_cropped) > A_in_analy_threshold;
filter_A_in_inv_pseudo(indicator) = 1 ./ A_in_analy_cropped(indicator);
filter_A_in_inv_pseudo(~indicator) = abs(A_in_analy_cropped(~indicator)) ./ (A_in_analy_cropped(~indicator) * A_in_analy_threshold);

figure(1);
subplot(2,2,1);
plot((1:N_samples_k), abs(filter_A_in_inv_pseudo));
subplot(2,2,2);
plot((1:N_samples_k), unwrap(angle(filter_A_in_inv_pseudo)));
subplot(2,2,3);
plot((1:N_samples_k), abs(filter_A_in_inv_pseudo .* A_in_analy_cropped));
subplot(2,2,4);
plot((1:N_samples_k), unwrap(angle(filter_A_in_inv_pseudo .* A_in_analy_cropped)));

data_RF_phasors_cropped = data_RF_phasors(indicator_k, :, :);
clear data_RF_phasors;

%compute frequency-continuous angular spectrum (correct phase shift caused by pos_rx_x(1) ~= 0)
data_RF_phasors_cropped_ang_spec = fft(data_RF_phasors_cropped, N_samples, 2) * samples_delta_x / sqrt(N_samples);
data_RF_phasors_cropped_ang_spec = fftshift(data_RF_phasors_cropped_ang_spec, 2) .* repmat(exp(-1j * axis_k_x * pos_rx_x(1)), [N_samples_k, 1, N_theta]);
clear data_RF_phasors_cropped;

%% create image using filtered backprojection for each frequency and each angle

%setup image geometry
N_image_axis = [factor_interp * N_samples, numel(image_pos_z)];   %number of lattice points on each axis
M_image_axis = (N_image_axis - 1) / 2;

image_delta_x = samples_delta_x / factor_interp;                  %spatial sampling interval in lateral direction
image_pos_x = (-M_image_axis(1):M_image_axis(1)) * image_delta_x; %lateral positions

% allocate memory for results
gamma_kappa_recon_theta = zeros( N_image_axis( 2 ), N_image_axis( 1 ), N_theta );

%iterate over transmit angles
for index_theta = 1:N_theta

    fprintf('index_theta = %d of %d\n', index_theta, N_theta);

    %TODO: window data in a smart way to utilize all measurements
    %TODO: compare to added DAS
    if 0
    %assertion: N_theta >= 2
        if index_theta == 1
         
            theta_lower_bound = 0;
            theta_upper_bound = theta_incident(index_theta + 1);
     
        elseif index_theta == N_theta
         
            theta_lower_bound = theta_incident(index_theta - 1);
            theta_upper_bound = pi;
        else
        
            theta_lower_bound = theta_incident(index_theta - 1);
            theta_upper_bound = theta_incident(index_theta + 1);
        end
    else
          
        theta_lower_bound = 0;
        theta_upper_bound = pi;
    end
    
    %modify angular spectrum to account for steering angle (element that fires first)
    % check steering direction to find the element that fires first
    if e_theta(1, index_theta) >= 0
        % steering direction: first quadrant
        pos_rx_x_ref = pos_rx_x(1);
    else
        % steering direction: second quadrant
        % assertion: e_theta(1, index_theta) < 0
        pos_rx_x_ref = pos_rx_x(N_elements);
    end

    %iterate over frequencies
    for index_k = 1:N_samples_k
        
        tic;
        fprintf('index_k = %d of %d...', index_k, N_samples_k);
    
        %compute filter kernel
        k_x_over_k_0 = axis_k_x ./ axis_k_norm(index_k);
        sqrt_k_x_over_k_0 = sqrt(1 - (k_x_over_k_0).^2);
    
        filter_1 = e_theta(1, index_theta) * k_x_over_k_0 + e_theta(2, index_theta) * sqrt_k_x_over_k_0 + 1;
        
        phase_pw = exp(1j * axis_k_norm(index_k) * (e_theta(1, index_theta) * repmat(image_pos_x - pos_rx_x_ref, [N_image_axis(2), 1]) + e_theta(2, index_theta) * repmat(image_pos_z', [1, N_image_axis(1)])) );
        filter_2 = exp(1j * axis_k_norm(index_k) * repmat(image_pos_z',  [1, N_samples]) .* repmat(sqrt_k_x_over_k_0, [N_image_axis(2), 1]));
        
        %compute valid spatial frequencies for backpropagation
        indicator_all = abs(axis_k_x) < axis_k_norm(index_k);
                
        indicator_unique_1 = axis_k_x < (axis_k_norm(index_k) * cos(theta_lower_bound));
        indicator_unique_2 = axis_k_x > (axis_k_norm(index_k) * cos(theta_upper_bound));
         
        indicator_unique = indicator_unique_1 & indicator_unique_2;
        indicator_multiple = indicator_all & ~indicator_unique;
        
        %indicator_lower = indicator_all & ~indicator_2;
        %indicator_upper = indicator_all & ~indicator_1;
        
        if sum(indicator_unique) > 0
            
           % commit backpropagation
           data_RF_phasors_cropped_ang_spec_filtered = zeros(1, N_samples);
           data_RF_phasors_cropped_ang_spec_filtered(indicator_unique) = data_RF_phasors_cropped_ang_spec(index_k, indicator_unique, index_theta) .* filter_1(indicator_unique);
           data_RF_phasors_cropped_ang_spec_filtered(indicator_multiple) = data_RF_phasors_cropped_ang_spec(index_k, indicator_multiple, index_theta) .* filter_1(indicator_multiple) / N_theta;
           
           % apply pseudo inverse filter
           data_RF_phasors_cropped_ang_spec_filtered = data_RF_phasors_cropped_ang_spec_filtered * filter_A_in_inv_pseudo(index_k);
           
           gamma_kappa_recon_dft_shift = zeros(N_image_axis(2), N_samples);
           gamma_kappa_recon_dft_shift(:, indicator_all) = repmat(data_RF_phasors_cropped_ang_spec_filtered(1,indicator_all), [N_image_axis(2), 1]) .* filter_2(:, indicator_all) / axis_k_norm(index_k);
           gamma_kappa_recon_dft = ifftshift(gamma_kappa_recon_dft_shift, 2);
           
           %interpolate backpropagated image in frequency domain
           gamma_kappa_recon_interp_dft = zeros(N_image_axis(2), N_image_axis(1));
           
           gamma_kappa_recon_interp_dft(:, 1:(M_samples + 1)) = gamma_kappa_recon_dft(:, 1:(M_samples + 1));   %zero insertion
           gamma_kappa_recon_interp_dft(:, (N_image_axis(1) - M_samples + 1):end) = gamma_kappa_recon_dft(:, (M_samples + 2):end);

           temp = fftshift(ifft(gamma_kappa_recon_interp_dft, [], 2), 2) .* phase_pw;
           gamma_kappa_recon_theta(:,:,index_theta) = gamma_kappa_recon_theta(:,:,index_theta) + temp;
         
        end %if sum(indicator_unique) > 0
        
        time = toc;
        fprintf('(%f ms)\n', time * 1000);
    end %for index_k = 1:N_samples_k
    
end %for theta_act = 1:N_theta

gamma_kappa_recon = sum( gamma_kappa_recon_theta, 3 ) / N_theta;

% compensate received power ratio in compound image
[ image_pos_X, image_pos_Z ] = meshgrid( image_pos_x, image_pos_z );

relative_power = (atan((pos_rx_x(N_elements) - image_pos_X) ./ image_pos_Z) - atan((pos_rx_x(1) - image_pos_X) ./ image_pos_Z)) / pi;
gamma_kappa_recon_compensated = gamma_kappa_recon ./ sqrt(relative_power);


