
% Example of use of the new Field II program running under Matlab
%
% This example shows how a linear array B-mode system scans an image
%
% This script assumes that the field_init procedure has been called
%
% Example by Joergen Arendt Jensen, Version 2.0, March 22, 2011.
% Generate the transducer apertures for send and receive

f0=3e6; % Transducer center frequency [Hz]
fs=100e6; % Sampling frequency [Hz]
c=1540; % Speed of sound [m/s]
lambda=c/f0; % Wave length [m]
width=lambda; % Width of element
element_height=5/1000; % Height of element [m]
kerf=width/20; % Kerf [m]
focus=[0 0 50]/1000; % Fixed focal point [m] {This sets the fixed focal point at 50 mm depth when creating the transducer aperture{
N_elements=192; % Number of elements in the transducer
N_active=64; % Active elements in the transducer
% Set the sampling frequency
set_sampling(fs);

% Generate aperture for emission
emit_aperture = xdc_linear_array (N_elements, width, element_height, kerf, 1, 5, focus);
% Set the impulse response and excitation of the emit aperture
impulse_response=sin(2*pi*f0*(0:1/fs:2/f0));
N = max(size(impulse_response));
impulse_response = impulse_response .* (0.5 - 0.5*cos(2*pi*(0:N-1)/(N-1)));
xdc_impulse (emit_aperture, impulse_response);
excitation=sin(2*pi*f0*(0:1/fs:2/f0));
xdc_excitation (emit_aperture, excitation);

% Generate aperture for reception
receive_aperture = xdc_linear_array (N_elements, width, element_height, kerf, 1, 5, focus);
% Set the impulse response for the receive aperture
xdc_impulse (receive_aperture, impulse_response);

% Load the computer phantom
[phantom_positions, phantom_amplitudes] = cyst_phantom(5000);

%Phanton would go here
[phantom_positions, phantom_amplitudes] = cyst_phantom(5000);

% Add a specular reflector (strong point target)

% դիր: center at 0 mm lateral, 35 mm depth
specular_position = [0 0 50]/1000;   % [x y z] in meters

% very strong amplitude (much higher than background)
specular_amplitude = 100;

% append to phantom
phantom_positions = [phantom_positions; specular_position];
phantom_amplitudes = [phantom_amplitudes; specular_amplitude];
% Do linear array imaging
no_lines=N_elements-N_active+1; % Number of A-lines in image
dx=width; % Increment for image
z_focus=50/1000;
% Pre-allocate some storage
image_data=zeros(1,no_lines);

for i=1:no_lines
i;
% Find position for imaging
x=(i-1-no_lines/2)*dx;

z_focus=50/1000;  % {this sets the dynamic focus to 50 mm during scanning.}

% Set the focus for this direction
xdc_center_focus (emit_aperture, [x 0 0]);
xdc_focus (emit_aperture, 0, [x 0 z_focus]);
xdc_center_focus (receive_aperture, [x 0 0]);
xdc_focus (receive_aperture, 0, [x 0 z_focus]);
% Set the active elements using the apodization
hamming_win = (0.54 - 0.46*cos(2*pi*(0:N_active-1)/(N_active-1)));
apo=[zeros(1, i-1) hamming_win zeros(1, N_elements-N_active-i+1)];
xdc_apodization (emit_aperture, 0, apo);
xdc_apodization (receive_aperture, 0, apo);
% Calculate the received response
[v, t1]=calc_scat(emit_aperture, receive_aperture, phantom_positions, phantom_amplitudes);
% Store the result
image_data(1:max(size(v)),i)=v;
times(i) = t1;
end
% Free space for apertures
xdc_free (emit_aperture)
xdc_free (receive_aperture)
% Adjust the data in time and display it as
% a gray scale image
min_sample=min(times)*fs;
for i=1:no_lines
sig = [zeros(round(times(i)*fs-min_sample),1); image_data(:,i)];
N = length(sig);
X = fft(sig);
h = zeros(N,1);
if mod(N,2) == 0
    h([1, N/2+1]) = 1;
    h(2:N/2) = 2;
else
    h(1) = 1;
    h(2:(N+1)/2) = 2;
end
rf_env = abs(ifft(X.*h));
env(1:size(rf_env,1),i)=rf_env;
end
% make logarithmic compression to a 60 dB dynamic range
% with proper units on the axis
env_dB=20*log10(env);
env_dB=env_dB-max(max(env_dB));
env_gray=127*(env_dB+60)/60;
depth=((0:size(env,1)-1)+min_sample)/fs*c/2;
x=((1:no_lines)-no_lines/2)*dx;
image(x*1000, depth*1000, env_gray)
xlabel('Lateral distance [mm]')
ylabel('Depth [mm]')
axis('image')
colormap(gray(128))
title('Image of cyst phantom (60 dB dynamic range)')