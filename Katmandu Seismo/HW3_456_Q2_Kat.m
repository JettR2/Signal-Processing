
% Reading in the Variables For use Throughout


% Scanning the TXT file
opts = detectImportOptions('KATNP.accel.e.txt');
preview('KATNP.accel.e.txt',opts)

% defining the time and acceleration vectors
[Time, Accel] = readvars('KATNP.accel.e.txt');
whos Time Accel

%%

% Download the KATNP acceleration file from Canvas. This is the east-west component of an
% accelerogram inside the city of Kathmandu during the 2015 M7.8 Gorkha, Nepal earthquake.
% Plot the waveform and plot the amplitude spectrum using semilogx(). Are there any obvious
% spectral peaks? What do you think they mean?

% Quesiton 1 
Q1Time = Time;
Q1Accel = Accel;

% Plotting Acceleration IN TIME DOMAIN
figure();
subplot(2,1,1);
plot(Q1Time,Q1Accel);
title('Q1 - Raw Data : Time vs Acceleration');
xlabel('Time (seconds)');
ylabel('E-W Acceleration (m/s^2)');
grid on;

% Fourier Transform
Q1FT_Accel = fft(Q1Accel);
% Amplitude spectrum
Q1Amplitude_Spectrum_Accel = abs(Q1FT_Accel);
% To Preserve the Magnitude of the Acceleration
Q1Amplitude_Spectrum_Accel = Q1Amplitude_Spectrum_Accel./1000;

% Parameters to build the frequency domain
Q1Sample_interval = 0.01;
Q1Total_samples = length(Q1Amplitude_Spectrum_Accel);
% Frequency dimension
Q1Sampling_rate = 1 / Q1Sample_interval; % 1/0.01s = 100 Hz
Q1Frequencies = (0:(Q1Total_samples-1)) * (Q1Sampling_rate / Q1Total_samples);

subplot(2,1,2);
semilogx(Q1Frequencies,Q1Amplitude_Spectrum_Accel);
title('Q1 - Amplitude Spectrum');
xlabel('Frequency (Hz)');
ylabel('Amplitude');
grid on;

%%


% Question 2

% Apply a taper of your choosing and plot the spectra with and 
% without the taper. What, if any, are the differences?

Q2Time = Time;
Q2Accel = Accel;


% The Spectra


% Acceleration Data

% Fourier Transform
Q2FT_Accel = fft(Q2Accel);

% Amplitude spectrum
Q2Amplitude_Spectrum_Accel = abs(Q2FT_Accel);
% To Preserve the Magnitude of the Acceleration
Q2Amplitude_Spectrum_Accel = Q2Amplitude_Spectrum_Accel./1000;

% Freq. Domain 


% Parameters to build the frequency domain
Q2Sample_interval = 0.001;
Q2Total_samples = length(Q2Amplitude_Spectrum_Accel);
% Frequency dimension
Q2Sampling_rate = 1 / Q2Sample_interval; % 1/0.001s = 1000 Hz
Q2Frequencies = (0:(Q2Total_samples-1)) * (Q2Sampling_rate / Q2Total_samples);


% The taper

% Hanning Window
Q2Amplitude = max(Accel); % So The Taper is about 1:1 in dimensions with the Data
Q2Hanning_freq = 1/130; % Data ends at 130 seconds
Q2Hanning_win = Q2Amplitude .* (1-cos(2*pi*Q2Hanning_freq*Q2Time)) ./2;

%Fourier Transform and Amp. Spec.
Q2FT_Hann = fft(Q2Hanning_win);
Q2Amp_Spec_Hann = abs(Q2FT_Hann);
Q2Amp_Spec_Hann = Q2Amp_Spec_Hann ./1000;


% The Tapered Data

Q2Tapered_Accel = Q2Accel .* Q2Hanning_win;
Q2FT_Tapered_Accel = fft(Q2Tapered_Accel);

% Amp Spec of Tapered Data
Q2Amp_Spec_Tapered_Accel = abs(Q2FT_Tapered_Accel);
Q2Amp_Spec_Tapered_Accel = Q2Amp_Spec_Tapered_Accel ./ 1000;


% PLotting in the Time Domain


figure();
subplot(2,1,1);
plot(Q2Time, Q2Accel,'DisplayName','Accel');
title('Q2 - Amp. Spec');
xlabel('');
ylabel('Acceleration (m/s^2)');
grid on;
hold on;

% Hanning Window
subplot(2,1,1);
plot(Q2Time,Q2Hanning_win,'DisplayName','Hann Window');
title('Q2 - Acceleration Data and Taper');
xlabel('Time (s)');
ylabel('Acceleration (m/s^2)');
legend();
grid on;
hold off;

% Tapered Acceleration
subplot(2,1,2);
plot(Q2Time,Q2Tapered_Accel,'DisplayName','Tapered Accel');
title('Q2 - Accel. Tapered');
xlabel('Time (s)');
ylabel('Acceleration (m/s^2)');
ylim([-1.5 1.5]);
legend();
grid on;
hold off;


% PLotting in the Frequency Domain

% Just the Amp Spec of Accel
figure();
subplot(3,1,2);
semilogx(Q2Frequencies,Q2Amplitude_Spectrum_Accel);
title('Q2 Acceleration');
ylabel('Amplitude');
xlabel('Frequency (Hz)');

% Acceleration and Hann Window in the same plot
subplot(3,1,1);
semilogx(Q2Frequencies, Q2Amplitude_Spectrum_Accel,'DisplayName','Accel');
title('Q2 - Accel and Taper');
xlabel('Frequency (Hz)');
ylabel('Amplitude');
xlim([0 5]);
grid on;
hold on;

% Hanning Window
subplot(3,1,1);
semilogx(Q2Frequencies,Q2Amp_Spec_Hann,'DisplayName','Amp.Spec Han Window');
title('Q2 - Amplitude Spectrum : Accel and Taper');
xlabel('Frequency (Hz)');
ylabel('Amplitude');
xlim([0 5]);
legend();
grid on;
hold off;

% Tapered Acceleration on the bottom
subplot(3,1,3);
semilogx(Q2Frequencies,Q2Amp_Spec_Tapered_Accel,'DisplayName','Tapered Accel');
title('Q2 -Hann. Win. Tapered Accel.');
xlabel('Frequency (Hz)');
ylabel('Amplitude');
ylim([0 1.2]);
legend();
grid on;
hold off;


%%


% Question 2 iii



% Define parameters
sample_interval = 0.01; % 0.01 seconds per sample
total_samples = 300; % 300 samples
frequency = 1; % Frequency of the sine wave in Hz
amplitude = 1;

% Create a vector of time samples
time = 0:sample_interval:(total_samples-1)*sample_interval;
% Frequency dimension
sampling_rate = 1 / sample_interval; % 1/0.01s = 100 Hz
frequencies = (0:(total_samples-1)) * (sampling_rate / total_samples);
% where (0:(total_samples-1)) is a vector from start to end. 

% Calculate the sine wave
sine_wave = amplitude * sin(2 * pi * frequency * time);

Integrated_Sine = cumtrapz(time,sine_wave);

FT_Sine = -abs(fft(sine_wave));
Division = FT_Sine ./ frequency;
Back_to_time = ifft(Division);

% Plot the sine wave
figure();
subplot(3,1,1);
plot(time, sine_wave);
title('Q1 - Sine Wave');
xlabel('Time (seconds)');
ylabel('Amplitude');
grid on;

subplot(3,1,2);
plot(time,Integrated_Sine);
title('Integrated Sine (-cos)')
xlabel('Time');

subplot(3,1,3);
plot(time,Back_to_time);


%%


% Question 2 i

% Download the KATNP acceleration file from Canvas. This is the east-west 
% component of an accelerogram inside the city of Kathmandu during the 
% 2015 M7.8 Gorkha, Nepal earthquake. Plot the waveform and plot the 
% amplitude spectrum using semilogx(). Are there any obvious spectral 
% peaks? What do you think they mean?

% write up

% There is an obnoxiously large peak around 10^-1 or 1 Hz range. I think
% this is most likely the actual ground motion, while the higher frequency
% action could be smaller tremors and other phenomena that are on a smaller
% scale. Shaking that occurs regardless of ground motion. Think of it, as
% you are standing on a platform, and that platform shakes pretty fast,
% enough that you want to bend your knees and brace, while the plateform
% moves all together, left and right, the large spike is that left and
% right motion.


% Question 2 ii

% Apply a taper of your choosing and plot the spectra with and without 
% the taper. What, if any, are the differences?

% write up

% Applying the Hann window taper looks like a bandpass filter while 
% looking in the time domain. Paying attention to the end and larger 
% peak magnitudes are relativley preserved compared to the samller
% frequencies, especially near the end. It would seem that the hann 
% window has performed as a low pass filter.
% In the ferquency domain, the main peak in the acceleration amplitude
% spectrum is shortened as well in the tapered acceleration. Off the main
% lobe, the lower frequencies are slightly decreased as well.



