%


% Question 1 i 

% Create a vector of 300 time samples with a sample interval of 0.01 seconds and use this to
% calculate a sine wave with a period of 1 second. Calculate the FFT
% of the sine wave and plot its amplitude against frequency. Use the xlim([0 5])
% command to zoom in on frequencies between 0 and 5 Hz. Use this as a "base map" onto which
% to plot other spectra you calculate.

% Define parameters
Q1_sample_interval = 0.01; % 0.01 seconds per sample
Q1_total_samples = 300; % 300 samples
Q1_frequency = 1; % Frequency of the sine wave in Hz
Q1_amplitude = 1;

% Create a vector of time samples
Q1_time = 0:Q1_sample_interval:(Q1_total_samples-1)*Q1_sample_interval;
% from 0 to (300-1) * 0.01 with an inteval of whatever the sample
% interval is

% Calculate the sine wave
Q1_sine_wave = Q1_amplitude * sin(2 * pi * Q1_frequency * Q1_time);

% Plot the sine wave
figure();
subplot(2,1,1);
plot(Q1_time, Q1_sine_wave);
title('Q1 - Sine Wave');
xlabel('Time (seconds)');
ylabel('Amplitude');
grid on;


% FFT of sine wave
Q1_Transform_SineWave = fft(Q1_sine_wave);
% Amplitude spectrum or Abs value of FT
Q1_Amp_Spec_SineWave = abs(Q1_Transform_SineWave);

% Frequency dimension
Q1_sampling_rate = 1 / Q1_sample_interval; % 1/0.01s = 100 Hz
Q1_frequencies = (0:(Q1_total_samples-1)) * (Q1_sampling_rate / Q1_total_samples);
% where (0:(total_samples-1)) is a vector from start to end. 
% and :  * (sampling_rate / total_samples) , is scaling the vector by the
% ratio of sampling_rate/total_samples , this is 100Hz/300 samples , which
% is essentially 1Hz/3samples

subplot(2,1,2);
plot(Q1_frequencies,Q1_Amp_Spec_SineWave,'DisplayName','300 Samples');
title('Q1 - Amplitude Spectrum');
xlabel('Frequency');
ylabel('Amplitude');
xlim([0 5]); % zooming in on frequencies between 0Hz and 5Hz
grid on;

%%

% Question 1ii

% Calculate and plot the amplitude of the FFT of the first 275
% samples of the sine wave. Explain why it differs from the 
% amplitude spectrum of the 300-sample sine wave, perhaps with
% the help of a sketch in the time domain. 
% (Note: Be sure to look at the width as well as the height 
% of the 1 Hz peak)

Q2_sample_interval = 0.01; % 0.01 seconds per sample
Q2_total_samples = 300; % 300 samples
Q2_frequency = 1; % Frequency of the sine wave in Hz
Q2_amplitude = 1;

% Create a vector of time samples
Q2_time = 0:Q2_sample_interval:(Q2_total_samples-1)*Q2_sample_interval;

Q2_sine_wave = Q2_amplitude * sin(2 * pi * Q2_frequency * Q2_time);
Q2_FFT_SineWave = fft(Q2_sine_wave);
Q2_Amp_Spec_SineWave = abs(Q2_FFT_SineWave);

% FFT of the first 275 samples
Q2_subset_sine_wave = Q2_sine_wave(1:275);
Q2_Transform_Subset_SineWave = fft(Q2_subset_sine_wave);
Q2_Amp_Spec_Subset_SineWave = abs(Q2_Transform_Subset_SineWave);

% Frequency dimension for the subset
Q2_sampling_rate = 1 / Q2_sample_interval;
Q2_frequencies_subset = (0:(length(Q2_subset_sine_wave)-1)) * (Q2_sampling_rate / length(Q2_subset_sine_wave));

Q2_sampling_rate = 1 / Q2_sample_interval; % 1/0.01s = 100 Hz
Q2_frequencies = (0:(Q2_total_samples-1)) * (Q2_sampling_rate / Q2_total_samples);

% Hold on to overlay the plots
figure();
plot(Q2_frequencies,Q2_Amp_Spec_SineWave,'DisplayName','300 Samples');
title('Q2 - Amplitude Spectrum');
xlabel('Frequency');
ylabel('Amplitude');
xlim([0 5]); % zooming in on frequencies between 0Hz and 5Hz
grid on;
hold on;

% Plotting Amplitude Spectrum of the first 275 samples
plot(Q2_frequencies_subset, Q2_Amp_Spec_Subset_SineWave, 'DisplayName','275 Samples');
title('Q2 - Amplitude Spectrum (With Subset)');
xlabel('Frequency (Hz)');
ylabel('Amplitude');
xlim([0 5]); % Zoom in on frequencies between 0Hz and 5Hz
legend();

% Release hold
hold off;


%%


% Question 1 iii


% Pad the 300-sample sine wave with 1748 zeros and plot the
% amplitude of the FFT against frequency on the plot. 
% Can you explain what you see? – It may help to 
% look at the amplitude spectrum of a boxcar function with 
% 2048 samples of which the first 300 are set to one and the 
% remainder to zero

Q3_sample_interval = 0.01; % 0.01 seconds per sample
Q3_total_samples = 300; % 300 samples
Q3_frequency = 1; % Frequency of the sine wave in Hz
Q3_amplitude = 1;
Q3_time = 0:Q3_sample_interval:(Q3_total_samples-1)*Q3_sample_interval;

%Sine wave and amplitude spec
Q3_sine_wave = Q3_amplitude * sin(2 * pi * Q3_frequency * Q3_time);
Q3_FFT_SineWave = fft(Q3_sine_wave);
Q3_Amp_Spec_SineWave = abs(Q3_FFT_SineWave);


% Pad the 300-sample sine wave with 1748 zeros
Q3_padded_sine_wave = [Q3_sine_wave, zeros(1, 1748)];

% FFT of the padded sine wave
Q3_Transform_Padded_SineWave = fft(Q3_padded_sine_wave);
Q3_Amp_Spec_Padded_SineWave = abs(Q3_Transform_Padded_SineWave);

Q3_sampling_rate = 1 / Q3_sample_interval; % 1/0.01s = 100 Hz

Q3_frequencies = (0:(Q3_total_samples-1)) * (Q3_sampling_rate / Q3_total_samples);

% Frequency dimension for the padded signal
Q3_total_samples_padded = length(Q3_padded_sine_wave);
Q3_frequencies_padded = (0:(Q3_total_samples_padded-1)) * (Q3_sampling_rate / Q3_total_samples_padded);

% Origional 300 Sample 

subplot(2,1,1);
plot(Q3_frequencies,Q3_Amp_Spec_SineWave,'DisplayName','300 Samples');
title('Amplitude Spectrum');
xlabel('Frequency');
ylabel('Amplitude');
xlim([0 5]); % zooming in on frequencies between 0Hz and 5Hz
grid on;

hold on;

% Plotting Amplitude Spectrum of the padded sine wave
plot(Q3_frequencies_padded, Q3_Amp_Spec_Padded_SineWave, 'DisplayName', 'Padded Sine Wave');
title('Q3 - Amplitude Spectrum (Padded)');
xlabel('Frequency (Hz)');
ylabel('Amplitude');
grid on;
legend;

hold off;

% Boxcar Parameters
Q3_total_samples = 2048;
Q3_boxcar_length = 300;

% Generate the boxcar function
Q3_boxcar = [ones(1, Q3_boxcar_length), zeros(1, Q3_total_samples - Q3_boxcar_length)];

% Compute the discrete Fourier transform (DFT)
Q3_FT_boxcar = fft(Q3_boxcar);

% Compute the amplitude spectrum
Q3_Boxcar_amp_spec = abs(Q3_FT_boxcar) / Q3_total_samples;

% Frequency axis
Q3_sampling_rate = 1; % Unit sampling rate (arbitrary for this example)
Q3_boxcar_frequencies = (0:Q3_total_samples-1) * Q3_sampling_rate / Q3_total_samples;

% Plot the amplitude spectrum
subplot(2,1,2);
plot(Q3_boxcar_frequencies, Q3_Boxcar_amp_spec);
title('Q3 - Amplitude Spectrum of Boxcar Function');
xlabel('Frequency');
ylabel('Amplitude');
xlim([0, 0.5]); % Limit the x-axis to the positive frequencies (Nyquist limit)
grid on;



%%


% Question 1 iv



% Hanning taper is 1 cycle of a cosine function to which 
% 1 has been added before dividing by two. Calculate a 
% 300-sample Hanning taper with the command scipy.signal. 
% hanning(300). In the time domain, plot the sine wave, 
% the Hanning taper and their product.


% Define parameters
Q4_sample_interval = 0.01;
Q4_total_samples = 300;
Q4_frequency = 1; 
Q4_amplitude = 1;

Q4_time = 0:Q4_sample_interval:(Q4_total_samples-1)*Q4_sample_interval;

% The Sine wave
Q4_sine_wave = Q4_amplitude * sin(2 * pi * Q4_frequency * Q4_time);

% This one actually works for some reason while the hanning function
% produces a weird looking product. I just recreated the Hanning window
% explicitly. This is what  I use for the rest of the Homework.

% Hanning Window
Q4_hanning_freq = 1/3;
Q4_Hanning_win = (1-cos(2*pi*Q4_hanning_freq*Q4_time)) ./2;

% Plot the sine wave
figure();
subplot(2,1,1);
plot(Q4_time, Q4_sine_wave,'DisplayName','Sine wave');
title('Q4 - Sine Wave');
xlabel(' ');
ylabel('Amplitude');
grid on;
hold on;

% Hanning Window
subplot(2,1,1);
plot(Q4_time,Q4_Hanning_win,'DisplayName','Explicit Han Window');
title('Q4 - Explicit Hanning Window on Sinewave');
xlabel('Time');
ylabel('Amplitude');
legend();
grid on;
hold off;

Q4_Hanning_Sine = Q4_sine_wave .* Q4_Hanning_win;

subplot(2,1,2);
plot(Q4_time,Q4_Hanning_Sine);
title('Q4 - Tapered Sine using explicit Han window');
xlabel('Time (seconds)');
ylabel('amplitude');
grid on;



%%


% Question 1 V

% Define parameters
Q5_sample_interval = 0.01;
Q5_total_samples = 300;
Q5_frequency = 1; 
Q5_amplitude = 1;
Q5_time = 0:Q5_sample_interval:(Q5_total_samples-1)*Q5_sample_interval;

Q5_sampling_rate = 1 / Q5_sample_interval; % 1/0.01s = 100 Hz
Q5_frequencies = (0:(Q5_total_samples-1)) * (Q5_sampling_rate / Q5_total_samples);

% Sine wave
Q5_sine_wave = Q5_amplitude * sin(2 * pi * Q5_frequency * Q5_time);

% FFT of sine wave
Q5_Transform_SineWave = fft(Q5_sine_wave);

% Amplitude spectrum or Abs value of FT
Q5_Amp_Spec_SineWave = abs(Q5_Transform_SineWave);


% Making the Tapered Sine wave 

% Hanning Window
Q5_hanning_freq = 1/3;
Q5_Hanning_win = (1-cos(2*pi*Q5_hanning_freq*Q5_time)) ./2;

% Tapered Sine wave
Q5_Hanning_Sine = Q5_sine_wave .* Q5_Hanning_win;
Q5_FT_Hanning_Sine = fft(Q5_Hanning_Sine);
Q5_Amp_Spec_Hanning_Sine = abs(Q5_FT_Hanning_Sine);

% Plot the amplitude spectra of the tapered 300-sample sine wave. 
% How does it compare with the spectrum of the untapered sine wave?

% Plotting top and bottom
figure();
subplot(3,1,1);
plot(Q5_frequencies,Q5_Amp_Spec_SineWave);
title('Q5 - Amplitude Spectrum : Sinewave');
xlabel('frequency');
ylabel('Amplitude');
xlim([0 5])
grid on;

subplot(3,1,2);
plot(Q5_frequencies,Q5_Amp_Spec_Hanning_Sine);
title('Q5 - Amplitude Spectrum : Tapered Sine')
xlabel('frequency');
ylabel('Amplitude');
xlim([0 5])
grid on;


% One on the other
subplot(3,1,3);
plot(Q5_frequencies,Q5_Amp_Spec_SineWave,'DisplayName','Sinewave');
title('Q5 - Amplitude Spectrum : Sinewave');
xlabel('frequency');
ylabel('Amplitude');
xlim([0 5])
grid on;
hold on;

plot(Q5_frequencies,Q5_Amp_Spec_Hanning_Sine,'Display','Tapered Sine');
title('Q5 - Amplitude Spectrum')
xlabel('frequency');
ylabel('Amplitude');
xlim([0 5])
grid on;
legend();
hold off;


%%
 
% Question 1 vi

% Amplitude spectrum of Taper
% Amp Spec of Sine wave is already done earlier
% Define parameters
Q6_sample_interval = 0.01;
Q6_total_samples = 300;
Q6_sampling_rate = 1 / Q6_sample_interval;
Q6_time = 0:Q6_sample_interval:(Q6_total_samples-1)*Q6_sample_interval;
Q6_frequencies = (0:(Q6_total_samples-1)) * (Q6_sampling_rate / Q6_total_samples);

% Hanning Window
Q6_hanning_freq = 1/3;
Q6_Hanning_win = (1-cos(2*pi*Q6_hanning_freq*Q6_time)) ./2;

% amplitude spectrum calculation
Q6_FT_Han_win = fft(Q6_Hanning_win);
Q6_Amp_Spec_Han_Win = abs(Q6_FT_Han_win);


% Hanning Window
figure();
plot(Q6_frequencies,Q6_Amp_Spec_Han_Win,'DisplayName','Han Window');
title('Q6 - Amplitude Spectrum : Hann Window');
xlabel('Frequency');
ylabel('Amplitude');
legend();
xlim([0 5]);
grid on;
hold off;


%%


% Question vii



% Plot the amplitude spectra of first 275 samples of sine 
% wave after tapering with a 275-sample Hanning window. 
% How does it compare with the spectrum of the untapered 275-sample
% sine wave?

% Define parameters
Q7_sample_interval = 0.01;
Q7_total_samples = 300;
Q7_frequency = 1; 
Q7_amplitude = 1;
Q7_time = 0:Q7_sample_interval:(Q7_total_samples-1)*Q7_sample_interval;
Q7_sampling_rate = 1 / Q7_sample_interval; % 1/0.01s = 100 Hz

%sine wave
Q7_sine_wave = Q7_amplitude * sin(2 * pi * Q7_frequency * Q7_time);

Q7_subset_sine_wave = Q7_sine_wave(1:275);

%Q7_frequencies = (0:(Q7_total_samples-1)) * (Q7_sampling_rate / Q7_total_samples);
Q7_frequencies_subset = (0:(length(Q7_subset_sine_wave)-1)) * (Q7_sampling_rate / length(Q7_subset_sine_wave));

% Hanning Window
Q7_hanning_freq = 1/3;
Q7_Hanning_win = (1-cos(2*pi*Q7_hanning_freq*Q7_time)) ./2;

Q7_Subset_Han_Win = Q7_Hanning_win(1:275);

% Tapered Sine

Q7_Sub_Tapered_Sine = Q7_subset_sine_wave .* Q7_Subset_Han_Win ;

% Creating the Amplitude Spectrum of the Tapered Sine
Q7_FT2_Tapered_Sine = fft(Q7_Sub_Tapered_Sine);
Q7_AS_Tapered_Sine = abs(Q7_FT2_Tapered_Sine);

%Creating the Amplitude spectrum of the Subset Sine
Q7_FT_Subset_SineWave = fft(Q7_subset_sine_wave);
Q7_Amp_Spec_Subset_SineWave = abs(Q7_FT_Subset_SineWave);



%plotting Amp Specturm of Tapered sine with 275 samples
figure();
plot(Q7_frequencies_subset,Q7_AS_Tapered_Sine,'DisplayName','Tapered Sine');
xlim([0 5]);
hold on;

% Plotting Amplitude Spectrum of the first 275 samples
plot(Q7_frequencies_subset, Q7_Amp_Spec_Subset_SineWave, 'DisplayName','Untapered Sine');
title('Q7 - Amplitude Spectrum : Subset 275 samples');
xlabel('Frequency (Hz)');
ylabel('Amplitude');
xlim([0 5]); % Zoom in on frequencies between 0Hz and 5Hz
legend();
hold off;

%%


% Question 1 viii



% Finally, taper the 300-sample sine wave with the 300-sample 
% Hanning window, pad with 1748 zeros and plot the amplitude 
% spectra. What is the effect of padding with zeros?

% Parametes
Q8_sample_interval = 0.01; 
Q8_total_samples = 300; 
Q8_frequency = 1; 
Q8_amplitude = 1;
Q8_sampling_rate = 1 / Q8_sample_interval;
Q8_time = 0:Q8_sample_interval:(Q8_total_samples-1)*Q8_sample_interval;
Q8_frequencies = (0:(Q8_total_samples-1)) * (Q8_sampling_rate / Q8_total_samples);

% Sine wave
Q8_sine_wave = Q8_amplitude * sin(2 * pi * Q8_frequency * Q8_time);
% hanning window
Q8_hanning_freq = 1/3;
Q8_Taper = (1-cos(2*pi*Q8_hanning_freq*Q8_time)) ./2;


% Unpadded Tapered Sine
Another_Tapered_Sine = Q8_sine_wave .* Q8_Taper;

% Padded Taper
Another_Padded_Taper = [Another_Tapered_Sine,  zeros(1,1748)];

% Amp Spec of Unpadded Tapered Sine
Q8_Amp_Spec_Another_Tapered_Sine = abs(fft(Another_Tapered_Sine));



% Padding the Tapered Sine
A_Very_Long_Name_For_Another_Tapered_Sine_With_Zeroes_This_Time = [Another_Tapered_Sine, zeros(1,1748)];

% Amplitude Spectrum of Padded Tapered Sine
FFT_AVLNFATSWZTT = fft(A_Very_Long_Name_For_Another_Tapered_Sine_With_Zeroes_This_Time);
Amp_Spec_AVLNFATSWZTT = abs(FFT_AVLNFATSWZTT);


% Padded Time domain
Q8_Total_samples_padded = length(Another_Padded_Taper);
Q8_time_padded = 0:Q8_sample_interval:(Q8_Total_samples_padded-1)*Q8_sample_interval;

% Frequency dimension for the padded signal
Q8_total_samples_padded = length(A_Very_Long_Name_For_Another_Tapered_Sine_With_Zeroes_This_Time);
Q8_frequencies_padded = (0:(Q8_total_samples_padded-1)) * (1 / (Q8_sample_interval * Q8_total_samples_padded));



% Time Domain of both padded and unpadded sines
%Tapered
figure();
subplot(2,1,1);
plot(Q8_time,Another_Tapered_Sine,'DisplayName','Tapered Sine')
title('Q8 - Amp. Spec. of Another Time Tapered Sine. 1-5 Hz')
%xlabel('Time (seconds)');
%ylabel('Magnitude');
xlim([0 5]);
grid on;

% Padded and tapered
subplot(2,1,2);
plot(Q8_time_padded, A_Very_Long_Name_For_Another_Tapered_Sine_With_Zeroes_This_Time,'DisplayName','Padded Tapered Sine');
title('Q8 - Another Tapered Sine wave, this time padded with zeros. 1-5Hz');
xlabel('Time (seconds)');
ylabel('Magnitude');
xlim([0 5]);
grid on;

% Amplitude Spectrum of Both padded and unpadded tapered sines
%unpadded
figure();
subplot(3,1,1);
plot(Q8_frequencies,Q8_Amp_Spec_Another_Tapered_Sine,'DisplayName','Tapered Sine')
title('Q8 - Amp. Spec. of Another Tapered Sine. Whole Window')
%xlabel('Time (seconds)');
%ylabel('Amplitude');
grid on;

%Padded
subplot(3,1,2);
plot(Q8_frequencies_padded, Amp_Spec_AVLNFATSWZTT,'DisplayName','Padded Tapered Sine');
title('Q8 - Amp. spec. Tapered Sine wave, padded zeros. Whole Window');
%('Frequency (Hz)');
ylabel('Amplitude');
grid on;

% Closer look at both
subplot(3,1,3);
plot(Q8_frequencies,Q8_Amp_Spec_Another_Tapered_Sine,'DisplayName','Tapered Sine')
%xlabel('Frequency (Hz)');
%ylabel('Amplitude');
xlim([0 2]);
grid on;
hold on;

%Padded
plot(Q8_frequencies_padded, Amp_Spec_AVLNFATSWZTT,'DisplayName','Padded Tapered Sine');
title('Q8 - Closer look at both, together');
xlabel('Frequency (Hz)');
%ylabel('Amplitude');
xlim([0 2]);
grid on;
hold on;


%%


% Question 1 i 

% Create a vector of 300 time samples with a sample interval of 0.01 
% seconds and use this to calculate a sine wave with a period of 1 second. 
% Calculate the FFT (see notes below for useful links) of the sine wave 
% and plot its amplitude against frequency. Use the xlim([0 5])
% command to zoom in on frequencies between 0 and 5 Hz. Use this as a 
% "base map" onto which to plot other spectra you calculate

% Write up

% n / A


% Question 1 ii

% Calculate and plot the amplitude of the FFT of the first 275 samples 
% of the sine wave. Explain why it differs from the amplitude spectrum of 
% the 300-sample sine wave, perhaps with the help of a sketch in the time 
% domain. (Note: Be sure to look at the width as well as the height
% of the 1 Hz peak)

% Write up

% As we have seen before, there is an inverse relationship between the size 
% of the window in the time domain and the width of the amplitude spectum.
% So, when we compare the amplitude spectum between the 300 samples vs 275
% samples, the spectrum of 300 samples makes a more narrow spectrum, and
% the 275 samples makes a more wide spectrum. This follows the uncertaintly
% principle, the closer we look, the less accurate we can really see. The
% smaller the window, the less accurate or certain we can be of the
% frequencies that make it up.


% Question 1 iii

% Pad the 300-sample sine wave with 1748 zeros and plot the amplitude of 
% the FFT against frequency on the plot. Can you explain what you see? 
% – It may help to look at the amplitude spectrum of a boxcar function 
% with 2048 samples of which the first 300 are set to one and the remainder 
% to zero.

% Write Up

% Padding with zeros makes the window a bit wider and curved, also not
% reachin all the way down to 0 amplitude on the main lobe. Now many side
% lobes are present, decaying towards zero with increasing frequency. This
% is because we have increased are window in the time domain by a huge
% amount, even though all zeros, zero is still data. Frequency data arises
% in places where we are fairly certain we shouldn't see ferquency.


% Question 1 iv

% Hanning taper is 1 cycle of a cosine function to which 1 has been added 
% before dividing by two. Calculate a 300-sample Hanning taper with the
% command scipy.signal.hanning(300). In the time domain, plot the sine 
% wave, the Hanning taper and their product.

% write up

% n/a


% Question 1 v

% Plot the amplitude spectra of the tapered 300-sample sine wave. 
% How does it compare with the spectrum of the untapered sine wave?

% write up

% The tapered sine wave is now shorter in amplitude as well as a much wider
% range of frequencies. This makes sense , because as we saw in question 4,
% the tapered sine wave has a less uniform wave, comprised of longer and
% shorter frequency wavelengths, that comprise smaller weights of the whole
% function than the fewer, if only one frequency that comprises the
% untapered sine wave.


% Question 1 vi

% Plot the amplitude spectra of the Hanning window and use it to explain 
% what you see in part 5. It is also interesting to use the semilogy 
% command to make the y-axis logarithmic.

% write up

% The hann window's amplitude spectrum is symmetric about 0, but ignoring
% the negative amplitudes, we can see that the amplitude spectrum is a
% fairly shapr spike, from 150, down to 0 over a range of frequencies 0 to
% approximately 0.667 Hz.


% Question 1 vii

% Plot the amplitude spectra of first 275 samples of sine wave after 
% tapering with a 275-sample Hanning window. How does it compare with the 
% spectrum of the untapered 275-sample sine wave?

% write up

% Both are less narrow because they are the spectrum over a smaller time
% interval, so their spikes will be more wide. The Taper is shorter, as we
% saw in the last question and the amplitude at either corner frequency is
% smaller than that of the untapered sine.


% Question 1 viii

% Finally, taper the 300-sample sine wave with the 300-sample Hanning 
% window, pad with 1748 zeros and plot the amplitude spectra. 
% What is the effect of padding with zeros?

% write up

% When we pad with zeros, in the time domain, we get a main lobe and side 
% lobes that decay to zero as frequency increases, looking at the whole
% spectrum, we see similar spikes near zero and near 100, and they are 
% identicle between the unpadded and padded tapered function. However a
% closer inspection shoes us they are not excatly the same. Both padded and
% unpadded functions contain nearly the same wavelenths, except the padded
% function has a larger window, creating some increased ambiguity in
% exactly what frequencies gather amplitude in the frequency domain.

