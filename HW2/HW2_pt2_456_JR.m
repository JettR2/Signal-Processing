
% Question 2 ii


% Plotting the Raw Mauna Loa Weekly CO2 data
% as well as plotting the Fourier Amplitude
% spectrum of the Mauna Loa Data


% Loading in the file from matlab files
filename = 'maunaloa_weekly.csv';
data = readtable(filename, 'ReadVariableNames', false, 'Delimiter', ',');
variableNames = {'Date', 'Values'};

% Assign variable names to the table
data.Properties.VariableNames = variableNames;

% Convert the dates to datetime format
dates = datetime(data.Date, 'Format', 'dd/MM/yyyy');
Values = data.Values;


% Plot the data
figure;
subplot(2, 1, 1);
plot(dates, Values, 'LineWidth', 2);
title('Maunaloa Weekly Data');
xlabel('Time');
ylabel('Floating Point Values');



% Fourier analysis plot
subplot(2, 1, 2);

% Time interval between data points
dt = days(7); % Weekly data

% Perform FFT
N = length(Values);
Fs = 1/days(dt); % Sampling frequency
frequencies = Fs*(0:(N/2))/N;
fft_values = fft(Values);
amplitudes = 2/N * abs(fft_values(1:N/2+1));

% Plot the amplitude spectrum with log scale on x-axis
semilogx(frequencies, amplitudes, 'LineWidth', 2);
title('Fourier Amplitude Spectrum');
xlabel('Frequency (Hz)');
ylabel('Amplitude');


%%

% Question 2 iii


% Plotting The Orgional Data as well as smoothing
% the data using the N-point moving average. 
% Also Plotting the Fourier Amplitude spectrum of the 
% Smoothed data now.


% Load your data
filename = 'maunaloa_weekly.csv';
data = readtable(filename, 'ReadVariableNames', false, 'Delimiter', ',');
variableNames = {'Date', 'Values'};
data.Properties.VariableNames = variableNames;

% Convert the dates to datetime format
dates = datetime(data.Date, 'Format', 'dd/MM/yyyy');
Values = data.Values;


% Plot the original data
figure;
subplot(3, 1, 1);
plot(dates, Values, 'LineWidth', 2);
title('Original Maunaloa Weekly Data');
xlabel('Time');
ylabel('Floating Point Values');

% Perform running mean smoothing
N = 48; % Can change N to smooth more or less
smoothed_values = smoothdata(Values, 'movmean', N);

% Plot the smoothed data
subplot(3, 1, 2);
plot(dates, smoothed_values, 'LineWidth', 2);
title('Smoothed Maunaloa Weekly Data');
xlabel('Time');
ylabel('Floating Point Values');

% Fourier analysis plot for smoothed data
subplot(3, 1, 3);

% Time interval between data points
dt = days(7); % Weekly data

% Perform FFT
N = length(smoothed_values);
Fs = 1/days(dt); % Sampling frequency
frequencies = Fs * (0:(N/2))/N;
fft_values_smoothed = fft(smoothed_values);
amplitudes_smoothed = 2/N * abs(fft_values_smoothed(1:N/2+1));

% Plot the amplitude spectrum with log scale on x-axis
semilogx(frequencies, amplitudes_smoothed, 'LineWidth', 2);
title('Fourier Amplitude Spectrum (Smoothed Data)');
xlabel('Frequency (Hz)');
ylabel('Amplitude');


%%

% Question 2 iv


% Now we are making two separate figures. 

% Both will include the origional data, as well 
% as a filtered version.

% The first figure we plot uses the Gaussian Distribution 
% as the kernel and fine the amplitude spectrum.

% The second, we use a boxcar function, who's have width is 
% defined as equal to tau, the RMS deviation of the gaussian.


% Loading in the data
filename = 'maunaloa_weekly.csv';
data = readtable(filename, 'ReadVariableNames', false, 'Delimiter', ',');
variableNames = {'Date', 'Values'};
data.Properties.VariableNames = variableNames;

% Make sure it's in the right time format
dates = datetime(data.Date, 'Format', 'dd/MM/yyyy');
Values = data.Values;

% Creating a time vector
t = linspace(-5, 5, length(Values));
% Defining Tau
tau = 1;

% Generate the Gaussian
g_w = (1 / sqrt(2 * pi * tau^2)) * exp(-(t.^2) / (2 * tau^2));
g_w = g_w / sum(g_w); % Normalize to sum to unity

% Convolve the original data with the Gaussian
smoothed_values_gaussian = conv(Values, g_w, 'same');

figure();
% Plot the smoothed data using Gaussian convolution
subplot(3, 1, 2);
plot(dates, smoothed_values_gaussian, 'LineWidth', 2);
title('Smoothed Maunaloa Weekly Data (Gaussian Convolution)');
xlabel('Time');
ylabel('CO2 measurements');

% Plot the original data
subplot(3,1,1);
plot(dates, Values, 'LineWidth', 2);
title('Original Maunaloa Weekly Data');
xlabel('Time');
ylabel('CO2 measurments');

% Time interval between data points
dt = days(7); % Weekly data

% Perform FFT
N = length(smoothed_values_gaussian);
Fs = 1/days(dt); % Sampling frequency
frequencies = Fs * (0:(N/2))/N;
fft_values_smoothed_boxcar = fft(smoothed_values_gaussian);
amplitudes_smoothed_boxcar = 2/N * abs(fft_values_smoothed_boxcar(1:N/2+1));

% Plot the amplitude spectrum with log scale on x-axis
% Fourier analysis plot for smoothed data
subplot(3, 1, 3);
semilogx(frequencies, amplitudes_smoothed_boxcar, 'LineWidth', 2);
title('Fourier Amplitude Spectrum (Gaussian Convolution)');
xlabel('Frequency (Hz)');
ylabel('Amplitude');



% Plotting Data, smoothed with boxcar




% tau = boxcar half width = RMS deviation
tau = 0.5;

% boxcar function
boxcar = rectangularPulse(-tau,tau,t);

% convolve the origional data with the Boxcar function
smoothed_values_boxcar = conv(Values, boxcar, 'same');


figure();
% Plot the original data
subplot(3,1,1);
plot(dates, Values, 'LineWidth', 2);
title('Original Maunaloa Weekly Data');
xlabel('Time');
ylabel('CO2 measurments');

% PLotting smoothed data with boxcar
subplot(3, 1, 2);
plot(dates, smoothed_values_boxcar, 'LineWidth', 2);
title('Smoothed Maunaloa Weekly Data (Boxcar)');
xlabel('Time');
ylabel('CO2 Measurments');


% Time interval between data points
dt = days(7); % Weekly data

% Perform FFT
N = length(smoothed_values_boxcar);
Fs = 1/days(dt); % Sampling frequency
frequencies = Fs * (0:(N/2))/N;
fft_values_smoothed_boxcar = fft(smoothed_values_boxcar);
amplitudes_smoothed_boxcar = 2/N * abs(fft_values_smoothed_boxcar(1:N/2+1));

% Plot the amplitude spectrum with log scale on x-axis
% Fourier analysis plot for smoothed data
subplot(3, 1, 3);
semilogx(frequencies, amplitudes_smoothed_boxcar, 'LineWidth', 2);
title('Fourier Amplitude Spectrum (Boxcar Convolution)');
xlabel('Frequency (Hz)');
ylabel('Amplitude');
