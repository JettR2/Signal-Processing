%

% Question 1 i

t = linspace(-5, 5, 1000);  % Time vector from -5 to 5 with 1000 data points
tau_1 = 1;  % Should be the same as RMS deviation

% Compute Gaussian distribution
g_w = (1 / sqrt(2 * pi * tau_1^2)) * exp(-(t.^2) / (2 * tau_1^2));


% Plot the Gaussian distribution
figure;
plot(t, g_w, 'LineWidth', 2);
title('Gaussian Distribution');
xlabel('Time');
ylabel('Amplitude');

%% 

% Question 1 ii


t = linspace(-5, 5, 1000);
tau_2 = tau_1;

g_w = (1 / sqrt(2 * pi * tau_2^2)) * exp(-(t.^2) / (2 * tau_2^2));


RMS_Deviation = sqrt((trapz((t.^2) .* g_w))/(trapz(g_w)));

text_printRMS = ['Root mean square:', num2str(RMS_Deviation)];
disp(text_printRMS);

text_RMS_is_tau = ['The value of tau, ',num2str(tau_2),', is equal to ',num2str(RMS_Deviation), ', our calculated RMS deviation of the gaussian!'];
disp(text_RMS_is_tau);

%% 
f = linspace(-5, 5 , 1000); 
tau_3 = tau_1;
omega = 2*pi*f;

% The Fourier Transform of the gaussian is
G_w = exp(-((omega.^2).*(tau_3.^2))./2);

% the calculated Amplitude Spectrum is the same
A_f = G_w;

figure;
plot(f, A_f, 'LineWidth', 2);
title('Amplitude Spectrum A(w)');
xlabel('Frequencey');
ylabel('Amplitude');

%%

% The Time series and frequency series have an inversely
% proportional relationship. When Tau gets large, the 
% gaussian gets more wide and the frequency series get 
% more narrow. 
%
