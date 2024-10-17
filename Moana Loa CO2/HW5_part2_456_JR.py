# -*- coding: utf-8 -*-
"""
Created on Sun Feb 25 14:26:40 2024

@author: jettr
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft
from scipy.signal import butter , filtfilt , freqz , lfilter

#%%


data = pd.read_csv('C:/Users/Jett R/Dropbox (University of Oregon)/23-24/Winter/Erth 456 Sig Proc/HW5/maunaloa_weekly.csv', header = None, names=['Dates' , 'CO2'])
df = pd.DataFrame(data)

#%%


# Question 2 a 



dt = 7 # One sample every 7 days , once a week

CO2 = df['CO2'] # CO2 Concentrations
CO2_array = np.array(CO2)
N = len(CO2)

dates = np.arange(0,N*dt,dt) # Days since first sample


# Time Series Plotting
plt.figure(figsize = (10 , 5))
plt.subplot(1,2,1)
plt.plot(dates,CO2_array)
plt.xticks(dates[::521],  rotation='vertical') # edit tick marks
plt.margins(0.05)
plt.xlabel('Time Since first sample')
plt.ylabel('CO2 ppm')
plt.title('Keeling Curve')


# Perform Fourier analysis
Sample_Rate = dt
N = len(CO2_array)
Fs = 1 / Sample_Rate
frequencies = np.linspace(0, Fs / 2, N // 2)
FFT_CO2 = fft(CO2_array) # FFT
Amplitudes = (2/N) * np.abs(FFT_CO2[0:N//2])

plt.subplot(1,2,2)
plt.semilogx(frequencies,Amplitudes)
plt.xlabel('Frequency (Hz)')
plt.ylabel('Amplitude')
plt.title('Amplitude Spectrum of Keeling Curve')
plt.grid(True)
plt.ylim(0 , 30)

plt.suptitle('Keeling and Amplitude Spectra' , fontsize = 15)
plt.tight_layout()


#%%

# Question 2 b
 



# Fourth Order High - Pass butterworth Filter 
# Removing multi-decadal trends

# Input function = CO2_array
# Time variable = dates

order = 4

Wn = 0.2 # Critical Frequency. For Digital filter, normalized from 0 to 1, where 
# 1 is the Nyquist frequency and fs is unspecified.

High_b , High_a = butter(order , Wn , 'highpass')

High_filtered_signal = lfilter(High_b,High_a,CO2_array)

# Whole Filtered Signal
plt.figure(figsize = (10 ,5 ))

plt.subplot(1,2,1)
plt.plot(dates , High_filtered_signal)
plt.xlabel('Time (Weeks) ~ 60 Years')
plt.ylabel('CO2')
plt.title('The Whole Biscuit')
plt.grid(True)



# Zoomed in filtered Siganl
plt.subplot(1,2,2)
plt.plot(dates , High_filtered_signal)
plt.xlabel('Days since first record')
plt.ylabel('CO2')
plt.title('A crumb')
plt.xlim(364,728) # 371 days (after 1st measure) to 728 days. 52weeks*samples/yr * 7days/week * sample = 364 days / year
plt.ylim(-10 , 10)
plt.grid(True)

plt.suptitle(f'Q 2 b - High Buttered Keel: Order {order} & Corner freq - {Wn}' , fontsize = 15)
plt.tight_layout()

#%%


# Question 2 c
 



# Fourth Order Low - Pass butterworth Filter 
# leaving multi-decadal trends

# Input function = CO2_array
# Time variable = dates

order = 4

Wn = 0.3 # Critical Frequency. For Digital filter, normalized from 0 to 1, where 
# 1 is the Nyquist frequency and fs is unspecified.

Low_b , Low_a = butter(order , Wn , 'lowpass') # Low Pass this time

Low_filtered_signal = lfilter(Low_b,Low_a,CO2_array)

# Whole Filtered Signal
plt.figure(figsize = (10 ,5 ))

plt.subplot(1,2,1)
plt.plot(dates , Low_filtered_signal)
plt.xlabel('Time (Weeks) ~ 60 Years')
plt.ylabel('CO2')
plt.title('The Whole Biscuit')
plt.grid(True)



# Zoomed in filtered Siganl
plt.subplot(1,2,2)
plt.plot(dates , Low_filtered_signal)
plt.xlabel('Days since first record')
plt.ylabel('CO2')
plt.title('A crumb')
plt.xlim(0,5000) # 371 days (after 1st measure) to 728 days. 52weeks*samples/yr * 7days/week * sample = 364 days / year
#plt.ylim(-10 , 10)
plt.grid(True)

plt.suptitle(f'Q 2 c - Low Buttered Keel: Order {order} & Corner freq - {Wn}' , fontsize = 15)
plt.tight_layout()