# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 20:03:08 2024

@author: jettr
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import butter , filtfilt , freqz

#%%

# Question 1

# Define the time vector
Q1time = np.arange(128)

# Delta Function
Q1Delta_Func = np.zeros(128)
Q1Delta_Func[0] = 1

# Plot Delta Function
plt.plot(Q1time, Q1Delta_Func)
plt.xlabel('Time')
plt.ylabel('Amplitude')
plt.title(' Q1 Delta Function at t=0')
plt.show()

#%%

# Question 1 a

Q2time = np.arange(128)
Q2Delta_Func = np.zeros(128)
Q2Delta_Func[0] = 1

# Nyquist Frequency
Q2nyquist_freq = 0.5

# Order of the butterworth Filter
Q2orders = [2,4,8,16]

fig = plt.figure()

# Loop Filter for different orders

for i, Q2order in enumerate(Q2orders, 1):
    
    Q2b,Q2a = butter(Q2order, Q2nyquist_freq, 'low') # Low Pass Butter
    
    Q2filtered_signal = filtfilt(Q2b,Q2a, Q2Delta_Func)
    
    #plt.subplot(len(Q2orders),1,i)
    plt.plot(Q2time,Q2filtered_signal, label = f'Order - {Q2order}')
    plt.xlim([0 , 20])
    
fig.text(0.5,0.001,'Time', ha = 'center' , rotation = 'horizontal')
fig.text(0.001, 0.5, 'Magnitude', ha='center', va='center', rotation='vertical')
plt.suptitle('Q1a Butterworth Filter - Orders : (2,4,8,16)' )
plt.legend()
plt.tight_layout()
plt.show()

#%%

# Question 1 b

# Generate delta function time series
Q3_time = np.arange(128)  # Time vector
Q3_DeltaFun = np.zeros(128)  # Delta function at time 0
Q3_DeltaFun[0] = 1

# Define Nyquist frequency
Q3_Nyquist_freq = 0.5

# Orders of the Butterworth filters
Q3_orders = [2, 4, 8, 16]

# Ampllitude Respone figure
fig = plt.figure()

# Apply filters and plot results
for i, Q3_order in enumerate(Q3_orders, 1):
    # Create Butterworth filter
    b, a = butter(Q3_order, Q3_Nyquist_freq, 'low')
    
    # Frequency response of the filter
    w, H = freqz(b, a, worN=Q3_time.size)
    
    # Plot phase response
    #plt.subplot(len(Q3_orders), 1, i)
    plt.plot(w/np.pi, np.angle(H) , label = f'Order - {Q3_order}')
    
fig.text(0.5,0.001,'Normalized Frequency (radians/sample)', ha = 'center' , rotation = 'horizontal')
fig.text(0.001, 0.5, 'Amplitude', ha='center', va='center', rotation='vertical')
fig.suptitle('Q 1b Amplitude Response')
plt.legend()
plt.tight_layout()
    

# Second Figure for the Phase Response

fig2 = plt.figure()
    
for i, Q3_order in enumerate(Q3_orders, 1):
    # Create Butterworth filter
    b, a = butter(Q3_order, Q3_Nyquist_freq, 'low')
    
    # Frequency response of the filter
    w, H = freqz(b, a, worN=Q3_time.size)
    
    # Plot amplitude response
    #plt.subplot(len(Q3_orders), 1, i)
    plt.plot(w/np.pi, np.abs(H), label = f'Order - {Q3_order}')
    plt.ylim([0, 1.2])  # Limit y-axis to [0, 1.2] for better visualization

# Adjust figure
fig2.text(0.5,0.001,'Normalized Frequency (radians/sample)', ha = 'center' , rotation = 'horizontal')
fig2.text(0.001, 0.5, 'Phase (Radians)', ha='center', va='center', rotation='vertical')
fig2.suptitle('Q 1b Phase Response')
plt.legend()
plt.tight_layout()
plt.show()

#%%

# Question 1 d

# Generate delta function time series
Q4_time = np.linspace(0, 127, 128)  # Time vector
Q4_DeltaFun = np.zeros(128)  # Delta function at time 0
Q4_DeltaFun[0] = 1

# Define Nyquist frequency
Q4_Nyquist_freq = 0.8

# Orders of the Butterworth filters
Q4_orders = [2, 4, 8, 16]

# Initialize figure
fig4 = plt.figure()

# Apply filters and plot results
for i, Q4_order in enumerate(Q4_orders, 1):
    # Create Butterworth filter
    b, a = butter(Q4_order, Q4_Nyquist_freq, 'high')
    
    # Apply filter
    filtered_signal = filtfilt(b, a, Q4_DeltaFun)
    
    # Plot filtered signal
    #plt.subplot(len(Q4_orders), 1, i)
    plt.plot(Q4_time, filtered_signal , label = f'Order - {Q4_order}')
    plt.xlim(0 , 30)

# Adjust figure
fig4.suptitle('Q 1d High-Pass Butterworth/Delta Function : Order 2 - 16')
fig4.text(0.5,0.001,'Time', ha = 'center' , rotation = 'horizontal')
fig4.text(0.001, 0.5, 'Magnitude', ha='center', va='center', rotation='vertical')
plt.tight_layout()
plt.legend()
plt.show()


#%%

# Question 1 e


# Generate delta function time series
Q5_time = np.linspace(0, 127, 128)  # Time vector
Q5_DeltaFun = np.zeros(128)  # Delta function at time 0
Q5_DeltaFun[0] = 1

# Define frequency range for the notch filter
Wn = [0.4, 0.6]

# Orders of the Butterworth filters
Q5orders = [2, 4, 8, 16]

# Initialize figure
fig5 = plt.figure()

# Apply filters and plot results
for i, Q5order in enumerate(Q5orders, 1):
    # Create notch filter
    b, a = butter(Q5order, Wn, 'bandstop')
    
    # Apply filter
    filtered_signal = filtfilt(b, a, Q5_DeltaFun)
    
    # Plot filtered signal
    #plt.subplot(len(Q5orders), 1, i)
    plt.plot(Q5_time, filtered_signal , label = f'Order {Q5order}')
    plt.xlim(0,30)

# Adjust figure
fig5.suptitle('Q 1e Notch Butterworth/Delta Function')
fig5.text(0.5,0.001,'Time', ha = 'center' , rotation = 'horizontal')
fig5.text(0.001, 0.5, 'Magnitude', ha='center', va='center', rotation='vertical')
plt.tight_layout()
plt.legend()
plt.show()

#%%


