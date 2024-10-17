# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 20:03:08 2024

@author: jettr
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import butter , filtfilt , freqz , lfilter

#%%


# Question 1

# Define the time vector
Q1time = np.arange(128)

# Delta Function
Q1Delta_Func = np.zeros(128)
Q1Delta_Func[0] = 1

# Plot Delta Function
plt.figure(figsize = (8,8))
plt.plot(Q1time, Q1Delta_Func)
plt.xlabel('Time (Seconds)')
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

fig = plt.figure(figsize = (8,8))

# Loop Filter for different orders

for i, Q2order in enumerate(Q2orders, 1):
    
    Q2b,Q2a = butter(Q2order, Q2nyquist_freq, 'low') # Low Pass Butter
    
    Q2filtered_signal = lfilter(Q2b,Q2a, Q2Delta_Func)
    
    #plt.subplot(len(Q2orders),1,i)
    plt.plot(Q2time,Q2filtered_signal, label = f'Order - {Q2order}')
    plt.xlim([0 , 20])
    
fig.text(0.5,0.001,'Time (Seconds)', ha = 'center' , rotation = 'horizontal')
fig.text(0.001, 0.5, 'Magnitude', ha='center', va='center', rotation='vertical')
plt.suptitle('Q1a Butterworth Filter - Orders - 2,4,8,16' )
plt.legend()
plt.tight_layout()
plt.show()


#%%

# Question 1 b

# Plot the Amplitude and Phase spectra of the delta function filtered with a butterworth filter

# Generate delta function time series
Q3_time = np.arange(128)  # Time vector
Q3_DeltaFun = np.zeros(128)  # Delta function at time 0
Q3_DeltaFun[0] = 1

# Define Nyquist frequency
Q3_Nyquist_freq = 0.5

# Orders of the Butterworth filters
Q3_orders = [2, 4, 8, 16]

# Phase Respone figure
fig = plt.figure(figsize = (8,8))

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
fig.text(0.001, 0.5, 'Phase', ha='center', va='center', rotation='vertical')
fig.suptitle('Q 1b Phase Response')
plt.legend()
plt.tight_layout()
    

# Second Figure for the Amplitude Response

fig2 = plt.figure(figsize = (8,8))
    
for i, Q3_order in enumerate(Q3_orders, 1):
    # Create Butterworth filter
    b, a = butter(Q3_order, Q3_Nyquist_freq, 'low')
    
    # Frequency response of the filter
    w, H = freqz(b, a, worN=Q3_time.size)
    
    # Plot amplitude response
    plt.plot(w/np.pi, np.abs(H), label = f'Order - {Q3_order}')
    plt.ylim([0, 1.2])  # Limit y-axis to [0, 1.2] for better visualization

# Adjust figure
fig2.text(0.5,0.001,'Frequency', ha = 'center' , rotation = 'horizontal')
fig2.text(0.001, 0.5, 'Amplitude', ha='center', va='center', rotation='vertical')
fig2.suptitle('Q 1b Amplitude Response')
plt.legend()
plt.tight_layout()
plt.show()


#%%
# =============================================================================
# 
#  Question 1 c 
# 
# The Phase Response, in the normalized frequency domain, from low to high order, are smooth, then become sharp and have
# 'teeth' . Order 2 has no teeth, smooth from 0 phase to -pi phase as freq 1.0 . Order two then has one
# tooth, at about 0.5 freq. This graph starts at 0, decays to -pi, concave down , spikes up to +pi , 
# then returns back to 0, concave up at the end. Order 8 acts similarly, with twice as many spikes , 
# then order 16 has twice as many spikes as order 8, however the end of the graph is very disorderly.
# 
# 
# The Amplitude response, in the frequency domain, all look about like the corner of a box. They start at 
# amplitude 1 at time 0, stay contant, then drop to 0 around time 0.4. The time at which they begin to
# drop gets closer to 0.5 at higher orders. The lower orders are much smoother, more like hills.
# They begin their descend much earlier, like order 2 , starting to decay at 0.2 time.
# 
# 
# =============================================================================
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
fig4 = plt.figure(figsize = (8,8))

# Apply filters and plot results
for i, Q4_order in enumerate(Q4_orders, 1):
    # Create Butterworth filter
    b, a = butter(Q4_order, Q4_Nyquist_freq, 'high')
    
    # Apply filter
    Q2filtered_signal = lfilter(b, a, Q4_DeltaFun)
    
    # Plot filtered signal
    #plt.subplot(len(Q4_orders), 1, i)
    plt.plot(Q4_time, Q2filtered_signal , label = f'Order - {Q4_order}')
    plt.xlim(0 , 30)

# Adjust figure
fig4.suptitle('Q 1d High-Pass Butterworth/Delta Function : Order 2 - 16')
fig4.text(0.5,0.001,'Time (Seconds)', ha = 'center' , rotation = 'horizontal')
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
Q5Wn = [0.4, 0.6]

# Orders of the Butterworth filters
Q5orders = [2, 4, 8, 16]

# Initialize figure
fig5 = plt.figure(figsize = (8 , 8))

# Apply filters and plot results
for i, Q5order in enumerate(Q5orders, 1):
    
    # Create notch filter
    b, a = butter(Q5order, Q5Wn, 'bandstop')
    
    # Apply filter
    Q5filtered_signal = lfilter(b, a, Q5_DeltaFun)
    
    # Plot filtered signal
    plt.plot(Q5_time, Q5filtered_signal , label = f'Order {Q5order}')
    plt.xlim(0,30)
    
    
    # Adjust figure
fig5.suptitle('Q 1e Notch Butterworth/Delta Function')
fig5.text(0.5,0.001,'Time (seconds)', ha = 'center' , rotation = 'horizontal')
fig5.text(0.001, 0.5, 'Magnitude', ha='center', va='center', rotation='vertical')
plt.tight_layout()
plt.legend()
plt.show()

fig52 = plt.figure(figsize = (18,5))

# Apply filters and plot results
for i, Q5order in enumerate(Q5orders, 1):
    
    # Create notch filter
    b, a = butter(Q5order, Q5Wn, 'bandstop')
    
    # Apply filter
    Q5filtered_signal = lfilter(b, a, Q5_DeltaFun)
    
    # Plot filtered signal
    plt.subplot(1, len(Q5orders), i)
    plt.plot(Q5_time, Q5filtered_signal)
    plt.title(f'Order {Q5order}')
    plt.xlim(0,30)

    # Adjust figure
fig52.suptitle('Q 1e Notch Butterworth/Delta Function')
fig52.text(0.5,0.001,'Time (seconds)', ha = 'center' , rotation = 'horizontal')
fig52.text(0.001, 0.5, 'Magnitude', ha='center', va='center', rotation='vertical')
plt.tight_layout()
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.25, hspace=None)
plt.show()



#%%


# Question 1 f


# Generate delta function time series
Q6_time = np.arange(-63, 65)  # Time vector
Q6_DeltaFun = np.zeros(Q6_time.shape)  # Delta function at time 0
Q6_DeltaFun[63] = 1  # Delta function at t = 0

# Nyquist frequency
Nyquist_freq = 0.5

# Orders of the Butterworth filters
Q6orders = [2, 4, 8, 16]

# Initialize figure
fig6 = plt.figure(figsize = (12,5))

# Apply filters and plot results
for i, Q6order in enumerate(Q6orders, 1):
    
    # Create Butterworth filter
    b, a = butter(Q6order, Nyquist_freq, 'low')
    
    # Apply zero-phase filter
    Q6filtered_signal = filtfilt(b, a, Q6_DeltaFun)
    
    # Plot filtered signal
    plt.plot(Q6_time, Q6filtered_signal, label = f'Order {Q6order}')
    plt.xlim(-20 , 20)
   

# Adjust figure
fig6.suptitle('Q 1f Zero-Phase Low-Pass Butterworth')
fig6.text(0.5,0.001,'Time (seconds)', ha = 'center' , rotation = 'horizontal')
fig6.text(0.001, 0.5, 'Magnitude', ha='center', va='center', rotation='vertical')
plt.legend()
plt.tight_layout()
plt.show()

