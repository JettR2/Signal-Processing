#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import scipy 
from scipy.signal import resample_poly, butter, lfilter,freqs,filtfilt,sosfilt,sosfiltfilt
from scipy.fft import fft

#get_ipython().run_line_magic('matplotlib', 'notebook')


# In[2]:


#read in csv file for Pitcarin buoy
iturup = pd.read_csv(r'C:/Users/Jett R/Dropbox (University of Oregon)/23-24/Winter/Erth 456 Sig Proc/456 Final/Buoys/NE_Tokyo.csv')
#convert month column '3' to '03' format
iturup['MM'] = iturup['MM'].apply(lambda x: '{:02d}'.format(x))

#drop values of height = 9999
iturup = iturup[iturup['HEIGHT'] != 9999]


# In[3]:


#converting columns into single datetime object
# assigning columns
year = '#YY'  
month = 'MM'  
day = 'DD' 
hour = 'hh'  
minute = 'mm' 
second = 'ss' 

# Adding a new datetime column for year, month, day, hour, minute, and second
iturup['Date'] = pd.to_datetime(iturup[year].astype(str) + 
                                       iturup[month].astype(str).str.zfill(2) + 
                                       iturup[day].astype(str).str.zfill(2) +
                                       iturup[hour].astype(str).str.zfill(2) +
                                       iturup[minute].astype(str).str.zfill(2) +
                                       iturup[second].astype(str).str.zfill(2), 
                                       format='%Y%m%d%H%M%S')
iturup.set_index(iturup['Date'], inplace=True)
time = iturup['Date']


# In[4]:


# Extract sample frequency "T" and "height" columns
T = iturup['T']
height = iturup['HEIGHT']

height_avg = np.mean(height)
iturup['height_norm'] = height - height_avg


# In[6]:


#separate data by sampling frequency
iturup1 = iturup[iturup['T'] == 1].copy()
iturup2 = iturup[iturup['T'] == 2].copy()
iturup3 = iturup[iturup['T'] == 3].copy()

#reset the index for each DataFrame
iturup1.reset_index(drop=True, inplace=True)
iturup2.reset_index(drop=True, inplace=True)
iturup3.reset_index(drop=True, inplace=True)

iturup1.set_index(iturup1['Date'], inplace=True)
time1 = iturup1['Date']
iturup1.drop(columns=['Date'], inplace=True)

iturup2.set_index(iturup2['Date'], inplace=True)
time2 = iturup2['Date']
iturup2.drop(columns=['Date'], inplace=True)

iturup3.set_index(iturup3['Date'], inplace=True)
time3 = iturup3['Date']
iturup3.drop(columns=['Date'], inplace=True)

iturup.drop(columns=['Date'], inplace=True)


# In[7]:


#plot data
cmap = ListedColormap(plt.cm.viridis(np.linspace(0, 1, 3)))
plt.figure()
plt.plot(time, iturup['height_norm'], linewidth=1)
plt.axvline(pd.Timestamp('2011-03-11 05:46:24'), c = 'red')
plt.scatter(time, iturup['height_norm'], c = iturup['T'], cmap = cmap, s = 5)
plt.colorbar()
plt.grid()
plt.show()


# In[8]:


# Create subplots
fig, axs = plt.subplots(3, 1, sharex = True, sharey = True)

# Scatter plot for pitcarin1
axs[0].scatter(iturup1.index, iturup1['height_norm'], s=1)
axs[0].set_title('15min sampling')
axs[0].grid()

# Scatter plot for pitcarin2
axs[1].plot(iturup2.index, iturup2['height_norm'])
axs[1].set_title('1min sampling')
axs[1].grid()

# Scatter plot for pitcarin3
axs[2].scatter(iturup3.index, iturup3['height_norm'], s=1)
axs[2].set_title('15s sampling')
axs[2].grid()

# Adjust layout
plt.tight_layout()

# Show plot
plt.show()


# In[9]:


# Resample the 'HEIGHT' column to 1 minute frequency
resamp = iturup['height_norm'].resample('15S').mean()

# Optionally, fill missing values using linear interpolation
resamp = resamp.interpolate(method='linear')

# Create a DataFrame from the resampled Series
resamp_df = pd.DataFrame(resamp)

# Assign column names to the DataFrame
resamp_df.columns = ['height_resampled']


# In[10]:


plt.figure()
plt.plot(resamp_df.index, resamp_df['height_resampled'], linewidth = 0.8)
plt.grid()
plt.xlabel('Date, Time')
plt.ylabel('Height - Avg Height (m)')
plt.axvline(pd.Timestamp('2011-03-11 05:46:24'), c = 'red')
plt.title('Resampled Data to 15s')
            


# In[11]:


heights = resamp_df['height_resampled']
N = len(resamp_df['height_resampled'])
sample_freq = 4 #4 samples per min
dt = 1 / sample_freq

iturup_fft = np.fft.fftshift(np.fft.fft(heights))
freq = np.fft.fftshift(np.fft.fftfreq(N, dt))
amp = np.abs(iturup_fft)

plt.figure()
plt.semilogx(freq, amp)
plt.grid(which = 'both')
plt.xlabel('Frequency')
plt.ylabel('Amplitude')


# In[20]:


# Define the filter parameters
poles = 4  # Filter order
fc =  0.0001 # Corner frequency in Hz
fs = 1/15 #1 sample every 15 sec 

# Calculate the normalized corner frequency
fnyquist = 0.5 * fs
normalized_corner_freq = fc / fnyquist

# Design the Butterworth highpass filter
b, a = butter(poles, normalized_corner_freq, btype='high', analog=False)


#get the frequency response
w, h = scipy.signal.freqz(b, a, worN=4096)

#convert from angular to linear freq
f = (w) / (np.pi)
f = f*fnyquist

iturup_filt = lfilter(b, a, heights)
iturup_filt2 = filtfilt(b, a, heights)

plt.figure()
plt.plot(resamp_df.index,heights,label='unfiltered')
plt.grid()

plt.figure()
plt.plot(resamp_df.index,iturup_filt,label=r'one pass, fc = %s' % fc, linewidth = 0.5)
plt.plot(resamp_df.index,iturup_filt2,label=r'two pass, fc = %s' % fc, linewidth = 0.5, c = 'green')
plt.axvline(pd.Timestamp('2011-03-11 05:46:24'), c = 'red')
plt.legend()
plt.grid()


# In[ ]:




