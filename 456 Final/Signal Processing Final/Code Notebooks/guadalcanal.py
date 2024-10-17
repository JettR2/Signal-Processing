#!/usr/bin/env python
# coding: utf-8

# In[5]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import scipy 
from scipy.signal import resample_poly, butter, lfilter,freqs,filtfilt,sosfilt,sosfiltfilt
from scipy.fft import fft

get_ipython().run_line_magic('matplotlib', 'notebook')


# In[6]:


#read in csv file for Pitcarin buoy
guadalcanal = pd.read_csv(r'C:\Users\kelsa\Dropbox\456Final\Buoys\guadalcanal.csv')
#convert month column '3' to '03' format
guadalcanal['MM'] = guadalcanal['MM'].apply(lambda x: '{:02d}'.format(x))

#drop values of height = 9999
guadalcanal = guadalcanal[guadalcanal['HEIGHT'] != 9999]


# In[8]:


#converting columns into single datetime object
# assigning columns
year = '#YY'  
month = 'MM'  
day = 'DD' 
hour = 'hh'  
minute = 'mm' 
second = 'ss' 

# Adding a new datetime column for year, month, day, hour, minute, and second
guadalcanal['Date'] = pd.to_datetime(guadalcanal[year].astype(str) + 
                                       guadalcanal[month].astype(str).str.zfill(2) + 
                                       guadalcanal[day].astype(str).str.zfill(2) +
                                       guadalcanal[hour].astype(str).str.zfill(2) +
                                       guadalcanal[minute].astype(str).str.zfill(2) +
                                       guadalcanal[second].astype(str).str.zfill(2), 
                                       format='%Y%m%d%H%M%S')
guadalcanal.set_index(guadalcanal['Date'], inplace=True)
time = guadalcanal['Date']


# In[14]:


# Extract sample frequency "T" and "height" columns
T = guadalcanal['T']
height = guadalcanal['HEIGHT']

height_avg = np.mean(height)
guadalcanal['height_norm'] = height - height_avg


# In[10]:


#separate data by sampling frequency
guadalcanal1 = guadalcanal[guadalcanal['T'] == 1].copy()
guadalcanal2 = guadalcanal[guadalcanal['T'] == 2].copy()
guadalcanal3 = guadalcanal[guadalcanal['T'] == 3].copy()

#reset the index for each DataFrame
guadalcanal1.reset_index(drop=True, inplace=True)
guadalcanal2.reset_index(drop=True, inplace=True)
guadalcanal3.reset_index(drop=True, inplace=True)


guadalcanal1.set_index(guadalcanal1['Date'], inplace=True)
time1 = guadalcanal1['Date']
guadalcanal1.drop(columns=['Date'], inplace=True)

guadalcanal2.set_index(guadalcanal2['Date'], inplace=True)
time2 = guadalcanal2['Date']
guadalcanal2.drop(columns=['Date'], inplace=True)

guadalcanal3.set_index(guadalcanal3['Date'], inplace=True)
time3 = guadalcanal3['Date']
guadalcanal3.drop(columns=['Date'], inplace=True)

guadalcanal.drop(columns=['Date'], inplace=True)


# In[17]:


#plot data
cmap = ListedColormap(plt.cm.viridis(np.linspace(0, 1, 3)))
plt.figure()
plt.plot(time, height, linewidth=1)
plt.axvline(pd.Timestamp('2011-03-11 05:46:24'), c = 'red')
plt.grid()
plt.show()


# In[22]:


# Create subplots
fig, axs = plt.subplots(3, 1, sharex = True, sharey = True)

# Scatter plot for pitcarin1
axs[0].scatter(guadalcanal1['datetime'], guadalcanal1['HEIGHT'], s=1)
axs[0].set_title('15min sampling')
axs[0].grid()

# Scatter plot for pitcarin2
axs[1].scatter(guadalcanal2['datetime'], guadalcanal2['HEIGHT'], s=1)
axs[1].set_title('1min sampling')
axs[1].grid()

# Scatter plot for pitcarin3
axs[2].scatter(guadalcanal3['datetime'], guadalcanal3['HEIGHT'], s=1)
axs[2].set_title('15s sampling')
axs[2].grid()

# Adjust layout
plt.tight_layout()

# Show plot
plt.show()


# In[15]:


# Resample the 'HEIGHT' column to 1 minute frequency
resamp = guadalcanal['height_norm'].resample('15S').mean()

# Optionally, fill missing values using linear interpolation
resamp = resamp.interpolate(method='linear')

# Create a DataFrame from the resampled Series
resamp_df = pd.DataFrame(resamp)

# Assign column names to the DataFrame
resamp_df.columns = ['height_resampled']


# In[16]:


plt.figure()
plt.plot(resamp_df.index, resamp_df['height_resampled'], linewidth = 0.8)
plt.grid()
plt.xlabel('Date, Time')
plt.ylabel('Height - Avg Height (m)')
plt.axvline(pd.Timestamp('2011-03-11 05:46:24'), c = 'red')
plt.title('Resampled Data to 15s')
            


# In[19]:


heights = resamp_df['height_resampled']
N = len(resamp_df['height_resampled'])
sample_freq = 4 #4 samples per min
dt = 1 / sample_freq

guadalcanal_fft = np.fft.fftshift(np.fft.fft(heights))
freq = np.fft.fftshift(np.fft.fftfreq(N, dt))
amp = np.abs(guadalcanal_fft)

plt.figure()
plt.semilogx(freq, amp)
plt.grid(which = 'both')
plt.xlabel('Frequency')
plt.ylabel('Amplitude')


# In[25]:


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

guadalcanal_filt = lfilter(b, a, heights)
guadalcanal_filt2 = filtfilt(b, a, heights)

plt.figure()
plt.plot(resamp_df.index,heights,label='unfiltered')
plt.grid()

plt.figure()
plt.plot(resamp_df.index,guadalcanal_filt,label=r'one pass, fc = %s' % fc, linewidth = 0.5)
plt.plot(resamp_df.index, guadalcanal_filt2,label=r'two pass, fc = %s' % fc, linewidth = 0.5, c = 'green')
plt.axvline(pd.Timestamp('2011-03-11 05:46:24'), c = 'red')
plt.legend()
plt.grid()


# In[ ]:




