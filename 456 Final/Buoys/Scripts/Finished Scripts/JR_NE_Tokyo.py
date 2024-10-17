# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 11:24:01 2024

@author: Jett R
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import scipy 
from scipy.signal import resample_poly

#get_ipython().run_line_magic('matplotlib', 'notebook')


# In[3]:


#read in csv file for Pitcarin buoy
Buoy = pd.read_csv(r'C:/Users/Jett R/Dropbox (University of Oregon)/23-24/Winter/Erth 456 Sig Proc/456 Final/Buoys/NE_Tokyo.csv')
#convert month column '3' to '03' format
Buoy['MM'] = Buoy['MM'].apply(lambda x: '{:02d}'.format(x))

#drop values of height = 9999
Buoy = Buoy[Buoy['HEIGHT'] != 9999]


# In[4]:


#converting columns into single datetime object
# assigning columns
year = '#YY'  
month = 'MM'  
day = 'DD' 
hour = 'hh'  
minute = 'mm' 
second = 'ss' 

# Adding a new datetime column for year, month, day, hour, minute, and second
Buoy['datetime'] = pd.to_datetime(Buoy[year].astype(str) + 
                                       Buoy[month].astype(str).str.zfill(2) + 
                                       Buoy[day].astype(str).str.zfill(2) +
                                       Buoy[hour].astype(str).str.zfill(2) +
                                       Buoy[minute].astype(str).str.zfill(2) +
                                       Buoy[second].astype(str).str.zfill(2), 
                                       format='%Y%m%d%H%M%S')
print(Buoy.head())
time = Buoy['datetime']


# In[5]:


# Extract sample frequency "T" and "height" columns
T = Buoy['T']
height = Buoy['HEIGHT']

# Calculate the average wave height
average_wave_height = height.mean()

# Normalize the wave height to the average wave height
normalized_height = height - average_wave_height


# In[6]:


#separate data by sampling frequency
Buoy1 = Buoy[Buoy['T'] == 1].copy()
Buoy2 = Buoy[Buoy['T'] == 2].copy()
Buoy3 = Buoy[Buoy['T'] == 3].copy()

#reset the index for each DataFrame
Buoy1.reset_index(drop=True, inplace=True)
Buoy2.reset_index(drop=True, inplace=True)
Buoy3.reset_index(drop=True, inplace=True)


# In[7]:


#plot data
cmap = ListedColormap(plt.cm.viridis(np.linspace(0, 1, 3)))
plt.figure(figsize = (10,8))
plt.plot(time, normalized_height, linewidth=1)
plt.title('Buoy Tides')
plt.axvline(pd.Timestamp('2011-03-11 05:46:24'), color='red', linestyle='--', label='Earthquake')
plt.ylabel('Normalized Wave Height')
plt.xlabel('Time')
plt.grid()
plt.show()


# In[10]:


# Create subplots
fig, axs = plt.subplots(3, 1, figsize=(12, 5), sharex=True, sharey=True)

# Scatter plot for Buoy1
axs[0].scatter(Buoy1['datetime'], normalized_height[Buoy['T'] == 1], s=1)
axs[0].set_title('15min sampling')
axs[0].axvline(pd.Timestamp('2011-03-11 05:46:24'), color='red', linestyle='--', label='Earthquake')
axs[0].grid()

# Scatter plot for Buoy2
axs[1].scatter(Buoy2['datetime'], normalized_height[Buoy['T'] == 2], s=1)
axs[1].set_title('1min sampling')
axs[1].axvline(pd.Timestamp('2011-03-11 05:46:24'), color='red', linestyle='--', label='Earthquake')
axs[1].grid()

# Scatter plot for Buoy3
axs[2].scatter(Buoy3['datetime'], normalized_height[Buoy['T'] == 3], s=1)
axs[2].set_title('15s sampling')
axs[2].axvline(pd.Timestamp('2011-03-11 05:46:24'), color='red', linestyle='--', label='Earthquake')
axs[2].grid()

# Adjust layout
plt.tight_layout()
plt.legend()

# Show plot
plt.show()



