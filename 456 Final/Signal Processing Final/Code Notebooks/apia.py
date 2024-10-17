#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import scipy 
from scipy.signal import resample_poly

get_ipython().run_line_magic('matplotlib', 'notebook')


# In[3]:


#read in csv file for Pitcarin buoy
apia = pd.read_csv(r'C:\Users\kelsa\Dropbox\456Final\Buoys\apia.csv')
#convert month column '3' to '03' format
apia['MM'] = apia['MM'].apply(lambda x: '{:02d}'.format(x))

#drop values of height = 9999
apia = apia[apia['HEIGHT'] != 9999]


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
apia['datetime'] = pd.to_datetime(apia[year].astype(str) + 
                                       apia[month].astype(str).str.zfill(2) + 
                                       apia[day].astype(str).str.zfill(2) +
                                       apia[hour].astype(str).str.zfill(2) +
                                       apia[minute].astype(str).str.zfill(2) +
                                       apia[second].astype(str).str.zfill(2), 
                                       format='%Y%m%d%H%M%S')
print(apia.head())
time = apia['datetime']


# In[5]:


# Extract sample frequency "T" and "height" columns
T = apia['T']
height = apia['HEIGHT']


# In[6]:


#separate data by sampling frequency
apia1 = apia[apia['T'] == 1].copy()
apia2 = apia[apia['T'] == 2].copy()
apia3 = apia[apia['T'] == 3].copy()

#reset the index for each DataFrame
apia1.reset_index(drop=True, inplace=True)
apia2.reset_index(drop=True, inplace=True)
apia3.reset_index(drop=True, inplace=True)


# In[7]:


#plot data
cmap = ListedColormap(plt.cm.viridis(np.linspace(0, 1, 3)))
plt.figure()
plt.plot(time, height, linewidth=1)
plt.grid()
plt.show()


# In[10]:


# Create subplots
fig, axs = plt.subplots(3, 1, sharex = True, sharey = True)

# Scatter plot for pitcarin1
axs[0].scatter(apia1['datetime'], apia1['HEIGHT'], s=1)
axs[0].set_title('15min sampling')
axs[0].grid()

# Scatter plot for pitcarin2
axs[1].scatter(apia2['datetime'], apia2['HEIGHT'], s=1)
axs[1].set_title('1min sampling')
axs[1].grid()

# Scatter plot for pitcarin3
axs[2].scatter(apia3['datetime'], apia3['HEIGHT'], s=1)
axs[2].set_title('15s sampling')
axs[2].grid()

# Adjust layout
plt.tight_layout()

# Show plot
plt.show()


# In[ ]:





# In[ ]:




