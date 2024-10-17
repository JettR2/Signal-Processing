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


# In[2]:


#read in csv file for Pitcarin buoy
pitcairn = pd.read_csv(r'C:\Users\hneuman\Downloads\Pitcarin_SPA_Buoy.csv')
#convert month column '3' to '03' format
pitcairn['MM'] = pitcairn['MM'].apply(lambda x: '{:02d}'.format(x))

#drop values of height = 9999
pitcairn = pitcairn[pitcairn['HEIGHT'] != 9999]


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
pitcairn['datetime'] = pd.to_datetime(pitcairn[year].astype(str) + 
                                       pitcairn[month].astype(str).str.zfill(2) + 
                                       pitcairn[day].astype(str).str.zfill(2) +
                                       pitcairn[hour].astype(str).str.zfill(2) +
                                       pitcairn[minute].astype(str).str.zfill(2) +
                                       pitcairn[second].astype(str).str.zfill(2), 
                                       format='%Y%m%d%H%M%S')
print(pitcairn.head())
time = pitcairn['datetime']


# In[31]:


# Extract sample frequency "T" and "height" columns
T = pitcairn['T']
height = pitcairn['HEIGHT']

avg_height = np.mean(pitcairn['HEIGHT'])

pitcairn['avg_height'] = np.mean(pitcairn['HEIGHT'])
pitcairn['height_norm'] = height - avg_height


# In[32]:


#separate data by sampling frequency
pitcairn1 = pitcairn[pitcairn['T'] == 1].copy()
pitcairn2 = pitcairn[pitcairn['T'] == 2].copy()
pitcairn3 = pitcairn[pitcairn['T'] == 3].copy()

#reset the index for each DataFrame
pitcairn1.reset_index(drop=True, inplace=True)
pitcairn2.reset_index(drop=True, inplace=True)
pitcairn3.reset_index(drop=True, inplace=True)


# In[33]:


#plot data
cmap = ListedColormap(plt.cm.viridis(np.linspace(0, 1, 3)))
plt.figure()
plt.plot(time, pitcairn['height_norm'], linewidth=1)
plt.scatter(time,  pitcairn['height_norm'], c = pitcairn['T'], cmap = cmap, s = 4)
plt.axvline(pd.Timestamp('2011-03-11 05:46:24'), c = 'red')
plt.grid()
plt.colorbar()
plt.show()


# In[34]:


# Create subplots
fig, axs = plt.subplots(3, 1, sharex = True, sharey = True)

# Scatter plot for pitcarin1
axs[0].scatter(pitcairn1['datetime'], pitcairn1['height_norm'], s = 3)
axs[0].set_title('15min sampling')
axs[0].axvline(pd.Timestamp('2011-03-11 05:46:24'), c = 'red')
axs[0].grid()

# Scatter plot for pitcarin2
axs[1].scatter(pitcairn2['datetime'], pitcairn2['height_norm'], s = 3)
#axs[1].plot(pitcairn2['datetime'], pitcairn2['height_norm'])
axs[1].axvline(pd.Timestamp('2011-03-11 05:46:24'), c = 'red')
axs[1].set_title('1min sampling')
axs[1].grid()

# Scatter plot for pitcarin3
axs[2].scatter(pitcairn3['datetime'], pitcairn3['height_norm'], s = 3)
axs[2].set_title('15s sampling')
axs[2].axvline(pd.Timestamp('2011-03-11 05:46:24'), c = 'red')
axs[2].grid()

# Adjust layout
plt.tight_layout()

# Show plot
plt.show()


# In[ ]:





# In[ ]:




