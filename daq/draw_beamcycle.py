#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import matplotlib.pyplot as plt


# # Data reading

# In[ ]:


rdata = pd.read_table("/home/yong/workspace/data/2.6/digi_select/P_2.6_2019_Aug_31_13_44_43_scalor.csv"
                      ,sep=r',',skipinitialspace=True)


# # Data cleaning

# ## Insert a new column 'Datetime' from the time information columns and set it as index
# 1. Convert column title 'Time_us to 'second' and 'Time_us' to 'microsecond'
# 2. Insert a datetime column with title 'Datatime' whose value is computed based on 'second' and 'microsecond'
# 3. Set 'Datatime' column as the index of the table

# In[2]:


rdata.rename(columns={'Time_s':'second',
                     'Time_us':'microsecond'},inplace= True)

rdata['Datetime']=pd.to_datetime(rdata.second*1e6+rdata.microsecond, unit= 'us')
rdata = rdata.set_index('Datetime')


# ## Insert records in between storage cycle with 1s second interval
# 1. During injection cycle (about 30s interval), DAQ is not recording data thus scalor not updated.
# 2. New records are inserted to have a smooth time series without big jumps.
# 3. Missing records are inserted with the copied value of the last record.

# In[3]:


# reindex with new 1s interval between beam cycle
nindex = []
l = rdata.index.size # lenth of the table
for index, item in list(enumerate(rdata.index)):
    nindex.append(item)
    if(index == l-1):
        break
    
    delta = (rdata.index[index+1]-item).total_seconds()
    if(delta>5): # injection cycle beginning found
        time = item.floor('s') # get rid of us
        for i in range(1,int(delta//5)):
            nindex.append(time + pd.Timedelta(i*5, unit='s'))

rdata = rdata.reindex(nindex, method='ffill', copy=False)


# # Data analysis
# ## Calculate the rates

# In[4]:


# calculate the difference with previous event, including time and counts
diff = rdata.diff()

# caluclate the rates and fill NA cells with zeros
tdiff=diff.second+1e-6*diff.microsecond

rate=diff.truediv(tdiff,axis='index')

rate=rate.fillna(0).add_suffix('_Rate')

rate=rate.drop(columns=['second_Rate','microsecond_Rate'])

# compute the ratio of DAQ efficiency
rate['DAQ_Efficiency']= rate['Trigger_Rate']/rate['CommonOR_Rate']*100
rate['Fwd_InOut_Ratio']= rate['Fwd_Inside_Rate']/rate['Fwd_Outside_Rate']
rate['Fwd_UpDown_Ratio']= rate['Fwd_Up_Rate']/rate['Fwd_Down_Rate']
rate=rate.fillna(0)

# resample every 2 seconds, and use the center of the time window as label
rate_reduced= rate.resample('2s',label= 'left', loffset= '1s').mean().fillna(0)


# # Data visualization

# In[144]:


list_ver = ['Fwd_Up_Rate', 'Fwd_Down_Rate']
list_hor = ['Fwd_Inside_Rate', 'Fwd_Outside_Rate']
cv_size = (14,6) # figure size


# ## Drawing the time serious
# ### Method1

# In[142]:


rate_reduced['Fwd_Up_Rate'].plot(figsize=cv_size)
ax1 = rate_reduced['Fwd_Down_Rate'].plot()
ax1.set_ylabel('Rate (1/s)')
ax2 = rate_reduced['Fwd_UpDown_Ratio'].plot(secondary_y=True)
ax2.set_ylabel('Ratio (Up/Down)')


# ### Method 2

# In[145]:


rate_reduced[list_hor].plot(figsize=cv_size)
rate_reduced.plot(y=list_ver, figsize=cv_size)


# ## Drawing histogram

# In[146]:


ax = rate_reduced[list_ver].plot.hist(bins=100, alpha=0.5, figsize=(9,6))
ax.set_ylabel('Events')
ax.set_xlabel('Rate')


# In[44]:


plt.savefig('/home/yong/Desktop/test.png', dpi=300)
plt.show()


# ## Use a different plot style

# In[49]:


import matplotlib.style as style
style.available


# In[171]:


style.use('bmh')
rate_reduced[list_ver].plot(figsize=cv_size)
style.use('fivethirtyeight')
rate_reduced[list_hor].plot(figsize=cv_size)


# ## Save figure into disk

# In[177]:


style.use('seaborn-white')

plt.figure()
ax1=rate_reduced.CommonOR_Rate['2019-08-31 12:00:00':'2019-08-31 12:15:00'].plot(figsize=(10,6),colormap='prism')
ax1.set_ylabel('Input Trigger Rate (1/s)',color='r',fontsize='x-large',fontweight='bold')
ax1.tick_params(axis='y',labelcolor='r', labelsize='large')
ax2=rate_reduced['DAQ_Efficiency']['2019-08-31 12:00:00':'2019-08-31 12:15:00'].plot(secondary_y=True,color='blue')
ax2.set_ylabel('DAQ Efficiency (%)',color='b', fontsize='x-large',fontweight='bold')
ax2.tick_params(axis='y',labelcolor='b',labelsize='large')

plt.savefig('/home/yong/Desktop/test.png',dpi=300)


# ## Subplots

# In[185]:


style.use('bmh')
list_all = list_hor + list_ver + ['DAQ_Efficiency']
ax1=rate_reduced[list_all]['2019-08-31 12:00:00':'2019-08-31 12:20:00'].plot(subplots=True,figsize=cv_size)


# In[ ]:




