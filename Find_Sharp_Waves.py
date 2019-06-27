
#extract and analyze sharp waves during rest

import VideoUtils as vu
import EEGUtils as eeg
import GeneralUtils as gu
from matplotlib import pyplot as plt
import matplotlib.lines as lines
import numpy as np
import sys
from scipy.signal import hilbert, chirp
import peakutils

def removeBlips(ts,index,time_thresh):
  """zero out short blips < time_thresh in indicator function"""
  #make an indicator function that has zeros where the index and ones everywhere else 
  mirror  = np.ones(len(ts))
  mirror[index] = 0

  #mirror index is 1 where above index is 0
  mirror_index = np.where(mirror == 1)[0]
  deblipped = np.abs(min_duration_indicator(ts, mirror_index,time_thresh) -1)

  return deblipped


def returnPeakIndices(eegvals, m, t):
     ind_peaks = peakutils.indexes(eegvals, thres=t, min_dist=m)
     ind_troughs = peakutils.indexes(eegvals*-1, thres=t, min_dist=m)
     return ind_peaks, ind_troughs

def min_duration_indicator(ts, index, thresh):
  """Returns groups of an indicator function that are  > thresh."""
  indicator = np.zeros(len(ts))
  indicator[index] = 1
 
  grouped = gu.group_consecutives(index)
  for g in grouped:
#    print(ts[g[-1]]-ts[g[0]])
    if (ts[g[-1]]-ts[g[0]] < thresh):
      indicator[g] = 0
 
  return indicator


fs = 512/.257552

eegfile = sys.argv[1]
eegdata, eegtimestamps = eeg.readEEG(eegfile)
#ts = np.linspace(eegtimestamps[0], eegtimestamps[-1], len(eegdata))
#print(1000000*len(eegdata)/(ts[-1]-ts[0]))

filtered = eeg.butter_bandpass_filter(eegdata, 140, 240, fs, order=4)

npzfile = np.load("./RawData/NOTMOVING.npz")
starts=npzfile['arr_0'].astype(int)
stops=npzfile['arr_1'].astype(int)

#find closest eegtimestamps to tstart, tstop; starts and stops * 100
EEGstart = [eeg.takeClosest(eegtimestamps, start*100) for start in starts]
EEGstop = [eeg.takeClosest(eegtimestamps, stop*100) for stop in stops] 

#there are 512 eeg samples per timestamp
startidx = [np.where(eegtimestamps==start)[0][0]*512 for start in EEGstart]
stopidx =  [np.where(eegtimestamps==stop)[0][0]*512 for stop in EEGstop]
raw_epochs = [eegdata[start:stop] for start, stop in zip(startidx,stopidx)]
#for ep in raw_epochs:
#  plt.plot(ep)
#  plt.show()

#the sum function below flattens the list of list
epoch = sum([list(filtered[start:stop]) for start, stop in zip(startidx,stopidx)], [])
raw_epoch = sum([list(eegdata[start:stop]) for start, stop in zip(startidx,stopidx)], [])

ts = np.linspace(0, 1000*len(epoch)/(1988), len(epoch))
#uncomment below to show all the epochs together
#plt.plot(ts,big_e)
#plt.show()

analytic_signal = hilbert(epoch)
amp_env = np.abs(analytic_signal)
   
rip_detection_thresh = np.mean(amp_env) + (np.std(amp_env) * 3.5)
rip_duration_thresh = np.mean(amp_env) + (np.std(amp_env) * 1.0)

'''
plt.plot(amp_env)
plt.plot(big_e, 'k')
plt.plot([0, len(amp_env)], [upper_rip_thresh, upper_rip_thresh], 'c')
plt.plot([0, len(amp_env)], [lower_rip_thresh, lower_rip_thresh], 'r')
plt.title("Amplitude Envelope")
plt.show()
'''

#the section below detects putative ripples using the detection threshold and
#finds their duration using the duration threshold

detection_index = np.where(amp_env > rip_detection_thresh)
duration_index = np.where(amp_env > rip_duration_thresh)

#deblipping eliminates short off conditions when the signal is high
deblipped_detection = removeBlips(ts, detection_index, 50)
deblipped_duration  = removeBlips(ts, duration_index, 50)

#grouping keeps indicators of at a minimum duration (35ms)
detection_indicator = min_duration_indicator(ts,np.where(deblipped_detection==1)[0], 35)
duration_indicator = min_duration_indicator(ts, np.where(deblipped_duration==1)[0], 85)

#rip_label indentifies each separate ripple; the height of the square reps the rip number
#like an indicator function with labels
rip_labels = np.zeros(len(ts), dtype=np.uint16)

#group duration_indicators 
grouped_rips = gu.group_consecutives(np.where(duration_indicator==1)[0])

for i,g in enumerate(grouped_rips):
  rip_labels[g] = i
   
durations = duration_indicator * rip_labels
detections = detection_indicator * rip_labels

unique_durations = np.unique(durations)
unique_detections = np.unique(detections)

rips = np.intersect1d(unique_durations, unique_detections)
print(len(rips))
'''
#remove detections that are not in rips
for rip in rips:
   rip_indx = np.where(durations==rip)
   #print(rip_indx[0][0])
   #print(rip_indx[0][-1])
   rip_start = rip_indx[0][0]
   rip_stop = rip_indx[0][-1]
   zeroed_ts = ts[rip_start:rip_stop]-ts[rip_start]
   
   plt.plot(zeroed_ts, epoch[rip_start:rip_stop])
   plt.plot(zeroed_ts, raw_epoch[rip_start:rip_stop], 'r')
   plt.show()
   print('*******')  


'''
#print(rips)

#remove durations with no detections
'''
plt.plot(ts,rips,'r')
plt.plot(ts,detections, 'k')
plt.show()


vals = np.unique(rips)

centered_rips = np.empty([len(vals), 4000])
centered_rips[:] = np.nan

for i,val in enumerate(vals):
  print("i: {}, val: {}".format(i,val))
  start = grouped_below[int(val)][0]
  stop = grouped_below[int(val)][-1]
#  peak = np.max(
  rip = filtered[start:stop]
  ripts = np.arange(len(rip))
  peaks, troughs = returnPeakIndices(rip ,10, .25)
  print(peaks)
  print(troughs)
  center_trough = troughs[np.round(len(troughs)/2)]
  print(center_trough)
  c_to_end = rip[center_trough:]
  c_to_begin = rip[:center_trough]

  centered_rips[i,2000:2000+len(c_to_end)] = c_to_end
  centered_rips[i,2000-len(c_to_begin):2000] = c_to_begin

 
  print("rip: {},  duration -- start: {}, stop: {}".format(int(val), start, stop))
  #plt.plot(ripts, rip)
  #plt.plot(ripts[troughs], rip[troughs], 'r.')
   
# plt.plot(theta[start-100:stop+100], 'r')
  #plt.show()
  #plt.plot(centered_rips[i,:])
  #plt.show()
#if plotit == 'plot':
#  plt.plot(filtered)
#  plt.plot(big_e)
#  plt.plot(rips*np.max(big_e), 'k') 
#  plt.plot(rips*upper_rip_thresh,'k')
#  plt.show()   
 
#outfile = eegfile.replace('.Ncs', '_RIPS')
#np.savez(outfile, outfiltered, outrips, outupper, outlower, outts)

rip_mean = np.nanmean(centered_rips,0)

plt.plot(rip_mean)
plt.show()

'''
