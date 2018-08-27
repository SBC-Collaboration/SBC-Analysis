from SBCcode.DataHandling.ReadBinary import ReadBlock as RB
import numpy as np
from SBCcode.DataHandling.GetSBCEvent import GetEvent as GE
import os
import scipy.signal as sig
# import matplotlib
# from matplotlib.widgets import TextBox

# Load necessary analysis binary files
RECON_DIRECTORY = "/pnfs/coupp/persistent/grid_output/SBC-17/output/"
RECON_DIRECTORY2 = "/pnfs/coupp/persistent/grid_output/SBC-17-T0Test3/output/"

coordinate_data = RB(os.path.join(RECON_DIRECTORY,"SimpleXYZ_all.bin"))
# acoustic_data = RB(os.path.join(RECON_DIRECTORY,"AcousticAnalysis_all.bin"))
acoustic_data = RB(os.path.join(RECON_DIRECTORY2,"AcousticAnalysis_all.bin"))
timing_data = RB(os.path.join(RECON_DIRECTORY2,"TimingAnalysis_all.bin"))
history_data = RB(os.path.join(RECON_DIRECTORY, "HistoryAnalysis_all.bin"))
event_data = RB(os.path.join(RECON_DIRECTORY, "EventAnalysis_all.bin"))
align_data = RB(os.path.join(RECON_DIRECTORY,"PMTfastDAQalignment_all.bin"))
human_get_bub_data = RB(os.path.join(RECON_DIRECTORY, "HumanGetBub_all.bin"))



# Create a boolean array to cut for only Californium runs
cali_run = (event_data['Event_run_type'] == 10252) | (event_data['Event_run_type'] == 20252)

# Save the event pressures for transducers 4 and 6 and create a boolean array to cut for low
# pressure events. No floor pressure is specified because most 10252 runs have pressures
# at -3 for transducer 4, which was broken
pressures = history_data['EventPressure'][0:,3]
pressures5 = history_data['EventPressure'][0:,5]
less_than_50 = (pressures<50)

# Combining Californium run boolean and pressure boolean arrays. Comment out first line and
# uncomment second line to use all events, not cutting for run type or pressure. May also use this cut to show specific
# events.

interesting_events = cali_run*less_than_50
# interesting_events=(event_data['runid'][:,0]==20171006) * (event_data['runid'][:,1]==3) * (event_data['ev']==38) #+\
# interesting_events = (event_data['runid'][:,0]==20171007) * (event_data['runid'][:,1]==6) * (event_data['ev']==3)
# interesting_events = np.ones(event_data['runid'].shape[0],dtype=bool)
#interesting_events = (pressures<30)



# Create a starting list of all of the runids and event numbers. eventsoi is a list of tuples
# of length two. The first element is a string with the date of the run followed by an underscore
# and the run number. The second element is the event number and is an integer. eventsoi is later
# turned into a numPy array to allow for boolean indexing.

eventsoi = []
last_eventoi = (None, None)
first_bubble = True
for run_id, ev in zip(timing_data["runid"], timing_data["ev"]):
    current_event = (str(run_id[0]) + "_" + str(run_id[1]), ev)
    if current_event != last_eventoi:
        first_bubble = True
    if first_bubble:
        first_bubble = False
        last_eventoi = current_event
        eventsoi.append(current_event)



# Create a cut to only include data taken when the camera was off. Create a cut that
# only includes good data where the pmt lag is not equal to -1.

cam_cut=timing_data['CAMstate'][:,0]==1
good_lag = timing_data['PMTmatch_lag'][:,1] != -1

# Create a boolean array for events where the realspace coordinate is known. Other elements are
# numPy nan types.

x = np.isfinite(coordinate_data['bubZ'])

# Merge camera cut, interesting events, good data and events with known realspace coordinates into one boolean array.
# Then save the real space coordinates of the data indexed by our final cut. Then save align times for PMT triggers.

final_cut=cam_cut * interesting_events * good_lag * x
# xcoords = coordinate_data['bubX'][final_cut]
# ycoords = coordinate_data['bubY'][final_cut]
zcoords = coordinate_data['bubZ'][final_cut]
align_times= (align_data['PMT_trigt0_sec']+align_data['PMT_trigt0_frac'])[final_cut]


# Load the pmt lag times calculated in TimingAnalysis_all.bin and acoustic t0s from AcousticAnalysis_all.bin,
# both indexed by our final cut. Calculate the aligned pmt t0 of the matched trigger from the two arrays.
pmt_lags = timing_data['PMTmatch_lag'][final_cut]
bubble_t0s = acoustic_data['bubble_t0'][final_cut]
pmt_match_t0 = bubble_t0s + pmt_lags


# Condense various array names for use in class definition
event_pressure3 = pressures[final_cut]
event_pressure5=pressures5[final_cut]
# bubble_t0s1=bubble_t0s[:,0]
# bubble_t0s2=bubble_t0s[:,1]
event_run_type=event_data['Event_run_type'][final_cut]

yoi = zcoords
eventsoi=np.array(eventsoi)[final_cut]
# num_human_bubbles = 0
# last_event= (None,None)
# for run,ev,nbub in zip(human_get_bub_data['runid'],human_get_bub_data['ev'],human_get_bub_data['nbubimage']):
#     runid = str(run[0])+'_'+str(run[1])
#     current_event = (runid, str(ev))
#     if current_event == last_event:
#         last_event = current_event
#         continue
#     else:
#         last_event = current_event
#     if current_event in eventsoi:
#         if nbub > 0:
#             num_human_bubbles += 1
#             continue
#         else: continue
#     else: continue



RAW_DIRECTORY = "/bluearc/storage/SBC-17-data/"

import matplotlib.pyplot as plt

bub1 = acoustic_data['bubble_t0'][final_cut][:,0]
bub2 = acoustic_data['bubble_t0'][final_cut][:,1]

finite_t0s = np.isfinite(acoustic_data['bubble_t0'])[final_cut]
# nans=np.sum(np.isnan(new_acoustic_data['bubble_t0'])[final_cut])
pmt_lags = pmt_lags[finite_t0s[:,0]*finite_t0s[:,1]]
pmt_lag1 = pmt_lags[:,0]
pmt_lag2 = pmt_lags[:,1]

xs1 = yoi[finite_t0s[:,0]*finite_t0s[:,1]]
ys1 = pmt_lag1
xs2 = yoi[finite_t0s[:,0]*finite_t0s[:,1]]
ys2 = pmt_lag2
eventsoi = eventsoi[finite_t0s[:,0]*finite_t0s[:,1]]
align_times = align_times[finite_t0s[:,0]*finite_t0s[:,1]]
bub1= bub1[finite_t0s[:,0]*finite_t0s[:,1]]
bub2 = bub2[finite_t0s[:,0]*finite_t0s[:,1]]
# for i in np.arange(eventsoi.shape[0]):
#     # try:
#         # print(os.path.join("/pnfs/coupp/persistent/grid_output/SBC-17-T0Test2/output",eventsoi[i,0],"AcousticAnalysis_"+ eventsoi[i,0]+".bin"))
#     if os.path.isfile(os.path.join("/pnfs/coupp/persistent/grid_output/SBC-17-T0Test3/output", eventsoi[i, 0],
#                             "AcousticAnalysis_" + eventsoi[i, 0] + ".bin")):
#         d = RB(os.path.join("/pnfs/coupp/persistent/grid_output/SBC-17-T0Test3/output", eventsoi[i, 0],
#                             "AcousticAnalysis_" + eventsoi[i, 0] + ".bin"))
#         if int(eventsoi[i,1]) <= d['bubble_t0'].shape[0]:
#             new_t0s = d['bubble_t0'](int([eventsoi[i, 1]]))
#         else:
#             print("RunID {} and event {} failed.".format(eventsoi[i, 0], eventsoi[i, 1]))
#             new_lags_cut[i] = False
#             continue
#         if np.isnan(new_t0s[0]) or np.isnan(new_t0s[1]):
#             print("RunID {} and event {} has nan value.".format(eventsoi[i, 0], str(eventsoi[i, 1])))
#             new_lags_cut[i] = False
#             new_bub1 = np.append(new_bub1, -1)
#             new_bub2 = np.append(new_bub2, -1)
#             nans +=1
#             continue
#         ys1[i] = ys1[i] + bubble_t0s[i, 0] - new_t0s[0]
#         ys2[i] = ys2[i] + bubble_t0s[i, 1] - new_t0s[1]
#         new_bub1 = np.append(new_bub1, new_t0s[0])
#         new_bub2 = np.append(new_bub2, new_t0s[1])
#     # except:
#     else:
#         print("Loading RunID {} failed.".format(eventsoi[i, 0]))
#         new_lags_cut[i] = False
#         new_bub1 = np.append(new_bub1, -1)
#         new_bub2 = np.append(new_bub2, -1)
# xs1 = xs1[new_lags_cut]
# ys1 = ys1[new_lags_cut]
# xs2 = xs2[new_lags_cut]
# ys2 = ys2[new_lags_cut]
# new_bub1 = new_bub1[new_lags_cut]
# new_bub2 = new_bub2[new_lags_cut]


# Class definition for interactive matplotlib pylot
class PointBrowser(object):
    """
    Click on a point to select and highlight it -- the data that
    generated the point will be shown in the lower axes.  Use the 'n'
    and 'p' keys to browse through the next and previous points
    """

    def __init__(self):
        # Define some class variables, including the last selected point index to index the eventsoi list
        self.lastind = 0

        # self.text = axes[0,0].text(0.05, 0.95, 'selected: none',
        #                     transform=axes[0,0].transAxes, va='top')
        self.selected1, = axes[0,0].plot([xs1[0]], [ys1[0]], 'o', ms=12, alpha=0.4,
                                 color='yellow', visible=False)
        self.selected2, = axes[0, 1].plot([xs2[0]], [ys2[0]], 'o', ms=12, alpha=0.4,
                                         color='yellow', visible=False)
        self.event_id = None

        self.ev_number = None
        self.new_event = True
    def onpress(self, action):
        # Method definition for changing selected point via keyboard commands. Press 'n' to switch to next chronological
        # event point, 'p' for previous. 'r' should move to next pmt trigger, 'l' should move to previous. Updates
        # class variable methods accordingly

        if self.lastind is None:
            return
        if action.key not in ('n', 'p','r','l'):
            return
        if action.key == 'n':
            inc = 1
            self.lastind += inc
            self.new_event=True
        elif action.key == 'p':
            inc = -1
            self.lastind += inc
            self.new_event=True
        elif action.key == 'r':
            self.pmt_index +=1
            self.new_event=False
        else:
            self.pmt_index -= 1
            self.new_event=False
        self.lastind = np.clip(self.lastind, 0, len(xs1) - 1)
        self.event_id = eventsoi[self.lastind][0]
        self.ev_number = eventsoi[self.lastind][1]
        self.bubble_time1=bub1[self.lastind]
        self.bubble_time2=bub2[self.lastind]
        self.run_type = event_run_type[self.lastind]
        self.pressure3 = event_pressure3[self.lastind]
        self.pressure5 = event_pressure5[self.lastind]
        self.update()

    def onpick(self, action):
        # Method definition for clicking event time lag point to display the appropriate Piezo traces in the second
        # and third row of plots.

        # Checks that mouse click is in one of the plots in the first row
        if (action.artist != line) and (action.artist != line2):
            return True

        N = len(action.ind)
        if not N:
            return True

        # the click locations
        x = action.mouseevent.xdata
        y = action.mouseevent.ydata

        # finds the closest event point to the mouse click and save that events index
        distances1 = np.hypot(x - xs1[action.ind], y - ys1[action.ind])
        distances2 = np.hypot(x - xs2[action.ind], y - ys2[action.ind])
        dist1min= distances1.argmin()
        dist2min = distances2.argmin()
        if dist1min < dist2min:
            indmin = dist1min
        else:
            indmin = dist2min
        dataind = action.ind[indmin]
        # changes the lastind variable to the index of the point clicked on, and updates the other variables accordingly
        self.lastind = dataind
        self.event_id = eventsoi[dataind][0]
        self.ev_number = eventsoi[dataind][1]
        self.bubble_time1=bub1[dataind]
        self.bubble_time2=bub2[dataind]
        self.run_type = event_run_type[dataind]
        self.pressure3 = event_pressure3[self.lastind]
        self.pressure5 = event_pressure5[self.lastind]
        self.new_event=True
        self.update()

    def update(self):
        # Method to update figure with new plots of selected event data.


        if self.lastind is None:
            return

        dataind = self.lastind

        # clear all plots except position vs time lag
        axes[1,0].cla()
        axes[1,1].cla()
        axes[2,0].cla()
        axes[2,1].cla()
        axes[3,0].cla()
        axes[3,1].cla()

        axes[1, 0].set_title(
            'Piezo 1 trace. Red line is acoustic t0')
        axes[1, 1].set_title(
            'Piezo 2 trace. Red line is acoustic t0')
        axes[2, 0].set_title(
            'Zoomed in on t0. Blue lines are close PMT triggers, purple is selected')
        axes[2, 1].set_title(
            'Zoomed in on t0. Blue lines are close PMT triggers, purple is selected')
        axes[3, 0].set_title(
            'Selected PMT trace. Press R for next trigger in time, L for previous')
        axes[3, 1].set_title(
            'Selected PMT trace. Press R for next trigger in time, L for previous')

        # Get the fastDAQ data of the selected event point. Saves the piezo traces.

        # Get the pmt trace data of the selected event point. Find the pmt triggers around the acoustic t0 and plot them
        # as lines in the zoomed in Piezo trace
        if self.new_event==True:
            self.ev_data = GE(os.path.join(RAW_DIRECTORY, self.event_id), self.ev_number, "fastDAQ")
            self.piezo1 = self.ev_data["fastDAQ"]["Piezo1"]
            self.piezo2 = self.ev_data["fastDAQ"]["Piezo2"]
            self.timebase = self.ev_data["fastDAQ"]["time"]

            self.pmt_data = GE(os.path.join(RAW_DIRECTORY, self.event_id), self.ev_number, "PMTtraces")
            self.pmt_times = (self.pmt_data['PMTtraces']['t0_sec'] + self.pmt_data['PMTtraces']['t0_frac'])-align_times[self.lastind]

            self.close_pmt_times1 = (self.pmt_times - self.bubble_time1 < .001) & (self.pmt_times - self.bubble_time1 > -.001)
            self.close_pmt_times2 = (self.pmt_times - self.bubble_time2 < .001) & (self.pmt_times - self.bubble_time2 > -.001)
            self.pmt_index=0

        for pmt in self.pmt_times[self.close_pmt_times1]:
            axes[2,0].axvline(x=pmt,color='b')
        for pmt in self.pmt_times[self.close_pmt_times2]:
            axes[2,1].axvline(x=pmt,color='b')


        # Plot the selected PMT trace in the bottom row plots. Assuming the matched PMT trigger is the same for both
        # Piezo traces, and the acoustic t0s are not significantly different between the Piezos, the plots will have
        # the same trace. If the index variable for the triggers close to the acoustic t0 is greater than the number
        # of plotted triggers, the try block will fail and the except block will reset the index to zero and start over.
        # Whichever pmt trigger is being read out will be yellow on the zoomed in Piezo trace plot.


        try:
            trigger_indexes1 = np.where(self.close_pmt_times1 == True)[0]
            traces1 = self.pmt_data['PMTtraces']['traces'][trigger_indexes1[self.pmt_index],0,:]
            dt1 = self.pmt_data['PMTtraces']['dt'][trigger_indexes1[self.pmt_index],0]
            v_scale1 = self.pmt_data['PMTtraces']['v_scale'][trigger_indexes1[self.pmt_index],0]
            v_offset1 = self.pmt_data['PMTtraces']['v_offset'][trigger_indexes1[self.pmt_index],0]
            xd1=np.arange(self.pmt_data['PMTtraces']['traces'].shape[2])*dt1
            yd_fine1 = traces1 * v_scale1 \
                           + v_offset1
            axes[2,0].axvline(x=self.pmt_times[self.close_pmt_times1][self.pmt_index],color='m')
            axes[3,0].plot(xd1,yd_fine1)

            trigger_indexes2 = np.where(self.close_pmt_times2 == True)[0]
            traces2 = self.pmt_data['PMTtraces']['traces'][trigger_indexes2[self.pmt_index], 0, :]
            dt2 = self.pmt_data['PMTtraces']['dt'][trigger_indexes2[self.pmt_index], 0]
            v_scale2 = self.pmt_data['PMTtraces']['v_scale'][trigger_indexes2[self.pmt_index], 0]
            v_offset2 = self.pmt_data['PMTtraces']['v_offset'][trigger_indexes2[self.pmt_index], 0]
            xd2 = np.arange(self.pmt_data['PMTtraces']['traces'].shape[2]) * dt2
            yd_fine2 = traces2 * v_scale2 \
                       + v_offset2
            axes[2,1].axvline(x=self.pmt_times[self.close_pmt_times2][self.pmt_index],color='m')
            axes[3, 1].plot(xd2, yd_fine2)
        except IndexError:
            try:
                self.pmt_index=0
                trigger_indexes1 = np.where(self.close_pmt_times1 == True)[0]
                traces1 = self.pmt_data['PMTtraces']['traces'][trigger_indexes1[self.pmt_index], 0, :]
                dt1 = self.pmt_data['PMTtraces']['dt'][trigger_indexes1[self.pmt_index], 0]
                v_scale1 = self.pmt_data['PMTtraces']['v_scale'][trigger_indexes1[self.pmt_index], 0]
                v_offset1 = self.pmt_data['PMTtraces']['v_offset'][trigger_indexes1[self.pmt_index], 0]
                xd1 = np.arange(self.pmt_data['PMTtraces']['traces'].shape[2]) * dt1
                yd_fine1 = traces1 * v_scale1 \
                           + v_offset1
                axes[2, 0].axvline(x=self.pmt_times[self.close_pmt_times1][self.pmt_index], color='m')
                axes[3, 0].plot(xd1, yd_fine1)

                trigger_indexes2 = np.where(self.close_pmt_times2 == True)[0]
                traces2 = self.pmt_data['PMTtraces']['traces'][trigger_indexes2[self.pmt_index], 0, :]
                dt2 = self.pmt_data['PMTtraces']['dt'][trigger_indexes2[self.pmt_index], 0]
                v_scale2 = self.pmt_data['PMTtraces']['v_scale'][trigger_indexes2[self.pmt_index], 0]
                v_offset2 = self.pmt_data['PMTtraces']['v_offset'][trigger_indexes2[self.pmt_index], 0]
                xd2 = np.arange(self.pmt_data['PMTtraces']['traces'].shape[2]) * dt2
                yd_fine2 = traces2 * v_scale2 \
                           + v_offset2
                axes[2, 1].axvline(x=self.pmt_times[self.close_pmt_times2][self.pmt_index], color='mnn')
                axes[3, 1].plot(xd2, yd_fine2)
            except IndexError: print('hello')



        # Plot the acoustic t0 on the zoomed in Piezo trace. Then plot the piezo traces and the filtered traces, and plot
        # the filtered traces zoomed in around the acoustic t0.
        # axes[1,0].axvline(x=self.bubble_time1, color="r")
        # axes[1,1].axvline(x=self.bubble_time2, color="r")
        axes[1,0].axvline(x=self.bubble_time1, color="r")
        axes[1,1].axvline(x=self.bubble_time2, color="r")
        axes[1,0].plot(self.timebase,self.piezo1)
        filter_d = np.exp(-.1)
        axes[1,0].plot(self.timebase, sig.lfilter(np.float64([1])-filter_d,np.float64([1, -filter_d]),self.piezo1),'g')
        axes[1,1].plot(self.timebase, self.piezo2)
        filter_d = np.exp(-.1)
        axes[1,1].plot(self.timebase, sig.lfilter(np.float64([1]) - filter_d, np.float64([1, -filter_d]), self.piezo2), 'g')



        # Shows the event run id, event number, run type, and fourth and sixth pressure tranducer measurements
        axes[2,0].text(0.05, 0.9, 'event_id={}\nevent={}'.format(self.event_id, self.ev_number),
                       transform=axes[2,0].transAxes, va='top')
        axes[2, 1].text(0.05, 0.9,
                        'run_type={}\npressure3={}\npressure5={}'.format(self.run_type, self.pressure3, self.pressure5),
                        transform=axes[2, 1].transAxes, va='top')

        axes[2,0].plot(self.timebase, sig.lfilter(np.float64([1])-filter_d,np.float64([1, -filter_d]),self.piezo1),'g')
        # axes[2,0].set_xlim(self.bubble_time1-.0005,self.bubble_time1+.0005)
        axes[2, 0].set_xlim(self.bubble_time1 - .0005, self.bubble_time1 + .0005)
        # axes[2,0].axvline(x=self.bubble_time1, color="r")
        axes[2,0].axvline(x=self.bubble_time1, color="r")

        axes[2, 1].plot(self.timebase, sig.lfilter(np.float64([1]) - filter_d, np.float64([1, -filter_d]), self.piezo2), 'g')
        # axes[2, 1].set_xlim(self.bubble_time2 - .0005, self.bubble_time2 + .0005)
        axes[2, 1].set_xlim(self.bubble_time2 - .0005, self.bubble_time2 + .0005)
        # axes[2, 1].axvline(x=self.bubble_time2, color="r")
        axes[2,1].axvline(x=self.bubble_time2, color="r")



        # Make the yellow marker for the new selected data points visible

        # ax2.set_ylim(-3, 3)
        self.selected1.set_visible(True)
        self.selected2.set_visible(True)

        self.selected1.set_data(xs1[dataind], ys1[dataind])
        self.selected2.set_data(xs2[dataind], ys2[dataind])



        # self.text.set_text('selected: %d' % dataind)
        fig.canvas.draw()


if __name__ == '__main__':
    # # import matplotlib
    # # matplotlib.rcParams['backend'] = 'TkAgg'
    # import matplotlib.pyplot as plt
    # new_lags_cut = np.ones(eventsoi.shape[0],dtype=bool)
    # xs1 = yoi
    # ys1 = pmt_lag1
    # xs2 = yoi
    # ys2 = pmt_lag2
    # new_bub1 = np.array([])
    # new_bub2=np.array([])
    # for i in np.arange(eventsoi.shape[0]):
    #     try:
    #         # print(os.path.join("/pnfs/coupp/persistent/grid_output/SBC-17-T0Test2/output",eventsoi[i,0],"AcousticAnalysis_"+ eventsoi[i,0]+".bin"))
    #         d = RB(os.path.join("/pnfs/coupp/persistent/grid_output/SBC-17-T0Test2/output",eventsoi[i,0],"AcousticAnalysis_"+ eventsoi[i,0]+".bin"))
    #         new_t0s = d['bubble_t0'][eventsoi[i,1]]
    #         ys1[i] = ys1[i] + bubble_t0s[i,0] - new_t0s[0]
    #         ys2[i] = ys2[i] + bubble_t0s[i,1] - new_t0s[1]
    #         new_bub1 = np.append(new_bub1,new_t0s[0])
    #         new_bub2 = np.append(new_bub2,new_t0s[1])
    #     except:
    #         print("RunID {} and event {} failed.".format(eventsoi[i,0],str(eventsoi[i,1])))
    #         new_lags_cut[i] = False
    #         new_bub1 = np.append(new_bub1,-1)
    #         new_bub2 = np.append(new_bub2,-1)
    # xs1 = xs1[new_lags_cut]
    # ys1 = ys1[new_lags_cut]
    # xs2 = xs2[new_lags_cut]
    # ys2 = ys2[new_lags_cut]
    # new_bub1=new_bub1[new_lags_cut]
    # new_bub2=new_bub2[new_lags_cut]
    # eventsoi=eventsoi[new_lags_cut]
    fig, axes = plt.subplots(4, 2)
    axes[0,0].set_title('Click on point to plot Piezo 1 trace \n press N for next point or P for previous')
    axes[0,1].set_title('Click on point to plot Piezo 2 trace \n press N for next point or P for previous')
    old_f = plt.gcf()
    plt.ioff()
    plt.figure()
    plt.plot(xs2, 1e6*ys2, "bo", picker=5)
    #x1, y1 = (-2.84015, -9.9795e1)
    #x2, y2 = (0.0237464, -3.85555e1)
    x1, y1= (-2.56452, -170.646)
    x2, y2 = (-0.247984, -195.688)
    slope = (y2-y1)/(x2-x1)/1e6
    print(slope)
    soundspeed = 1/(slope*100)
    print(soundspeed)
    xar = np.array([x1, x2])
    yar = np.array([y1, y2])
    plt.plot(xar, yar, color="darkorange", linewidth=4)

    plt.ylim(-.0005e6,.0002e6)
    plt.xlabel("Z-position (cm)")
    plt.ylabel("t0 lag between PMT and acoustic analysis (us)")
    plt.show()

    # plt.figure(old_f.number)
    line, = axes[0,0].plot(xs1, ys1, 'bo', picker=5)  # 5 points tolerance
    line2, = axes[0,1].plot(xs2, ys2, 'bo', picker=5)


    # axbox = plt.axes([0.1, 0.05, 0.8, 0.075])nnn
    # text_box = TextBox(axbox, 'Evaluate', initial='')
    axes[0,0].set_ylim(-.0005,.0002)
    axes[0,1].set_ylim(-.0005,.0002)
    browser = PointBrowser()
    axes[1, 0].set_title(
        'Piezo 1 trace. Red line is acoustic t0')
    axes[1, 1].set_title(
        'Piezo 2 trace. Red line is acoustic t0')
    axes[2, 0].set_title(
        'Zoomed in on t0. Blue lines are close PMT triggers, purple is selected')
    axes[2, 1].set_title(
        'Zoomed in on t0. Blue lines are close PMT triggers, purple is selected')
    axes[3, 0].set_title(
        'Selected PMT trace. Press R for next trigger in time, L for previous')
    axes[3, 1].set_title(
        'Selected PMT trace. Press R for next trigger in time, L for previous')
    fig.canvas.mpl_connect('pick_event', browser.onpick)
    fig.canvas.mpl_connect('key_press_event', browser.onpress)
    figManager = plt.get_current_fig_manager()
    plt.show()