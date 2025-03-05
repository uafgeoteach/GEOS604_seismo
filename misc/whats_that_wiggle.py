"""
Download and view seismic waveforms used in GEOS604: Seismology to train
students in seeing waveforms, identifying key features, and diagnosing
station and or source characteristics.

Things I want to train students to see, be able to do:
    - Main phases: P, S, Love, Rayleigh
    - Earthquake distance with P-minus-S
    - Surface wave dispersion
    - Love and Rayleigh wave horizontal polarization
    - Rayleigh wave ellipticity
    - Surface waves from shallow vs. deep events
    - Multi-orbit surface waves
    - Data glitches

To Do List:
    - Earthquake triangulation with multiple P-minus-S times
    - Earthquake location with P-phases only
    - Polarization of P and S waves
    - Shear wave splitting

..rubric::
    
    To run this you can run the code with python with the first argument being
    the KEY you want

    $ python whats_that_wiggle.py LE18_MDM

    Some events will provide arrival times on the waveforms, you just need 
    to add a second argument (can be anything)

    $ python whats_that_wiggle.py LE18_MDM 1
"""
import os
import sys
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import obspy

from obspy import UTCDateTime, read
from obspy.taup import TauPyModel
from obspy.clients.fdsn import Client
from obspy.geodetics import locations2degrees


# CONSTANT SET
DATA_DIR = "./"
OVERWRITE = False# True
PHASE_LIST = ["P", "S", "Pn", "Sn", "pP", "PP", "SS", "sS"]
WIGGLES = {
    # KEY    STARTTIME              ENDTIME                CODE           OUTPUT FILE        
    # Leilani Estates to Murphy Dome
    "LE18_MDM": ["2018-05-04T22:40:00", "2018-05-04T23:00:00", "AK.MDM..BH?", "V",  (1, None)],  
    "LE18_MDM_2": ["2018-05-04T22:32:54", "2018-05-05T03:00:00", "AK.MDM..BH?", "V",  (0.05, None)],  
    "LE18_NEA2": ["2018-05-04T22:40:00", "2018-05-04T23:00:00", "AK.NEA2..BH?", "V",  (1, None)],  
    # Leilani Estates to Homer
    "LE18_HOM": ["2018-05-04T22:38:00", "2018-05-04T22:58:00", "AK.HOM..BH?", "V",  (1, None)],  
    # Mineral Springs
    "MS11_CCB_1": ["2011-08-23T17:58:00", "2011-08-23T18:25:00", "AK.CCB..BH?", "V",  (1, None)],  
    "MS11_CCB_2": ["2011-08-23T17:58:00", "2011-08-23T18:25:00", "AK.CCB..BH?", "V",  (0.01, 0.1)],  
    # Denali Fault
    "DF02_COLA": ["2002-11-03T22:12:41", "2002-11-03T22:32:41", "IU.COLA.10.BH?", "V", (None, None)],
    "DF02_KDAK": ["2002-11-03T22:12:41", "2002-11-03T22:32:41", "II.KDAK.00.BH?", "V", (1, None)],
    # NK6
    "NK6_COLA": ["2017-09-03T03:38:00", "2017-09-03T04:07:00", "IU.COLA.00.BH?", "V", (1, None)],
    # Japan to Cola (deep and shallow)
    "BI15_COLA": ["2015-05-30T11:30:00", "2015-05-30T12:20:00", "IU.COLA.00.BH?", "V", (1, None)],
    "NT24_COLA": ["2024-01-01T07:15:00", "2024-01-01T08:05:00", "IU.COLA.00.BH?", "V", (1, None)],
            }

EVENTS = {
    # KEY        ORIGINTIME             LAT_EV  LON_EV    LAT_STA LON_STA   DEPTH
    "LE18_MDM": ["2018-05-04T22:32:54", 19.318, -155.000,  64.960, -148.231, 5.8],
    "LE18_HOM": ["2018-05-04T22:32:54", 19.318, -155.000,  59.657, -151.652, 5.8],
    "DF02_KDAK": ["2002-11-03 22:12:41", 63.514, -147.453, 56.782, -152.384, 4.2],
        }


def get_data(starttime, endtime, code, output="VEL", buffer=2 * 60):
    """
    Download waveform data from IRIS
    """
    fid = f"{starttime}_{code}"
    full_path = os.path.join(DATA_DIR, fid)

    if os.path.exists(full_path) and not OVERWRITE:
        st = read(full_path)
    else:
        c = Client("IRIS")

        starttime -= buffer
        endtime += buffer
        network, station, location, channel = code.split(".")

        st = c.get_waveforms(
                network=network, station=station, location=location, 
                channel=channel, starttime=starttime, endtime=endtime,
                attach_response=True
                )

        inv = c.get_stations(
                network=network, station=station, location=location, 
                channel=channel, starttime=starttime, endtime=endtime,
                level="channel"
                )

        # Basic processing prior to instrument response
        st.detrend("linear")
        st.filter("highpass", freq=1/120)
        st.remove_response(output=output)
        st.rotate("->ZNE", inventory=inv)

    return st


def process_data(st, filt):
    """
    Filter the data
    """
    low, high = filt

    # Remove any lingering trends before filter
    st.detrend("polynomial", order=2)
    st.taper(0.05)

    if high and low:
        st.filter("bandpass", freqmin=low, freqmax=high)
    elif high:
        st.filter("highpass", freq=high)
    elif low:
        st.filter("lowpass", freq=low)
    else:
        pass

    return st


def get_arrival_times(event):
    """
    Get TauP traveltimes
    """
    origintime, ev_lat, ev_lon, sta_lat, sta_lon, depth = EVENTS[event]
    dist_deg = locations2degrees(ev_lat, ev_lon, sta_lat, sta_lon)
    print(f"{dist_deg}deg", end="...")
    model = TauPyModel(model="PREM")
    arrivals = model.get_travel_times(source_depth_in_km=depth, 
                                      distance_in_degree=dist_deg,
                                      phase_list=PHASE_LIST)
    arrival_times = {}
    origintime = UTCDateTime(origintime)
    # Get actual date time for arrival, rather than relative arrival
    for arrival in arrivals:
        arrival_times[arrival.name] = origintime + arrival.time

    return arrival_times
    

def plot_waveforms(st):
    """
    Generic Obspy Stream waveform plotter for use in GEOS419 Time Series Lab

    :type st: obspy.core.stream.Stream
    :param st: Stream object to plot
    :type fmt: str
    :param fmt: x-axis label formatter, 'absolute' for timestamps, 'relative' 
        for starttime t=0
    """
    f, axs = plt.subplots(len(st), figsize=(24, 10), sharex=True, sharey=True)
    plt.subplots_adjust(hspace=0.025)

    # X-axis formatting for datetime objects
    axs[-1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M:%S"))  # %H:%M:%S.%f
    axs[-1].tick_params(axis="x", labelrotation=20)

    # Set regular timing intervals rather than whatever is given by starttime
    # If traces are longer than 1 hour, change how the ticks are made
    if st[0].times().max() / 60 / 60 > 1:
        xlocator = mdates.MinuteLocator(byminute=np.arange(0, 61, 10))
        axs[-1].xaxis.set_major_locator(xlocator)
        xlocator = mdates.MinuteLocator(byminute=np.arange(0, 61, 1))
        axs[-1].xaxis.set_minor_locator(xlocator)
    else:
        xlocator = mdates.MinuteLocator(byminute=np.arange(0, 61, 1))
        axs[-1].xaxis.set_major_locator(xlocator)
        xlocator = mdates.SecondLocator(bysecond=[15, 30, 45])
        axs[-1].xaxis.set_minor_locator(xlocator)

    # Make copies so we can sort the waveforms in-place
    st_plot = st.copy()
    st_plot.sort()
    time = st_plot[0].times("matplotlib")  # assuming all have the same time axis

    # Loop through each trace and plot
    for i, tr in enumerate(st_plot):
        # Trace normalize
        axs[i].plot(time, tr.data, label=tr.get_id(), c="k", 
                    lw=0.75, zorder=10)
        axs[i].legend(loc="upper left", frameon=False)
        # Remove y-axis but add a line at y=0
        # axs[i].get_yaxis().set_visible(False)
        axs[i].grid(visible=True, which="major", linestyle="--", alpha=0.5)
        axs[i].grid(visible=True, which="minor", linestyle=":", alpha=0.5)

    # General Plot Aesthetics
    axs[-1].set_xlabel("Time [s]")

    # !!! Not generic if we change output
    axs[int((len(st)-1)/2)].set_ylabel("Velocity [m/s]")

    plt.xlim([time.min(), time.max()])
    plt.tight_layout()

    return f, axs


if __name__ == "__main__":
    user_input = sys.argv[1]
    try:
        plot_arrivals = bool(sys.argv[2])
    except IndexError:
        plot_arrivals = False

    assert(user_input in WIGGLES), f"No key {user_input} in data dictionary"

    if not os.path.exists(DATA_DIR):
        os.makedirs(DATA_DIR)

    starttime, endtime, code, output, filt = WIGGLES[user_input]
    starttime = UTCDateTime(starttime)
    endtime = UTCDateTime(endtime)

    output = {"D": "DISP", "V": "VEL", "A":"ACC"}[output]

    print("get", end="...")
    st = get_data(starttime, endtime, code, output)
    st.write(os.path.join(DATA_DIR, user_input), format="MSEED")
    print("proc", end="...")
    st = process_data(st, filt)
    st.trim(starttime, endtime)
    print("plot", end="...")
    f, axs = plot_waveforms(st)
    if plot_arrivals:
        print("arr", end="...")
        arrival_times = get_arrival_times(user_input)
        for phase, time in arrival_times.items():
            for ax in axs:
                ymin, ymax = ax.get_ylim()
                ax.axvline(time.matplotlib_date, c="r", alpha=0.5)
                ax.text(time.matplotlib_date, ymax * 2/3, s=phase, c="r")

    # Figure out how to make the title
    low, high = filt
    if low and high:
        title = f"{low}-{high}Hz Bandpass Filter"
    elif low:
        title = f"{low}Hz Lowpass Filter"
    elif high:
        title = f"{high}Hz Highpass Filter"
    else:
        title = "Raw Waveform"

    axs[0].set_title(title)
    plt.savefig(os.path.join(DATA_DIR, f"{user_input}.png"))
    plt.show()

