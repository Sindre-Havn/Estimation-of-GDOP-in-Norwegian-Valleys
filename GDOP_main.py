"""
  ____   ____     ___    ____    _____   ____      _                       
 / ___| |  _ \   / _ \  |  _ \  | ____| |  _ \    | |__    _   _           
| |  _  | | | | | | | | | |_) | |  _|   | |_) |   | '_ \  | | | |          
| |_| | | |_| | | |_| | |  __/  | |___  |  _ <    | |_) | | |_| |          
 \____| |____/   \___/  |_|     |_____| |_| \_\  _|_.__/   \__, |          
/ ___|  (_)  _ __     __| |  _ __    ___    | | | |   __ _ |___/ __  _ __  
\___ \  | | | '_ \   / _` | | '__|  / _ \   | |_| |  / _` | \ \ / / | '_ \ 
 ___) | | | | | | | | (_| | | |    |  __/   |  _  | | (_| |  \ V /  | | | |
|____/  |_| |_| |_|  \__,_| |_|     \___|   |_| |_|  \__,_|   \_/   |_| |_|


Main sources:
* https://s-taka.org/en/making-skyplot-with-python/
* https://trepo.tuni.fi/bitstream/handle/10024/132157/TampierJaraFelipe.pdf?sequence=2&isAllowed=y

Nice to know:
* Dont confuse geo_angle (clockwise, 0 degree at north) with
  'math' angle (cunterclockwise, 0 degree at east).

"""
from skyfield.api import Topos, load, wgs84
from skyfield.sgp4lib import EarthSatellite
import numpy as np
from matplotlib import pyplot as plt
from datetime import datetime, timedelta
import pandas as pd
import json
import os
from common import EPHEMERIS_FOLDER, ARCGIS_DATA_FOLDER, SKYLINE_GRAPHS_FOLDER, DOP_RESULTS_FOLDER
from pathlib import Path
import pytz

# Select GNSS constellations
USE_GPS     = True
USE_GALILEO = True
USE_GLONASS = False
USE_BEIDOU  = False

PLOT_SAT = True
USE_WDOP = True # Not finished
UERE = {'gps':1.9, 'galileo':1.8, 'glonass':2.8, 'beidou':1.7}

# Select time-stamps for calculating GDOP
UTC_DIFFERENCE = timedelta(hours=2) # Since Norway is UTC+2. This should be improved to solve for winter/summer time.
TIME_START = datetime.now() - UTC_DIFFERENCE # Normalize to UTC time
# TIME_START = datetime(year=2024, month=10, day=10, hour=11) # Set custom time
TIME_DURATION = timedelta(
                      days=0,
                      hours=2,
                      minutes=0,
                      seconds=0 )
TIME_STOP = TIME_START + TIME_DURATION 
TIME_STEP = timedelta( hours=0, minutes=5, seconds=0)

def load_sat_const_ephem(gnss) -> dict:
    """Takes in a gnss name and load the corresponding ephemeris .json file."""
    folder = os.listdir(EPHEMERIS_FOLDER)
    for filename in folder:
        if filename[:len(gnss)] == gnss:
            print(filename)
            location = Path(EPHEMERIS_FOLDER / filename)
            with open(location, 'r') as file:
                    data = json.load(file)
                    return data

def load_GNSS_data(USE_GPS, USE_GALILEO, USE_GLONASS, USE_BEIDOU):
    if not(USE_GPS or USE_GALILEO or USE_GLONASS or USE_BEIDOU): # If all USE_cases are False
        raise ValueError('The function have to take in atleast one GNSS constellation, zero constellations given.')
    sat_const_ephem = dict() # Jeg synes det er veldig gøy at han ved siden av hører på en kjemi 101 spilleliste med breaking bad sanger
    if USE_GPS:
        gps_ephem = load_sat_const_ephem('gps')
        sat_const_ephem['gps'] = gps_ephem
    if USE_GALILEO:
        galileo_ephem = load_sat_const_ephem('galileo')
        sat_const_ephem['galileo'] = galileo_ephem
    if USE_GLONASS:
        glonass_ephem = load_sat_const_ephem('glonass')
        sat_const_ephem['glonass'] = glonass_ephem
    if USE_BEIDOU:
        beidou_ephem = load_sat_const_ephem('beidou')
        sat_const_ephem['beidou'] = beidou_ephem
    return sat_const_ephem

def generate_time_interval(start, end, step):
    sample_times = pd.date_range(start, end, freq=step)
    print(sample_times)
    return sample_times

def configure_plot(ax):
    FONTSIZE = 8
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.set_ylim(0, 90)  # Altitude 0° (horizon) to 90° (zenith)
    ax.set_rmax(90)
    ax.set_rticks([0, 30, 60, 90])  # Altitude rings
    ax.set_rlabel_position(-22.5)  # Move radial labels
    ax.set_rgrids([0, 30, 60, 90], labels=["90", "60", "30", "0"], fontsize=FONTSIZE)
    ax.grid(True)
    ax.set_thetalim([0, 2*np.pi])
    ax.set_thetagrids(np.rad2deg(np.linspace(0, 2*np.pi, 5)[1:]), 
                      labels=["E", "S", "W", "N"], fontsize=FONTSIZE)
    return ax

def geodeg2deg(deg):
    """Rotates the angle so it rotates counterclockwise
    from 0 degree in east,
    instead of clockwise with 0 degree in north."""
    ang = 90 - deg
    if ang < 0: ang += 360
    return ang

def geo_polar_to_cartesian(deg, radius):
    math_rad = np.deg2rad(geodeg2deg(deg))
    x = np.cos(math_rad)*radius
    y = np.sin(math_rad)*radius
    return x, y


def area(x1, y1, x2, y2, x3, y3):
    return abs((x1 * (y2 - y3) + x2 * (y3 - y1) 
                + x3 * (y1 - y2)) / 2.0)

def is_inside_triangle(x1, y1, x2, y2, x3, y3, x, y):
    """Used to approximate the zentih angle between two bearing angles,
    example between 203-204 degree, with 360x1 degree resolution.
    Return True if point (x,y) is inside the triangle 
    of the points (x1, y1), (x2, y2), (x3, y3)."""
    # Calculate area of triangle ABC
    A = area (x1, y1, x2, y2, x3, y3)
    # Calculate area of triangle PBC 
    A1 = area (x, y, x2, y2, x3, y3)
    # Calculate area of triangle PAC 
    A2 = area (x1, y1, x, y, x3, y3)
    # Calculate area of triangle PAB 
    A3 = area (x1, y1, x2, y2, x, y)
    # Check if sum of A1, A2 and A3 
    # is same as A
    diff = A - (A1 + A2 + A3)
    # Due to rounding error, check if the areas
    # is close enough to equal.
    if(round(diff, 10) == 0):
        return True
    else:
        return False

def polar2cart(r, theta, phi):
    return [
         r * np.sin(theta) * np.cos(phi),
         r * np.sin(theta) * np.sin(phi),
         r * np.cos(theta)
    ]


def calc_dop(gnss_sats_rel_pos):
        satelite_count_in_first_gnss = len(list(gnss_sats_rel_pos.values())[0])
        if satelite_count_in_first_gnss < 4:
            return None, None, None, None, None, None

        all_sats = []
        for const_sats in gnss_sats_rel_pos.values():
            for pos in list(const_sats):
                all_sats.append(pos)
        H: np.array = []
        # psd - distance from observer to sat
        # Calc vis sats is list of sat names
        for pos in all_sats:
            psd = pos[3]
            # Row is normalized vector, with absolute length equal 1.
            row = [-pos[0] / psd, -pos[1] / psd, -pos[2] / psd, 1]
            H.append(row)
        H = np.array(H)
        m = H.T @ H
        Q = np.linalg.inv(m)
        T = np.trace(Q)
        EDOP = np.sqrt(Q[0][0]) # East DOP
        NDOP = np.sqrt(Q[1][1]) # North DOP
        HDOP = np.sqrt(EDOP**2 + NDOP**2) # Horizontal DOP
        VDOP = Q[2][2]                    # Vertical DOP
        TDOP = Q[3][3]                    # Time DOP
        GDOP = np.sqrt(HDOP**2 + VDOP**2 + TDOP**2) # Geometric DOP
        return EDOP, NDOP, HDOP, VDOP, TDOP, GDOP


def calc_wdop(gnss_sats_rel_pos):
    satelite_count_in_first_gnss = len(list(gnss_sats_rel_pos.values())[0])
    if satelite_count_in_first_gnss < 4:
        return None, None, None, None, None, None

    all_sats = []
    for const_sats in gnss_sats_rel_pos.values():
        for pos in list(const_sats):
            all_sats.append(pos)
    H: np.array = []
    # psd - distance from observer to sat
    # Calc vis sats is list of sat names
    for pos in all_sats:
        psd = pos[3]
        # Row is normalized vector, with absolute length equal 1.
        row = [-pos[0] / psd, -pos[1] / psd, -pos[2] / psd, 1]
        H.append(row)
    H = np.array(H)
    # Create weights matrix W
    gnss_count = gnss_sats_rel_pos.keys()
    sat_UERE_vec = []
    W_side_length = 0
    for gnss in gnss_sats_rel_pos.keys():
        sat_count = len(gnss_sats_rel_pos[gnss])
        W_side_length += sat_count
        for sat in range(sat_count):
            sat_UERE_vec.append(1/ UERE[gnss]**2)
    W = np.zeros((W_side_length,W_side_length))
    for i in range(len(sat_UERE_vec)):
        W[i,i] = sat_UERE_vec[i]

    m = H.T @ W @ H
    Q = np.linalg.inv(m)
    T = np.trace(Q)
    EDOP = np.sqrt(Q[0][0]) # East DOP
    NDOP = np.sqrt(Q[1][1]) # North DOP
    HDOP = np.sqrt(EDOP**2 + NDOP**2) # Horizontal DOP
    VDOP = Q[2][2]                    # Vertical DOP
    TDOP = Q[3][3]                    # Time DOP
    GDOP = np.sqrt(HDOP**2 + VDOP**2 + TDOP**2) # Geometric DOP
    return EDOP, NDOP, HDOP, VDOP, TDOP, GDOP


def load_skylines(obs_point_count):
    angle2zenith_dicts = []
    for i in range(obs_point_count):
        try: # Bug in arcgis code fails to generate ~50% of skylines. Failed skyline is not written to files.
            location = Path(SKYLINE_GRAPHS_FOLDER / f"angles_table{i}.csv")
            angles_table_df = pd.read_csv(location, sep=';', decimal=",")
            angles_table_df['HOR_AN_GEO'] = angles_table_df['HOR_AN_GEO'].apply(np.round).apply(int)
            skyline_dict = angles_table_df.set_index('HOR_AN_GEO')['ZENITH_ANG'].to_dict()
            angle2zenith_dicts.append(skyline_dict)
        except FileNotFoundError:
            angle2zenith_dicts.append(None)
    return angle2zenith_dicts

def load_obs_points():
    location = Path(ARCGIS_DATA_FOLDER / 'observation_points.csv')
    point_df = pd.read_csv(location, sep=';', decimal=",")
    print(point_df)
    # Note 'POINT_X', 'POINT_Y', 'POINT_Z' should actually be long, lat, alt
    observation_points = point_df[['LONG', 'LAT', 'ALT']].to_dict('index')
    obs_point_count = len(point_df.index)
    return observation_points, obs_point_count

def is_visible(r, circle_sector):
    angle_lower = int(np.floor(circle_sector))
    angle_upper = int(np.ceil(circle_sector))
    zenith1 = skyline[angle_lower]
    zenith2 = skyline[angle_lower]
    sect_vec_lower = geo_polar_to_cartesian(angle_lower, zenith1)
    sect_vec_upper = geo_polar_to_cartesian(angle_upper, zenith2)
    sat_vec = geo_polar_to_cartesian(circle_sector, r)
    inside = is_inside_triangle(0, 0,
                           sect_vec_lower[0], sect_vec_lower[1],
                           sect_vec_upper[0], sect_vec_upper[1],
                           sat_vec[0], sat_vec[1])
    return inside

def plot_sat(sat_is_visible, ax, gnss, theta, r, idx, sat):
    letter_prefix = {'gps'   :'G',
                    'galileo':'E',
                    'glonass':'R',
                    'beidou' :'B'}
    identifier = sat[-4:-1].split(' ')[-1]
    name = letter_prefix[gnss]+identifier if idx == 0 else ''
    FONTSIZE = 8
    if sat_is_visible:
        ax.plot(theta, r, 'o', c=(0.2+0.8/len(times)*idx, 0.0, 0.0, 0.5))
        ax.text(theta, r, name, fontsize=FONTSIZE, ha='right', va='bottom')
    else:
        ax.plot(theta, r, 'o', c=(0.1, 0.1, 0.1, 0.5))
        ax.text(theta, r, name, fontsize=FONTSIZE, ha='right', va='bottom')

def calc_DOP_trough_time(observer, skyline, gnss_ephem, times):
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'},
                           figsize=(6, 6))
    ax = configure_plot(ax)
    skyline_circle_angle = np.deg2rad(list(skyline.keys()))
    plt.fill_between(skyline_circle_angle, 0, skyline.values(), alpha=0.2)
    idx = 0
    ts = load.timescale()
    for time in times:
        gnss_sats_rel_pos = dict()
        for gnss in gnss_ephem.keys():
            const_ephem = gnss_ephem[gnss]
            sat_objects = {sat: EarthSatellite(const_ephem[sat][0], const_ephem[sat][1], sat, ts) for sat in const_ephem}
            t = ts.utc(time.year, time.month, time.day, time.hour, time.minute, time.second)
            # Initialize satellite objects
            gnss_sats_rel_pos[gnss] = []
            for sat in sat_objects:
                difference = sat_objects[sat] - observer
                topocentric = difference.at(t)
                alt, az, distance = topocentric.altaz()
    
                theta = np.radians(az.degrees)
                r = 90 - alt.degrees
                circle_sector = az.degrees
                sat_is_visible = is_visible(r, circle_sector)
                if sat_is_visible:
                    x,y,z = polar2cart(distance.m, theta, r)
                    gnss_sats_rel_pos[gnss].append((x,y,z,distance.m))
                    #print(sat, r, theta, zenith1, zenith2, angle_lower, angle_upper)
                if PLOT_SAT: plot_sat(sat_is_visible, ax, gnss, theta, r, idx, sat)
        idx += 1
        if USE_WDOP:
            EDOP, NDOP, HDOP, VDOP, TDOP, GDOP = calc_wdop(gnss_sats_rel_pos)
        else:
            EDOP, NDOP, HDOP, VDOP, TDOP, GDOP = calc_dop(gnss_sats_rel_pos)
        append_at_last_idx = len(GDOP_dfs[time].index)
        GDOP_dfs[time].loc[append_at_last_idx] = [EDOP, NDOP, HDOP, VDOP, TDOP, GDOP] 
        print(t, GDOP)
    plt.show()

def append_None_to_DOP_rows():
    for key in GDOP_dfs:
        append_at_last_idx = len(GDOP_dfs[key].index)
        GDOP_dfs[key].loc[append_at_last_idx] = [None, None, None, None, None, None]

if __name__ == '__main__':
    old_results = os.listdir(DOP_RESULTS_FOLDER)
    if len(old_results) > 0:
        input('There are old files in results folder that will be DELETED,\npress [enter] to proceed anyway.\n>:')
        for file in old_results:
            location = Path(DOP_RESULTS_FOLDER / file)
            os.remove(location)

    obs_points_dict, obs_point_count = load_obs_points()
    angle2zenith_dicts = load_skylines(obs_point_count)
    gnss_ephem = load_GNSS_data(USE_GPS, USE_GALILEO, USE_GLONASS, USE_BEIDOU)
    times = generate_time_interval(TIME_START, TIME_STOP, TIME_STEP)
    global GDOP_dfs
    GDOP_dfs = dict()
    for time in times:
        GDOP_dfs[time] = pd.DataFrame(columns=['EDOP', 'NDOP', 'HDOP', 'VDOP', 'TDOP', 'GDOP'])
    
    for i in range(obs_point_count): # iterate trough observation points
        obs_pt = obs_points_dict[i]
        skyline = angle2zenith_dicts[i]
        if skyline is None: # Bug in arcgis code fails to generate ~50% of skylines.
            append_None_to_DOP_rows()
            continue
        observer = Topos(latitude_degrees=obs_pt['LAT'], longitude_degrees=obs_pt['LONG'], elevation_m=obs_pt['ALT'])
        calc_DOP_trough_time(observer, skyline, gnss_ephem, times)
    
    for time in times:
        time_str = time.strftime('%Y.%m.%d-%H.%M')
        filename = 'WDOP_'+time_str+'.csv' if USE_WDOP else 'DOP_'+time_str+'.csv'
        location = Path(DOP_RESULTS_FOLDER / filename)
        GDOP_dfs[time].to_csv(location)