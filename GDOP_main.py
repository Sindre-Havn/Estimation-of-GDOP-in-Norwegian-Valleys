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
import urllib.request
import json
import os
from common import EPHEMERIS_FOLDER, ARCGIS_DATA_FOLDER, SKYLINE_GRAPHS_FOLDER
from pathlib import Path

# Select GNSS constellations
USE_GPS     = True
USE_GALILEO = True
USE_GLONASS = False
USE_BEIDOU  = False

# Select time-stamps for calculating GDOP
TIME_START = datetime.now()
TIME_DURATION = timedelta(
                      days=0,
                      hours=1,
                      minutes=0,
                      seconds=0 )
TIME_STOP = TIME_START + TIME_DURATION
TIME_STEP = timedelta( hours=0, minutes=5, seconds=0 )


def merge(dict1, dict2):
    dict1.update(dict2)

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
    sat_const_ephem = dict()
    if USE_GPS:
        gps_ephem = load_sat_const_ephem('gps')
        merge(sat_const_ephem, gps_ephem)
    if USE_GALILEO:
        galileo_ephem = load_sat_const_ephem('galileo')
        merge(sat_const_ephem, galileo_ephem)
    if USE_GLONASS:
        glonass_ephem = load_sat_const_ephem('glonass')
        merge(sat_const_ephem, glonass_ephem)
    if USE_BEIDOU:
        beidou_ephem = load_sat_const_ephem('beidou')
        merge(sat_const_ephem, beidou_ephem)
    return sat_const_ephem

def generate_time_interval(start, end, step):
    sample_times = pd.date_range(start, end, freq=step)
    return sample_times

def configure_plot(ax):
    FILENAME = 'skyplot'
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

def is_inside_triangle_sector(x1, y1, x2, y2, x3, y3, x, y):
    """Return True if point (x,y) is inside the triangle 
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
    #print('DIFF:',diff)
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


def get_single_gdop(sat_rel_pos_list):
        if len(sat_rel_pos_list) < 4:
            return None

        mat: np.array = []
        # d(0) - satellite x-pos
        # psd - distance from observer to sat
        # Calc vis sats is list of sat names
        for pos in sat_rel_pos_list:
            psd = pos[3]
            row = [-pos[0] / psd, -pos[1] / psd, -pos[2] / psd, 1]
            mat.append(row)

        m = np.matmul(np.transpose(mat), mat)
        Q = np.linalg.inv(m)
        T = np.trace(Q)
        print('HDOP:', np.sqrt(Q[0][0]**2 + Q[1][1]**2), 'VDOP:', Q[2][2], 'TDOP:', Q[3][3])
        G = np.sqrt(T)
        return G

def load_skylines(skyline_count):
    angle2zenith_dicts = []
    for i in range(1,skyline_count+1):
        try:
            location = Path(SKYLINE_GRAPHS_FOLDER / f"angles_table{i}.csv")
            angles_table_df = pd.read_csv(location, sep=';', decimal=",") # ../
            angles_table_df['HOR_AN_GEO'] = angles_table_df['HOR_AN_GEO'].apply(np.round).apply(int)
            skyline_dict = angles_table_df.set_index('HOR_AN_GEO')['ZENITH_ANG'].to_dict()
            angle2zenith_dicts.append(skyline_dict)
        except FileNotFoundError:
            angle2zenith_dicts.append(None)
    return angle2zenith_dicts

def load_obs_points():
    location = Path(ARCGIS_DATA_FOLDER / 'observation_points.csv')
    point_df = pd.read_csv(location, sep=';', decimal=",") # ../
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
    inside = is_inside_triangle_sector(0, 0,
                           sect_vec_lower[0], sect_vec_lower[1],
                           sect_vec_upper[0], sect_vec_upper[1],
                           sat_vec[0], sat_vec[1])
    return inside

def plot_sat(sat_is_visible, ax, name, theta, r, idx):
    FONTSIZE = 8
    if sat_is_visible:
        ax.plot(theta, r, 'o', c=(0.5+0.5/len(times)*idx, 0.0, 0.0, 0.5))
        ax.text(theta, r, name, fontsize=FONTSIZE, ha='right', va='bottom')
    else:
        ax.plot(theta, r, 'o', c=(0.1, 0.1, 0.1, 0.5))
        ax.text(theta, r, name, fontsize=FONTSIZE, ha='right', va='bottom')

def calc_DOP_trough_time(observer, skyline, sat_ephem, times):
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'},
                           figsize=(6, 6))
    ax = configure_plot(ax)
    skyline_circle_angle = np.deg2rad(list(skyline.keys()))
    plt.fill_between(skyline_circle_angle, 0, skyline.values(), alpha=0.2)
    
    idx = 0
    ts = load.timescale()
    sat_objects = {sat: EarthSatellite(sat_ephem[sat][0], sat_ephem[sat][1], sat, ts) for sat in sat_ephem}
    for time in times:
        # Load timescale and observer position
        t = ts.utc(time.year, time.month, time.day, time.hour, time.minute, time.second)
        # Initialize satellite objects
        visible_sats_rel_pos = []
        for sat in sat_objects:
            difference = sat_objects[sat] - observer
            topocentric = difference.at(t)
            alt, az, distance = topocentric.altaz()

            theta = np.radians(az.degrees)
            name = f'{sat}' if idx in (0, len(times)) else ''
            r = 90 - alt.degrees
            circle_sector = az.degrees
            sat_is_visible = is_visible(r, circle_sector)
            if sat_is_visible:
                x,y,z = polar2cart(distance.m, theta, r)
                visible_sats_rel_pos.append((x,y,z,distance.m))
                #print(sat, r, theta, zenith1, zenith2, angle_lower, angle_upper)
            plot_sat(sat_is_visible, ax, name, theta, r, idx)
        #print('FIX', observer_pos, visible_sats)
        GDOP = get_single_gdop(visible_sats_rel_pos)
        print(t, GDOP)
        idx += 1
    plt.show()


if __name__ == '__main__':
    obs_points_dict, obs_point_count = load_obs_points()
    angle2zenith_dicts = load_skylines(obs_point_count)
    sat_ephem = load_GNSS_data(USE_GPS, USE_GALILEO, USE_GLONASS, USE_BEIDOU)
    times = generate_time_interval(TIME_START, TIME_STOP, TIME_STEP)
    global GDOP_list
    GDOP_list = []
    
    for i in range(obs_point_count):
        obs_pt = obs_points_dict[i]
        skyline = angle2zenith_dicts[i]
        if skyline is None:
            GDOP_list.append(None)
            continue
        observer = Topos(latitude_degrees=obs_pt['LAT'], longitude_degrees=obs_pt['LONG'], elevation_m=obs_pt['ALT'])
        calc_DOP_trough_time(observer, skyline, sat_ephem, times)