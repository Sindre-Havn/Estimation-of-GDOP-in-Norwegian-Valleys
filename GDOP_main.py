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


def isEven(n):
  # if n&1 == 0, then num is even
  if n & 1:
    return False
  # if n&1 == 1, then num is odd
  else:
    return True

def merge(dict1, dict2):
    dict1.update(dict2)

def get_sat_const_ephem(url, prefix_str):
    """Takes inn a url from https://celestrak.org/NORAD/elements/
    and converts the response to a dictionary with "prefix_str"(+number)
    as key, and TLE data from url as values."""
    with urllib.request.urlopen(url) as response:
        res = response.read()
        formated = res.split(b'\r\n')
        sat_const_ephem = dict()
        for i in range(0,len(formated)-2,3):
            key = prefix_str+formated[i].decode("utf-8").rstrip()[-4:-1]
            value = [ formated[i+1].decode("utf-8"), formated[i+2].decode("utf-8") ]
            sat_const_ephem[key] = value
        return sat_const_ephem

def get_GNSS_data(USE_GPS, USE_GALILEO, USE_GLONASS, USE_BEIDOU):
    if not(USE_GPS or USE_GALILEO or USE_GLONASS or USE_BEIDOU): # If all USE_cases are False
        raise ValueError('The function have to take in atleast one GNSS constellation, zero constellations given.')
    sat_const_ephem = dict()
    if USE_GPS:               
        gps_url = 'https://celestrak.org/NORAD/elements/gp.php?GROUP=gps-ops&FORMAT=tle'
        gps_prefix = 'G' # for 'GPS'
        gps_ephem = get_sat_const_ephem(gps_url, gps_prefix)
        merge(sat_const_ephem, gps_ephem)
    if USE_GALILEO:
        galileo_url = 'https://celestrak.org/NORAD/elements/gp.php?GROUP=galileo&FORMAT=tle'
        galileo_prefix = 'E' # for 'Europe'
        galileo_ephem = get_sat_const_ephem(galileo_url, galileo_prefix)
        merge(sat_const_ephem, galileo_ephem)
    if USE_GLONASS:
        glonass_url = 'https://celestrak.org/NORAD/elements/gp.php?GROUP=glo-ops&FORMAT=tle'
        glonass_prefix = 'R' # for 'Russia'
        glonass_ephem = get_sat_const_ephem(glonass_url, glonass_prefix)
        merge(sat_const_ephem, glonass_ephem)
    if USE_BEIDOU:
        beidou_url = 'https://celestrak.org/NORAD/elements/gp.php?GROUP=beidou&FORMAT=tle'
        beidou_prefix = 'B' # for 'Beidou'
        beidou_ephem = get_sat_const_ephem(beidou_url, beidou_prefix)
        merge(sat_const_ephem, beidou_ephem)
    return sat_const_ephem

def generate_time_interval(start, end, step):
    sample_times = pd.date_range(start, end, freq=step)
    #sample_times = np.array(sample_times)
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

def geo_polar_to_kartesian(deg, radius):
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

def load_skylines(rel_path, skyline_count):
    angle2zenith_dicts = []
    for i in range(1,skyline_count+1):
        angles_table_df = pd.read_csv(f'test_data/skyline_graphs/angles_table{i}.csv', sep=';', decimal=",") # ../
        angles_table_df['HOR_AN_GEO'] = angles_table_df['HOR_AN_GEO'].apply(np.round).apply(int)
        skyline_dict = angles_table_df.set_index('HOR_AN_GEO')['ZENITH_ANG'].to_dict()
        angle2zenith_dicts.append(skyline_dict)
    return angle2zenith_dicts

def load_obs_points(rel_path, skyline_count):
    point_df = pd.read_csv(f'test_data/new_E136_S8_17_pts_3D.csv', sep=';', decimal=",") # ../
    # Note 'POINT_X', 'POINT_Y', 'POINT_Z' should actually be long, lat, alt
    observation_points = point_df[['POINT_X', 'POINT_Y', 'POINT_Z']].to_dict('index')
    return observation_points

def plot_trough_time(observer, skyline, sat_ephem, times):
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'},
                           figsize=(6, 6))
    ax = configure_plot(ax)
    skyline_circle_angle = np.deg2rad(list(skyline.keys()))
    plt.fill_between(skyline_circle_angle, 0, skyline.values(), alpha=0.2)
    FONTSIZE = 8
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

            r = 90 - alt.degrees
            circle_sector = az.degrees
            angle_lower = int(np.floor(circle_sector))
            angle_upper = int(np.ceil(circle_sector))
            zenith1 = skyline[angle_lower]
            zenith2 = skyline[angle_lower]
            sect_vec_lower = geo_polar_to_kartesian(angle_lower, zenith1)
            sect_vec_upper = geo_polar_to_kartesian(angle_upper, zenith2)
            sat_vec = geo_polar_to_kartesian(circle_sector, r)
            inside = is_inside_triangle_sector(0, 0,
                                   sect_vec_lower[0], sect_vec_lower[1],
                                   sect_vec_upper[0], sect_vec_upper[1],
                                   sat_vec[0], sat_vec[1])
            theta = np.radians(az.degrees)
            name = f'{sat}' if idx in (0, len(times)) else ''
            if inside:
                x,y,z = polar2cart(distance.m, theta, r)
                visible_sats_rel_pos.append((x,y,z,distance.m))
                #print(sat, r, theta, zenith1, zenith2, angle_lower, angle_upper)
                ax.plot(theta, r, 'o', c=(0.5+0.5/len(times)*idx, 0.0, 0.0, 0.5))
                ax.text(theta, r, name, fontsize=FONTSIZE, ha='right', va='bottom')
            else:
                ax.plot(theta, r, 'o', c=(0.1, 0.1, 0.1, 0.5))
                ax.text(theta, r, name, fontsize=FONTSIZE, ha='right', va='bottom')
        #print('FIX', observer_pos, visible_sats)
        GDOP = get_single_gdop(visible_sats_rel_pos)
        print(t, GDOP)
        idx += 1
    plt.show()


if __name__ == '__main__':
    skyline_path = None
    skyline_count = 20
    angle2zenith_dicts = load_skylines(skyline_path, skyline_count)
    obs_pts_path = None
    obs_points_dict = load_obs_points(obs_pts_path, skyline_count)
    
    for i in range(20):#skyline_count):
        obs_pt = obs_points_dict[i]
        skyline = angle2zenith_dicts[i]
        sat_ephem = get_GNSS_data(USE_GPS, USE_GALILEO, USE_GLONASS, USE_BEIDOU)
        observer = Topos(latitude_degrees=obs_pt['POINT_Y'], longitude_degrees=obs_pt['POINT_X'], elevation_m=obs_pt['POINT_Z'])
        times = generate_time_interval(TIME_START, TIME_STOP, TIME_STEP)
        plot_trough_time(observer, skyline, sat_ephem, times)

    #lat, lon, alt = 62.553035, 7.675825, 4.9081
    #print(lla2ecef(lat, lon, alt))

