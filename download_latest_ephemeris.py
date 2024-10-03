from common import EPHEMERIS_FOLDER
import datetime
import urllib
import urllib.request
import json
import os

GNSS_URLS = {
    'gps'     : 'https://celestrak.org/NORAD/elements/gp.php?GROUP=gps-ops&FORMAT=tle',
    'galileo' : 'https://celestrak.org/NORAD/elements/gp.php?GROUP=galileo&FORMAT=tle',
    'glonass' : 'https://celestrak.org/NORAD/elements/gp.php?GROUP=glo-ops&FORMAT=tle',
    'beidou'  : 'https://celestrak.org/NORAD/elements/gp.php?GROUP=beidou&FORMAT=tle'
}

def get_sat_const_ephem(url) -> dict:
    """Takes inn a url from https://celestrak.org/NORAD/elements/
    and converts the response to a dictionary of and TLE formated ephemeris data."""
    with urllib.request.urlopen(url) as response:
        res = response.read()
        formated = res.split(b'\r\n')
        strings = [line.decode("utf-8") for line in formated]
        sat_const_ephem = dict()
        for i in range(0,len(strings)-2,3):
            key = strings[i].rstrip()
            value = [ strings[i+1], strings[i+2] ]
            sat_const_ephem[key] = value
        return sat_const_ephem

def main():
    current_time = datetime.datetime.now().strftime('%Y.%m.%d-%H.%M')
    os.chdir(EPHEMERIS_FOLDER)
    old_files_in_folder = os.listdir()
    for gnss in GNSS_URLS:
        gnss_ephem = get_sat_const_ephem(GNSS_URLS[gnss])
        new_filename = f'{gnss}_{current_time}.json'
        if new_filename in old_files_in_folder:
            print('The ephemeris data has already been renewed this minute.')
            return
        with open(new_filename, 'w') as file:
            file.write(
                json.dumps(gnss_ephem, indent=4)
            )
    for old_file in old_files_in_folder:
        os.remove(old_file)


if __name__ == '__main__':
    main()


