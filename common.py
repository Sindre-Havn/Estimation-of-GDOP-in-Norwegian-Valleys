from pathlib import Path
import os

BASE_FOLDER = Path(__file__).resolve().parent
EPHEMERIS_FOLDER = Path(BASE_FOLDER / "gnss_ephemeris_data")
ARCGIS_DATA_FOLDER = Path(BASE_FOLDER / "arcgis_data")
SKYLINE_GRAPHS_FOLDER = Path(ARCGIS_DATA_FOLDER / "skyline_graphs")
DOP_OUTPUT_CSV = Path(BASE_FOLDER / "DOP.csv")

if not os.path.isdir(EPHEMERIS_FOLDER):
    os.mkdir(EPHEMERIS_FOLDER)