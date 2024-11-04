# Logger Details
# -------------------------------------------------
LOGFILE_DIR = "./logs"

# LOG_FORMAT = "%(asctime)s [%(levelname)s] File %(module)s: line %(lineno)d: %(message)s"
# LOG_FORMAT = "%(asctime)s [%(levelname)s] File %(pathname)s: line %(lineno)d:\n%(message)s"
# LOG_FORMAT = "%(asctime)s [%(levelname)s] File %(pathname)s: line %(lineno)d: %(message)s"
LOG_FORMAT = "%(asctime)s [%(levelname)s] %(filename)s: line %(lineno)d: %(message)s"

LOG_SEPARATOR = f"\n\n\n{'-'*50}\n\n\n"



# Data Details
# -------------------------------------------------
REQUIRED_DATA_COLUMNS = [
    'COLLECTION_DATE',
    'STATION_ID',
    'ANALYTE_NAME',
    'RESULT',
    'RESULT_UNITS',
]

REQUIRED_CONSTRUCTION_DATA_COLUMNS = [
    'STATION_ID',
    'AQUIFER',
    # 'WELL_USE',
    'STATION_USE',
    'LATITUDE',
    'LONGITUDE',
    'GROUND_ELEVATION',
    'TOTAL_DEPTH',
]

COLLECTION_DATE_FORMAT = "%m/%d/%y"
COLLECTION_TIME_FORMAT = "%H:%M"



# Coordinates information
# -------------------------------------------------
DEFAULT_SOURCE_COORDINATES = (436642.70, 3681927.09)




# Requirement Messages
# -------------------------------------------------
PYLENM_DATA_REQUIREMENTS = f"""
PYLENM DATA REQUIREMENTS:
The imported data needs to meet ALL of the following conditions to have a successful import:
    1) Data should be a pandas dataframe.
    2) Data must have these column names: {REQUIRED_DATA_COLUMNS}
"""

PYLENM_CONSTRUCTION_DATA_REQUIREMENTS = f"""
PYLENM CONSTRUCTION DATA REQUIREMENTS:
The imported data needs to meet ALL of the following conditions to have a successful import:
    1) Data should be a pandas dataframe.
    2) Data must have these column names: {REQUIRED_CONSTRUCTION_DATA_COLUMNS}
"""

