from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent

# Material data
CLADDING_INDEX_DATA = Path(
   ...
)

WG_CORE_INDEX_DATA = Path(
    ...
)

DATA_FILE_LOCATION=  Path(...)
# Material names

WG_CORE_MATERIAL = "core"

WG_CLADDING_MATERIAL = "Cladding"



#optical parameters


WAVELENGTH = 1.55e-6 #central
BANDWIDTH = 0.1e-6
SPEED_OF_LIGHT = 299792458.0



#waveguide geometry

WG_WIDTH = 3.2e-6
WG_HEIGHT = 3.2e-6
CLADDING_THICKNESS = 4*WG_WIDTH
STRAIGHT_SECTION_LENGTH = 20e-6

#taper parameters
TAPER_WIDTH= 8e-6
TAPER_LENGTH= 220e-6

#substrate rules
SUBSTRATE_MATERIAL= "Si (Silicon) - Palik"
SUBSTRATE_TOP = 0e-6
SUBSTRATE_THICKNESS = 10e-6

# Bend_radius
MIN_BEND_RADIUS = 1000e-6



