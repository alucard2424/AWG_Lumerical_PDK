import ansys.lumerical.core as lumapi
import sys
from pathlib import Path
import tempfile
import shutil
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pdk import create_material_from_data, SM_straight_waveguide_geometry, slab_geometry, add_mesh_override, add_frequency_domain_monitor, FPR,draw_port_marker
from pdk.config import DATA_FILE_LOCATION, MIN_BEND_RADIUS, STRAIGHT_SECTION_LENGTH, SUBSTRATE_TOP, CLADDING_THICKNESS
from pdk.config import (
    WG_WIDTH,
    WG_HEIGHT,
    WG_CORE_MATERIAL,
    WG_CLADDING_MATERIAL,
    CLADDING_INDEX_DATA,
    WG_CORE_INDEX_DATA,
)

NARRAY=64
Rin= 500e-6
Rout =1000e-6
output_angle = 80


with lumapi.FDTD() as solver:
    solver.eval("redrawoff;")
    print(f"This is Lumerical FDTD "f"{solver.version()}")
    create_material_from_data(WG_CORE_INDEX_DATA,WG_CORE_MATERIAL, solver, mesh_order=1)
    create_material_from_data( CLADDING_INDEX_DATA, WG_CLADDING_MATERIAL, solver, mesh_order=2)
    
    #test a single FPR
    
    fpr = FPR(
            solver=solver,
            input_aperture_spacing=9e-6,
            output_aperture_spacing=9e-6,
            Nin=1,
            Narray=NARRAY,
            Radius_input_arc=Rin,
            Radius_output_arc=Rout,
            FPR_length=100e-6,
            input_opening_angle=output_angle,
            output_opening_angle=output_angle,
            arc_points=100,
            group_name="FPR_IN",
            rotation_angle=0

    )
    array_ports = fpr["positions"]
    array_angles = fpr["angles"]
    insertion_port_in = fpr["insertion positions"]
    insertion_port_in_angles = fpr["insertion angles"]
    print(array_ports)
    print()
    print(array_ports[NARRAY-1][1])
    print(array_ports[0][1])


    solver.addfdtd()
    xbuffer=50e-6
    ybuffer= 50e-6
    solver.set("dimension", "2D")
    solver.set("x min",insertion_port_in[0][0]-xbuffer)
    solver.set("x max",array_ports[0][0]+xbuffer)


    y_arc_min = -Rout * np.sin(np.deg2rad(output_angle) / 2)
    y_arc_max =  Rout * np.sin(np.deg2rad(output_angle) / 2)
    solver.set("y min", y_arc_min - ybuffer)
    solver.set("y max", y_arc_max + ybuffer)
    solver.set("z",(CLADDING_THICKNESS+SUBSTRATE_TOP)/2)
    #solver.set("z span", 5e-6)
    solver.set("mesh accuracy", 1)
    solver.set("auto shutoff min", 1e-3)
    solver.set("auto shutoff max", 100)
    solver.set("background material", WG_CLADDING_MATERIAL)
    solver.set("simulation time", 8000e-15)
    



    solver.addmode()
    solver.set("injection axis","x")
    solver.set("x",insertion_port_in[0][0]+0.1e-6)
    solver.set("y",insertion_port_in[0][1])
    solver.set("y span",WG_WIDTH+6e-6)
    solver.set("z",(CLADDING_THICKNESS+SUBSTRATE_TOP)/2)
    solver.set("z span",WG_WIDTH+6e-6)
    solver.set("theta", 0)
    solver.set("wavelength start", 1.54e-6)
    solver.set("wavelength stop", 1.56e-6)
    solver.set("mode selection", 2) #TE mode

    
    
    input("Press any key to close session")
    