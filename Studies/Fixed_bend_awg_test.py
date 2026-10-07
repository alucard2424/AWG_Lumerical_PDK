import ansys.lumerical.core as lumapi
import sys
from pathlib import Path
import numpy as np



sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pdk import create_material_from_data, get_index_data, FPR, draw_port_marker, fixed_bend_radius_route_wg, optimize_awg_routes_after_port0_with_fixed_bend_radius, draw_awg_routes_fixed_radius#, AWG_full_geometry, AWG_array_geometry, aperture
from pdk.config import DATA_FILE_LOCATION, MIN_BEND_RADIUS, STRAIGHT_SECTION_LENGTH, SUBSTRATE_TOP, CLADDING_THICKNESS


from pdk.config import (
    WG_WIDTH,
    WG_HEIGHT,
    WG_CORE_MATERIAL,
    WG_CLADDING_MATERIAL,
    CLADDING_INDEX_DATA,
    WG_CORE_INDEX_DATA,
)


retrieval_location = DATA_FILE_LOCATION / "Mode_results.csv"

slab_data = get_index_data("slab",retrieval_location)
waveguide_data = get_index_data("waveguide",   retrieval_location)
neff_slab = float(slab_data["Effective Index"].iloc[0])
ng_slab   = float(slab_data["Group Index"].iloc[0])

neff = float(waveguide_data["Effective Index"].iloc[0])
ng   = float(waveguide_data["Group Index"].iloc[0])

#

#design parameters
N_receivers=8
fpr_rotation_angle= 60 #degs
NARRAY = 100
LAMBDA_CENTER = 1550e-9
CHANNEL_SPACING = 2e-9
ORDER = 45
RA = 1200e-6
RMIN = 1000e-6
ARRAY_WAVEGUIDE_SPACING = 2e-6
NEFF = neff
NG = ng
NEFF_SLAB = neff_slab
def calculate_awg_lengths(
        Narray,
        order,
        wavelength,
        neff,
        base_length
):

    delta_L = order * wavelength / neff

    indices = np.arange(Narray)

    lengths = (
        base_length
        + indices * delta_L
    )

    return lengths, delta_L
delta_L = calculate_awg_lengths(NARRAY, ORDER, LAMBDA_CENTER, NEFF, STRAIGHT_SECTION_LENGTH)[1]


with lumapi.FDTD() as solver:
    solver.eval("redrawoff;")
    print(f"This is Lumerical FDTD " f"{solver.version()}" )
    create_material_from_data(WG_CORE_INDEX_DATA,WG_CORE_MATERIAL, solver,mesh_order=1)
    create_material_from_data(CLADDING_INDEX_DATA,WG_CLADDING_MATERIAL, solver, mesh_order=2)

    
    fpr_in = FPR(
        solver=solver,
        input_aperture_spacing=9e-6,
        output_aperture_spacing=9e-6,
        Nin=1,
        Narray=NARRAY,
        Radius_input_arc=500e-6,
        Radius_output_arc=1000e-6,
        FPR_length=100e-6,
        input_opening_angle=80,
        output_opening_angle=80,
        arc_points=100,
        group_name="FPR_IN",
        rotation_angle=fpr_rotation_angle

)
    fpr_out = FPR(
            solver=solver,
            input_aperture_spacing=9e-6,
            output_aperture_spacing=9e-6,
            Nin=N_receivers,
            Narray=NARRAY,
            Radius_input_arc=500e-6,
            Radius_output_arc=1000e-6,
            FPR_length=100e-6,
            input_opening_angle=80,
            output_opening_angle=80,
            arc_points=100,
            group_name="FPR_out",
            rotation_angle=fpr_rotation_angle,
            output_FPR=True,
            FPR_separation=3000e-6

    )
    input_ports = fpr_in["positions"]
    input_angles = fpr_in["angles"]
    output_ports = fpr_out["positions"]
    output_angles = fpr_out["angles"]
    insertion_port_in = fpr_in["insertion positions"]
    insertion_port_in_angles = fpr_in["insertion angles"]
    receiver_ports = fpr_out["insertion positions"]
    receiver_port_angles = fpr_out["insertion angles"]
    # print(f"first Input port: {input_ports[0]}, angle: {input_angles[0]}")
    # print(f"last Output port: {output_ports[-1]}, angle: {output_angles[-1]}")
    # print(f"Insertion port: {insertion_port_in[0]}, angle: {insertion_port_in_angles[0]}")
    # draw_port_marker(solver,insertion_port_in[0],"insertion_port_marker")
    # draw_port_marker(solver,receiver_ports[N_receivers-1],"receiver_port_marker")
    # draw_port_marker(solver,input_ports[0],"input_port_1_marker")
    # draw_port_marker(solver,output_ports[0],"output_port_1_marker" )
   
    left_ports = input_ports
    left_angles = input_angles

    right_ports = output_ports
    right_angles = output_angles
    route_0 = fixed_bend_radius_route_wg(
        solver=solver,
        input_port=input_ports[0],
        input_angle=input_angles[0],
        output_port=output_ports[0],
        output_angle=output_angles[0],
        bend_radius=MIN_BEND_RADIUS,
        delay_length=0.0,
        base_straight_length=STRAIGHT_SECTION_LENGTH,
        WG_WIDTH=WG_WIDTH,
        WG_CORE_MATERIAL=WG_CORE_MATERIAL,
        SUBSTRATE_TOP=SUBSTRATE_TOP,
        CLADDING_THICKNESS=CLADDING_THICKNESS,
        WG_HEIGHT=WG_HEIGHT,
        target_input_angle=0.0,
        target_output_angle=np.pi,
    )
    
    
    
    L_base = route_0["total_length"]

    routes= optimize_awg_routes_after_port0_with_fixed_bend_radius(
        input_ports=input_ports,
        output_ports=output_ports,
        input_angles=input_angles,
        output_angles=output_angles,
        port0_geometry=route_0,
        delta_L=delta_L,
        WG_WIDTH=WG_WIDTH,
        min_bend_radius=MIN_BEND_RADIUS,
        min_wg_spacing=ARRAY_WAVEGUIDE_SPACING,
        
    )
    solver.eval("redrawoff;")
    draw_awg_routes_fixed_radius(solver,routes,input_ports,input_angles,output_ports,output_angles,WG_WIDTH,WG_CORE_MATERIAL,SUBSTRATE_TOP, CLADDING_THICKNESS,WG_HEIGHT, npoints=200, name_prefix="AWG")
    

    

    
    
    
    
    solver.addfdtd()
    final_route = routes[-1]
    end_point = final_route["input_bend_end"]
    y_max= end_point[1]+10e-6
    xbuffer=150e-6
    ybuffer= 10e-6
    solver.set("dimension", "2D")
    solver.set("x min",insertion_port_in[0][0]-xbuffer)
    solver.set("x max",receiver_ports[N_receivers-1][0]+xbuffer)
    solver.set("y min",receiver_ports[0][1]-ybuffer)
    solver.set("y max",y_max+ybuffer)
    solver.set("z",(CLADDING_THICKNESS+SUBSTRATE_TOP)/2)
    #solver.set("z span", 5e-6)
    solver.set("mesh accuracy", 1)
    solver.set("auto shutoff min", 1e-3)
    solver.set("auto shutoff max", 100)
    solver.set("background material", WG_CLADDING_MATERIAL)
    solver.set("simulation time", 40000e-15)



    solver.addmode()
    solver.set("injection axis","x")
    solver.set("x",insertion_port_in[0][0]+0.1e-6)
    solver.set("y",insertion_port_in[0][1]+0.1e-6)
    solver.set("y span",WG_WIDTH+6e-6)
    solver.set("z",(CLADDING_THICKNESS+SUBSTRATE_TOP)/2)
    solver.set("z span",WG_WIDTH+6e-6)
    solver.set("theta", fpr_rotation_angle)
    solver.set("wavelength start", 1.54e-6)
    solver.set("wavelength stop", 1.56e-6)
    solver.set("mode selection", 2) #TE mode



    


    

    print("\nAWG array construction complete.")
    input("Press any key to close session")
    