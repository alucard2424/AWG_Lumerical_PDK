import numpy as np
from .config import WAVELENGTH, WG_HEIGHT, WG_WIDTH, SUBSTRATE_TOP, CLADDING_THICKNESS, BANDWIDTH
from .simulation_tools import add_frequency_domain_monitor



def slab_mode_solver(solver,
                      TE_fraction=1.0
                     ):

    solver.addfde()
    solver.set("solver type", "1D Z:X Prop")
    solver.set("x", 0)
    solver.set("y", 0)
    solver.set("z", SUBSTRATE_TOP + CLADDING_THICKNESS/2)
    solver.set("z span", CLADDING_THICKNESS )

    solver.setanalysis("wavelength", WAVELENGTH)
    solver.setanalysis("search", "near n")

    solver.findmodes()

    for i in range(1, int(solver.nummodes()) + 1):

        solver.selectmode(i)

        te_fraction = solver.getdata(f"mode{i}", "TE polarization fraction")

        if np.isclose(te_fraction, TE_fraction, atol=1e-3):

            neff = np.real(solver.getdata(f"mode{i}","neff")[0, 0])

            ng = np.real(solver.getdata(f"mode{i}","ng" )[0, 0])

            return {
                "mode": i,
                "te_fraction": te_fraction,
                "neff": neff,
                "ng": ng,
            }

    raise RuntimeError("No mode found for required TE fraction.")



def waveguide_mode_solver(
        solver,
        TE_fraction
        ):
    #add smallest effective area mode check to make sure no higher order model with te_frac=1 is chosen.
    solver.addfde()
    
    solver.set("solver type", "2D X Normal")

    solver.set("x", 0)
    solver.set("y", 0)
    solver.set('y span', CLADDING_THICKNESS )
    solver.set("z", SUBSTRATE_TOP + CLADDING_THICKNESS/2)
    solver.set("z span", CLADDING_THICKNESS )

    solver.setanalysis("wavelength", WAVELENGTH)
    solver.setanalysis("search", "near n")
    # solver.set("y min bc", 'PML')
    # solver.set("y max bc", 'PML')
    # solver.set("z min bc", 'PML')
    # solver.set("z max bc", 'PML')
    solver.set("y min bc", 'Anti-Symmetric')

    solver.findmodes()
    for i in range(1, int(solver.nummodes()) + 1):
    
            solver.selectmode(i)
    
            te_fraction = solver.getdata(f"mode{i}", "TE polarization fraction")
    
            if np.isclose(te_fraction, TE_fraction, atol=1e-3):
    
                neff = np.real(solver.getdata(f"mode{i}","neff")[0, 0])
    
                ng = np.real(solver.getdata(f"mode{i}","ng" )[0, 0])
    
                return {
                    "mode": i,
                    "te_fraction": te_fraction,
                    "neff": neff,
                    "ng": ng,
                }
    
    raise RuntimeError("No mode found for required TE fraction.")


def grating_coupler_solver(solver,
                           SM_WG_length,
                           slab_width,
                           Monitor_position =50e-6,
                           simulation_time=2000e-15):

    Monitor_position_adjusted = SM_WG_length/2+Monitor_position
    simulation_buffer = 4e-6
    mode_source_position = SM_WG_length/2-2e-6


     
    solver.addvarfdtd() # need to fix this and create a wrapper to do the job for me
    solver.set('x min', mode_source_position-simulation_buffer)
    solver.set('x max', Monitor_position_adjusted+simulation_buffer)
    solver.set('y', 0e-6)
    solver.set("y span", slab_width)
    solver.set("z min", SUBSTRATE_TOP)
    solver.set("z max", CLADDING_THICKNESS)
    solver.set('simulation time', simulation_time)
    solver.set('bandwidth', 'broadband')
    #slab mode position
    solver.set("x0", 10e-6)
    solver.set("y0", 0e-6)
    # test points
    test_points = np.array([10e-6, 10e-6])
    solver.set("number of test points", len(test_points))
    solver.set("test points", test_points)


    solver.addmodesource()
    solver.set('injection axis', "x")
    solver.set('x', mode_source_position)
    solver.set('y', 0e-6)
    solver.set("y span", CLADDING_THICKNESS+2e-6)
    solver.set('wavelength start', WAVELENGTH-BANDWIDTH/2)
    solver.set('wavelength stop', WAVELENGTH+BANDWIDTH/2)

    add_frequency_domain_monitor(solver, 'linear_y', monitor_name= "Field Profile", x= Monitor_position,
                                     x_span=0, y=0, y_span=slab_width, z=CLADDING_THICKNESS/2, z_span=0)