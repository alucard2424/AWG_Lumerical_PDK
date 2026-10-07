from .config import SUBSTRATE_TOP


def add_mesh_override(
        solver,
        mesh_name,
        structure_name,
        override_stepsize_x = 0.05e-6,
        override_stepsize_y = 0.05e-6,
        override_stepsize_z = 0.05e-6,
):
    solver.addmesh()
    solver.set("name", mesh_name)
    solver.set('based on a structure', 1)
    solver.set('structure', structure_name)
    solver.set('dx', override_stepsize_x)
    solver.set('dy', override_stepsize_y)
    solver.set('dz', override_stepsize_z)




def add_frequency_domain_monitor(
    solver,
    monitor_type="2d_xy",
    monitor_name="field_profile",
    x=0,
    x_span=5e-6,
    y=0, y_span=5e-6,
    z=0,
    z_span=5e-6,
):
    """
    Add a frequency-domain monitor.

    monitor_type:
        linear_x
        linear_y
        linear_z
        2d_xy
        2d_xz
        2d_yz
        3d
    """

    monitor_types = {
        "linear_x": 2,
        "linear_y": 3,
        "linear_z": 4,

        "2d_xy": 5,
        "2d_xz": 6,
        "2d_yz": 7,

        "3d": 8, }

    if monitor_type not in monitor_types:
        raise ValueError(
            f"Unknown monitor type '{monitor_type}'. "
            f"Choose from: {list(monitor_types.keys())}" )

    solver.adddftmonitor()
    solver.set("name", monitor_name)
    solver.set("monitor type", monitor_types[monitor_type])
    z = SUBSTRATE_TOP+z
    
    # Line monitors
   

    if monitor_type == "linear_x":
        solver.set("y", y)
        solver.set("z", z)
        solver.set("x span", x_span)

    elif monitor_type == "linear_y":
        solver.set("x", x)
        solver.set("z", z)
        solver.set("y span", y_span)

    elif monitor_type == "linear_z":
        solver.set("x", x)
        solver.set("y", y)
        solver.set("z span", z_span)

   
    #2D Monitors
   

    elif monitor_type == "2d_xy":
        solver.set("z", z)
        solver.set("x span", x_span)
        solver.set("y span", y_span)

    elif monitor_type == "2d_xz":
        solver.set("y", y)
        solver.set("x span", x_span)
        solver.set("z span", z_span)

    elif monitor_type == "2d_yz":
        solver.set("x", x)
        solver.set("y span", y_span)
        solver.set("z span", z_span)
    #3D Monitor
    

    elif monitor_type == "3d":
        solver.set("x", x)
        solver.set("y", y)
        solver.set("z", z)
        solver.set("x span", x_span)
        solver.set("y span", y_span)
        solver.set("z span", z_span)
    


