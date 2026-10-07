from .materials import (
    create_material_from_data,
    create_material_from_Sellmeier_Coefficients,
    create_material_from_index,
)


from . import config


from .geometry import (
    SM_straight_waveguide_geometry,
    add_substrate,
    slab_geometry,
)

from .modes import(
    slab_mode_solver,
    waveguide_mode_solver,
    grating_coupler_solver
)

from .database import (
    store_index_data,
    get_index_data,
)


from .simulation_tools import (
    add_mesh_override,
    add_frequency_domain_monitor,
)

from .AWG_components import(
    FPR,
    draw_port_marker,
    fixed_bend_radius_route_wg,
    # optimize_awg_routes_after_port0,
    # draw_optimized_awg_route,
)
from .fixed_bend_central_array import(
    optimize_awg_routes_after_port0_with_fixed_bend_radius,
    draw_awg_routes_fixed_radius,
)