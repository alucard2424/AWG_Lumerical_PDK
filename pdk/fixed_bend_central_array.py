import numpy as np
from scipy.constants import c
import math
from .config import WG_CLADDING_MATERIAL, WG_CORE_MATERIAL, WG_WIDTH, WG_HEIGHT, MIN_BEND_RADIUS, SUBSTRATE_TOP, CLADDING_THICKNESS, TAPER_LENGTH, TAPER_WIDTH, STRAIGHT_SECTION_LENGTH
from scipy.constants import c


def transform_ports(points, angles,
                    rotation_deg=0.0,
                    translation=(0.0, 0.0),
                    flip_x=False):
    """

    angles are in radians.

    rotation_deg:
        rotation about z-axis (rotation 1)

    flip_x:
        equivalent to a 180-degree rotation about y-axis
        in the 2D layout plane: (x, y) -> (-x, y)
    """

    points = np.asarray(points, dtype=float)
    angles = np.asarray(angles, dtype=float)

    #rotation about z
    alpha = np.deg2rad(rotation_deg)
    R = np.array([
        [np.cos(alpha), -np.sin(alpha)],
        [np.sin(alpha),  np.cos(alpha)] ])
    points_global = points @ R.T
    angles_global = angles + alpha
    # rotation 2 flip if output fpr, in the XY plane: x->-x,y -> y
    if flip_x:
        points_global[:, 0] *= -1
        angles_global = np.pi - angles_global

    #translate group
    points_global += np.asarray(translation)
    # Normalize angles
    angles_global = ( angles_global + np.pi) % (2*np.pi) - np.pi
    return points_global, angles_global

def draw_port_marker(solver, point, name, size=10e-6):

    x,y =point
    h = size /2

    vertices = np.array([[x-h, y-h],[x+h, y-h],
        [x+h, y+h],
        [x-h, y+h]
    ])

    solver.addpoly()
    solver.set("name", name)
    solver.set("material", WG_CORE_MATERIAL)
    solver.set("z", SUBSTRATE_TOP + CLADDING_THICKNESS / 2)
    solver.set("z span", WG_HEIGHT)
    solver.set("vertices", vertices)









########################################



def fixed_bend_radius_route_wg(
    solver, input_port, input_angle, output_port, output_angle,
    bend_radius, delay_length=0.0, base_straight_length=0.0,
    WG_WIDTH=WG_WIDTH, WG_CORE_MATERIAL=None, SUBSTRATE_TOP=SUBSTRATE_TOP,
    CLADDING_THICKNESS=CLADDING_THICKNESS,  WG_HEIGHT=WG_HEIGHT,
    target_input_angle=0.0, target_output_angle=np.pi,
    npoints=100, tolerance=1e-12, max_iterations=100,
):
    """Complete fixed-radius two-bend route with exact target length.

    The same straight length is used on both ports. It is solved so that
    actual centerline length equals the zero-delay length plus delay_length.
    """
    input_port = np.asarray(input_port, dtype=float)
    output_port = np.asarray(output_port, dtype=float)
    if bend_radius <= WG_WIDTH / 2:
        raise ValueError("bend_radius must be larger than WG_WIDTH / 2")
    if delay_length < 0:
        raise ValueError("delay_length must be non-negative")

    def wrap(a):
        return (a + np.pi) % (2*np.pi) - np.pi

    def u(a):
        return np.array([np.cos(a), np.sin(a)])

    def bend_geometry(start, start_angle, end_angle):
        ba = wrap(end_angle - start_angle)
        if abs(ba) < 1e-12:
            return start.copy(), 0.0, None, ba
        sign = np.sign(ba)
        left = np.array([-np.sin(start_angle), np.cos(start_angle)])
        center = start + sign * bend_radius * left
        r0 = start - center
        phi0 = np.arctan2(r0[1], r0[0])
        phi1 = phi0 + ba
        end = center + bend_radius * u(phi1)
        return end, bend_radius * abs(ba), center, ba

    def geometry(port_straight):
        inp_end = input_port + port_straight * u(input_angle)
        in_bend_end, in_bl, in_center, in_ba = bend_geometry(
            inp_end, input_angle, target_input_angle)

        out_dir = -u(output_angle)
        out_inner = output_port - port_straight * out_dir
        out_bend_start_angle = np.arctan2(
            (out_inner - output_port)[1], (out_inner - output_port)[0])
        out_bend_end, out_bl, out_center, out_ba = bend_geometry(
            out_inner, out_bend_start_angle, target_output_angle)

        middle_vec = out_bend_end - in_bend_end
        middle_len = np.linalg.norm(middle_vec)
        middle_dir = (middle_vec / middle_len if middle_len > 1e-15
                      else u(target_input_angle))
        total = 2*port_straight + in_bl + middle_len + out_bl
        return dict(
            input_straight_end=inp_end, input_bend_end=in_bend_end,
            output_straight_end=out_inner, output_bend_end=out_bend_end,
            input_bend_length=in_bl, output_bend_length=out_bl,
            middle_straight_length=middle_len, middle_direction=middle_dir,
            input_bend_center=in_center, output_bend_center=out_center,
            input_bend_angle=in_ba, output_bend_angle=out_ba,
            port_straight_length=port_straight, total_length=total)

    zero = geometry(base_straight_length)
    target = zero['total_length'] + delay_length

    def f(s):
        return geometry(s)['total_length'] - target

    lo = 0.0
    hi = max(base_straight_length + delay_length, base_straight_length, 1e-9)
    while f(hi) < 0:
        hi *= 2.0
        if hi > 1e3:
            raise RuntimeError('Could not bracket exact-length solution')
    for _ in range(max_iterations):
        mid = 0.5*(lo + hi)
        if f(mid) < 0:
            lo = mid
        else:
            hi = mid
        if abs(f(mid)) <= tolerance:
            break
    g = geometry(0.5*(lo + hi))
    if abs(g['total_length'] - target) > tolerance:
        raise RuntimeError('Exact length solve did not converge')

    def draw_straight(name, start, direction, length):
        direction = np.asarray(direction, float)
        direction /= np.linalg.norm(direction)
        center = start + 0.5*length*direction
        if solver is not None:
            solver.addrect(); solver.set('name', name)
            solver.set('material', WG_CORE_MATERIAL)
            solver.set('z', SUBSTRATE_TOP + CLADDING_THICKNESS/2)
            solver.set('z span', WG_HEIGHT)
            solver.set('x span', length); solver.set('y span', WG_WIDTH)
            solver.set('x', center[0]); solver.set('y', center[1])
            solver.set('first axis', 'z')
            solver.set('rotation 1', np.degrees(np.arctan2(direction[1], direction[0])))

    def draw_bend(name, start, start_angle, end_angle):
        end, length, center, ba = bend_geometry(start, start_angle, end_angle)
        if solver is not None and length > 0:
            ro = bend_radius + WG_WIDTH/2; ri = bend_radius - WG_WIDTH/2
            phi0 = np.arctan2((start-center)[1], (start-center)[0])
            phi = np.linspace(phi0, phi0+ba, npoints)
            outer = np.column_stack((center[0]+ro*np.cos(phi), center[1]+ro*np.sin(phi)))
            inner = np.column_stack((center[0]+ri*np.cos(phi[::-1]), center[1]+ri*np.sin(phi[::-1])))
            solver.addpoly(); solver.set('name', name)
            solver.set('material', WG_CORE_MATERIAL)
            solver.set('z', SUBSTRATE_TOP + CLADDING_THICKNESS/2)
            solver.set('z span', WG_HEIGHT); solver.set('vertices', np.vstack((outer, inner)))
        return end

    s = g['port_straight_length']
    draw_straight('input_straight_wg', input_port, u(input_angle), s)
    draw_bend('input_bend_to_horizontal', g['input_straight_end'], input_angle, target_input_angle)
    draw_straight('output_straight_wg', output_port, u(output_angle), s)
    draw_bend('output_bend_to_horizontal', g['output_straight_end'], np.arctan2((u(output_angle))[1], (u(output_angle))[0]), target_output_angle)
    draw_straight('middle_straight_wg', g['input_bend_end'], g['middle_direction'], g['middle_straight_length'])

    g.update(zero_delay_total_length=zero['total_length'], target_total_length=target,
             delay_length=delay_length, length_error=g['total_length']-target)
    return g






def optimize_awg_routes_after_port0_with_fixed_bend_radius(
    input_ports,
    output_ports,
    input_angles,
    output_angles,
    port0_geometry,
    delta_L,
    min_bend_radius,
    min_wg_spacing,
    base_straight_length=None,
    WG_WIDTH=WG_WIDTH,

    # Optimization bounds
    fixed_bend_radius=MIN_BEND_RADIUS,
    min_straight_length=0e-6,
    min_middle_length=0e-6,
    verbose=True,

):

    input_ports = np.asarray(input_ports, dtype=float)
    output_ports = np.asarray(output_ports, dtype=float)
    input_angles = np.asarray(input_angles, dtype=float)
    output_angles = np.asarray(output_angles, dtype=float)

    N = len(input_ports)

    if not (len(output_ports) == N and len(input_angles) == N and len(output_angles) == N):
            raise ValueError("All port arrays must have the same length.")   
    if N < 2:
        raise ValueError("Need at least port 0 and one additional port.")
    if delta_L <= 0:
        raise ValueError("delta_L must be positive." )

    if WG_WIDTH <= 0:
        raise ValueError("WG_WIDTH must be positive.")

    if min_wg_spacing < 0:
        raise ValueError("min_wg_spacing must be >= 0.")

    if min_bend_radius <= WG_WIDTH / 2.0:
        raise ValueError("min_bend_radius is too small.")

    if min_middle_length < 0:
        raise ValueError("min_middle_length must be >= 0.")

    if min_straight_length < 0:
        raise ValueError("min_straight_length must be >= 0.")
    if "total_length" in port0_geometry:
            L_base = float(port0_geometry["total_length"])
    
    elif "target_total_length" in port0_geometry:
        L_base = float(port0_geometry["target_total_length"])

    elif "zero_delay_total_length" in port0_geometry:
        L_base = float(port0_geometry["zero_delay_total_length"])

    else:
        raise ValueError("port0_geometry must contain ""'total_length', 'target_total_length', " "or 'zero_delay_total_length'.")

    
    def wrap(angle):
        return ((angle + np.pi)% (2.0 * np.pi)- np.pi)

    def unit(angle):
        return np.array([np.cos(angle),np.sin(angle),])

    
    #bend geomtry
    

    def bend_geometry(
        start,
        start_angle,
        end_angle,
        R,
    ):

        start = np.asarray(start, dtype=float,)
        theta = wrap(end_angle - start_angle)

        if abs(theta) < 1e-14:
            return (start.copy(), None, 0.0, 0.0,)

        sign = np.sign(theta)
        left = np.array([-np.sin(start_angle), np.cos(start_angle),])
        center = (start+ sign * R * left )
        r0 = start - center
        phi0 = np.arctan2(r0[1], r0[0],)
        phi1 = phi0 + theta
        end = ( center+ R * np.array([np.cos(phi1),np.sin(phi1),]))
        return ( end,center,theta, R * abs(theta),)


    R = fixed_bend_radius
    s0 = base_straight_length if base_straight_length is not None else 0.0

    def geometry(i, s, R):
        ip, ia = input_ports[i], input_angles[i]
        op, oa = output_ports[i], output_angles[i]

        inp_end = ip + s * unit(ia)
        in_end, in_c, in_ba, in_bl = bend_geometry(inp_end, ia, 0.0, R)

        out_dir = -unit(oa)
        out_inner = op - s * out_dir
        out_start_angle = np.arctan2((out_inner - op)[1], (out_inner - op)[0])
        out_end, out_c, out_ba, out_bl = bend_geometry(out_inner, out_start_angle, np.pi, R)

        mid_vec = out_end - in_end
        mid_len = np.linalg.norm(mid_vec)
        mid_dir = mid_vec / mid_len if mid_len > 1e-15 else unit(0.0)

        total = 2 * s + in_bl + mid_len + out_bl
        return dict(
            R=R, S=s,
            input_straight_end=inp_end,
            input_bend_end=in_end,
            input_bend_center=in_c,
            input_bend_angle=in_ba,
            input_bend_length=in_bl,
            output_straight_end=out_inner,
            output_bend_end=out_end,
            output_bend_center=out_c,
            output_bend_angle=out_ba,
            output_bend_length=out_bl,
            middle_direction=mid_dir,
            middle_length=mid_len,          
            middle_straight_length=mid_len, 
            port_straight_length=s,
            total_length=total,
        )

    def solve_channel(i, target_length, R):
        def f(s):
            return geometry(i, s, R)['total_length'] - target_length

        lo, hi = 0.0, max(s0, 1e-9)
        tries = 0
        while f(hi) < 0:
            hi *= 2.0
            tries += 1
            if tries > 60:
                raise RuntimeError(f"Could not bracket length solution for channel {i}")
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if f(mid) < 0:
                lo = mid
            else:
                hi = mid
            if abs(f(mid)) < 1e-12:
                break
        g = geometry(i, 0.5 * (lo + hi), R)
        g['target_total_length'] = target_length
        g['length_error'] = g['total_length']-target_length
        return g

    routes = [None] * N
    routes[0] = dict(port0_geometry)
    routes[0].setdefault('R', fixed_bend_radius)

    if 'S' not in routes[0]:
        routes[0]['S'] = routes[0].get('port_straight_length', s0)

    routes[0].setdefault('middle_length',routes[0].get('middle_straight_length', 0.0))

    routes[0].setdefault('middle_straight_length',routes[0]['middle_length'] )

    routes[0].setdefault('target_total_length', L_base)
    for i in range(1, N):
        target_length_i = L_base + i * delta_L
        g_i = solve_channel(i, target_length_i, R)

        if g_i['port_straight_length'] < min_straight_length:
            raise ValueError(
                f"Channel {i}: solved straight length "
                f"{g_i['port_straight_length']:.3e} m is below "
                f"min_straight_length. Increase base_straight_length, "
                f"delta_L, or bend_radius.")
        if g_i['middle_straight_length'] < min_middle_length:
            raise ValueError(
                f"Channel {i}: middle straight length "
                f"{g_i['middle_straight_length']:.3e} m is below "
                f"min_middle_length." )
        routes[i] = g_i

    # --- crude adjacent-channel spacing sanity check ---
    # Assumes the "onto the bus" bend-end points line up transversely
    # (i.e. the middle straights run parallel). Adjust the axis/points
    # you compare if your star-coupler geometry differs.
    for i in range(1, N):
        y_prev = routes[i - 1]['input_bend_end'][1]
        y_curr = routes[i]['input_bend_end'][1]
        gap = abs(y_curr - y_prev)
        if gap < (min_wg_spacing + WG_WIDTH):
            if verbose:
                print(
                    f"Warning: channels {i-1} and {i} appear closer "
                    f"({gap*1e6:.3f} µm) than min_wg_spacing + WG_WIDTH "
                    f"({(min_wg_spacing+WG_WIDTH)*1e6:.3f} µm)."
                )

    return routes




def draw_awg_routes_fixed_radius(
    solver,
    routes,
    input_ports,
    input_angles,
    output_ports,
    output_angles,
    WG_WIDTH=0.5e-6,
    WG_CORE_MATERIAL=None,
    SUBSTRATE_TOP=SUBSTRATE_TOP,
    CLADDING_THICKNESS=CLADDING_THICKNESS,
    WG_HEIGHT=WG_HEIGHT,
    npoints=100,
    name_prefix="awg",
):
    """
    Draw every optimized AWG central-array route into Lumerical.

    Physical route (per channel):

        input_port
            |
            | input straight
            v
           p1
            )
           )
          p2
           ---------------- p3
                            )
                           )
                          p4
                           |
                           | output straight
                           v
                       output_port

    routes[i] is expected to be a dict with keys:
        R, S,
        input_straight_end, input_bend_end,
        input_bend_center, input_bend_angle,
        output_straight_end, output_bend_end,
        output_bend_center, output_bend_angle,
        middle_direction, middle_length (or middle_straight_length)

    input_ports / output_ports / input_angles / output_angles are the
    same per-channel arrays passed to the optimizer.
    """

    input_ports = np.asarray(input_ports, dtype=float)
    output_ports = np.asarray(output_ports, dtype=float)
    input_angles = np.asarray(input_angles, dtype=float)
    output_angles = np.asarray(output_angles, dtype=float)

    N = len(routes)
    if not (len(input_ports) == N and len(output_ports) == N
            and len(input_angles) == N and len(output_angles) == N):
        raise ValueError("routes and port arrays must all have length N.")

    def unit(angle):
        return np.array([np.cos(angle), np.sin(angle)])

    
    #straights
    def draw_straight(name, start, direction, length):
        direction = np.asarray(direction, dtype=float)
        norm = np.linalg.norm(direction)
        if norm < 1e-15:
            raise RuntimeError(f"{name}: zero-length direction.")
        direction = direction / norm
        center = start + 0.5 * length * direction
        solver.addrect()
        solver.set('name', name)
        solver.set('material', WG_CORE_MATERIAL)
        solver.set('z', SUBSTRATE_TOP + CLADDING_THICKNESS / 2)
        solver.set('z span', WG_HEIGHT)
        solver.set('x span', length)
        solver.set('y span', WG_WIDTH)
        solver.set('x', center[0])
        solver.set('y', center[1])
        solver.set('first axis', 'z')
        solver.set('rotation 1', np.degrees(np.arctan2(direction[1], direction[0])))

    def draw_bend(name, arc_start_point, center, bend_angle, R):
        if abs(bend_angle) < 1e-12:
            return

        center = np.asarray(center, dtype=float)
        arc_start_point = np.asarray(arc_start_point, dtype=float)
        phi0 = np.arctan2(
            arc_start_point[1] - center[1],
            arc_start_point[0] - center[0],
        )
        phi = np.linspace(phi0, phi0 + bend_angle, npoints)

        ro = R + WG_WIDTH / 2
        ri = R - WG_WIDTH / 2
        if ri <= 0:
            raise ValueError(f"{name}: inner bend radius <= 0.")

        outer = np.column_stack((center[0] + ro * np.cos(phi), center[1] + ro * np.sin(phi)))
        inner = np.column_stack((center[0] + ri * np.cos(phi[::-1]),center[1] + ri * np.sin(phi[::-1])))
        vertices = np.vstack((outer, inner))
        solver.addpoly()
        solver.set('name', name)
        solver.set('material', WG_CORE_MATERIAL)
        solver.set('z', SUBSTRATE_TOP + CLADDING_THICKNESS / 2)
        solver.set('z span', WG_HEIGHT)
        solver.set('vertices', vertices)

    #loop over channels
   
    for i in range(N):
        g = routes[i]

        R = g["R"]
        S = g["S"]
        if R <= WG_WIDTH / 2:
            raise ValueError(f"Port {i}: bend radius too small.")
        if S < 0:
            raise ValueError(f"Port {i}: negative straight length.")
        p1 = np.asarray(g["input_straight_end"], dtype=float)
        p2 = np.asarray(g["input_bend_end"], dtype=float)
        p3 = np.asarray(g["output_bend_end"], dtype=float)
        p4 = np.asarray(g["output_straight_end"], dtype=float)

        middle_direction = np.asarray(g["middle_direction"], dtype=float)
        middle_length = g.get("middle_length", g.get("middle_straight_length"))

        port_name = f"port_{i:02d}"

        #input straight
        draw_straight(f"{name_prefix}_{port_name}_input_straight", input_ports[i], unit(input_angles[i]), S, )

        #input bend: p1 -> p2, arc starts at p1
        draw_bend( f"{name_prefix}_{port_name}_input_bend", p1, g["input_bend_center"], g["input_bend_angle"], R,)

        #middle straight: p2 -> p3
        draw_straight( f"{name_prefix}_{port_name}_middle_straight",p2, middle_direction, middle_length,)

        #output bend:optimizer traces it p4 -> p3, so arc starts at p4
        draw_bend(f"{name_prefix}_{port_name}_output_bend", p4, g["output_bend_center"], g["output_bend_angle"], R,)

        #output straight
        draw_straight(f"{name_prefix}_{port_name}_output_straight", output_ports[i], unit(output_angles[i]), S,)

        
        #sanity check vs. expected straight endpoints
        
        expected_p1 = input_ports[i] + S * unit(input_angles[i])
        expected_p4 = output_ports[i] + S * unit(output_angles[i])

        p1_error = np.linalg.norm(p1 - expected_p1)
        p4_error = np.linalg.norm(p4 - expected_p4)

        if p1_error > 1e-9:
            print(f"WARNING Port {i}: input endpoint mismatch = {p1_error*1e6:.6f} um")
        if p4_error > 1e-9:
            print(f"WARNING Port {i}: output endpoint mismatch = {p4_error*1e6:.6f} um")

    return routes