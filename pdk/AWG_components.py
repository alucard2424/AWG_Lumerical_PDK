from tracemalloc import start

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
    Apply the same transformation used by the FPR structure group.

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

def FPR(
        solver,
        input_aperture_spacing,
        output_aperture_spacing,
        Nin,
        Narray,
        Radius_input_arc,
        Radius_output_arc,
        FPR_length,
        input_opening_angle,
        output_opening_angle,
        arc_points=500,
        group_name = "FPR",
        rotation_angle=0.0, #degrees
        output_FPR=False,
        FPR_separation=0.0
):

    input_arc_taper_name = "input_taper_channel_"
    array_taper_name = "Input_FPR_array_taper_"
    if output_FPR:
        input_arc_taper_name = "output_taper_channel_"
        array_taper_name = "Output_FPR_array_taper_"
    solver.addstructuregroup()
    solver.set("name", group_name)
    Nout= Narray
    input_angle  = input_opening_angle*np.pi/180
    output_angle = output_opening_angle*np.pi/180
    WG_buffer= 0#0.1e-6
    Npts = 2*arc_points + 2
    V = np.zeros((Npts, 2))
    n=0

    for i in range(0, arc_points):
        theta = -input_angle/2 +i*input_angle/(arc_points-1)
        V[n,0]= -Radius_input_arc*np.cos(theta)
        V[n,1]= Radius_input_arc*np.sin(theta)
        n=n+1

    #connect top of input arc to top of output arc
    V[n,0] = FPR_length + Radius_output_arc*np.cos(output_angle/2)
    V[n,1] = Radius_output_arc*np.sin(output_angle/2)
    n = n + 1
    for i in range (0, arc_points):
        theta = output_angle/2 - i*output_angle/(arc_points-1)

        V[n,0] = FPR_length + Radius_output_arc*np.cos(theta)
        V[n,1] = Radius_output_arc*np.sin(theta)
        n = n + 1

    #connect bottom of output arc back to input arc
    V[n,0] =-Radius_input_arc*np.cos(input_angle/2)
    V[n,1] =-Radius_input_arc*np.sin(input_angle/2)

    #FPR Geometry
    

    solver.addpoly()
    solver.set("name","Polygon")
    solver.set("material",WG_CORE_MATERIAL)
    solver.set("override mesh order from material database", 1)
    solver.set("mesh order", 1)
    if output_FPR:
        solver.set(  "z",-(SUBSTRATE_TOP + CLADDING_THICKNESS / 2))
    else:
        solver.set(  "z",+(SUBSTRATE_TOP + CLADDING_THICKNESS / 2))
    solver.set("z span", WG_HEIGHT)
    solver.set("vertices",V)
    solver.select("Polygon")
    solver.addtogroup(group_name)
    

    #input wgs
    total_input_span = (Nin-1)*input_aperture_spacing
    insertion_ports=np.zeros((Narray, 2))
    insertion_angles= np.zeros(Narray)
    for i in range(1, Nin+1):

        arc_pos = -total_input_span/2 + (i-1)*input_aperture_spacing

        theta = arc_pos/Radius_input_arc

        x0 = -Radius_input_arc*np.cos(theta)+WG_buffer
        y0 =  Radius_input_arc*np.sin(theta)
        x_port = x0 - TAPER_LENGTH*np.cos(theta)
        y_port = y0 + TAPER_LENGTH*np.sin(theta)

        insertion_ports[i-1] = [x_port, y_port]

        # Direction of the waveguide/taper at the external port
        insertion_angles[i-1] = -theta
        solver.addobject("linear_taper")

        solver.set("name",input_arc_taper_name+str(i))
        solver.set("material",WG_CORE_MATERIAL)
       
        if output_FPR:
            solver.set(  "z",-(SUBSTRATE_TOP + CLADDING_THICKNESS / 2))
        else:
            solver.set(  "z",+(SUBSTRATE_TOP + CLADDING_THICKNESS / 2))
        solver.set("thickness", WG_HEIGHT)
        
        # wide side faces FPR
        solver.set("width_r", TAPER_WIDTH )
        solver.set("width_l", WG_WIDTH)
        
        solver.set("len", TAPER_LENGTH)
        solver.set("angle_side", 90)
        
        # position so taper connects waveguide to slab
        solver.set("x",x0 - 0.5*(TAPER_LENGTH-WG_buffer)*np.cos(theta))
        solver.set("y",y0 + 0.5*TAPER_LENGTH*np.sin(theta))
        
        solver.set("first axis","z")
        solver.set("rotation 1",-theta*180/np.pi)
        #solver.select("input_taper_"+str(i))
        #solver.addtogroup(group_name)

    #output wgs
    total_output_span = (Nout-1)*output_aperture_spacing
    array_ports = np.zeros((Narray, 2))
    array_angles = np.zeros(Narray)
    for i in range(1, Nout+1):

        arc_pos = -total_output_span/2 + (i-1)*output_aperture_spacing
        theta = arc_pos/Radius_output_arc
        x0 = FPR_length + Radius_output_arc*np.cos(theta)-WG_buffer
        y0 = Radius_output_arc*np.sin(theta)
        y_facet = Radius_output_arc*np.sin(theta)
        y_output_end = y_facet + TAPER_LENGTH*np.sin(theta)

        x_port = ( x0+ TAPER_LENGTH*np.cos(theta))

        y_port = (y_facet+TAPER_LENGTH*np.sin(theta))

        # Store local port
        array_ports[i-1] = [x_port, y_port]
        array_angles[i-1] = theta

        solver.addobject("linear_taper")
        solver.set("name",array_taper_name+str(i))
        solver.set("material",WG_CORE_MATERIAL)
       
        if output_FPR:
            solver.set(  "z",-(SUBSTRATE_TOP + CLADDING_THICKNESS / 2))
        else:
            solver.set(  "z",+(SUBSTRATE_TOP + CLADDING_THICKNESS / 2))
        solver.set("thickness", WG_HEIGHT)
        # wide side faces FPR
        solver.set("width_l", TAPER_WIDTH)
        solver.set("width_r", WG_WIDTH)
        solver.set("len", TAPER_LENGTH)
        solver.set("angle_side", 90)
        solver.set("x",x0 + 0.5*(TAPER_LENGTH-WG_buffer)*np.cos(theta))
        solver.set("y",y0 + 0.5*TAPER_LENGTH*np.sin(theta))
        
        solver.set("first axis","z")
        solver.set("rotation 1",theta*180/np.pi)
        #solver.select("output_taper_"+str(i))
        #solver.addtogroup(group_name)
    translation = np.array([0.0, 0.0])
    flip_x = False
    solver.select(group_name)
    solver.set('first axis', 'z')
    solver.set("rotation 1", rotation_angle )
    if output_FPR:
        x = ( FPR_separation+ 2 * (FPR_length+ Radius_output_arc * np.cos(output_angle/2))+ 2 * TAPER_LENGTH)
        translation = np.array([x, 0.0])
        flip_x = True
        solver.set("x", x)
        solver.set('second axis', 'y')
        solver.set("rotation 2", 180)#flip the FPR to face the other way   
    # Transform insertion ports
    insertion_ports_global, insertion_angles_global = transform_ports(
        insertion_ports,
        insertion_angles,
        rotation_deg=rotation_angle,
        translation=translation,
        flip_x=flip_x
    ) 
    array_ports_global, array_angles_global = transform_ports(
        array_ports,
        array_angles,
        rotation_deg=rotation_angle,
        translation=translation,
        flip_x=flip_x
    )
    return {
        "positions": array_ports_global,
        "angles": array_angles_global,
        "insertion positions": insertion_ports_global,
        "insertion angles": insertion_angles_global
    }







########################################



def fixed_bend_radius_route_wg(
    solver, input_port, input_angle, output_port, output_angle,
    bend_radius, delay_length=0.0, base_straight_length=0.0,
    WG_WIDTH=0.5e-6, WG_CORE_MATERIAL=None, SUBSTRATE_TOP=0.0,
    CLADDING_THICKNESS=0.0, WG_HEIGHT=0.22e-6,
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






