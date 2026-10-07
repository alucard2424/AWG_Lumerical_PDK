from .config import SUBSTRATE_TOP, SUBSTRATE_THICKNESS, SUBSTRATE_MATERIAL, CLADDING_THICKNESS, WG_HEIGHT, WG_CORE_MATERIAL, WG_WIDTH, WG_CLADDING_MATERIAL, TAPER_LENGTH, TAPER_WIDTH, STRAIGHT_SECTION_LENGTH
import numpy as np

def SM_straight_waveguide_geometry(wg_width, wg_height, wg_span, solver, 
                                   cladding_material_name=WG_CLADDING_MATERIAL, 
                                   core_material_name=WG_CORE_MATERIAL, 
                                   WG_name='SM Waveguide'):
    solver.addstructuregroup()
    
    solver.set("name", WG_name)
    solver.addrect()
    solver.set("name", 'Cladding')
    solver.addtogroup(WG_name)
    solver.set("x", 0e-6)
    solver.set("x span", wg_span)
    solver.set("y", 0e-6)
    solver.set("y span", CLADDING_THICKNESS)
    solver.set("z min", SUBSTRATE_TOP)
    solver.set("z max", CLADDING_THICKNESS)
    solver.set("material", cladding_material_name)
    solver.set('override mesh order from material database', 1)
    solver.set('mesh order', 2)
    solver.set('render type', "wireframe")
    
        

    solver.addrect()
    solver.set("name", 'Waveguide')
    solver.addtogroup(WG_name)
    solver.set("x", 0e-6)
    solver.set("x span", wg_span)
    solver.set("y", 0e-6)
    solver.set("y span", wg_width)
    solver.set("z", SUBSTRATE_TOP + CLADDING_THICKNESS/2 )
    solver.set("z span", wg_height)   
    solver.set("material", core_material_name)
    solver.set('override mesh order from material database', 1)
    solver.set('mesh order', 1)



def add_substrate(
    x_span,
    y_span,
    solver,
    substrate_name="Substrate"
):
    solver.addrect()

    solver.set("name", substrate_name)

    solver.set("x", 0)
    solver.set("x span", x_span)

    solver.set("y", 0)
    solver.set("y span", y_span)

    
    solver.set("z min", SUBSTRATE_TOP - SUBSTRATE_THICKNESS)
    solver.set("z max", SUBSTRATE_TOP)

    solver.set("material", SUBSTRATE_MATERIAL)

    solver.set("override mesh order from material database", 1)
    solver.set("mesh order", 5)


def slab_geometry(SM_wg_height,
                length, 
                width, solver, 
                cladding_material_name=WG_CLADDING_MATERIAL, 
                core_material_name=WG_CORE_MATERIAL, 
):
    solver.addstructuregroup()
        
    solver.set("name", "Slab")
    solver.addrect()
    solver.set("name", 'Cladding')
    solver.addtogroup('Slab')
    solver.set("x", 0e-6)
    solver.set("x span", length)
    solver.set("y", 0e-6)
    solver.set("y span", CLADDING_THICKNESS+width)
    solver.set("z min", SUBSTRATE_TOP)
    solver.set("z max", CLADDING_THICKNESS)
    solver.set("material", cladding_material_name)
    solver.set('override mesh order from material database', 1)
    solver.set('mesh order', 2)
    solver.set('render type', "wireframe")
    
        

    solver.addrect()
    solver.set("name", 'Waveguide')
    solver.addtogroup('Slab')
    solver.set("x", 0e-6)
    solver.set("x span", length)
    solver.set("y", 0e-6)
    solver.set("y span", width)
    solver.set("z", SUBSTRATE_TOP + CLADDING_THICKNESS/2 )
    solver.set("z span", SM_wg_height)   
    solver.set("material", core_material_name)
    solver.set('override mesh order from material database', 1)
    solver.set('mesh order', 1)
    


    

