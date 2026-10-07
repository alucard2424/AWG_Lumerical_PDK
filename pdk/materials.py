import pandas as pd
import numpy as np
from .config import SPEED_OF_LIGHT

def create_material_from_data(file_path, material_name, solver, mesh_order):
        """
        This function reads a material data file, processes the data, and creates a material in Lumerical FDTD.
        Args:
        - file_path (str): Path to the material data file (in txt format).
        - material_name (str): Name of the material to create in the FDTD simulation.
        """


        # Step 1: Read the data from the file

        df = pd.read_csv(file_path, delimiter="\t", header=None)
        df.columns = ["Wavelength (m)", "Re(index)", "Im(index)"]

        # Print the first few rows of the data
        print("Data from file:")
        print(df.head())
        
        
        # Step 2: Convert wavelength to frequency and calculate permittivity
        wavelengths = df.iloc[:, 0].values * 10**-6  
        frequency = SPEED_OF_LIGHT / wavelengths  # Frequency in Hz
        # Refractive index (n) from the real part of the refractive index (n = Re(index))
        n = df.iloc[:, 1].values  # Real part of the refractive index
        permittivity = n**2  # Permittivity is the square of the refractive index
        sampledData = np.column_stack([frequency, permittivity])
        Material=solver.addmaterial("Sampled 3D data")
        solver.setmaterial(Material,"name",material_name)
        
        
        try:
            #solver.setmaterial(material_name, "Sampled 3D data", sampledData)
            solver.setmaterial(material_name, "sampled data", sampledData)
            solver.setmaterial(material_name,"Mesh order",mesh_order)
            print(f"Material '{material_name}' created successfully!")
        except Exception as e:
             print(f"Error setting material: {e}")



def create_material_from_Sellmeier_Coefficients (Coefficients_array, Material_name, solver):
        try:
            Material=solver.addmaterial("Sellmeier")
            solver.setmaterial(Material,"name",Material_name)
            solver.setmaterial(Material_name,"B1",Coefficients_array[0])
            solver.setmaterial(Material_name,"C1",Coefficients_array[1])
            solver.setmaterial(Material_name,"B2",Coefficients_array[2])
            solver.setmaterial(Material_name,"C2",Coefficients_array[3])
            solver.setmaterial(Material_name,"B3",Coefficients_array[4])
            solver.setmaterial(Material_name,"C3",Coefficients_array[5])
            print(f"Sellmeier Material '{Material_name}' created successfully!")

        except Exception as e:
            print(f"Error while creating Sellmeier material '{Material_name}': {e}")
    

def create_material_from_index(real_refractive_index, imaginary_refractive_index, material_name, solver):
    try:
        # Step 1: Take the refractive index data as input argument
        n=real_refractive_index
        k=imaginary_refractive_index
        # Step 2: Add material to FDTD environment
        solver.setmaterial(solver.addmaterial("(n,k) Material"), "name", material_name)
        # Step 3: Assign index data to material in FDTD environment
        solver.setmaterial(material_name, {"Refractive Index": n, "Imaginary Refractive Index": k})
        print(f"Material '{material_name}' created successfully!")

    except Exception as e:
        print(f"Error while creating material '{material_name}': {e}")