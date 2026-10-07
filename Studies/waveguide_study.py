import ansys.lumerical.core as lumapi
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pdk import create_material_from_data, waveguide_mode_solver, SM_straight_waveguide_geometry, store_index_data, get_index_data
from pdk.config import WG_WIDTH, WG_HEIGHT, WG_CORE_MATERIAL, WG_CLADDING_MATERIAL, CLADDING_INDEX_DATA, WG_CORE_INDEX_DATA, DATA_FILE_LOCATION

data_storage_file_path = DATA_FILE_LOCATION 
data_name = "Mode_results.csv"
with lumapi.MODE() as solver:
    print(f"This is Lumerical solver {solver.version()}")
    create_material_from_data(WG_CORE_INDEX_DATA, WG_CORE_MATERIAL,solver, mesh_order=1)
    create_material_from_data(CLADDING_INDEX_DATA, WG_CLADDING_MATERIAL,solver, mesh_order=2)
    SM_straight_waveguide_geometry(WG_WIDTH, WG_HEIGHT, 100e-6,solver, WG_CLADDING_MATERIAL, WG_CORE_MATERIAL, "SM_WG")
    solver.addmesh()
    solver.set("name", "WG mesh")
    solver.set('based on a structure', 1)
    solver.set('structure', 'SM_WG')
    solver.set('dx', 0.05e-6)
    solver.set('dy', 0.05e-6)
    solver.set('dz', 0.05e-6)
    mode_data = waveguide_mode_solver(solver=solver, TE_fraction=1.0)
    store_index_data(
        study_name="waveguide",
        mode=mode_data["mode"],
        te_fraction=mode_data["te_fraction"],
        effective_index=mode_data["neff"],
        group_index=mode_data["ng"],
        file_path=data_storage_file_path,
        data_name= data_name
    )
    # Retrieve results
    retrieval_location = data_storage_file_path / data_name
    WG_data = get_index_data(
        "waveguide",
        retrieval_location
    )

    print("\nWaveguide data:")
    print(WG_data)
    input("Press any key to close session")