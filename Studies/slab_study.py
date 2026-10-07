import ansys.lumerical.core as lumapi
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pdk import create_material_from_data, slab_mode_solver, slab_geometry, store_index_data, get_index_data
from pdk.config import WG_WIDTH, WG_HEIGHT, WG_CORE_MATERIAL, WG_CLADDING_MATERIAL, CLADDING_INDEX_DATA, WG_CORE_INDEX_DATA, DATA_FILE_LOCATION

data_storage_file_path = DATA_FILE_LOCATION
data_name = "Mode_results.csv"
with lumapi.MODE() as solver:
    print(f"This is Lumerical MODE {solver.version()}")
    create_material_from_data(WG_CORE_INDEX_DATA, WG_CORE_MATERIAL,solver)
    create_material_from_data(CLADDING_INDEX_DATA, WG_CLADDING_MATERIAL,solver)
    slab_geometry(WG_HEIGHT, 400e-6,200e-6, solver, WG_CLADDING_MATERIAL, WG_CORE_MATERIAL)
    mode_data =slab_mode_solver(solver=solver, TE_fraction=1.0)
    store_index_data(
        study_name="slab",
        mode=mode_data["mode"],
        te_fraction=mode_data["te_fraction"],
        effective_index=mode_data["neff"],
        group_index=mode_data["ng"],
        file_path=data_storage_file_path,
        data_name= data_name
    )
    retrieval_location = data_storage_file_path / data_name
    # Retrieve results
    slab_data = get_index_data(
        "slab",
        retrieval_location
    )

    print("\nSlab data:")
    print(slab_data)
    input("Press any key to close session")


