from pathlib import Path
import pandas as pd
import os


def store_index_data(
    study_name,
    mode,
    te_fraction,
    effective_index,
    group_index,
    file_path,
    data_name
):
    """Store mode-index results for a study."""

    file_path = Path(file_path)
    file_path.mkdir(parents=True, exist_ok=True)

    data = pd.DataFrame({
        "Study": [study_name],
        "Mode": [mode],
        "TE Fraction": [te_fraction],
        "Effective Index": [effective_index],
        "Group Index": [group_index]
    })

    results_file = file_path / data_name

    if results_file.exists():
        df_existing = pd.read_csv(results_file)

        # Remove the old entry for THIS study only
        df_existing = df_existing[
            df_existing["Study"] != study_name
        ]

        # Add the updated study entry
        data = pd.concat(
            [df_existing, data],
            ignore_index=True
        )

    # Write to the actual results file
    data.to_csv(results_file, index=False)

def get_index_data(study_name, file_path):
    """Retrieve index data belonging to a specific study."""

    data=pd.read_csv(file_path)

    return data[data["Study"] ==study_name]