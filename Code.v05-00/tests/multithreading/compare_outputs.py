"""
Script to compare the outputs of two APCEMM runs to ensure data variables are
bitwise identical. Uses Python PEP 723 for script dependency declaration, works
nicely with uv. Run this with 'uv run compare_outputs.py dir1/path1/ dir2/path2/

Returns 0 if files in directories are identical.
"""
# /// script
# requires-python = ">=3.12"
# dependencies = ["xarray", "netCDF4"]
# ///

import sys
from pathlib import Path

import xarray as xr

ALL_IDENTICAL = 0
NOT_IDENTICAL = 1
BAD_USAGE = 2


def compare_identical(path1: Path, path2: Path) -> int:
    """
    Compare the ts_aerosol.nc files in two APCEMM output directories.
    Use xarray.equals() to ensure that each variable is bitwise identical
    without enforcing that attrs must be identical

    Parameters
    ----------
    path1 : Path
        Path to APCEMM output dir
    path2 : Path
        Path to APCEMM output dir

    Returns
    -------
    int
        0 if all equal
        1 if at least an output is different
        2 bad call to function
    """

    if path1.resolve() == path2.resolve():
        print(
            "Do not compare the same output directory (paths are the same)", flush=True
        )
        return BAD_USAGE

    # Fetch all ts_aerosol files
    files1 = sorted(path1.glob("ts*.nc"))
    files2 = sorted(path2.glob("ts*.nc"))

    if len(files1) != len(files2):
        print("Not same number of files", flush=True)
        return NOT_IDENTICAL

    if len(files1) == 0:
        print("Found 0 files to compare", flush=True)
        return NOT_IDENTICAL

    print(f"Found {len(files1)} files to compare", flush=True)

    with (
        xr.open_dataset(path1 / "epm-output.nc", decode_times=False) as ds1,
        xr.open_dataset(path2 / "epm-output.nc", decode_times=False) as ds2,
    ):
        if not ds1.equals(ds2):
            print("EPM outputs are different", flush=True)
            return NOT_IDENTICAL

    all_equal = True

    for f1, f2 in zip(files1, files2, strict=True):
        with (
            xr.open_dataset(f1, decode_times=False) as ds1,
            xr.open_dataset(f2, decode_times=False) as ds2,
        ):
            is_equal = ds1.equals(ds2)
            all_equal = all_equal and is_equal

            if not is_equal:
                print(f"{f1} and {f2} are different!", flush=True)
                break

    print(f"{all_equal = }", flush=True)

    if not all_equal:
        return NOT_IDENTICAL

    return ALL_IDENTICAL


if __name__ == "__main__":
    # Check we have both paths as arguments
    if len(sys.argv) != 3:
        print("Script usage: compare_outputs.py dir1/path1/ dir2/path2/", flush=True)
        sys.exit(BAD_USAGE)

    path1 = Path(sys.argv[1])
    path2 = Path(sys.argv[2])

    # Propagate return code of function to script return code
    sys.exit(compare_identical(path1, path2))
