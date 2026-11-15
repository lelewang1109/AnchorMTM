"""Extract the January 2017 northern-China subset from ChinaVis CAQRA CSVs."""

import argparse
from pathlib import Path

import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv_directory", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    import pandas as pd

    files = sorted(args.csv_directory.glob("CN-Reanalysis*.csv"))
    if len(files) != 744:
        parser.error(f"expected 744 hourly CSVs, found {len(files)}")
    times = np.array([path.stem[-10:] for path in files])
    expected = np.array(
        [f"201701{day:02d}{hour:02d}" for day in range(1, 32) for hour in range(24)]
    )
    np.testing.assert_array_equal(times, expected)

    def load(path):
        return pd.read_csv(
            path, usecols=[0, 6, 7, 11, 12], skipinitialspace=True, encoding="latin1"
        ).to_numpy(float)

    first = load(files[0])
    lat, lon = first[:, 3], first[:, 4]
    mask = (109.5 <= lon) & (lon <= 120.5) & (31.5 <= lat) & (lat <= 41.5)
    coordinates = np.column_stack((lon[mask], lat[mask]))
    pm25, wind = [], []
    for path in files:
        frame = load(path)
        np.testing.assert_array_equal(frame[:, 3:], first[:, 3:])
        pm25.append(frame[mask, 0])
        wind.append(frame[mask, 1:3])
    pm25 = np.asarray(pm25)
    if pm25.size == 0 or not np.isfinite(pm25).all() or np.any(pm25 < 0):
        raise ValueError("The subset must contain finite, nonnegative concentrations")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.output, pm25=pm25, wind=np.asarray(wind), xy=coordinates, times=times
    )


if __name__ == "__main__":
    main()
