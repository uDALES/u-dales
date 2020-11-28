""" Plot outputs for the current examples.
"""

from typing import List
from pathlib import Path
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


def plot_static(cases=List, output_dir=Path) -> None:
    for case in cases:
        print(f'Plotting case {case}...')
        ds = xr.open_dataset(str(output_dir / case / f'fielddump.{case}.nc'))
        _, ax = plt.subplots(1,2, figsize=(12,4))
        ds['u'].isel(yt=-1).isel(time=-1).plot(ax=ax[0])
        ds['w'].isel(yt=-1).isel(time=-1).plot(ax=ax[1])
        plt.savefig(output_dir / case / f'{case}.png')


def plot_animation(case=str, output_dir=Path) -> None:
    # Modified from
    # https://brushingupscience.com/2016/06/21/matplotlib-animations-the-easy-way/

    def animate(i=int) -> None:
        print(f'Drawing frame {i}')
        img.set_data(ds['u'].isel(yt=-1).isel(time=i))

    print(f'Plotting gif animaton for case {case}')
    ds = xr.open_dataset(output_dir / case / f'fielddump.{case}.nc')
    ds = ds.isel(time=slice(10, None)) # Remove the first timesteps
    fig, ax = plt.subplots(figsize=(4, 4))
    img = ax.imshow(ds['u'].isel(yt=64).isel(time=1))

    times = list(ds['u'].time)
    anim = FuncAnimation(
        fig, animate, interval=100, frames=len(times))

    ax.set_ylabel('z in m')
    ax.set_xlabel('x in m')
    plt.draw()
    plt.axis('off')
    plt.tight_layout()
    anim.save(output_dir / case / f'{case}.gif', writer='imagemagick')

if __name__ == "__main__":
    PROJ_DIR = Path(__file__).parents[2]
    output_dir = PROJ_DIR / 'outputs'
    cases = ['001', '101', '102', '201', '502']
    plot_static(cases=cases, output_dir=output_dir)
    #plot_animation(case='201_extended', output_dir=output_dir)
