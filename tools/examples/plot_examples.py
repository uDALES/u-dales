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
        img.set_data(ds['u'].isel(yt=-1).isel(time=i))
        ax.set_title('Time: {:3.0f} s'.format(times[i].values))

    print(f'Plotting gif animaton for case {case}')
    ds = xr.open_dataset(output_dir / case / f'fielddump.{case}.nc')
    fig, ax = plt.subplots(figsize=(8, 4))
    img = ax.imshow(ds['u'].isel(yt=-1).isel(time=1))

    times = list(ds['u'].time)
    anim = FuncAnimation(
        fig, animate, interval=100, frames=len(times))

    ax.set_ylabel('z in m')
    ax.set_xlabel('x in m')
    plt.draw()
    anim.save(output_dir / case / f'{case}.gif', writer='imagemagick')

if __name__ == "__main__":
    cases = ['001', '101', '102', '201', '502']
    output_dir = Path(__file__).parents[2] / 'outputs'
    print(output_dir.resolve())
    plot_static(cases=cases, output_dir=output_dir)
    plot_animation(case='502', output_dir=output_dir)
