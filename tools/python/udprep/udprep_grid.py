from __future__ import annotations

from typing import Any, Dict, List
import numpy as np
import sys
from scipy.optimize import fsolve
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for headless systems
import matplotlib.pyplot as plt

# Configure matplotlib to use LaTeX rendering
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Computer Modern Roman']

from .udprep import Section, SectionSpec

DEFAULTS: Dict[str, Any] = Section.load_defaults_json().get("grid", {})
FIELDS: List[str] = list(DEFAULTS.keys())


class GridSection(Section):
    def run_all(self) -> None:
        """Run grid generation preprocessing steps in the standard order."""
        steps = [
            ("generate_xygrid", self.generate_xygrid),
            ("generate_zgrid", self.generate_zgrid),
        ]
        self.run_steps("grid", steps)

    def generate_xygrid(self) -> None:
        """Create staggered x/y grids for preprocessing.
        
        Generates cell-centered (xt, yt) and face-centered (xm, ym) grids
        based on domain size and grid spacing.
        """
        # Cell-centered grids (staggered, starting at 0.5*dx)
        self.xt = np.arange(0.5 * self.dx, self.xlen, self.dx)
        self.yt = np.arange(0.5 * self.dy, self.ylen, self.dy)
        
        # Cell-face-centered grids (starting at 0)
        self.xm = np.arange(0, self.xlen + self.dx, self.dx)
        self.ym = np.arange(0, self.ylen + self.dy, self.dy)

    def generate_zgrid(self) -> None:
        """Create vertical grid.
        
        Generates cell-centered (zt) and cell-faces (zm) vertical grids
        along with cell spacing (dzt). Uses uniform or stretched grids
        depending on configuration flags.
        """
        if not self.lzstretch:
            # Uniform grid spacing
            self.zt = np.arange(0.5 * self.dz, self.zsize, self.dz)
            self.zm = np.arange(0, self.zsize + self.dz, self.dz)
            self.dzt = self.zm[1:] - self.zm[:-1]
        else:
            # Stretched grid - call appropriate stretching method
            if self.lstretchexp:
                self.stretch_exp()
            elif self.lstretchexpcheck:
                self.stretch_exp_check()
            elif self.lstretchtanh:
                self.stretch_tanh()
            elif self.lstretch2tanh:
                self.stretch_2tanh()
            else:
                print("="*67, file=sys.stderr)
                print("ERROR: Invalid stretch configuration", file=sys.stderr)
                print("="*67, file=sys.stderr)
                raise ValueError("Invalid stretch")
            
            # Plot and save dzt variation
            fig, ax = plt.subplots(figsize=(6, 4))
            ax.plot(self.dzt, linewidth=1.5, color='#2C3E50')
            ax.set_title(r'\textbf{Vertical Grid Spacing Variation}', fontsize=12)
            ax.set_xlabel(r'Grid Index $k$', fontsize=11)
            ax.set_ylabel(r'Grid Spacing $\Delta z$ [m]', fontsize=11)
            ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)
            ax.tick_params(labelsize=10)
            ax.set_xlim(0, len(self.dzt) - 1)
            fig.savefig('dz_variation.png', dpi=300, bbox_inches='tight')
            plt.close(fig)

    def stretch_exp(self) -> None:
        """Generate exponential stretch grid.
        
        Creates a vertical grid with linear spacing in the lower part and 
        exponentially stretched spacing in the upper part.
        """
        il = round(self.hlin / self.dzlin)
        ir = self.ktot - il
        
        # Initialize arrays
        self.zt = np.zeros(self.ktot)
        self.dzt = np.zeros(self.ktot)
        self.zm = np.zeros(self.ktot + 1)
        
        # Linear spacing in lower part
        self.zt[:il] = np.arange(0.5 * self.dzlin, self.hlin, self.dzlin)
        self.zm[:il + 1] = np.arange(0, self.hlin + self.dzlin, self.dzlin)
        
        # Exponential stretching in upper part
        gf = self.stretchconst
        
        while True:
            # Exponential spacing
            indices = np.arange(0, ir + 1)
            self.zm[il:] = self.zm[il] + (self.zsize - self.zm[il]) * (np.exp(gf * indices / ir) - 1) / (np.exp(gf) - 1)
            
            if (self.zm[il + 1] - self.zm[il]) < self.dzlin:
                gf -= 0.01  # Make sufficiently small steps to avoid an initial bump in dz
            else:
                if (self.zm[-1] - self.zm[-2]) > 3 * self.dzlin:
                    print("="*67, file=sys.stderr)
                    print("WARNING: Final grid spacing large - consider reducing domain height", file=sys.stderr)
                    print("="*67, file=sys.stderr)
                break
        
        # Calculate cell centers and spacings
        for i in range(self.ktot):
            self.zt[i] = (self.zm[i] + self.zm[i + 1]) / 2
            self.dzt[i] = self.zm[i + 1] - self.zm[i]

    def stretch_exp_check(self) -> None:
        """Generate exponential stretch grid after validating exponential stretch.
        
        Creates a vertical grid with exponential stretching using a root-finding
        procedure to ensure smooth grid transition and quality checks.
        """
        il = round(self.hlin / self.dzlin)
        ir = self.ktot - il
        z0 = il * self.dzlin  # hlin will be modified as z0
        
        # Initialize arrays
        self.zt = np.zeros(self.ktot)
        self.dzt = np.zeros(self.ktot)
        self.zm = np.zeros(self.ktot + 1)
        
        L = self.zsize - z0
        dxi = 1 / ir
        xi = np.arange(0, 1 + dxi, dxi)
        
        # Determine alpha using root finding
        # alpha / (exp(alpha)-1) = (dzlin*ir)/L
        def func(alpha):
            return alpha - (self.dzlin * ir) / L * (np.exp(alpha) - 1)
        
        alpha = fsolve(func, 1.0)[0]
        A = 1 / (np.exp(alpha) - 1)
        
        # Define grid functions
        zhat = A * (np.exp(alpha * xi) - 1)
        z = z0 + zhat * L
        
        # Create grid
        self.zm[:il + 1] = np.arange(0, z0 + self.dzlin, self.dzlin)  # Linear part
        self.zm[il:] = z  # Stretched part
        
        # Perform grid quality checks
        dz = np.diff(self.zm)
        stretch = dz[1:] / dz[:-1]
        
        if (np.min(stretch) < 0.95 or np.max(stretch) > 1.05):
            print("="*67, file=sys.stderr)
            print("WARNING: Generated grid is of bad quality", file=sys.stderr)
            print("Stretching factor = dz(n+1)/dz(n) should be between 0.95 and 1.05", file=sys.stderr)
            print(f"min value = {np.min(stretch):.3f}", file=sys.stderr)
            print(f"max value = {np.max(stretch):.3f}", file=sys.stderr)
            print("="*67, file=sys.stderr)
        
        # Warning if grid is refined near the top
        if alpha < 0:
            print("="*67, file=sys.stderr)
            print("WARNING: Possibly incorrect value for alpha", file=sys.stderr)
            print("The calculated value of alpha is less than zero, which implies", file=sys.stderr)
            print("you are refining the grid towards the domain top.", file=sys.stderr)
            print("="*67, file=sys.stderr)
        
        # Calculate cell centers and spacings
        self.zt = (self.zm[1:] + self.zm[:-1]) / 2.0
        self.dzt = self.zm[1:] - self.zm[:-1]

    def stretch_tanh(self) -> None:
        """Generate tanh-based stretch grid.
        
        Creates a vertical grid with linear spacing in the lower part and
        tanh-based stretched spacing in the upper part.
        """
        il = round(self.hlin / self.dzlin)
        ir = self.ktot - il
        
        # Initialize arrays
        self.zt = np.zeros(self.ktot)
        self.dzt = np.zeros(self.ktot)
        self.zm = np.zeros(self.ktot + 1)
        
        # Linear spacing in lower part
        self.zt[:il] = np.arange(0.5 * self.dzlin, self.hlin, self.dzlin)
        self.zm[:il + 1] = np.arange(0, self.hlin + self.dzlin, self.dzlin)
        
        # Tanh stretching in upper part
        gf = self.stretchconst
        
        while True:
            # Tanh-based spacing
            indices = np.arange(0, ir + 1)
            self.zm[il:] = self.zm[il] + (self.zsize - self.zm[il]) * (1 - np.tanh(gf * (1 - 2 * indices / (2 * ir))) / np.tanh(gf))
            
            if (self.zm[il + 1] - self.zm[il]) < self.dzlin:
                gf -= 0.01  # Make sufficiently small steps to avoid an initial bump in dz
            else:
                if (self.zm[-1] - self.zm[-2]) > 3 * self.dzlin:
                    print("="*67, file=sys.stderr)
                    print("WARNING: Final grid spacing large - consider reducing domain height", file=sys.stderr)
                    print("="*67, file=sys.stderr)
                break
        
        # Calculate cell centers and spacings
        for i in range(self.ktot):
            self.zt[i] = 0.5 * (self.zm[i] + self.zm[i + 1])
            self.dzt[i] = self.zm[i + 1] - self.zm[i]

    def stretch_2tanh(self) -> None:
        """Generate double-tanh stretch grid.
        
        Creates a vertical grid with linear spacing in the lower part and
        double-tanh stretched spacing in the upper part.
        """
        il = round(self.hlin / self.dzlin)
        ir = self.ktot - il
        
        # Initialize arrays
        self.zt = np.zeros(self.ktot)
        self.dzt = np.zeros(self.ktot)
        self.zm = np.zeros(self.ktot + 1)
        
        # Linear spacing in lower part
        self.zt[:il] = np.arange(0.5 * self.dzlin, self.hlin, self.dzlin)
        self.zm[:il + 1] = np.arange(0, self.hlin + self.dzlin, self.dzlin)
        
        # Double-tanh stretching in upper part
        gf = self.stretchconst
        
        while True:
            # Double-tanh based spacing
            indices = np.arange(0, ir + 1)
            self.zm[il:] = self.zm[il] + (self.zsize - self.zm[il]) / 2 * (1 - np.tanh(gf * (1 - 2 * indices / ir)) / np.tanh(gf))
            
            if (self.zm[il + 1] - self.zm[il]) < self.dzlin:
                gf -= 0.01  # Make sufficiently small steps to avoid an initial bump in dz
            else:
                if np.max(np.diff(self.zm)) > 3 * self.dzlin:
                    print("="*67, file=sys.stderr)
                    print("WARNING: Final grid spacing large - consider reducing domain height", file=sys.stderr)
                    print("="*67, file=sys.stderr)
                break
        
        # Calculate cell centers and spacings
        for i in range(self.ktot):
            self.zt[i] = (self.zm[i] + self.zm[i + 1]) / 2
            self.dzt[i] = self.zm[i + 1] - self.zm[i]


SPEC = SectionSpec(name="grid", fields=FIELDS, defaults=DEFAULTS, section_cls=GridSection)
