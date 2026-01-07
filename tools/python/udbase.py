"""
uDALES Post-Processing Module

Python implementation of the MATLAB udbase class for analyzing
uDALES simulation outputs.

Aug 2024, Maarten van Reeuwijk, Jingzi Huang. first version
Dec 2024, Maarten van Reeuwijk, Chris Wilson. added facet functionality.
Oct 2025, Maarten van Reeuwijk, Jingzi Huang, Dipanjan Majumdar. Major upgrade of facet functionality and tree handling.

Copyright (C) 2016- the uDALES Team.
"""

import numpy as np
import xarray as xr
from pathlib import Path
from typing import Optional, Union, Dict, Any, List
import warnings

# Import UDGeom from the udgeom package
try:
    from udgeom import UDGeom
except ImportError:
    # Fallback for old structure
    try:
        from .udgeom import UDGeom
    except ImportError:
        UDGeom = None
        warnings.warn("Could not import UDGeom. Geometry functionality will be limited.")


def _file_has_data(path: Path, skiprows: int = 0) -> bool:
    try:
        with path.open("r", encoding="ascii", errors="ignore") as f:
            for _ in range(skiprows):
                next(f, None)
            for line in f:
                stripped = line.strip()
                if stripped and not stripped.startswith("#"):
                    return True
    except OSError:
        return False
    return False


class UDBase:
    """
    Post-processing class for uDALES simulations.
    
    This class provides methods to load and analyze uDALES simulation outputs,
    including field data, statistics, facet data, and geometry.
    
    Parameters
    ----------
    expnr : int or str
        Experiment number
    path : str or Path, optional
        Path to experiment directory. Defaults to current directory.
    
    Examples
    --------
    >>> sim = UDBase(expnr=65, path='experiments/065')
    >>> u_field = sim.load_field('u')
    >>> stats = sim.load_stat_xyt('u')
    """
    
    def __init__(self, expnr: Union[int, str], path: Optional[Union[str, Path]] = None, load_geometry: bool = True):
        """Initialize UDBase instance."""
        
        # Store experiment number
        self.expnr = f"{int(expnr):03d}"
        
        # Set paths
        self.cpath = Path.cwd()
        self.path = Path(path) if path else self.cpath
        
        # File presence flags
        self._lfprof = True
        self._lfsolid = True
        self._lffacets = True
        self._lffactypes = True
        self._lffacetarea = True
        self._lffacet_sections = True
        self._lfgeom = True
        self._lftrees = True
        
        # Read namoptions file
        self._read_namoptions()
        
        # Load grid
        self._load_grid()
        
        # Load geometry if present (can be disabled for faster startup)
        if load_geometry:
            self._load_geometry()
        
        # Load solid masks if present
        self._load_solid_masks()
        
        # Load facet data if present
        self._load_facet_data()
        
        # Load tree data if present
        self._load_tree_data()
    
    def _read_namoptions(self):
        """
        Read namoptions file and store all parameters as attributes.
        
        Parses the namoptions.XXX file and converts values:
        - .true./.false. -> bool
        - Numbers -> int or float
        - Strings -> str (with quotes removed)
        """
        filepath = self.path / f"namoptions.{self.expnr}"
        
        if not filepath.exists():
            raise FileNotFoundError(f"namoptions.{self.expnr} not found in {self.path}")
        
        with open(filepath, 'r') as f:
            for line in f:
                # Skip comments and namelist headers
                if line.strip().startswith('&') or line.strip().startswith('!') or '=' not in line:
                    continue
                
                # Parse key=value pairs
                if '=' in line:
                    parts = line.split('=', 1)
                    key = parts[0].strip()
                    val = parts[1].split('!')[0].strip()  # Remove inline comments
                    
                    # Convert value types
                    if val.lower() == '.true.':
                        val = True
                    elif val.lower() == '.false.':
                        val = False
                    else:
                        # Try to parse as number
                        try:
                            if '.' in val or 'e' in val.lower():
                                val = float(val)
                            else:
                                val = int(val)
                        except ValueError:
                            # String value - remove quotes
                            val = val.strip("'\"")
                    
                    # Store as attribute
                    setattr(self, key, val)
        
        # Compute derived values
        if hasattr(self, 'xlen') and hasattr(self, 'itot'):
            self.dx = self.xlen / self.itot
        if hasattr(self, 'ylen') and hasattr(self, 'jtot'):
            self.dy = self.ylen / self.jtot
    
    def _load_grid(self):
        """
        Load grid coordinates.
        
        Attempts to load from prof.inp.XXX file. If not present,
        generates uniform grid based on domain size and number of cells.
        """
        prof_file = self.path / f"prof.inp.{self.expnr}"
        
        if prof_file.exists():
            try:
                # Load prof.inp - skip header lines
                data = np.loadtxt(prof_file, skiprows=1)
                
                # Cell centers (zt)
                self.zt = data[:, 0]
                
                # Cell edges (zm)
                zm_full = np.concatenate([[0], 0.5 * (self.zt[:-1] + self.zt[1:]), [self.zsize]])
                self.zm = zm_full[:-1]  # Match dimensions
                
                # Grid spacing
                self.dzt = np.diff(np.concatenate([self.zm, [self.zsize]]))
                self.dzm = np.concatenate([[2 * self.zt[0]], np.diff(self.zt)])
                
            except Exception as e:
                warnings.warn(f"Error loading prof.inp.{self.expnr}: {e}. Using uniform grid.")
                self._lfprof = False
                self._generate_uniform_zgrid()
        else:
            warnings.warn(f"prof.inp.{self.expnr} not found. Assuming equidistant grid.")
            self._lfprof = False
            self._generate_uniform_zgrid()
        
        # Generate x and y grids
        self.xm = np.arange(self.itot) * self.dx
        self.xt = self.xm + 0.5 * self.dx
        
        self.ym = np.arange(self.jtot) * self.dy
        self.yt = self.ym + 0.5 * self.dy
    
    def _generate_uniform_zgrid(self):
        """Generate uniform z-grid when prof.inp is not available."""
        if hasattr(self, 'zsize') and hasattr(self, 'ktot'):
            dz = self.zsize / self.ktot
            
            self.zm = np.arange(self.ktot) * dz
            self.zt = self.zm + 0.5 * dz
            
            self.dzm = np.full(self.ktot, dz)
            self.dzt = np.full(self.ktot, dz)
        else:
            warnings.warn("Cannot generate z-grid: zsize or ktot not found in namoptions")
    
    def _load_geometry(self):
        """Load STL geometry file if present."""
        if hasattr(self, 'stl_file'):
            stl_path = self.path / self.stl_file
            if stl_path.exists():
                try:
                    # Import udgeom module
                    from udgeom import UDGeom
                    self.geom = UDGeom(self.path)
                    self.geom.load(self.stl_file)
                except ImportError as e:
                    warnings.warn(f"Cannot load geometry: {e}. Install trimesh: pip install trimesh")
                    self.geom = None
                except Exception as e:
                    warnings.warn(f"Error loading geometry: {e}")
                    self.geom = None
            else:
                warnings.warn(f"STL file {self.stl_file} not found.")
                self.geom = None
        else:
            self.geom = None
    
    def _load_solid_masks(self):
        """Load solid point arrays for u, v, w, c grids."""
        for grid_type in ['u', 'v', 'w', 'c']:
            solid_file = self.path / f"solid_{grid_type}.txt"
            
            if solid_file.exists():
                try:
                    if not _file_has_data(solid_file, skiprows=1):
                        mask = np.zeros((self.itot, self.jtot, self.ktot), dtype=bool)
                        setattr(self, f'S{grid_type}', mask)
                        continue
                    # Read indices (1-based from Fortran)
                    indices = np.loadtxt(solid_file, skiprows=1, dtype=int)
                    if indices.size == 0:
                        mask = np.zeros((self.itot, self.jtot, self.ktot), dtype=bool)
                        setattr(self, f'S{grid_type}', mask)
                        continue

                    # Create boolean mask
                    mask = np.zeros((self.itot, self.jtot, self.ktot), dtype=bool)
                    
                    if indices.ndim == 1:
                        indices = indices.reshape(1, -1)
                    
                    # Convert to 0-based indexing
                    mask[indices[:, 0] - 1, indices[:, 1] - 1, indices[:, 2] - 1] = True
                    
                    # Store as attribute (Su, Sv, Sw, Sc)
                    setattr(self, f'S{grid_type}', mask)
                    
                except Exception as e:
                    warnings.warn(f"Error loading solid_{grid_type}.txt: {e}")
                    setattr(self, f'S{grid_type}', None)
                    self._lfsolid = False
            else:
                setattr(self, f'S{grid_type}', None)
    
    def _load_facet_data(self):
        """Load facet information if present."""
        self.facs = {}
        self.factypes = {}
        
        # Load facet areas
        facetarea_file = self.path / f"facetarea.inp.{self.expnr}"
        if facetarea_file.exists():
            try:
                self.facs['area'] = np.loadtxt(facetarea_file, skiprows=1)
            except Exception as e:
                warnings.warn(f"Error loading facetarea.inp.{self.expnr}: {e}")
                self._lffacetarea = False
        else:
            self._lffacetarea = False
        
        # Load facet types
        facets_file = self.path / f"facets.inp.{self.expnr}"
        if facets_file.exists():
            try:
                data = np.loadtxt(facets_file, skiprows=1)
                self.facs['typeid'] = data[:, 0].astype(int)
                self.facs['normals'] = data[:, 1:4]  # Surface normals
            except Exception as e:
                warnings.warn(f"Error loading facets.inp.{self.expnr}: {e}")
                self._lffacets = False
        else:
            self._lffacets = False
        
        # Load facet type properties
        factypes_file = self.path / f"factypes.inp.{self.expnr}"
        if factypes_file.exists():
            try:
                data = np.loadtxt(factypes_file, skiprows=3)
                
                if not hasattr(self, 'nfaclyrs'):
                    self.nfaclyrs = 3  # Default
                
                self.factypes['id'] = data[:, 0].astype(int)
                self.factypes['lGR'] = data[:, 1].astype(bool)
                self.factypes['z0'] = data[:, 2]
                self.factypes['z0h'] = data[:, 3]
                self.factypes['al'] = data[:, 4]  # Albedo
                self.factypes['em'] = data[:, 5]  # Emissivity
                self.factypes['d'] = data[:, 6:6 + self.nfaclyrs]  # Layer thicknesses
                self.factypes['C'] = data[:, 6 + self.nfaclyrs:6 + 2 * self.nfaclyrs]  # Heat capacity
                self.factypes['lam'] = data[:, 6 + 2 * self.nfaclyrs:6 + 3 * self.nfaclyrs + 1]  # Conductivity
                
                # Add names for common wall types
                wall_names = {
                    0: "Default dummy",
                    -1: "Asphalt floor",
                    -101: "Concrete bounding wall",
                    1: "Concrete",
                    2: "Brick",
                    3: "Stone",
                    4: "Painted wood",
                    11: "Green 1",
                    12: "Green 2"
                }
                
                self.factypes['name'] = [
                    wall_names.get(int(id_), f"Custom walltype {int(id_)}")
                    for id_ in self.factypes['id']
                ]
                
            except Exception as e:
                warnings.warn(f"Error loading factypes.inp.{self.expnr}: {e}")
                self._lffactypes = False
        else:
            self._lffactypes = False
        
        # Load facet sections
        self._load_facet_sections()
    
    def _load_facet_sections(self):
        """Load facet section information for u, v, w, c grids."""
        self.facsec = {}
        
        for grid_type in ['u', 'v', 'w', 'c']:
            facsec_file = self.path / f"facet_sections_{grid_type}.txt"
            fluid_boundary_file = self.path / f"fluid_boundary_{grid_type}.txt"
            
            if facsec_file.exists() and fluid_boundary_file.exists():
                try:
                    if not _file_has_data(facsec_file, skiprows=1) or not _file_has_data(fluid_boundary_file, skiprows=1):
                        self._lffacet_sections = False
                        continue
                    # Load facet section data
                    facsec_data = np.loadtxt(facsec_file, skiprows=1)
                    
                    # Load fluid boundary locations
                    fluid_boundary = np.loadtxt(fluid_boundary_file, skiprows=1, dtype=int)
                    if facsec_data.size == 0 or fluid_boundary.size == 0:
                        self._lffacet_sections = False
                        continue

                    if facsec_data.ndim == 1:
                        facsec_data = facsec_data.reshape(1, -1)
                    if fluid_boundary.ndim == 1:
                        fluid_boundary = fluid_boundary.reshape(1, -1)
                    
                    # Store in structure
                    self.facsec[grid_type] = {
                        'facid': facsec_data[:, 0].astype(int) - 1,
                        'area': facsec_data[:, 1],
                        'locs': fluid_boundary[facsec_data[:, 2].astype(int) - 1, :].astype(int) - 1,
                        'distance': facsec_data[:, 3]
                    }
                    
                except Exception as e:
                    warnings.warn(f"Error loading facet_sections_{grid_type}.txt: {e}")
                    self._lffacet_sections = False
            else:
                self._lffacet_sections = False
    
    def _load_tree_data(self):
        """Load tree bounding box information if present."""
        trees_file = self.path / f"trees.inp.{self.expnr}"
        
        if trees_file.exists():
            try:
                # Load tree data (skip first 2 header lines)
                self.trees = np.loadtxt(trees_file, skiprows=2, dtype=int)
                if self.trees.ndim == 1:
                    self.trees = self.trees.reshape(1, -1)
                self.trees = self.trees.astype(int) - 1
            except Exception as e:
                warnings.warn(f"Error loading trees.inp.{self.expnr}: {e}")
                self._lftrees = False
                self.trees = None
        else:
            self._lftrees = False
            self.trees = None
    
    def __repr__(self):
        """String representation of UDBase object."""
        info = [
            f"UDBase(expnr='{self.expnr}')",
            f"  path: {self.path}",
        ]
        
        if hasattr(self, 'itot'):
            info.append(f"  grid: {self.itot} x {self.jtot} x {self.ktot}")
        
        if hasattr(self, 'xlen'):
            info.append(f"  domain: {self.xlen} x {self.ylen} x {self.zsize}")
        
        if self.geom is not None:
            info.append(f"  geometry: loaded")
        
        if hasattr(self, 'nfcts'):
            info.append(f"  facets: {self.nfcts}")
        
        return "\n".join(info)
    
    # ===== Data Loading Methods =====
    
    def load_field(self, var: Optional[str] = None) -> Union[xr.Dataset, np.ndarray]:
        """
        Load 3D instantaneous field data.
        
        Parameters
        ----------
        var : str, optional
            Variable name to load. If None, displays available variables and returns full dataset.
        
        Returns
        -------
        xarray.Dataset or numpy.ndarray
            Full dataset if var is None (for browsing), otherwise numpy array of the requested variable.
        
        Examples
        --------
        >>> sim.load_field()  # Display available variables
        >>> u = sim.load_field('u')  # Returns numpy array
        >>> print(u.shape)
        """
        filename = self.path / f"fielddump.{self.expnr}.nc"
        return self._load_ncdata(filename, var)
    
    def load_stat_xyt(self, var: Optional[str] = None) -> Union[xr.Dataset, np.ndarray]:
        """
        Load slab-averaged (xy-plane) statistics.
        
        Parameters
        ----------
        var : str, optional
            Variable name to load. If None, displays available variables and returns full dataset.
        
        Returns
        -------
        xarray.Dataset or numpy.ndarray
            Full dataset if var is None (for browsing), otherwise numpy array of the requested variable.
        
        Examples
        --------
        >>> sim.load_stat_xyt()  # Display available variables
        >>> u_avg = sim.load_stat_xyt('u')  # Returns numpy array
        """
        filename = self.path / f"xytdump.{self.expnr}.nc"
        return self._load_ncdata(filename, var)
    
    def load_stat_t(self, var: Optional[str] = None) -> Union[xr.Dataset, np.ndarray]:
        """
        Load time-averaged 3D statistics.
        
        Parameters
        ----------
        var : str, optional
            Variable name to load. If None, displays available variables and returns full dataset.
        
        Returns
        -------
        xarray.Dataset or numpy.ndarray
            Full dataset if var is None (for browsing), otherwise numpy array of the requested variable.
        
        Examples
        --------
        >>> sim.load_stat_t()  # Display available variables
        >>> u_tavg = sim.load_stat_t('u')  # Returns numpy array
        """
        filename = self.path / f"tdump.{self.expnr}.nc"
        return self._load_ncdata(filename, var)
    
    def load_stat_tree(self, var: Optional[str] = None) -> Union[xr.Dataset, np.ndarray]:
        """
        Load time-averaged statistics of tree source terms.
        
        Retrieves tree drag, heat, and moisture source terms from the treedump file.
        
        Parameters
        ----------
        var : str, optional
            Variable name to load. If None, displays available variables and returns full dataset.
        
        Returns
        -------
        xarray.Dataset or numpy.ndarray
            Full dataset if var is None (for browsing), otherwise numpy array of the requested variable.
        
        Examples
        --------
        >>> sim.load_stat_tree()  # Display available variables
        >>> tree_drag = sim.load_stat_tree('tree_drag_u')  # Returns numpy array
        """
        filename = self.path / f"treedump.{self.expnr}.nc"
        return self._load_ncdata(filename, var)
    
    def load_slice(self, plane: str, var: Optional[str] = None) -> Union[xr.Dataset, np.ndarray]:
        """
        Load 2D slice data.
        
        Parameters
        ----------
        plane : str
            Slice plane: 'i', 'j', or 'k' for x, y, or z slices
        var : str, optional
            Variable name to load. If None, displays available variables and returns full dataset.
        
        Returns
        -------
        xarray.Dataset or numpy.ndarray
            Full dataset if var is None (for browsing), otherwise numpy array of the requested variable.
        
        Examples
        --------
        >>> sim.load_slice('k')  # Display available variables in horizontal slice
        >>> u_slice = sim.load_slice('k', 'u')  # Returns numpy array
        """
        if plane not in ['i', 'j', 'k']:
            raise ValueError("plane must be 'i', 'j', or 'k'")
        
        filename = self.path / f"{plane}slicedump.{self.expnr}.nc"
        return self._load_ncdata(filename, var)
    
    def _load_ncdata(self, filename: Path, var: Optional[str]) -> Union[xr.Dataset, np.ndarray]:
        """
        Helper method to load NetCDF data using xarray.
        
        Automatically reverses dimension order to match MATLAB conventions.
        When a specific variable is requested, returns it as a numpy array for
        memory efficiency with large datasets.
        
        NetCDF files use C-style (row-major) dimension ordering, while MATLAB
        uses Fortran-style (column-major) ordering. Due to this fundamental
        difference, MATLAB automatically reverses dimension order when reading
        NetCDF files. This method reverses dimensions to match MATLAB's behavior.
        
        Examples of transformations:
        - (time, zt) -> (zt, time)
        - (time, zt, yt, xm) -> (xm, yt, zt, time)
        - (time, fct) -> (fct, time)
        - (time, lyr, fct) -> (fct, lyr, time)
        
        Parameters
        ----------
        filename : Path
            Path to NetCDF file
        var : str, optional
            Variable to extract. If None, displays available variables.
        
        Returns
        -------
        xarray.Dataset or numpy.ndarray
            Full dataset if var is None (for browsing), numpy array if var is specified.
        """
        if not filename.exists():
            raise FileNotFoundError(f"File not found: {filename}")
        
        ds = xr.open_dataset(filename)
        
        # Transpose all data variables to match MATLAB's column-major convention
        # Reverse dimension order for all variables with 2+ dimensions
        transposed_vars = {}
        for var_name in ds.data_vars:
            data_var = ds[var_name]
            if len(data_var.dims) >= 2:
                transposed_vars[var_name] = data_var.transpose(*reversed(data_var.dims))
            else:
                transposed_vars[var_name] = data_var
        
        # Create new dataset with transposed variables
        ds_transposed = xr.Dataset(transposed_vars, coords=ds.coords, attrs=ds.attrs)
        
        if var is None:
            # Display file contents
            self._display_ncinfo(ds_transposed, filename.name)
            return ds_transposed
        else:
            if var not in ds_transposed:
                raise KeyError(f"Variable '{var}' not found in {filename.name}")
            
            # Return as numpy array for consistency with MATLAB and memory efficiency
            return ds_transposed[var].values
    
    def _display_ncinfo(self, ds: xr.Dataset, filename: str):
        """Display information about NetCDF dataset."""
        print(f"\nContents of {filename}:")
        print(f"{'Name':<20} {'Dimensions':<30} {'Shape':<20}")
        print("-" * 70)
        
        for var_name in sorted(ds.data_vars):
            var = ds[var_name]
            dims = ', '.join(var.dims)
            shape = ' x '.join(map(str, var.shape))
            print(f"{var_name:<20} {dims:<30} {shape:<20}")
    
    # ===== Facet Data Loading Methods =====
    
    def load_fac_momentum(self, var: Optional[str] = None) -> Union[xr.Dataset, np.ndarray]:
        """
        Load facet momentum data (pressure and shear stress).
        
        Parameters
        ----------
        var : str, optional
            Variable name to load. If None, displays available variables and returns full dataset.
        
        Returns
        -------
        xarray.Dataset or numpy.ndarray
            Full dataset if var is None (for browsing), otherwise numpy array of the requested variable.
        
        Examples
        --------
        >>> sim.load_fac_momentum()  # Display available variables
        >>> pressure = sim.load_fac_momentum('pres')  # Returns numpy array
        """
        filename = self.path / f"fac.{self.expnr}.nc"
        return self._load_ncdata(filename, var)
    
    def load_fac_eb(self, var: Optional[str] = None) -> Union[xr.Dataset, np.ndarray]:
        """
        Load facet surface energy balance data.
        
        Parameters
        ----------
        var : str, optional
            Variable name to load. If None, displays available variables and returns full dataset.
        
        Returns
        -------
        xarray.Dataset or numpy.ndarray
            Full dataset if var is None (for browsing), otherwise numpy array of the requested variable.
        
        Examples
        --------
        >>> sim.load_fac_eb()  # Display available variables
        >>> H = sim.load_fac_eb('hf')  # Sensible heat flux - Returns numpy array
        >>> K = sim.load_fac_eb('netsw')  # Net shortwave - Returns numpy array
        """
        filename = self.path / f"facEB.{self.expnr}.nc"
        return self._load_ncdata(filename, var)
    
    def load_fac_temperature(self, var: Optional[str] = None) -> Union[xr.Dataset, np.ndarray]:
        """
        Load facet temperature data.
        
        Parameters
        ----------
        var : str, optional
            Variable name to load. If None, displays available variables and returns full dataset.
        
        Returns
        -------
        xarray.Dataset or numpy.ndarray
            Full dataset if var is None (for browsing), otherwise numpy array of the requested variable.
        
        Examples
        --------
        >>> sim.load_fac_temperature()  # Display available variables
        >>> T = sim.load_fac_temperature('T')  # Temperature in layers - Returns numpy array
        >>> dTdz = sim.load_fac_temperature('dTdz')  # Temperature gradient - Returns numpy array
        """
        filename = self.path / f"facT.{self.expnr}.nc"
        return self._load_ncdata(filename, var)
    
    def load_seb(self) -> Dict[str, np.ndarray]:
        """
        Load all surface energy balance terms.
        
        Returns
        -------
        dict
            Dictionary containing all SEB terms:
            - Kstar: Net shortwave radiation
            - Lstar: Net longwave radiation
            - Lin: Incoming longwave
            - Lout: Outgoing longwave
            - H: Sensible heat flux
            - E: Latent heat flux
            - G: Ground heat flux
            - Tsurf: Surface temperature
            - t: Time array
        
        Examples
        --------
        >>> seb = sim.load_seb()
        >>> print(seb['Kstar'].shape)
        >>> seb_avg = sim.area_average_seb(seb)
        """
        # Load energy balance terms
        try:
            t = self.load_fac_eb('t')
        except FileNotFoundError:
            raise FileNotFoundError(
                f"Surface energy balance file facEB.{self.expnr}.nc not found in {self.path}"
            )
        K = self.load_fac_eb('netsw')
        Lin = self.load_fac_eb('LWin')
        Lout = self.load_fac_eb('LWout')
        H = self.load_fac_eb('hf')
        E = self.load_fac_eb('ef')
        
        # Load temperature data for ground heat flux calculation
        T = self.load_fac_temperature('T')
        dTdz = self.load_fac_temperature('dTdz')
        lam = self.assign_prop_to_fac('lam')
        
        # Calculate derived quantities
        L = Lin - Lout
        # Data now comes as (fct, lyr, time) from _load_ncdata, transposed to match MATLAB
        # lam is (n_facets, n_layers), dTdz is (n_facets, n_layers, n_time)
        G = -lam[:, 0, np.newaxis] * dTdz[:, 0, :]  # Ground heat flux: (n_facets, n_time)
        Tsurf = T[:, 0, :]  # Surface temperature: (n_facets, n_time)
        
        # All data is already in MATLAB convention: (n_facets, n_time)
        # _load_ncdata transposes from NetCDF (time, fct) to (fct, time)
        return {
            'Kstar': K,
            'Lstar': L,
            'Lin': Lin,
            'Lout': Lout,
            'H': -H,  # Sign convention
            'E': -E,  # Sign convention
            'G': G,
            'Tsurf': Tsurf,
            't': t
        }
    
    # ===== Facet Analysis Methods =====
    
    def assign_prop_to_fac(self, prop: str) -> np.ndarray:
        """
        Assign facet type property to individual facets.
        
        Parameters
        ----------
        prop : str
            Property name from factypes (e.g., 'al', 'em', 'z0', 'lam')
        
        Returns
        -------
        ndarray
            Property values for each facet
        
        Examples
        --------
        >>> albedo = sim.assign_prop_to_fac('al')
        >>> emissivity = sim.assign_prop_to_fac('em')
        """
        if not self._lffactypes:
            raise ValueError(f"factypes.inp.{self.expnr} required for this method")
        
        if not self._lffacets:
            raise ValueError(f"facets.inp.{self.expnr} required for this method")
        
        if prop not in self.factypes:
            raise KeyError(f"Property '{prop}' not found in factypes")
        
        # Get property values and type IDs
        prop_vals = self.factypes[prop]
        type_ids = self.factypes['id']
        facet_types = self.facs['typeid']
        
        # Map type IDs to property values
        if prop_vals.ndim == 1:
            # 1D property (e.g., albedo)
            result = np.full_like(facet_types, np.nan, dtype=float)
            for i, tid in enumerate(type_ids):
                result[facet_types == tid] = prop_vals[i]
        else:
            # 2D property (e.g., layer thicknesses)
            result = np.full((len(facet_types), prop_vals.shape[1]), np.nan)
            for i, tid in enumerate(type_ids):
                result[facet_types == tid, :] = prop_vals[i, :]
        
        return result
    
    def area_average_fac(self, var: np.ndarray, sel: Optional[np.ndarray] = None) -> np.ndarray:
        """
        Area-weighted average of facet data.
        
        Parameters
        ----------
        var : ndarray
            Facet variable to average. Can be:
            - 1D array with shape (n_facets,)
            - 2D array with shape (n_facets, n_times) or (n_times, n_facets)
            Method auto-detects which dimension is facets based on size.
        sel : ndarray, optional
            Boolean mask or indices to select subset of facets
        
        Returns
        -------
        ndarray
            Area-averaged values. Same shape as input except facet dimension is removed.
        
        Examples
        --------
        >>> H_avg = sim.area_average_fac(H)
        >>> # Average over specific facets
        >>> roof_facets = sim.facs['typeid'] == 1
        >>> H_avg_roof = sim.area_average_fac(H, roof_facets)
        """
        if not self._lffacetarea:
            raise ValueError(f"facetarea.inp.{self.expnr} required for this method")
        
        areas = self.facs['area']
        n_facets = len(areas)
        
        if sel is None:
            sel = slice(None)
            sel_areas = areas
        else:
            sel_areas = areas[sel]
        
        # Handle different array dimensions
        if var.ndim == 1:
            # Simple 1D case
            total_var = np.sum(var[sel] * sel_areas)
            total_area = np.sum(sel_areas)
            return total_var / total_area
        elif var.ndim == 2:
            # Determine which axis is facets
            if var.shape[0] == n_facets:
                # Shape (n_facets, n_other) - facets in first dimension
                if isinstance(sel, slice):
                    total_var = np.sum(var * areas[:, np.newaxis], axis=0)
                else:
                    total_var = np.sum(var[sel, :] * sel_areas[:, np.newaxis], axis=0)
                total_area = np.sum(sel_areas)
                return total_var / total_area
            elif var.shape[1] == n_facets:
                # Shape (n_other, n_facets) - facets in second dimension  
                if isinstance(sel, slice):
                    total_var = np.sum(var * areas[np.newaxis, :], axis=1)
                else:
                    total_var = np.sum(var[:, sel] * sel_areas[np.newaxis, :], axis=1)
                total_area = np.sum(sel_areas)
                return total_var / total_area
            else:
                raise ValueError(f"Neither dimension of var {var.shape} matches number of facets {n_facets}")
        else:
            raise ValueError(f"area_average_fac only supports 1D and 2D arrays, got shape {var.shape}")
    
    def area_average_seb(self, seb: Dict[str, np.ndarray]) -> Dict[str, np.ndarray]:
        """
        Area-average all surface energy balance terms.
        
        Parameters
        ----------
        seb : dict
            Dictionary from load_seb()
        
        Returns
        -------
        dict
            Area-averaged SEB terms
        
        Examples
        --------
        >>> seb = sim.load_seb()
        >>> seb_avg = sim.area_average_seb(seb)
        >>> import matplotlib.pyplot as plt
        >>> plt.plot(seb_avg['t'], seb_avg['Kstar'])
        """
        return {
            'Kstar': self.area_average_fac(seb['Kstar']),
            'Lstar': self.area_average_fac(seb['Lstar']),
            'Lin': self.area_average_fac(seb['Lin']),
            'Lout': self.area_average_fac(seb['Lout']),
            'H': self.area_average_fac(seb['H']),
            'E': self.area_average_fac(seb['E']),
            'G': self.area_average_fac(seb['G']),
            't': seb['t']
        }
    
    @staticmethod
    def time_average(var: np.ndarray, other: Optional[np.ndarray] = None):
        """
        Time-average variables.
        
        Parameters
        ----------
        var : ndarray
            Variable to average. Time must be the last dimension.
        other : ndarray, optional
            Second variable (same shape as ``var``) to compute covariance with.
        
        Returns
        -------
        tuple
            - (mean, variance) when only ``var`` is provided
            - (x_mean, y_mean, covariance) when ``other`` is provided
        
        Notes
        -----
        Variance/covariance are computed by delegating to ``merge_stat`` with
        zero instantaneous contributions.
        """
        X = np.asarray(var)
        n = X.shape[-1]
        
        if other is None:
            zeros = np.zeros_like(X)
            return UDBase.merge_stat(X, zeros, n)
        
        Y = np.asarray(other)
        zeros = np.zeros_like(X)
        return UDBase.merge_stat(X, Y, zeros, n)

    @staticmethod
    def merge_stat(X: np.ndarray, *args, Y: Optional[np.ndarray] = None,
                   XpXp: Optional[np.ndarray] = None,
                   XpYp: Optional[np.ndarray] = None):
        """
        Merge short-term statistics into longer windows.

        Supported call patterns:
          - merge_stat(X, n)
          - merge_stat(X, XpXp, n)                  
          - merge_stat(X, Y, XpYp, n)               
          - merge_stat(X, n, XpXp=array)            # single variable mean/variance
          - merge_stat(X, n, Y=array, XpYp=array)   # two variables mean/cov

        Parameters
        ----------
        X : ndarray
            First variable. Final dimension is time or short-term windows.
        args : tuple
            Positional parameters parsed according to the patterns above.
        Y : ndarray, optional
            Second variable (same trailing dimension as ``X``).
        XpXp : ndarray, optional
            Variance contribution for ``X`` aligned with ``X``.
        XpYp : ndarray, optional
            Covariance contribution for ``X`` and ``Y`` aligned with ``X``.

        Returns
        -------
        tuple | ndarray
            ``Xmean`` if only ``X`` provided;
            ``Xmean, var`` if ``XpXp`` given;
            ``Xmean, Ymean, cov`` if ``Y`` provided (and optionally ``XpYp``).

        Notes
        -----
        Discards the oldest samples that do not fill a complete window so the
        most recent data is retained. Variance/covariance combine the mean of
        short-window contributions with the variance/covariance of the short
        means inside each merged window.
        """
        # Parse positional arguments to support both MATLAB and Python styles
        n = None
        X = np.asarray(X)
        if len(args) == 1:
            n = int(args[0])
        elif len(args) == 2 and Y is None:
            # MATLAB style: merge_stat(X, XpXp, n)
            XpXp = np.asarray(args[0])
            n = int(args[1])
        elif len(args) == 3:
            # MATLAB style: merge_stat(X, Y, XpYp, n)
            Y = np.asarray(args[0])
            XpYp = np.asarray(args[1])
            n = int(args[2])
        else:
            raise ValueError("merge_stat expects 1, 2, or 3 positional arguments after X")

        if n <= 0:
            raise ValueError("n must be positive")
        if X.shape[-1] < n:
            raise ValueError("Not enough samples to form a single merged window")

        nwin = X.shape[-1] // n
        start = X.shape[-1] - nwin * n  # discard oldest incomplete window
        X_use = X[..., start:]
        X_group = X_use.reshape(*X.shape[:-1], nwin, n)
        Xmean = X_group.mean(axis=-1)

        if Y is None:
            if XpXp is None:
                return Xmean

            XpXp = np.asarray(XpXp)
            if XpXp.shape[-1] != X.shape[-1]:
                raise ValueError("XpXp must match X shape in the last dimension")
            XpXp_use = XpXp[..., start:]
            XpXp_group = XpXp_use.reshape(*XpXp.shape[:-1], nwin, n)
            within = XpXp_group.mean(axis=-1)
            between = ((X_group - Xmean[..., None]) ** 2).mean(axis=-1)
            return Xmean, within + between

        Y = np.asarray(Y)
        if Y.shape[-1] < n:
            raise ValueError("Y does not have enough samples to form a merged window")
        if Y.shape[-1] != X.shape[-1]:
            raise ValueError("X and Y must share the same length in the last dimension")

        Y_use = Y[..., start:]
        Y_group = Y_use.reshape(*Y.shape[:-1], nwin, n)
        Ymean = Y_group.mean(axis=-1)

        if XpYp is None:
            cov_within = ((X_group - Xmean[..., None]) * (Y_group - Ymean[..., None])).mean(axis=-1)
            between = 0.0
        else:
            XpYp = np.asarray(XpYp)
            if XpYp.shape[-1] != X.shape[-1]:
                raise ValueError("XpYp must match X and Y in the last dimension")
            XpYp_use = XpYp[..., start:]
            XpYp_group = XpYp_use.reshape(*XpYp.shape[:-1], nwin, n)
            cov_within = XpYp_group.mean(axis=-1)
            between = ((X_group - Xmean[..., None]) * (Y_group - Ymean[..., None])).mean(axis=-1)

        cov = cov_within + between
        return Xmean, Ymean, cov

    @staticmethod
    def coarsegrain_field(var: np.ndarray, Lflt: np.ndarray,
                          xm: np.ndarray, ym: np.ndarray) -> np.ndarray:
        """
        Apply 2D periodic box filters to a 3D field.

        Parameters
        ----------
        var : ndarray, shape (nx, ny, nz)
            Field data (time or height on the last axis).
        Lflt : array-like
            Filter lengths (meters). Multiple lengths are allowed.
        xm, ym : array-like
            Grid coordinates in meters; must be 1D and roughly uniform.

        Returns
        -------
        ndarray
            Filtered field with shape (nx, ny, nz, n_filters).
        """
        var = np.asarray(var)
        if var.ndim != 3:
            raise ValueError("var must be 3D with shape (nx, ny, nz)")

        xm = np.asarray(xm).ravel()
        ym = np.asarray(ym).ravel()
        if xm.size < 2 or ym.size < 2:
            raise ValueError("xm and ym must contain at least two points")

        dx = float(np.mean(np.diff(xm)))
        dy = float(np.mean(np.diff(ym)))
        if dx <= 0 or dy <= 0:
            raise ValueError("Grid spacings must be positive")

        Lflt_arr = np.atleast_1d(Lflt)
        nx, ny, nz = var.shape
        out = np.empty((nx, ny, nz, len(Lflt_arr)))

        # Build filters using half-width convention (matches MATLAB implementation)
        ii, jj = np.meshgrid(np.arange(nx), np.arange(ny), indexing="ij")
        di = np.minimum(ii, nx - ii)  # periodic distance in i
        dj = np.minimum(jj, ny - jj)  # periodic distance in j

        for i, L in enumerate(Lflt_arr):
            ngx = max(int(round((L / dx) / 2.0)), 1)
            ngy = max(int(round((L / dy) / 2.0)), 1)
            mask = (di <= ngx) & (dj <= ngy)
            kernel = mask.astype(float)
            kernel /= kernel.sum()
            k_hat = np.fft.fftn(kernel)

            for k in range(nz):
                v_hat = np.fft.fftn(var[:, :, k])
                out[:, :, k, i] = np.real(np.fft.ifftn(v_hat * k_hat))

        return out
    
    def plot_trees(self, show: bool = True):
        """
        Plot tree volumetric regions on top of the geometry using trimesh/plotly,
        matching the rendering pipeline used by plot_fac.
        
        Parameters
        ----------
        show : bool, default=True
            Display the plot immediately. If False, return the figure without showing.
        """
        if not self._lfgeom or self.geom is None:
            raise ValueError("Geometry (STL) file required for plot_trees()")
        if not self._lftrees or self.trees is None:
            raise ValueError("trees.inp file required for plot_trees()")
        
        try:
            import trimesh
        except ImportError:
            raise ImportError("trimesh is required for visualization. Install with: pip install trimesh")
        
        # Base geometry mesh (light gray)
        base_mesh = self.geom.stl.copy()
        base_color = np.array([220, 220, 220, 255], dtype=np.uint8)
        base_mesh.visual.face_colors = np.tile(base_color, (len(base_mesh.faces), 1))
        
        tree_meshes = []
        tree_color = np.array([34, 139, 34, 140], dtype=np.uint8)  # forest green, semi-transparent
        
        for i in range(self.trees.shape[0]):
            il, iu, jl, ju, kl, ku = self.trees[i, :]
            # Convert 1-based to 0-based
            il, iu, jl, ju, kl, ku = il - 1, iu - 1, jl - 1, ju - 1, kl - 1, ku - 1
            
            xmin, xmax = self.xm[il], self.xm[iu] + self.dx
            ymin, ymax = self.ym[jl], self.ym[ju] + self.dy
            zmin, zmax = self.zm[kl], self.zm[ku] + self.dzm[ku]
            
            extents = np.array([xmax - xmin, ymax - ymin, zmax - zmin])
            center = np.array([xmin + extents[0] / 2, ymin + extents[1] / 2, zmin + extents[2] / 2])
            
            box = trimesh.creation.box(extents=extents, transform=trimesh.transformations.translation_matrix(center))
            box.visual.face_colors = np.tile(tree_color, (len(box.faces), 1))
            tree_meshes.append(box)
        
        meshes = [base_mesh] + tree_meshes if tree_meshes else [base_mesh]
        
        # Build full edge list so all facet edges (front and back) are shown
        faces = self.geom.stl.faces
        edges = set()
        for tri in faces:
            e0 = tuple(sorted((tri[0], tri[1])))
            e1 = tuple(sorted((tri[1], tri[2])))
            e2 = tuple(sorted((tri[2], tri[0])))
            edges.update([e0, e1, e2])
        outline_edges = list(edges)
        
        fig = self._render_scene(meshes, show_outlines=True, custom_edges=outline_edges, show=show)
        if fig is not None:
            fig.update_layout(title=f'Geometry with Trees ({len(tree_meshes)} regions)')
        return fig
    
    def plot_fac(self, var: np.ndarray, building_ids: Optional[np.ndarray] = None, show: bool = True):
        """
        Plot facet data as a 3D surface.
        
        Parameters
        ----------
        var : ndarray, shape (n_faces,)
            Facet variable to plot (one value per facet)
        building_ids : array-like, optional
            Building IDs to plot. If None, plots all buildings.
        show : bool, default=True
            Display the plot immediately. If False, return the figure without showing.
        
        Returns
        -------
        fig : plotly.graph_objects.Figure or None
            Plotly figure object (if in notebook), None otherwise.
            Can be used to further customize the plot.
            
        Examples
        --------
        >>> # Plot net shortwave radiation for all buildings
        >>> fig = sim.plot_fac(K)
        >>> 
        >>> # Customize the returned figure
        >>> fig = sim.plot_fac(K)
        >>> fig.update_layout(title='Custom Title')
        >>> fig.show()
        >>> 
        >>> # Plot only for specific buildings
        >>> sim.plot_fac(K, building_ids=[1, 5, 10])
        """
        # Function only works when required data has been loaded
        if self.geom is None:
            raise ValueError("This method requires a geometry (STL) file.")
        
        # Validate input
        if len(var) != self.geom.n_faces:
            raise ValueError(f"Variable length ({len(var)}) must match number of facets ({self.geom.n_faces})")
        
        # Create colored mesh and render
        mesh = self._create_colored_mesh(var, building_ids)
        fig = self._render_scene(mesh, building_ids=building_ids, show=show)
        
        # Add building outlines
        self._add_building_outlines_to_scene(building_ids)
        
        return fig
    
    def _create_colored_mesh(self, var: np.ndarray, building_ids: Optional[np.ndarray] = None):
        """Create a colored trimesh object from facet data, optionally filtered by building IDs."""
        try:
            import trimesh
        except ImportError:
            raise ImportError("trimesh is required. Install with: pip install trimesh")
        
        import matplotlib.pyplot as plt
        import matplotlib.cm as cm
        
        # Get geometry data
        vertices = self.geom.stl.vertices
        faces = self.geom.stl.faces
        
        # Filter faces by building IDs if specified
        face_mask = None
        if building_ids is not None:
            # Get face to building mapping
            face_to_building = self.geom.get_face_to_building_map()
            building_ids = np.asarray(building_ids)
            
            # Create face mask for requested building IDs
            face_mask = np.isin(face_to_building, building_ids)
            
            # Apply face mask to select only specified buildings
            if np.any(face_mask):
                selected_faces = faces[face_mask]
                selected_var = var[face_mask]
            else:
                # No valid faces found - show warning and use all
                print('Warning: No valid faces found for the specified building IDs')
                selected_faces = faces
                selected_var = var
                face_mask = None
        else:
            selected_faces = faces
            selected_var = var
        
        # Create trimesh object with selected faces
        # Keep all original vertices but only include selected faces
        mesh = trimesh.Trimesh(vertices=vertices, faces=selected_faces, process=False)
        
        # Map variable values to face colors using colormap
        valid_mask = ~np.isnan(selected_var)
        if np.any(valid_mask):
            vmin = np.nanmin(selected_var[valid_mask])
            vmax = np.nanmax(selected_var[valid_mask])
            norm = plt.Normalize(vmin=vmin, vmax=vmax)
            cmap = cm.get_cmap('viridis')
            
            # Create face colors (RGBA)
            face_colors = np.ones((len(selected_faces), 4))  # Default white
            face_colors[valid_mask] = cmap(norm(selected_var[valid_mask]))
            mesh.visual.face_colors = face_colors
        
        return mesh
    
    def _render_scene(self, mesh, show_outlines: bool = True, angle_threshold: float = 45.0, building_ids: Optional[np.ndarray] = None, custom_edges: Optional[List[tuple]] = None, show: bool = True):
        """Render the mesh scene using trimesh/plotly. Returns figure handle if available."""
        try:
            import trimesh
        except ImportError:
            raise ImportError("trimesh is required. Install with: pip install trimesh")
        
        # Normalize mesh input to a list
        meshes = mesh if isinstance(mesh, (list, tuple)) else [mesh]
        
        # Create scene
        scene = trimesh.Scene()
        for m in meshes:
            scene.add_geometry(m)
        
        # Add outline edges using udgeom's outline calculation
        outline_edges = []
        if show_outlines:
            if custom_edges is not None:
                outline_edges = custom_edges
            else:
                outline_edges = self.geom._calculate_outline_edges(angle_threshold=angle_threshold)

            # Filter outline edges by building IDs if specified
            if building_ids is not None and len(outline_edges) > 0:
                face_to_building = self.geom.get_face_to_building_map()
                building_ids = np.asarray(building_ids)

                # Filter edges: keep only edges where both vertices belong to faces in selected buildings
                filtered_edges = []
                vertices = self.geom.stl.vertices
                faces = self.geom.stl.faces

                for edge in outline_edges:
                    # Find faces that use these vertices
                    v0, v1 = edge
                    # Check if any selected building's faces use this edge
                    faces_with_edge = np.where(
                        ((faces[:, 0] == v0) | (faces[:, 1] == v0) | (faces[:, 2] == v0)) &
                        ((faces[:, 0] == v1) | (faces[:, 1] == v1) | (faces[:, 2] == v1))
                    )[0]

                    # Keep edge if any of these faces belong to selected buildings
                    if len(faces_with_edge) > 0:
                        if np.any(np.isin(face_to_building[faces_with_edge], building_ids)):
                            filtered_edges.append(edge)

                outline_edges = filtered_edges

            if len(outline_edges) > 0:
                vertices = self.geom.stl.vertices
                entities = [trimesh.path.entities.Line([edge[0], edge[1]]) for edge in outline_edges]
                entity_colors = np.tile([0, 0, 0, 255], (len(entities), 1))
                path = trimesh.path.Path3D(entities=entities, vertices=vertices, colors=entity_colors)
                scene.add_geometry(path)
        
        # Check if running in Jupyter notebook
        try:
            from IPython.display import display
            in_notebook = True
        except ImportError:
            in_notebook = False
        
        if in_notebook:
            return self._render_plotly(meshes, outline_edges, show=show)
        else:
            if show:
                self._render_trimesh(scene, len(outline_edges))
            return None
    
    def _render_plotly(self, meshes, outline_edges, show: bool = True):
        """Render using plotly for notebook display. Returns the figure object."""
        try:
            import plotly.graph_objects as go
            import plotly.io as pio
            
            pio.renderers.default = 'notebook'
            traces = []
            for m in meshes:
                vertices = m.vertices
                faces = m.faces
                colors = m.visual.face_colors
                
                # Derive opacity from alpha channel if present; default to 1.0
                opacity = 1.0
                if colors.shape[1] == 4:
                    # Use mean alpha normalized to [0,1]
                    opacity = np.clip(np.mean(colors[:, 3]) / 255.0, 0.0, 1.0)
                
                traces.append(
                    go.Mesh3d(
                        x=vertices[:, 0], y=vertices[:, 1], z=vertices[:, 2],
                        i=faces[:, 0], j=faces[:, 1], k=faces[:, 2],
                        facecolor=[f'rgb({c[0]},{c[1]},{c[2]})' for c in colors[:, :3]],
                        opacity=opacity,
                        flatshading=True
                    )
                )
            
            fig = go.Figure(data=traces)
            
            # Add outline edges
            if len(outline_edges) > 0:
                edge_x, edge_y, edge_z = [], [], []
                z_offset = 0.1
                # Use vertices from the first mesh for outlines
                base_vertices = meshes[0].vertices
                for edge in outline_edges:
                    p0 = base_vertices[edge[0]]
                    p1 = base_vertices[edge[1]]
                    z0 = p0[2] + z_offset if abs(p0[2]) < 0.01 else p0[2]
                    z1 = p1[2] + z_offset if abs(p1[2]) < 0.01 else p1[2]
                    edge_x.extend([p0[0], p1[0], None])
                    edge_y.extend([p0[1], p1[1], None])
                    edge_z.extend([z0, z1, None])
                
                fig.add_trace(go.Scatter3d(
                    x=edge_x, y=edge_y, z=edge_z,
                    mode='lines', line=dict(color='black', width=2),
                    showlegend=False, hoverinfo='skip'
                ))
            
            fig.update_layout(
                scene=dict(
                    aspectmode='data',
                    xaxis_title='x (m)', yaxis_title='y (m)', zaxis_title='z (m)',
                    xaxis=dict(showgrid=False, showbackground=False),
                    yaxis=dict(showgrid=False, showbackground=False),
                    zaxis=dict(showgrid=False, showbackground=False),
                    camera=dict(
                        projection=dict(type='orthographic'),
                        eye=dict(x=-1.25, y=-1.25, z=1.25)
                    )
                ),
                showlegend=False
            )
            
            if show:
                fig.show()
            return fig
            
        except ImportError:
            print("Plotly not available. Install with: pip install plotly")
            return None
    
    def _render_trimesh(self, scene, num_outline_edges, show: bool = True):
        """Render using trimesh viewer."""
        if not show:
            return
        try:
            scene.show()
        except Exception as e:
            print(f"Could not open trimesh viewer: {e}")
            print("Install pyglet or pyrender: pip install pyglet")
    
    def _add_building_outlines_to_scene(self, building_ids=None):
        """Add building outlines to the current scene (placeholder for future implementation)."""
        # This method matches MATLAB's add_building_outlines() call structure
        # Currently handled within _render_scene via udgeom's _calculate_outline_edges
        pass
    
    def _add_building_outlines(self, ax, building_ids=None, angle_threshold: float = 45.0):
        """
        Helper method to add building outline edges to current 3D matplotlib plot.
        
        Used by matplotlib-based plotting methods like plot_fac_type.
        
        Parameters
        ----------
        ax : matplotlib Axes3D
            The 3D axes to add outlines to
        building_ids : list of int, optional
            Building IDs to outline. If None or empty, outline all buildings.
        angle_threshold : float, default=45.0
            Angle threshold in degrees for edge detection
        """
        if self.geom is None or not hasattr(self.geom, 'stl') or self.geom.stl is None:
            return
        
        # Use the geometry's built-in outline edge calculation method
        outline_edges = self.geom._calculate_outline_edges(angle_threshold=angle_threshold)
        
        if len(outline_edges) == 0:
            return
        
        # Get vertices
        vertices = self.geom.stl.vertices
        
        # Plot outline edges (including ground facet edges)
        z_offset = 0.1  # Offset to prevent z-fighting with ground plane
        for edge in outline_edges:
            p0 = vertices[edge[0]]
            p1 = vertices[edge[1]]
            
            # Add z-offset for edges at ground level to make them visible
            z0 = p0[2] + z_offset if abs(p0[2]) < 0.01 else p0[2]
            z1 = p1[2] + z_offset if abs(p1[2]) < 0.01 else p1[2]
            
            ax.plot([p0[0], p1[0]], 
                   [p0[1], p1[1]], 
                   [z0, z1], 
                   'k-', linewidth=2, alpha=1.0, zorder=10)
    
    def plot_fac_type(self, building_ids: Optional[np.ndarray] = None, 
                      show_outlines: bool = True, angle_threshold: float = 45.0, show: bool = True):
        """
        Plot the different surface types in the geometry using trimesh/Plotly.
        
        Parameters
        ----------
        building_ids : array-like, optional
            Building IDs to plot. If None, plots all buildings.
        show_outlines : bool, default=True
            Whether to show building outline edges
        angle_threshold : float, default=45.0
            Angle threshold in degrees for outline edge detection
        show : bool, default=True
            Display the plot immediately. If False, return the figure without showing.
            
        Returns
        -------
        fig : plotly.graph_objects.Figure or None
            Plotly figure object (if in notebook), None otherwise.
            Can be used to further customize the plot.
            
        Raises
        ------
        ValueError
            If required data (geometry, facets, factypes) is not loaded
        ImportError
            If trimesh is not installed
            
        Examples
        --------
        >>> sim = UDBase(101, 'experiments/101')
        >>> fig = sim.plot_fac_type()
        >>> 
        >>> # Customize the returned figure
        >>> fig = sim.plot_fac_type()
        >>> fig.update_layout(title='My Custom Title')
        >>> fig.show()
        """
        # Check required data
        if self.geom is None:
            raise ValueError("This method requires a geometry (STL) file. "
                           "Ensure stl_file is specified in namoptions.")
        
        if not hasattr(self, 'facs') or self.facs is None:
            raise ValueError("This method requires facet data. "
                           f"Ensure {self.ffacets}.{self.expnr} exists.")
        
        if not hasattr(self, 'factypes') or self.factypes is None:
            raise ValueError("This method requires facet type data. "
                           f"Ensure {self.ffactypes}.{self.expnr} exists.")
        
        try:
            import trimesh
        except ImportError:
            raise ImportError("trimesh is required for this visualization. "
                            "Install with: pip install trimesh")
        
        # Get data
        facids = self.facs['typeid']
        typeids = self.factypes['id']
        names = self.factypes['name']
        unique_ids = np.unique(facids)
        
        # Get default matplotlib color cycle for consistency
        import matplotlib.pyplot as plt
        prop_cycle = plt.rcParams['axes.prop_cycle']
        default_colors = prop_cycle.by_key()['color']
        
        # Get geometry data
        vertices = self.geom.stl.vertices
        faces = self.geom.stl.faces
        
        # Create face colors array (RGBA)
        face_colors = np.ones((len(faces), 4)) * 255  # Default white
        
        # Map each type to a color
        type_labels = []
        type_colors = []
        for idx, type_id in enumerate(unique_ids):
            type_mask = (facids == type_id)
            
            # Get matplotlib color (hex) and convert to RGB
            hex_color = default_colors[idx % len(default_colors)]
            # Convert hex to RGB (0-255)
            hex_color = hex_color.lstrip('#')
            rgb = tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4))
            
            # Set face colors for this type
            face_colors[type_mask, :3] = rgb
            face_colors[type_mask, 3] = 230  # Alpha
            
            # Get label
            name_idx = np.where(typeids == type_id)[0]
            if len(name_idx) > 0:
                label = names[name_idx[0]]
            else:
                label = f"Type {type_id}"
            
            type_labels.append(label)
            type_colors.append(rgb)
        
        # Create trimesh object
        mesh = trimesh.Trimesh(vertices=vertices, faces=faces, process=False)
        mesh.visual.face_colors = face_colors
        
        # Check if running in Jupyter notebook
        try:
            from IPython.display import display
            in_notebook = True
        except ImportError:
            in_notebook = False
        
        if in_notebook:
            print(f"Rendering {len(mesh.faces)} faces for notebook display...")
            
            try:
                import plotly.graph_objects as go
                import plotly.io as pio
                
                # Configure plotly for notebook display
                pio.renderers.default = 'notebook'
                
                # Create figure with separate mesh traces for each type (for legend)
                fig = go.Figure()
                
                for idx, type_id in enumerate(unique_ids):
                    type_mask = (facids == type_id)
                    type_face_indices = np.where(type_mask)[0]
                    
                    if len(type_face_indices) == 0:
                        continue
                    
                    # Get color
                    rgb = type_colors[idx]
                    color_str = f'rgb({rgb[0]},{rgb[1]},{rgb[2]})'
                    
                    # Get faces for this type
                    type_faces = faces[type_face_indices]
                    
                    # Add mesh trace for this type
                    fig.add_trace(go.Mesh3d(
                        x=vertices[:, 0],
                        y=vertices[:, 1],
                        z=vertices[:, 2],
                        i=type_faces[:, 0],
                        j=type_faces[:, 1],
                        k=type_faces[:, 2],
                        color=color_str,
                        opacity=1.0,
                        flatshading=True,
                        name=type_labels[idx],
                        showlegend=True
                    ))
                
                # Add outline edges if requested
                if show_outlines:
                    outline_edges = self.geom._calculate_outline_edges(angle_threshold=angle_threshold)
                    if len(outline_edges) > 0:
                        print(f"Added {len(outline_edges)} outline edges")
                        edge_x = []
                        edge_y = []
                        edge_z = []
                        z_offset = 0.1  # Offset to prevent z-fighting with ground plane
                        for edge in outline_edges:
                            p0 = vertices[edge[0]]
                            p1 = vertices[edge[1]]
                            # Add z-offset for edges at ground level
                            z0 = p0[2] + z_offset if abs(p0[2]) < 0.01 else p0[2]
                            z1 = p1[2] + z_offset if abs(p1[2]) < 0.01 else p1[2]
                            edge_x.extend([p0[0], p1[0], None])
                            edge_y.extend([p0[1], p1[1], None])
                            edge_z.extend([z0, z1, None])
                        
                        fig.add_trace(go.Scatter3d(
                            x=edge_x, y=edge_y, z=edge_z,
                            mode='lines',
                            line=dict(color='black', width=2),
                            name='Outlines',
                            showlegend=False,
                            hoverinfo='skip'
                        ))
                
                fig.update_layout(
                    scene=dict(
                        aspectmode='data',
                        xaxis_title='x (m)',
                        yaxis_title='y (m)',
                        zaxis_title='z (m)',
                        xaxis=dict(showgrid=False, showbackground=False),
                        yaxis=dict(showgrid=False, showbackground=False),
                        zaxis=dict(showgrid=False, showbackground=False),
                        camera=dict(
                            projection=dict(type='orthographic'),
                            eye=dict(x=-1.25, y=-1.25, z=1.25)
                        )
                    ),
                    title='Surface Types',
                    showlegend=True
                )
                
                if show:
                    fig.show()
                return fig
                
            except ImportError:
                print("Plotly not available. Falling back to static rendering.")
                print("Install plotly for interactive 3D: pip install plotly")
                # Fall back to showing in external window
                scene = trimesh.Scene(mesh)
                if show_outlines:
                    outline_edges = self.geom._calculate_outline_edges(angle_threshold=angle_threshold)
                    if len(outline_edges) > 0:
                        entities = [trimesh.path.entities.Line([edge[0], edge[1]]) for edge in outline_edges]
                        entity_colors = np.tile([0, 0, 0, 255], (len(entities), 1))
                        path = trimesh.path.Path3D(entities=entities, vertices=vertices, colors=entity_colors)
                        scene.add_geometry(path)
                if show:
                    try:
                        scene.show()
                    except:
                        print("Could not display. Try installing: pip install plotly or pyglet")
        else:
            # Show in external window
            print(f"Opening trimesh viewer with {len(mesh.faces)} faces...")
            scene = trimesh.Scene(mesh)
            if show_outlines:
                outline_edges = self.geom._calculate_outline_edges(angle_threshold=angle_threshold)
                if len(outline_edges) > 0:
                    print(f"Added {len(outline_edges)} outline edges")
                    entities = [trimesh.path.entities.Line([edge[0], edge[1]]) for edge in outline_edges]
                    entity_colors = np.tile([0, 0, 0, 255], (len(entities), 1))
                    path = trimesh.path.Path3D(entities=entities, vertices=vertices, colors=entity_colors)
                    scene.add_geometry(path)
            if show:
                try:
                    scene.show()
                except Exception as e:
                    print(f"Could not open trimesh viewer: {e}")
                    print("You may need to install pyglet or pyrender: pip install pyglet")
    
    def convert_fac_to_field(self, var: np.ndarray, facsec: Optional[Dict] = None,
                            dz: Optional[np.ndarray] = None) -> np.ndarray:
        """
        Convert a facet variable to a 3D field density.
        
        This method converts facet-based data to a 3D volumetric field by
        distributing facet areas into grid cells. Useful for computing volume
        integrals of surface quantities.
        
        Parameters
        ----------
        var : ndarray, shape (n_facets,)
            Facet variable to convert to field
        facsec : dict, optional
            Facet sections dictionary. If None, uses self.facsec['c']
        dz : ndarray, optional
            Vertical grid spacing. If None, uses self.dzt
            
        Returns
        -------
        fld : ndarray, shape (itot, jtot, ktot)
            3D field density [var_units / m]
            
        Raises
        ------
        ValueError
            If facet section data is not available
            
        Notes
        -----
        The conversion computes a density field where each cell contains:
        fld[i,j,k] = sum(var[facid] * area / (dx * dy * dz))
        
        This is used internally by calculate_frontal_properties() to compute
        frontal area densities from projected facet areas.
        
        Examples
        --------
        Convert surface temperature to a field:
        >>> Ts = sim.load_fac_temperature('Ts')
        >>> T_field = sim.convert_fac_to_field(Ts[:, 0])
        
        Use custom facet sections (e.g., for u-grid):
        >>> var = np.random.randn(sim.geom.n_faces)
        >>> fld = sim.convert_fac_to_field(var, facsec=sim.facsec['u'])
        """
        # Check that facet sections are available
        if not hasattr(self, 'facsec') or self.facsec is None:
            raise ValueError(
                "This method requires facet section data. "
                f"Ensure {self.ffacet_sections}_(u,v,w,c).{self.expnr} and "
                f"{self.ffluid_boundary}_(u,v,w,c).{self.expnr} files exist."
            )
        
        # Use defaults if not provided
        if facsec is None:
            facsec = self.facsec['c']
        
        if dz is None:
            dz = self.dzt
        
        # Initialize field
        fld = np.zeros((self.itot, self.jtot, self.ktot), dtype=np.float32)
        
        # Get facet section data
        facids = facsec['facid']
        areas = facsec['area']
        locs = facsec['locs']  # (i, j, k) locations (0-based)
        
        i_idx = locs[:, 0].astype(int)
        j_idx = locs[:, 1].astype(int)
        k_idx = locs[:, 2].astype(int)
        
        # Loop over all facet sections and create density field
        for m in range(len(areas)):
            facid = facids[m]
            i, j, k = i_idx[m], j_idx[m], k_idx[m]
            
            # Add contribution to cell
            cell_volume = self.dx * self.dy * dz[k]
            fld[i, j, k] += var[facid] * areas[m] / cell_volume
        
        return fld
    
    def convert_facvar_to_field(self, var: np.ndarray, facsec: Dict,
                                building_ids: Optional[np.ndarray] = None) -> np.ndarray:
        """
        Transfer a facet variable onto the grid.
        
        Matches MATLAB implementation: convert_facvar_to_field(obj, var, facsec, building_ids)
        
        Parameters
        ----------
        var : ndarray
            Facet variable (e.g., from load_fac_eb, load_fac_temperature)
        facsec : dict
            Facet section structure (e.g., obj.facsec['u'])
        building_ids : array-like, optional
            Array of building IDs to include. If None, all buildings included.
            
        Returns
        -------
        ndarray
            Variable on the grid (itot x jtot x ktot)
            
        Examples
        --------
        >>> # Convert all facets
        >>> fld = sim.convert_facvar_to_field(var, sim.facsec['c'])
        
        >>> # Convert only specific buildings
        >>> fld = sim.convert_facvar_to_field(var, sim.facsec['c'], [1, 5, 10])
        """
        # Normalize by converting ones, then divide actual values by normalization
        norm = self.convert_facflx_to_field(np.ones_like(var), facsec, self.dzt, building_ids)
        norm[norm == 0] = 1  # Avoid NaNs
        fld = self.convert_facflx_to_field(var, facsec, self.dzt, building_ids) / norm
        return fld
    
    def convert_facflx_to_field(self, var: np.ndarray, facsec: Dict, dz: np.ndarray,
                                building_ids: Optional[np.ndarray] = None) -> np.ndarray:
        """
        Convert a facet flux variable to a density in a 3D field.
        
        Matches MATLAB implementation: convert_facflx_to_field(obj, var, facsec, dz, building_ids)
        
        Parameters
        ----------
        var : ndarray
            Facet flux variable (e.g., from load_fac_eb)
        facsec : dict
            Facet section structure (e.g., obj.facsec['u'])
        dz : ndarray
            Vertical grid spacing at cell centers (obj.dzt)
        building_ids : array-like, optional
            Array of building IDs to include. If None, all buildings included.
            
        Returns
        -------
        ndarray
            3D field density (itot x jtot x ktot)
            
        Raises
        ------
        ValueError
            If facet section data is not loaded
            
        Examples
        --------
        >>> # Convert all facets
        >>> fld = sim.convert_facflx_to_field(var, sim.facsec['c'], sim.dzt)
        
        >>> # Convert only specific buildings
        >>> fld = sim.convert_facflx_to_field(var, sim.facsec['c'], sim.dzt, [1, 5, 10])
        """
        if not hasattr(self, 'facsec') or self.facsec is None:
            raise ValueError("This method requires facet section data. "
                           "Ensure facet_sections_*.txt and fluid_boundary_*.txt files exist.")
        
        # Initialize output field
        fld = np.zeros((self.itot, self.jtot, self.ktot), dtype=np.float32)
        
        # Get facet section data
        facids = facsec['facid']
        areas = facsec['area']
        locs = facsec['locs']  # (i, j, k) locations (0-based)
        
        # If building IDs specified, get face mask for filtering
        face_mask = None
        if building_ids is not None:
            # Get face to building mapping
            face_to_building = self.geom.get_face_to_building_map()
            building_ids = np.asarray(building_ids)
            
            # Create face mask for requested building IDs
            face_mask = np.isin(face_to_building, building_ids)
        
        # Loop over all facet sections and create density field
        for m in range(len(facids)):
            facid = int(facids[m])
            
            # If building filtering is active, check if this face should be included
            if face_mask is not None and not face_mask[facid]:
                continue
            
            # Get grid location (convert from 1-indexed to 0-indexed)
            i, j, k = int(locs[m, 0]), int(locs[m, 1]), int(locs[m, 2])
            
            # Add contribution to cell
            cell_volume = self.dx * self.dy * dz[k]
            fld[i, j, k] += var[facid] * areas[m] / cell_volume
        
        return fld
    
    def load_prof(self) -> np.ndarray:
        """
        Load information from prof.inp file.
        
        Matches MATLAB implementation: load_prof(obj)
        
        Returns
        -------
        ndarray
            Profile data from prof.inp file
            
        Raises
        ------
        FileNotFoundError
            If prof.inp file does not exist
            
        Examples
        --------
        >>> prof_data = sim.load_prof()
        """
        fname = f"{self.fprof}.{self.expnr}"
        fpath = self.path / fname
        
        if not fpath.exists():
            raise FileNotFoundError(f"Profile file not found: {fpath}")
        
        # Read data starting from line 3 (skip 2 header lines)
        data = np.loadtxt(fpath, skiprows=2)
        return data
    
    def load_facsec(self, var: str) -> Dict[str, np.ndarray]:
        """
        Load facet section information for a specific variable.
        
        Matches MATLAB implementation: load_facsec(obj, strvar)
        
        Parameters
        ----------
        var : str
            Variable name ('u', 'v', 'w', or 'c')
            
        Returns
        -------
        dict
            Dictionary with keys:
            - 'facid': Facet IDs (0-based)
            - 'area': Facet section areas
            - 'locs': Fluid boundary point locations (i, j, k), 0-based
            - 'distance': Distances from facet to fluid point
            
        Raises
        ------
        FileNotFoundError
            If facet section or fluid boundary files don't exist
            
        Examples
        --------
        >>> facsec_u = sim.load_facsec('u')
        >>> print(facsec_u.keys())
        """
        # Load facet section data
        fname_sec = self.path / f"{self.ffacet_sections}_{var}.txt"
        if not fname_sec.exists():
            raise FileNotFoundError(f"Facet section file not found: {fname_sec}")
        
        facsecs = np.loadtxt(fname_sec, skiprows=1)
        
        # Load fluid boundary points
        fname_fluid = self.path / f"{self.ffluid_boundary}_{var}.txt"
        if not fname_fluid.exists():
            raise FileNotFoundError(f"Fluid boundary file not found: {fname_fluid}")
        
        fluid_boundary = np.loadtxt(fname_fluid, skiprows=1)
        
        # Create structure matching MATLAB
        data = {
            'facid': facsecs[:, 0].astype(int) - 1,
            'area': facsecs[:, 1],
            'locs': fluid_boundary[facsecs[:, 2].astype(int) - 1, :].astype(int) - 1,
            'distance': facsecs[:, 3]
        }
        
        return data
    
    def calculate_frontal_properties(self) -> Dict[str, Any]:
        """
        Calculate skyline, frontal areas, and blockage ratios.
        
        Computes geometric properties of the urban canopy including:
        - Skyline profiles in x and y directions
        - Total frontal areas perpendicular to x and y
        - Blockage ratios (fraction of frontal area blocked)
        
        Returns
        -------
        res : dict
            Dictionary containing:
            - 'skylinex' : ndarray, shape (jtot, ktot)
                Binary indicator of blocked cells in x-direction (y-z plane)
            - 'skyliney' : ndarray, shape (itot, ktot)
                Binary indicator of blocked cells in y-direction (x-z plane)
            - 'Afx' : float
                Total frontal area perpendicular to x-direction [m²]
            - 'Afy' : float
                Total frontal area perpendicular to y-direction [m²]
            - 'brx' : float
                Blockage ratio in x-direction (dimensionless, 0-1)
            - 'bry' : float
                Blockage ratio in y-direction (dimensionless, 0-1)
                
        Raises
        ------
        ValueError
            If geometry or facet section data is not available
            
        Notes
        -----
        The frontal area is the projected area of all surfaces onto a plane
        perpendicular to the flow direction. The blockage ratio is the fraction
        of the domain cross-section that is blocked by buildings.
        
        These quantities are important for understanding drag and flow resistance
        in urban canopies.
        
        Examples
        --------
        >>> props = sim.calculate_frontal_properties()
        >>> print(f"Frontal area (x): {props['Afx']:.1f} m²")
        >>> print(f"Blockage ratio (x): {props['brx']:.3f}")
        >>> 
        >>> # Visualize skyline
        >>> import matplotlib.pyplot as plt
        >>> plt.imshow(props['skylinex'].T, origin='lower', aspect='auto')
        >>> plt.xlabel('j (grid points)')
        >>> plt.ylabel('k (grid points)')
        >>> plt.title('Skyline in x-direction')
        >>> plt.colorbar(label='Blocked (1) / Open (0)')
        >>> plt.show()
        """
        # Check required data
        if self.geom is None:
            raise ValueError("This method requires a geometry (STL) file. "
                           "Ensure stl_file is specified in namoptions.")
        
        if not hasattr(self, 'facsec') or self.facsec is None:
            raise ValueError(
                "This method requires facet section data. "
                f"Ensure {self.ffacet_sections}_(u,v,w,c).{self.expnr} and "
                f"{self.ffluid_boundary}_(u,v,w,c).{self.expnr} files exist."
            )
        
        # Get face normals
        norms = self.geom.face_normals
        
        # Create surface quantity for projected area in x and y directions
        # Project onto planes perpendicular to x and y axes
        # Only count outward-facing (negative dot product with axis)
        phix = -np.minimum(np.dot(norms, np.array([1, 0, 0])), 0)
        phiy = -np.minimum(np.dot(norms, np.array([0, 1, 0])), 0)
        
        # Convert to density fields using convert_facflx_to_field (matches MATLAB)
        rhoLx = self.convert_facflx_to_field(phix, self.facsec['c'], self.dzt)
        rhoLy = self.convert_facflx_to_field(phiy, self.facsec['c'], self.dzt)
        
        # Calculate indicator functions for blockage
        # Ibx[j,k] = 1 if any cell along x-direction at (j,k) is blocked
        Ibx = (np.sum(rhoLx, axis=0) > 0).astype(float)
        # Iby[i,k] = 1 if any cell along y-direction at (i,k) is blocked
        Iby = (np.sum(rhoLy, axis=1) > 0).astype(float)
        
        # Integrate to get frontal areas and blockage ratios
        Afx = 0.0
        Afy = 0.0
        brx = 0.0
        bry = 0.0
        
        for k in range(self.ktot):
            # Frontal areas (integrate density over volume)
            Afx += np.sum(rhoLx[:, :, k]) * self.dx * self.dy * self.dzt[k]
            Afy += np.sum(rhoLy[:, :, k]) * self.dx * self.dy * self.dzt[k]
            
            # Blockage ratios (integrate indicator over cross-section)
            brx += np.sum(Ibx[:, k]) * self.dy * self.dzt[k]
            bry += np.sum(Iby[:, k]) * self.dx * self.dzt[k]
        
        # Normalize blockage ratios by total cross-sectional area
        brx /= (self.ylen * self.zsize)
        bry /= (self.xlen * self.zsize)
        
        # Print results
        print(f"x-direction: frontal area = {Afx:8.1f} m², blockage ratio = {brx:8.3f}")
        print(f"y-direction: frontal area = {Afy:8.1f} m², blockage ratio = {bry:8.3f}")
        
        return {
            'skylinex': Ibx,
            'skyliney': Iby,
            'Afx': Afx,
            'Afy': Afy,
            'brx': brx,
            'bry': bry
        }
    
    def plot_building_ids(self, show: bool = True):
        """
        Plot building IDs from above (x,y view) with distinct colors.
        
        Matches MATLAB implementation: plot_building_ids(obj)
        
        Creates a top-view plot showing buildings in different colors with
        building IDs displayed at the center of gravity of each building.
        Buildings are numbered from left-bottom to right-top based on their
        centroid positions.

        Parameters
        ----------
        show : bool, default=True
            Display the plot immediately. If False, return the figure without showing.
            
        Raises
        ------
        ValueError
            If geometry is not loaded or has no buildings
        ImportError
            If matplotlib is not installed
        show : bool, default=True
            Display the plot immediately. If False, return the figure without showing.
            
        Examples
        --------
        Plot building IDs with default settings:
        >>> sim.plot_building_ids()
        
        See Also
        --------
        plot_2dmap : Plot 2D map with custom values and labels
        """
        if self.geom is None:
            raise ValueError("This method requires a geometry (STL) file. "
                           "Ensure stl_file is specified in namoptions.")
        
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            raise ImportError("matplotlib is required for visualization. "
                            "Install with: pip install matplotlib")
        
        # Get building outlines
        outlines = self.geom.calculate_outline2d()
        if not outlines:
            raise ValueError("No buildings found in geometry")
        
        num_buildings = len(outlines)
        
        # Create values (building indices) and labels
        values = np.arange(1, num_buildings + 1)
        labels = [str(i) for i in range(1, num_buildings + 1)]
        
        # Randomize colors to avoid adjacent buildings having similar colors
        color_order = np.random.permutation(num_buildings)
        color_values = color_order.copy()
        
        # Use plot_2dmap to create the visualization
        fig, ax = self.plot_2dmap(color_values, labels, show=show)
        pc = ax.collections[0]  # the PatchCollection from plot_2dmap
        fig.colorbar(pc, ax=ax, cmap='hsv')
        ax.set_title(f'Building Layout with IDs (Total: {num_buildings})')
        ax.set_xlabel('x [m]')
        ax.set_ylabel('y [m]')
        ax.set_aspect('equal')
        if show:
            plt.show()
        return fig, ax
    
    def plot_2dmap(self, val: Union[float, np.ndarray], 
                   labels: Optional[Union[str, list]] = None,
                   show: bool = True):
        """
        Plot a 2D map of buildings colored by a value per building.
        
        Matches MATLAB implementation: plot_2dmap(obj, val, labels)
        
        Creates a top-down view showing building outlines colored according
        to specified values, with optional text labels at building centroids.
        
        Parameters
        ----------
        val : float or ndarray
            Scalar value (applied to all buildings) or vector with one value
            per building. If vector, length must match number of buildings.
        labels : str, list of str, or None, optional
            Text labels to display at the centroid of each building.
            - If str: same label for all buildings
            - If list: must have length equal to number of buildings
            - If None: no labels displayed
        show : bool, default=True
            Display the plot immediately. If False, return the figure without showing.
            
        Raises
        ------
        ValueError
            If val length doesn't match number of buildings
            If labels length doesn't match number of buildings
            If geometry is not loaded
        ImportError
            If matplotlib is not installed
            
        Examples
        --------
        Plot all buildings with same color:
        >>> sim.plot_2dmap(1.0)
        
        Plot building heights with labels:
        >>> buildings = sim.geom.get_buildings()
        >>> heights = [max(b['triangulation'].Points[:, 2]) for b in buildings]
        >>> labels = [f"{i+1}" for i in range(len(buildings))]
        >>> sim.plot_2dmap(heights, labels)
        >>> plt.title('Maximum Building Height')
        
        See Also
        --------
        plot_building_ids : Plot buildings with ID labels
        """
        if self.geom is None or not hasattr(self.geom, 'stl') or self.geom.stl is None:
            raise ValueError("Geometry data not available. Cannot compute outlines.")
        
        try:
            import matplotlib.pyplot as plt
            from matplotlib.patches import Polygon as mplPolygon
            from matplotlib.collections import PatchCollection
        except ImportError:
            raise ImportError("matplotlib is required for visualization. "
                            "Install with: pip install matplotlib")
        
        # Get building outlines
        outlines = self.geom.calculate_outline2d()
        if not outlines:
            raise ValueError("No building outlines found in geometry")
        
        num_buildings = len(outlines)
        
        # Normalize values to per-building vector
        if np.isscalar(val):
            values = np.full(num_buildings, float(val))
        else:
            val_array = np.asarray(val)
            if len(val_array) != num_buildings:
                raise ValueError(f"Length of val ({len(val_array)}) must match "
                               f"number of buildings ({num_buildings})")
            values = val_array.astype(float)
        
        # Handle labels
        if labels is not None:
            if isinstance(labels, str):
                label_array = [labels] * num_buildings
            else:
                label_array = list(labels)
                if len(label_array) != num_buildings:
                    raise ValueError(f"Number of labels ({len(label_array)}) must match "
                                   f"number of buildings ({num_buildings})")
        else:
            label_array = None
        
        # Create figure and plot
        fig, ax = plt.subplots(figsize=(10, 8))
        
        patches = []
        colors = []
        
        for i, outline in enumerate(outlines):
            if np.isnan(values[i]):
                continue
            
            polygon = outline.get('polygon', None)
            centroid = outline.get('centroid', None)
            
            if polygon is None or len(polygon) == 0:
                continue
            
            # Create polygon patch (project to x-y plane)
            xy = polygon[:, :2]  # Take only x, y coordinates
            patch = mplPolygon(xy, closed=True)
            patches.append(patch)
            colors.append(values[i])
            
            # Add text label if provided
            if label_array is not None and centroid is not None:
                if not np.any(np.isnan(centroid[:2])):
                    ax.text(centroid[0], centroid[1], label_array[i],
                           ha='center', va='center',
                           fontsize=10, fontweight='bold',
                           color='black',
                           bbox=dict(boxstyle='round,pad=0.3',
                                   facecolor='white',
                                   edgecolor='none',
                                   alpha=0.7))
        
        # Create patch collection
        if patches:
            pc = PatchCollection(patches, cmap='viridis', edgecolor='black', linewidth=0.5)
            pc.set_array(np.array(colors))
            ax.add_collection(pc)
            plt.colorbar(pc, ax=ax)
        
        # Set axis properties
        ax.set_aspect('equal')
        ax.set_xlabel('x [m]')
        ax.set_ylabel('y [m]')
        ax.set_xlim(0, self.xlen)
        ax.set_ylim(0, self.ylen)
        ax.grid(True, alpha=0.3)
        
        if show:
            plt.show()
        return fig, ax
    
    def __str__(self):
        """User-friendly string representation."""
        return self.__repr__()


if __name__ == "__main__":
    # Basic test
    print("UDBase module loaded successfully")
    print("\nUsage:")
    print("  from udbase import UDBase")
    print("  sim = UDBase(expnr=65, path='experiments/065')")
