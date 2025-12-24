"""
uDALES Post-Processing Module

Python implementation of the MATLAB udbase class for analyzing
uDALES simulation outputs.

Copyright (C) 2024 the uDALES Team.
"""

import numpy as np
import xarray as xr
from pathlib import Path
from typing import Optional, Union, Dict, Any
import warnings


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
    
    def __init__(self, expnr: Union[int, str], path: Optional[Union[str, Path]] = None):
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
        
        # Read namoptions file
        self._read_namoptions()
        
        # Load grid
        self._load_grid()
        
        # Load geometry if present
        self._load_geometry()
        
        # Load solid masks if present
        self._load_solid_masks()
        
        # Load facet data if present
        self._load_facet_data()
    
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
                    # Read indices (1-based from Fortran)
                    indices = np.loadtxt(solid_file, skiprows=1, dtype=int)
                    
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
                    # Load facet section data
                    facsec_data = np.loadtxt(facsec_file, skiprows=1)
                    
                    # Load fluid boundary locations
                    fluid_boundary = np.loadtxt(fluid_boundary_file, skiprows=1, dtype=int)
                    
                    # Store in structure
                    self.facsec[grid_type] = {
                        'facid': facsec_data[:, 0].astype(int),
                        'area': facsec_data[:, 1],
                        'locs': fluid_boundary[facsec_data[:, 2].astype(int) - 1, :],  # Convert to 0-based
                        'distance': facsec_data[:, 3]
                    }
                    
                except Exception as e:
                    warnings.warn(f"Error loading facet_sections_{grid_type}.txt: {e}")
                    self._lffacet_sections = False
            else:
                self._lffacet_sections = False
    
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
    
    def load_field(self, var: Optional[str] = None) -> Union[xr.Dataset, xr.DataArray]:
        """
        Load 3D instantaneous field data.
        
        Parameters
        ----------
        var : str, optional
            Variable name to load. If None, displays available variables.
        
        Returns
        -------
        xarray.Dataset or xarray.DataArray
            Complete dataset if var is None, otherwise the requested variable.
        
        Examples
        --------
        >>> sim.load_field()  # Display available variables
        >>> u = sim.load_field('u')
        >>> print(u.shape)
        """
        filename = self.path / f"fielddump.{self.expnr}.nc"
        return self._load_ncdata(filename, var)
    
    def load_stat_xyt(self, var: Optional[str] = None) -> Union[xr.Dataset, xr.DataArray]:
        """
        Load slab-averaged (xy-plane) statistics.
        
        Parameters
        ----------
        var : str, optional
            Variable name to load. If None, displays available variables.
        
        Returns
        -------
        xarray.Dataset or xarray.DataArray
            Complete dataset if var is None, otherwise the requested variable.
        
        Examples
        --------
        >>> sim.load_stat_xyt()  # Display available variables
        >>> u_avg = sim.load_stat_xyt('u')
        """
        filename = self.path / f"xytdump.{self.expnr}.nc"
        return self._load_ncdata(filename, var)
    
    def load_stat_t(self, var: Optional[str] = None) -> Union[xr.Dataset, xr.DataArray]:
        """
        Load time-averaged 3D statistics.
        
        Parameters
        ----------
        var : str, optional
            Variable name to load. If None, displays available variables.
        
        Returns
        -------
        xarray.Dataset or xarray.DataArray
            Complete dataset if var is None, otherwise the requested variable.
        
        Examples
        --------
        >>> sim.load_stat_t()  # Display available variables
        >>> u_tavg = sim.load_stat_t('u')
        """
        filename = self.path / f"tdump.{self.expnr}.nc"
        print(filename.exists())
        return self._load_ncdata(filename, var)
    
    def load_slice(self, plane: str, var: Optional[str] = None) -> Union[xr.Dataset, xr.DataArray]:
        """
        Load 2D slice data.
        
        Parameters
        ----------
        plane : str
            Slice plane: 'i', 'j', or 'k' for x, y, or z slices
        var : str, optional
            Variable name to load. If None, displays available variables.
        
        Returns
        -------
        xarray.Dataset or xarray.DataArray
            Complete dataset if var is None, otherwise the requested variable.
        
        Examples
        --------
        >>> sim.load_slice('k')  # Display available variables in horizontal slice
        >>> u_slice = sim.load_slice('k', 'u')
        """
        if plane not in ['i', 'j', 'k']:
            raise ValueError("plane must be 'i', 'j', or 'k'")
        
        filename = self.path / f"{plane}slicedump.{self.expnr}.nc"
        return self._load_ncdata(filename, var)
    
    def _load_ncdata(self, filename: Path, var: Optional[str]) -> Union[xr.Dataset, xr.DataArray]:
        """
        Helper method to load NetCDF data using xarray.
        
        Automatically reverses dimension order to match MATLAB conventions.
        
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
            Variable to extract
        
        Returns
        -------
        xarray.Dataset or xarray.DataArray
            Data with dimensions reversed to match MATLAB conventions
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
            
            return ds_transposed[var]
    
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
    
    def load_fac_momentum(self, var: Optional[str] = None) -> Union[xr.Dataset, xr.DataArray]:
        """
        Load facet momentum data (pressure and shear stress).
        
        Parameters
        ----------
        var : str, optional
            Variable name to load. If None, displays available variables.
        
        Returns
        -------
        xarray.Dataset or xarray.DataArray
            Complete dataset if var is None, otherwise the requested variable.
        
        Examples
        --------
        >>> sim.load_fac_momentum()  # Display available variables
        >>> pressure = sim.load_fac_momentum('pf')
        """
        filename = self.path / f"fac.{self.expnr}.nc"
        return self._load_ncdata(filename, var)
    
    def load_fac_eb(self, var: Optional[str] = None) -> Union[xr.Dataset, xr.DataArray]:
        """
        Load facet surface energy balance data.
        
        Parameters
        ----------
        var : str, optional
            Variable name to load. If None, displays available variables.
        
        Returns
        -------
        xarray.Dataset or xarray.DataArray
            Complete dataset if var is None, otherwise the requested variable.
        
        Examples
        --------
        >>> sim.load_fac_eb()  # Display available variables
        >>> H = sim.load_fac_eb('hf')  # Sensible heat flux
        >>> K = sim.load_fac_eb('netsw')  # Net shortwave
        """
        filename = self.path / f"facEB.{self.expnr}.nc"
        return self._load_ncdata(filename, var)
    
    def load_fac_temperature(self, var: Optional[str] = None) -> Union[xr.Dataset, xr.DataArray]:
        """
        Load facet temperature data.
        
        Parameters
        ----------
        var : str, optional
            Variable name to load. If None, displays available variables.
        
        Returns
        -------
        xarray.Dataset or xarray.DataArray
            Complete dataset if var is None, otherwise the requested variable.
        
        Examples
        --------
        >>> sim.load_fac_temperature()  # Display available variables
        >>> T = sim.load_fac_temperature('T')  # Temperature in layers
        >>> dTdz = sim.load_fac_temperature('dTdz')  # Temperature gradient
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
        G = -lam[:, 0, np.newaxis] * dTdz[:, 0, :].values  # Ground heat flux: (n_facets, n_time)
        Tsurf = T[:, 0, :].values  # Surface temperature: (n_facets, n_time)
        
        # All data is already in MATLAB convention: (n_facets, n_time)
        # _load_ncdata transposes from NetCDF (time, fct) to (fct, time)
        return {
            'Kstar': K.values,
            'Lstar': L.values,
            'Lin': Lin.values,
            'Lout': Lout.values,
            'H': -H.values,  # Sign convention
            'E': -E.values,  # Sign convention
            'G': G,
            'Tsurf': Tsurf,
            't': t.values
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
    def time_average(var: np.ndarray, time: np.ndarray, 
                     tstart: Optional[float] = None, 
                     tstop: Optional[float] = None) -> np.ndarray:
        """
        Time-average a variable.
        
        Parameters
        ----------
        var : ndarray
            Variable to average. Time must be the last dimension.
        time : ndarray
            Time array
        tstart : float, optional
            Start time for averaging
        tstop : float, optional
            Stop time for averaging
        
        Returns
        -------
        ndarray
            Time-averaged values
        
        Examples
        --------
        >>> H_avg = UDBase.time_average(H, t, tstart=3600, tstop=7200)
        """
        if tstart is None:
            indstart = 0
        else:
            indstart = np.argmax(time >= tstart)
        
        if tstop is None:
            indstop = len(time)
        else:
            indstop = np.argmax(time >= tstop)
            if indstop == 0:
                indstop = len(time)
        
        return np.mean(var[..., indstart:indstop], axis=-1)

    @staticmethod
    def merge_stat(X: np.ndarray, n: int, Y: Optional[np.ndarray] = None,
                   XpXp: Optional[np.ndarray] = None,
                   XpYp: Optional[np.ndarray] = None):
        """
        Merge short-term statistics into longer windows.

        Parameters
        ----------
        X : ndarray
            First variable. Final dimension is time or short-term windows.
        n : int
            Number of consecutive samples/windows to merge.
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
        X = np.asarray(X)
        if X.shape[-1] < n:
            raise ValueError("Not enough samples to form a single merged window")

        nwin = X.shape[-1] // n
        start = X.shape[-1] - nwin * n
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

        for i, L in enumerate(Lflt_arr):
            wx = max(1, int(round(L / dx)))
            wy = max(1, int(round(L / dy)))

            kernel = np.zeros((nx, ny))
            kernel[:wx, :wy] = 1.0 / (wx * wy)
            kernel = np.roll(kernel, -wx // 2, axis=0)
            kernel = np.roll(kernel, -wy // 2, axis=1)
            k_hat = np.fft.fftn(kernel)

            for k in range(nz):
                v_hat = np.fft.fftn(var[:, :, k])
                out[:, :, k, i] = np.real(np.fft.ifftn(v_hat * k_hat))

        return out
    
    def plot_fac(self, var: np.ndarray, cmap: str = 'viridis', 
                 show_edges: bool = False, colorbar: bool = True,
                 title: Optional[str] = None, figsize: tuple = (10, 8),
                 vmin: Optional[float] = None, vmax: Optional[float] = None):
        """
        Plot a facet variable as a 3D colored surface.
        
        This method visualizes scalar data defined on each facet of the geometry
        by coloring the triangular faces according to the variable values.
        
        Parameters
        ----------
        var : ndarray, shape (n_faces,)
            Facet variable to plot (one value per facet)
        cmap : str, default='viridis'
            Matplotlib colormap name
        show_edges : bool, default=False
            If True, show triangle edges. Use False for cleaner visualization.
        colorbar : bool, default=True
            If True, display a colorbar
        title : str, optional
            Plot title. If None, no title is displayed.
        figsize : tuple, default=(10, 8)
            Figure size in inches (width, height)
        vmin, vmax : float, optional
            Min and max values for colormap. If None, uses data min/max.
            
        Raises
        ------
        ValueError
            If geometry is not loaded or if var has wrong dimensions
        ImportError
            If matplotlib is not installed
            
        Examples
        --------
        Plot net shortwave radiation on facets:
        >>> K = sim.load_fac_eb('K')
        >>> sim.plot_fac(K, title='Net Shortwave Radiation (W/m²)', cmap='hot')
        
        Plot momentum flux with custom range:
        >>> taux = sim.load_fac_momentum('taux')
        >>> sim.plot_fac(taux[:, 0], vmin=-1, vmax=1, title='Tau_x')
        """
        if self.geom is None:
            raise ValueError("This method requires a geometry (STL) file. "
                           "Ensure stl_file is specified in namoptions.")
        
        try:
            import matplotlib.pyplot as plt
            from mpl_toolkits.mplot3d.art3d import Poly3DCollection
        except ImportError:
            raise ImportError("matplotlib is required for visualization. "
                            "Install with: pip install matplotlib")
        
        # Validate input
        if len(var) != self.geom.n_faces:
            raise ValueError(f"Variable length ({len(var)}) must match number of facets ({self.geom.n_faces})")
        
        # Create figure
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111, projection='3d')
        
        # Get geometry data
        vertices = self.geom.stl.vertices
        faces = self.geom.stl.faces
        triangles = vertices[faces]
        
        # Create color-mapped collection
        edge_color = 'k' if show_edges else 'none'
        collection = Poly3DCollection(
            triangles,
            facecolors=plt.cm.get_cmap(cmap)(plt.Normalize(vmin=vmin, vmax=vmax)(var)),
            edgecolors=edge_color,
            linewidths=0.1 if show_edges else 0
        )
        ax.add_collection3d(collection)
        
        # Set axis limits
        bounds = self.geom.bounds
        ax.set_xlim(bounds[0, 0], bounds[1, 0])
        ax.set_ylim(bounds[0, 1], bounds[1, 1])
        ax.set_zlim(bounds[0, 2], bounds[1, 2])
        
        # Labels
        ax.set_xlabel('X (m)')
        ax.set_ylabel('Y (m)')
        ax.set_zlabel('Z (m)')
        
        if title:
            ax.set_title(title)
        
        # Equal aspect ratio
        max_range = np.array([
            bounds[1, 0] - bounds[0, 0],
            bounds[1, 1] - bounds[0, 1],
            bounds[1, 2] - bounds[0, 2]
        ]).max() / 2.0
        
        mid_x = (bounds[1, 0] + bounds[0, 0]) * 0.5
        mid_y = (bounds[1, 1] + bounds[0, 1]) * 0.5
        mid_z = (bounds[1, 2] + bounds[0, 2]) * 0.5
        
        ax.set_xlim(mid_x - max_range, mid_x + max_range)
        ax.set_ylim(mid_y - max_range, mid_y + max_range)
        ax.set_zlim(mid_z - max_range, mid_z + max_range)
        
        # Viewing angle
        ax.view_init(elev=30, azim=45)
        
        # Colorbar
        if colorbar:
            mappable = plt.cm.ScalarMappable(
                cmap=cmap,
                norm=plt.Normalize(vmin=vmin if vmin is not None else var.min(),
                                 vmax=vmax if vmax is not None else var.max())
            )
            mappable.set_array(var)
            plt.colorbar(mappable, ax=ax, shrink=0.5, aspect=5)
        
        plt.tight_layout()
        plt.show()
    
    def plot_fac_type(self, figsize: tuple = (12, 10), show_legend: bool = True):
        """
        Plot the different surface types in the geometry.
        
        This method visualizes the facet types defined in the simulation,
        coloring each surface type differently. Useful for verifying that
        wall properties are correctly assigned.
        
        Parameters
        ----------
        figsize : tuple, default=(12, 10)
            Figure size in inches (width, height)
        show_legend : bool, default=True
            If True, display a legend with surface type names
            
        Raises
        ------
        ValueError
            If required data (geometry, facets, factypes) is not loaded
        ImportError
            If matplotlib is not installed
            
        Examples
        --------
        >>> sim = UDBase(101, 'experiments/101')
        >>> sim.plot_fac_type()
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
            import matplotlib.pyplot as plt
            from mpl_toolkits.mplot3d.art3d import Poly3DCollection
        except ImportError:
            raise ImportError("matplotlib is required for visualization. "
                            "Install with: pip install matplotlib")
        
        # Get data
        facids = self.facs['typeid']
        typeids = self.factypes['id']
        names = self.factypes['name']
        unique_ids = np.unique(facids)
        
        # Get default matplotlib color cycle
        prop_cycle = plt.rcParams['axes.prop_cycle']
        default_colors = prop_cycle.by_key()['color']
        
        # Check if we have enough colors
        if len(unique_ids) > len(default_colors):
            warnings.warn(f"Too many surface types ({len(unique_ids)}) for unique colors. "
                        f"Only {len(default_colors)} available. Some colors will repeat.")
        
        # Create figure
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111, projection='3d')
        
        # Get geometry data
        vertices = self.geom.stl.vertices
        faces = self.geom.stl.faces
        
        # Plot each surface type
        labels = []
        for idx, type_id in enumerate(unique_ids):
            # Select facets of this type
            type_mask = (facids == type_id)
            type_faces = faces[type_mask]
            type_triangles = vertices[type_faces]
            
            # Get color (cycle through if needed)
            color = default_colors[idx % len(default_colors)]
            
            # Get name
            name_idx = np.where(typeids == type_id)[0]
            if len(name_idx) > 0:
                label = names[name_idx[0]]
            else:
                label = f"Type {type_id}"
            
            labels.append(label)
            
            # Create collection
            collection = Poly3DCollection(
                type_triangles,
                facecolors=color,
                edgecolors='none',
                alpha=0.9,
                label=label
            )
            ax.add_collection3d(collection)
        
        # Set axis limits
        bounds = self.geom.bounds
        ax.set_xlim(bounds[0, 0], bounds[1, 0])
        ax.set_ylim(bounds[0, 1], bounds[1, 1])
        ax.set_zlim(bounds[0, 2], bounds[1, 2])
        
        # Labels
        ax.set_xlabel('X (m)')
        ax.set_ylabel('Y (m)')
        ax.set_zlabel('Z (m)')
        ax.set_title('Surface Types')
        
        # Equal aspect ratio
        max_range = np.array([
            bounds[1, 0] - bounds[0, 0],
            bounds[1, 1] - bounds[0, 1],
            bounds[1, 2] - bounds[0, 2]
        ]).max() / 2.0
        
        mid_x = (bounds[1, 0] + bounds[0, 0]) * 0.5
        mid_y = (bounds[1, 1] + bounds[0, 1]) * 0.5
        mid_z = (bounds[1, 2] + bounds[0, 2]) * 0.5
        
        ax.set_xlim(mid_x - max_range, mid_x + max_range)
        ax.set_ylim(mid_y - max_range, mid_y + max_range)
        ax.set_zlim(mid_z - max_range, mid_z + max_range)
        
        # Viewing angle
        ax.view_init(elev=30, azim=45)
        
        # Legend
        if show_legend and len(labels) > 0:
            # Create custom legend handles
            from matplotlib.patches import Patch
            legend_elements = [
                Patch(facecolor=default_colors[i % len(default_colors)], 
                      edgecolor='none', label=labels[i])
                for i in range(len(labels))
            ]
            ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.05, 1))
        
        plt.tight_layout()
        plt.show()
    
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
        locs = facsec['locs']  # (i, j, k) locations (1-based from Fortran)
        
        # Convert locations to 0-based indexing
        i_idx = locs[:, 0] - 1
        j_idx = locs[:, 1] - 1
        k_idx = locs[:, 2] - 1
        
        # Loop over all facet sections and create density field
        for m in range(len(areas)):
            facid = facids[m] - 1  # Convert to 0-based
            i, j, k = i_idx[m], j_idx[m], k_idx[m]
            
            # Add contribution to cell
            cell_volume = self.dx * self.dy * dz[k]
            fld[i, j, k] += var[facid] * areas[m] / cell_volume
        
        return fld
    
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
        
        # Convert to density fields
        rhoLx = self.convert_fac_to_field(phix, self.facsec['c'], self.dzt)
        rhoLy = self.convert_fac_to_field(phiy, self.facsec['c'], self.dzt)
        
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
    
    def plot_building_ids(self, figsize: tuple = (10, 8), cmap: str = 'hsv'):
        """
        Plot building IDs from above (x,y view) with distinct colors.
        
        Creates a top-view plot showing buildings in different colors with
        building IDs displayed at the center of gravity of each building.
        Buildings are numbered from left-bottom to right-top based on their
        centroid positions.
        
        Parameters
        ----------
        figsize : tuple, default=(10, 8)
            Figure size in inches (width, height)
        cmap : str, default='hsv'
            Matplotlib colormap name for distinct colors
            
        Raises
        ------
        ValueError
            If geometry is not loaded or has no buildings
        ImportError
            If matplotlib is not installed
            
        Examples
        --------
        Plot building IDs with default settings:
        >>> sim.plot_building_ids()
        
        Use custom figure size and colormap:
        >>> sim.plot_building_ids(figsize=(12, 10), cmap='tab20')
        
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
        self.plot_2dmap(color_values, labels, figsize=figsize, cmap=cmap)
        
        plt.title(f'Building Layout with IDs (Total: {num_buildings})')
        plt.xlabel('x [m]')
        plt.ylabel('y [m]')
    
    def plot_2dmap(self, val: Union[float, np.ndarray], 
                   labels: Optional[Union[str, list]] = None,
                   figsize: tuple = (10, 8), cmap: str = 'viridis',
                   show_colorbar: bool = True):
        """
        Plot a 2D map of buildings colored by a value per building.
        
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
        figsize : tuple, default=(10, 8)
            Figure size in inches (width, height)
        cmap : str, default='viridis'
            Matplotlib colormap name
        show_colorbar : bool, default=True
            Whether to display a colorbar
            
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
        fig, ax = plt.subplots(figsize=figsize)
        
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
            pc = PatchCollection(patches, cmap=cmap, edgecolor='black', linewidth=0.5)
            pc.set_array(np.array(colors))
            ax.add_collection(pc)
            
            if show_colorbar:
                plt.colorbar(pc, ax=ax)
        
        # Set axis properties
        ax.set_aspect('equal')
        ax.set_xlabel('x [m]')
        ax.set_ylabel('y [m]')
        ax.set_xlim(0, self.xlen)
        ax.set_ylim(0, self.ylen)
        ax.grid(True, alpha=0.3)
    
    def __str__(self):
        """User-friendly string representation."""
        return self.__repr__()


if __name__ == "__main__":
    # Basic test
    print("UDBase module loaded successfully")
    print("\nUsage:")
    print("  from udbase import UDBase")
    print("  sim = UDBase(expnr=65, path='experiments/065')")