"""Shared exception hierarchy for the uDALES Python postprocessing tools.

Library code raises these (never ``sys.exit``), so callers — notebooks, scripts,
workflow managers — can catch failures without the host interpreter being killed.

    UDALESError                     base for everything raised by the tools
    ├── ConfigurationError          invalid/inconsistent case config or namelist
    ├── DataFormatError             an input file does not match the expected format
    ├── DependencyError             a required optional dependency is missing
    ├── GeometryError               invalid or unprocessable geometry
    └── RadiationError              failure in a radiation/shortwave/SEB computation
"""


class UDALESError(Exception):
    """Base class for all errors raised by the uDALES postprocessing tools."""


class ConfigurationError(UDALESError):
    """Invalid or inconsistent case configuration or namelist."""


class DataFormatError(UDALESError):
    """An input/data file does not match the expected format."""


class DependencyError(UDALESError):
    """A required optional dependency is not installed."""


class GeometryError(UDALESError):
    """Invalid or unprocessable geometry."""


class RadiationError(UDALESError):
    """Failure in a radiation, shortwave, or surface-energy-balance computation."""
