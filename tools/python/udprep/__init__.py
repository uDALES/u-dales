"""uDALES preprocessing package.

Public names are imported lazily (PEP 562) so that ``import udprep`` — or reaching
a pure helper such as ``udprep._radiation_compute`` — does not eagerly import the
numba-JIT direct-shortwave kernels. The kernels are only compiled when the
direct-shortwave path is actually used.
"""

__all__ = ["DirectShortwaveSolver", "UDPrep", "vegetation_block_to_veg"]

_LAZY = {
    "DirectShortwaveSolver": (".directshortwave", "DirectShortwaveSolver"),
    "UDPrep": (".udprep", "UDPrep"),
    "vegetation_block_to_veg": (".udprep_vegetation", "vegetation_block_to_veg"),
}


def __getattr__(name):
    try:
        module_name, attr = _LAZY[name]
    except KeyError:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}") from None
    import importlib

    value = getattr(importlib.import_module(module_name, __name__), attr)
    globals()[name] = value  # cache for subsequent access
    return value


def __dir__():
    return sorted(set(list(globals().keys()) + __all__))
