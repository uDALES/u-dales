# Parallelisation

Parallelisation of u-DALES 2.0 is achieved using the [2DECOMP&FFT library](https://github.com/2decomp-fft/2decomp-fft).
The 2DECOMP&FFT library decomposes Cartesian domains using a 2-D pencil-based decomposition, assigning a pencil to each MPI rank.
Operations on each pencil can then be performed in parallel.
Users should ensure the 2-D decomposition grid factors the mesh exactly to ensure load balance.
By doing so each pencil contains an equal number of grid points for an equal volume of work per MPI rank.

Two communication patterns implemented by 2DECOMP&FFT are used in u-DALES 2.0.

## Data transposes

Solving the Poisson problem using FFTs requires access to the full extent of an axis in each orientation.
The 2DECOMP&FFT library implements transpose operations from Z->Y, Y->X, X->Y and Y->Z pencil orientations allowing the solver to process each orientation in turn.

## Halo exchanges

The discretisation stencils in the transport equations require neighbour data.
At pencil boundaries this requires exchanging data between neighbouring pencils.
The 2DECOMP&FFT library provides halo exchange subroutines to support such data exchange.
