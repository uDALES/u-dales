# fielddump output-time diagnostics

`test_fielddump_diagnostics.py` checks the fielddump fields that the model does
not store but computes when the dump is written: the divergence of the
corrected velocity and its three terms (`di`, `ux`, `vy`, `wz`) and the
right-hand side of the Poisson equation and its three terms (`rs`, `du`, `dv`,
`dw`). It also exercises the pointer-only codes `pd`, `pu`, `pv`, `pw`
(pressure correction and predicted velocities of the last substep).

Case 300 (16x16x32, stretched grid, no buildings) runs for 2 s on a 2x2
decomposition, once with all of these codes in `fieldvars` and once without
any fielddump. The checks are:

- each derivative equals the finite difference of the field it is derived
  from, as written in the same record (`dudx` against `u`, `dpupdx` against
  `pup`, and so on), on the interior cells where both sides are available;
- each divergence equals the sum of its three written terms;
- every value is finite, including the first record;
- the restart files of the two runs are byte-identical, so requesting the
  diagnostics does not alter the solution.

The dump is written in single precision, so all tolerances are float32-scale
(relative 1e-5 of the field's maximum).

Runtime class: CI integration (fast, four ranks). Environment: `UDALES_BUILD`,
`UDALES_RUNTIME_MODULES`, `MPIEXEC`, as for the other integration tests.
