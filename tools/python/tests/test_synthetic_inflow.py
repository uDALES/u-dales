from __future__ import annotations

import struct
import tempfile
import unittest
from pathlib import Path

import numpy as np

from _common import PYTHON_DIR  # noqa: F401  (bootstraps sys.path via conftest)

from exceptions import DataFormatError
from udprep.synthetic_inflow import (
    LengthTimeScales,
    ReynoldsStressProfile,
    _autocorr_1d,
    _integral_scale_to_first_zero,
    assemble_driver_plane,
    centres_from_prof,
    edges_from_prof,
    integral_length_scale_periodic,
    integral_time_scale,
    jh_kh_for_advection,
    profile_from_driver_files,
    profile_from_tdump,
    read_length_time_scales_file,
    read_reynolds_stress_file,
    scales_from_driver_files,
    uniform_edges,
    uniform_zf_zh,
    write_length_time_scales_file,
    write_reynolds_stress_file,
)

# The real GMD 2024 indoor-outdoor prof.inp this fix targets: 128 levels,
# dz = 3.5 mm near the wall, exponentially stretched up to zsize = 0.96 m.
# Skipped if not present on the current machine (it lives on ephemeral
# scratch, not in the repo).
_GMD_PROF_INP = Path(
    "/rds/general/ephemeral/user/mvr/ephemeral/gmd-indoor/prod567/prof.inp.567"
)


def _make_profile(ktot=5, zsize=5.0):
    z = uniform_edges(zsize, ktot)
    umean = 0.2 * z
    r11 = 0.05 + 0.0 * z
    r22 = 0.03 + 0.0 * z
    r33 = 0.02 + 0.0 * z
    r21 = 0.3 * np.sqrt(r11 * r22)
    r31 = 0.2 * np.sqrt(r11 * r33)
    r32 = 0.1 * np.sqrt(r22 * r33)
    return ReynoldsStressProfile(z=z, umean=umean, R11=r11, R21=r21, R22=r22, R31=r31, R32=r32, R33=r33)


class TestReynoldsStressFileRoundtrip(unittest.TestCase):
    def test_write_read_roundtrip(self):
        profile = _make_profile(ktot=6, zsize=6.0)
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "Reynolds_stress_profiles_velocity.txt"
            write_reynolds_stress_file(path, profile)

            # header + ktot+1 data rows, 8 columns, parseable with skiprows=1
            raw = np.loadtxt(path, skiprows=1)
            self.assertEqual(raw.shape, (7, 8))

            back = read_reynolds_stress_file(path)
            np.testing.assert_allclose(back.z, profile.z)
            np.testing.assert_allclose(back.umean, profile.umean)
            np.testing.assert_allclose(back.R11, profile.R11)
            np.testing.assert_allclose(back.R21, profile.R21)
            np.testing.assert_allclose(back.R33, profile.R33)

    def test_header_is_single_line(self):
        profile = _make_profile()
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "f.txt"
            write_reynolds_stress_file(path, profile)
            with path.open() as fh:
                header = fh.readline()
            self.assertFalse(header.startswith(" "))
            self.assertEqual(header.count("\n"), 1)


class TestValidation(unittest.TestCase):
    def test_rejects_nonzero_bottom_edge(self):
        profile = _make_profile()
        profile.z = profile.z.copy()
        profile.z[0] = 0.1
        with self.assertRaises(DataFormatError):
            profile.validate()

    def test_rejects_non_monotone_z(self):
        profile = _make_profile()
        profile.z = profile.z.copy()
        profile.z[2] = profile.z[1]
        with self.assertRaises(DataFormatError):
            profile.validate()

    def test_rejects_negative_variance(self):
        profile = _make_profile()
        profile.R11 = profile.R11.copy()
        profile.R11[1] = -1.0
        with self.assertRaises(DataFormatError):
            profile.validate()

    def test_rejects_cauchy_schwarz_violation(self):
        profile = _make_profile()
        profile.R21 = profile.R21.copy()
        profile.R21[1] = 10.0 * np.sqrt(profile.R11[1] * profile.R22[1]) + 1.0
        with self.assertRaises(DataFormatError):
            profile.validate()

    def test_length_time_scale_rejects_zero_scale(self):
        z = uniform_edges(4.0, 4)
        scales = LengthTimeScales(
            z=z, nl_y=np.ones_like(z, dtype=int), nl_z=np.ones_like(z, dtype=int), t_scale=np.zeros_like(z)
        )
        with self.assertRaises(DataFormatError):
            scales.validate()

    def test_length_time_scale_roundtrip(self):
        z = uniform_edges(4.0, 4)
        scales = LengthTimeScales(
            z=z, nl_y=np.array([1, 2, 3, 4, 4]), nl_z=np.array([1, 1, 2, 2, 3]), t_scale=0.5 * np.ones_like(z)
        )
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "length_time_scales_u.txt"
            write_length_time_scales_file(path, scales)
            raw = np.loadtxt(path, skiprows=1)
            self.assertEqual(raw.shape, (5, 4))
            back = read_length_time_scales_file(path)
            np.testing.assert_array_equal(back.nl_y, scales.nl_y)
            np.testing.assert_array_equal(back.nl_z, scales.nl_z)
            np.testing.assert_allclose(back.t_scale, scales.t_scale)


class TestAdvectionHaloWidths(unittest.TestCase):
    def test_cd2_gives_jh1_kh1(self):
        self.assertEqual(jh_kh_for_advection(2), (1, 1))

    def test_kappa_gives_jh2_kh1(self):
        self.assertEqual(jh_kh_for_advection(7), (2, 1))

    def test_fluxlimiter_gives_jh2_kh2(self):
        self.assertEqual(jh_kh_for_advection(5), (2, 2))


class TestProfileFromTdump(unittest.TestCase):
    def test_reproduces_profile_at_iplane(self):
        try:
            from netCDF4 import Dataset
        except ImportError:
            self.skipTest("netCDF4 not installed")

        nx, ny, nzt, ntime = 6, 4, 8, 3
        zsize, ktot = 8.0, 8
        zt = (np.arange(nzt) + 0.5) * (zsize / nzt)
        zm = np.arange(nzt) * (zsize / nzt)
        iplane = 3

        rng = np.random.default_rng(0)

        def make_var(scale):
            # varies by x (so selecting iplane matters) and by z; identical
            # across y and time so the y/time averaging is a no-op we can
            # check against exactly.
            per_x = 1.0 + 0.01 * np.arange(nx)
            per_z = scale * (1.0 + 0.1 * np.arange(nzt))
            base = np.outer(per_x, per_z)  # (x, z)
            full = np.broadcast_to(base, (ntime, ny, nx, nzt)).copy()
            full += 0.0 * rng.standard_normal(full.shape)  # keep deterministic
            return full  # dims order: (time, yt, xt, zt) on purpose (not x-first)

        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "tdump.001.nc"
            with Dataset(path, "w") as ds:
                ds.createDimension("time", ntime)
                ds.createDimension("yt", ny)
                ds.createDimension("xt", nx)
                ds.createDimension("zt", nzt)
                ds.createDimension("zm", nzt)
                ds.createVariable("zt", "f8", ("zt",))[:] = zt
                ds.createVariable("zm", "f8", ("zm",))[:] = zm

                for name, dimz, scale in [
                    ("ut", "zt", 1.0),
                    ("upuptc", "zt", 0.10),
                    ("vpvptc", "zt", 0.08),
                    ("wpwptc", "zt", 0.06),
                    ("upvpt", "zt", 0.01),
                    ("upwpt", "zm", 0.02),
                    ("vpwpt", "zm", 0.015),
                ]:
                    v = ds.createVariable(name, "f8", ("time", "yt", "xt", dimz))
                    v[:] = make_var(scale)

            profile = profile_from_tdump(path, iplane=iplane, zsize=zsize, ktot=ktot)
            profile.validate()

            # ground truth at the chosen x index, interpolated the same way
            per_x_val = 1.0 + 0.01 * iplane
            z_target = uniform_edges(zsize, ktot)

            def expected(scale, z_source):
                vals = per_x_val * scale * (1.0 + 0.1 * np.arange(nzt))
                out = np.interp(z_target, z_source, vals)
                out[0] = 0.0
                slope = (vals[-1] - vals[-2]) / (z_source[-1] - z_source[-2])
                out[-1] = vals[-1] + slope * (z_target[-1] - z_source[-1])
                return out

            np.testing.assert_allclose(profile.z, z_target)
            np.testing.assert_allclose(profile.umean, expected(1.0, zt), rtol=1e-8)
            np.testing.assert_allclose(profile.R11, expected(0.10, zt), rtol=1e-8)
            np.testing.assert_allclose(profile.R33, expected(0.06, zt), rtol=1e-8)
            np.testing.assert_allclose(profile.R31, expected(0.02, zm), rtol=1e-8)


def _write_driver_binary(directory: Path, expnr: int, nprocy: int, jtot: int, ktot: int, jh: int, kh: int, t, u, v, w):
    """Write fake tdriver_000/{u,v,w}driver_<ipy> files in the solver's layout.

    ``u``, ``v``, ``w`` each have shape (nt, jtot, ktot) -- already the
    *physical* (halo-free, global-y) field; halos are zero-filled since
    :func:`read_driver_component` strips them before use.
    """
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    nt = len(t)
    (directory / f"tdriver_000.{expnr:03d}").write_bytes(
        b"".join(struct.pack("<d", float(tt)) for tt in t)
    )
    jtot_local = jtot // nprocy
    n_j = jtot_local + 2 * jh
    n_k = ktot + 2 * kh
    for ipy in range(nprocy):
        j0 = ipy * jtot_local
        for comp, field in (("u", u), ("v", v), ("w", w)):
            path = directory / f"{comp}driver_{ipy:03d}.{expnr:03d}"
            with path.open("wb") as f:
                for it in range(nt):
                    rec = np.zeros((n_j, n_k), dtype="<f8")
                    rec[jh : jh + jtot_local, kh : kh + ktot] = field[it, j0 : j0 + jtot_local, :]
                    f.write(rec.tobytes(order="F"))


class TestDriverFilePath(unittest.TestCase):
    def setUp(self):
        self.rng = np.random.default_rng(42)

        self.expnr = 123
        self.nprocy = 2
        self.jtot = 20
        self.ktot = 6
        self.jh, self.kh = 1, 1
        self.zsize = 6.0
        self.ylen = 20.0
        self.dt = 0.05
        self.nt = 4000
        self.dy = self.ylen / self.jtot
        self.dz = self.zsize / self.ktot
        self.T_true = 1.0  # s, AR(1) correlation time
        self.wavelength = 10.0  # m, sinusoidal y structure

        zf, zh = uniform_zf_zh(self.zsize, self.ktot)
        self.zf, self.zh = zf, zh

        # AR(1) process, rho(lag) = exp(-lag*dt/T) exactly at sampled lags
        phi = np.exp(-self.dt / self.T_true)
        s = np.empty(self.nt)
        s[0] = self.rng.standard_normal()
        noise = self.rng.standard_normal(self.nt - 1)
        for n in range(1, self.nt):
            s[n] = phi * s[n - 1] + np.sqrt(1.0 - phi**2) * noise[n - 1]
        s /= s.std()

        y = np.arange(self.jtot) * self.dy
        g = np.cos(2.0 * np.pi * y / self.wavelength)

        self.u_mean_true = 0.3 * zf  # m/s, zero at z=0
        self.v_mean_true = np.zeros_like(zf)
        self.w_mean_true = np.zeros_like(zh)

        sigma_u, sigma_v, sigma_w = 0.5, 0.3, 0.2
        st_g = np.outer(s, g)  # (nt, jtot)

        u = self.u_mean_true[None, None, :] + sigma_u * st_g[:, :, None] * np.ones((1, 1, self.ktot))
        v = self.v_mean_true[None, None, :] + sigma_v * st_g[:, :, None] * np.ones((1, 1, self.ktot))
        w = self.w_mean_true[None, None, :] + sigma_w * st_g[:, :, None] * np.ones((1, 1, self.ktot))

        self.sigma_u, self.sigma_v, self.sigma_w = sigma_u, sigma_v, sigma_w
        self.g_var = float(np.mean(g * g))  # ~0.5 for this wavelength/domain

        self.tmpdir = tempfile.TemporaryDirectory()
        self.directory = Path(self.tmpdir.name)
        t = np.arange(self.nt) * self.dt
        _write_driver_binary(
            self.directory, self.expnr, self.nprocy, self.jtot, self.ktot, self.jh, self.kh, t, u, v, w
        )

    def tearDown(self):
        self.tmpdir.cleanup()

    def test_assemble_driver_plane_shape_and_values(self):
        t, field = assemble_driver_plane(
            self.directory, self.expnr, "u", self.nprocy, self.jtot, self.ktot, self.jh, self.kh
        )
        self.assertEqual(field.shape, (self.nt, self.jtot, self.ktot))
        self.assertEqual(len(t), self.nt)
        np.testing.assert_allclose(t, np.arange(self.nt) * self.dt)

    def test_time_record_is_increasing(self):
        t, _ = assemble_driver_plane(
            self.directory, self.expnr, "u", self.nprocy, self.jtot, self.ktot, self.jh, self.kh
        )
        self.assertTrue(np.all(np.diff(t) > 0))

    def test_profile_recovers_mean_and_stresses(self):
        profile = profile_from_driver_files(
            self.directory, self.expnr, self.nprocy, self.jtot, self.ktot, self.zsize, self.jh, self.kh
        )
        profile.validate()

        # bottom edge is forced to zero by convention; check interior + top
        np.testing.assert_allclose(profile.umean[1:], 0.3 * profile.z[1:], atol=0.02)

        expected_R11 = self.sigma_u**2 * self.g_var
        expected_R22 = self.sigma_v**2 * self.g_var
        expected_R33 = self.sigma_w**2 * self.g_var
        expected_R21 = self.sigma_u * self.sigma_v * self.g_var

        np.testing.assert_allclose(profile.R11[1:], expected_R11, rtol=0.15)
        np.testing.assert_allclose(profile.R22[1:], expected_R22, rtol=0.15)
        np.testing.assert_allclose(profile.R33[1:], expected_R33, rtol=0.15)
        np.testing.assert_allclose(profile.R21[1:], expected_R21, rtol=0.15)

    def test_length_time_scales_recover_known_time_and_y_length_scale(self):
        scales_u, scales_v, scales_w = scales_from_driver_files(
            self.directory,
            self.expnr,
            self.nprocy,
            self.jtot,
            self.ktot,
            self.zsize,
            self.ylen,
            self.dt,
            self.jh,
            self.kh,
        )
        for scales in (scales_u, scales_v, scales_w):
            scales.validate()
            # recovered integral time scale should be within a generous
            # tolerance of the AR(1) process's true correlation time
            np.testing.assert_allclose(scales.t_scale, self.T_true, rtol=0.5)

        # analytic reference for the y integral length scale of a pure
        # periodic cosine of the same wavelength/domain/grid, run through the
        # same "integrate the autocorrelation to its first zero crossing"
        # helper used by the estimator -- the two should closely agree since
        # the AR(1) time-noise cancels exactly out of the *normalized* y
        # correlation (it is a common factor at every y).
        lag_m = np.arange(self.jtot // 2 + 1) * self.dy
        rho_analytic = np.cos(2.0 * np.pi * lag_m / self.wavelength)
        l_y_expected = _integral_scale_to_first_zero(rho_analytic, self.dy, max_lag=self.jtot // 2)
        nl_y_expected = max(int(l_y_expected / self.dy), 1)

        for scales in (scales_u, scales_v, scales_w):
            self.assertEqual(int(scales.nl_y[3]), nl_y_expected)

    def test_output_files_parse_back_with_loadtxt(self):
        profile = profile_from_driver_files(
            self.directory, self.expnr, self.nprocy, self.jtot, self.ktot, self.zsize, self.jh, self.kh
        )
        scales_u, scales_v, scales_w = scales_from_driver_files(
            self.directory,
            self.expnr,
            self.nprocy,
            self.jtot,
            self.ktot,
            self.zsize,
            self.ylen,
            self.dt,
            self.jh,
            self.kh,
        )
        out_dir = self.directory / "syntheticInflow_inputs"
        from udprep.synthetic_inflow import write_reynolds_stress_and_scales

        write_reynolds_stress_and_scales(out_dir, profile, scales_u, scales_v, scales_w)

        rs = np.loadtxt(out_dir / "Reynolds_stress_profiles_velocity.txt", skiprows=1)
        self.assertEqual(rs.shape, (self.ktot + 1, 8))
        for name in ("u", "v", "w"):
            sc = np.loadtxt(out_dir / f"length_time_scales_{name}.txt", skiprows=1)
            self.assertEqual(sc.shape, (self.ktot + 1, 4))


class TestEdgesFromProf(unittest.TestCase):
    def test_reproduces_gmd_grid(self):
        if not _GMD_PROF_INP.is_file():
            self.skipTest(f"{_GMD_PROF_INP} not available on this machine")

        zf = centres_from_prof(_GMD_PROF_INP)
        zh = edges_from_prof(_GMD_PROF_INP)

        self.assertEqual(len(zf), 128)
        self.assertEqual(len(zh), 129)
        # dz = 3.5mm near the wall: zf[0] = dz/2, zf[1] = 3*dz/2.
        self.assertAlmostEqual(zf[0], 0.00175, places=6)
        self.assertAlmostEqual(zf[1], 0.00525, places=6)
        # bottom edge is exactly 0, top edge is the file's own zsize (0.96),
        # to within round-off from the recursive zh(k+1)=zh(k)+2*(zf(k)-zh(k))
        # construction.
        self.assertEqual(zh[0], 0.0)
        self.assertAlmostEqual(zh[-1], 0.96, places=9)
        self.assertTrue(np.all(np.diff(zh) > 0.0))


class TestStretchedGridDriverPath(unittest.TestCase):
    """Reproduces the GMD-style stretched vertical grid (linear near the wall,
    then exponential stretching), on which the generator's input files must be
    placed on the target's real prof.inp edges, not on a uniform grid.
    """

    def setUp(self):
        self.expnr = 555
        self.nprocy = 2
        self.jtot = 20
        self.ktot = 8
        self.jh, self.kh = 1, 1
        self.zsize = 6.0
        self.ylen = 20.0
        self.dt = 0.1
        self.nt = 3

        # stretched edges: zh = zsize * (k/ktot)**1.5 (same recipe as the GMD
        # case's exponential stretching, simplified to a single power law for
        # a compact, deterministic test grid).
        k = np.arange(self.ktot + 1)
        self.zh = self.zsize * (k / self.ktot) ** 1.5
        self.zf = 0.5 * (self.zh[:-1] + self.zh[1:])  # native u/v cell centres
        self.zh_native = self.zh[:-1]  # native w edges

        self.u_true = lambda z: 0.3 * z  # known, exactly-linear mean profile

        # no turbulence: an exactly-known mean with zero variance keeps the
        # Reynolds-stress validation trivially satisfied and isolates the
        # umean check from interpolation noise unrelated to the z-grid bug.
        u_native = self.u_true(self.zf)
        u = np.broadcast_to(u_native[None, None, :], (self.nt, self.jtot, self.ktot)).copy()
        v = np.zeros_like(u)
        w = np.zeros_like(u)

        self.tmpdir = tempfile.TemporaryDirectory()
        self.directory = Path(self.tmpdir.name)
        t = np.arange(self.nt) * self.dt
        _write_driver_binary(
            self.directory, self.expnr, self.nprocy, self.jtot, self.ktot, self.jh, self.kh, t, u, v, w
        )

    def tearDown(self):
        self.tmpdir.cleanup()

    def test_zh_argument_is_required_and_recovers_true_profile(self):
        # Without zh, profile_from_driver_files resynthesises a *uniform*
        # ktot+1 edge grid from zsize/ktot and assumes the driver planes'
        # native data sit on that same uniform grid -- but the fixture above
        # recorded its values at the real *stretched* zf/zh positions, so the
        # two disagree and the recovered profile is wrong (this is exactly
        # the bug: the generator reads these files by index against its own
        # prof.inp-derived stretched grid, not a uniform one).
        stale_profile = profile_from_driver_files(
            self.directory, self.expnr, self.nprocy, self.jtot, self.ktot, self.zsize, self.jh, self.kh
        )
        with self.assertRaises(AssertionError):
            np.testing.assert_allclose(
                stale_profile.umean, self.u_true(stale_profile.z), atol=1e-6
            )

        # With the real zh supplied, both the native source grid (zf/zh from
        # zh) and the output target grid are the true stretched grid, so a
        # piecewise-linear interpolation of an exactly-linear function is
        # exact (up to floating point) everywhere, including the
        # linearly-extrapolated top edge.
        profile = profile_from_driver_files(
            self.directory,
            self.expnr,
            self.nprocy,
            self.jtot,
            self.ktot,
            self.zsize,
            self.jh,
            self.kh,
            zh=self.zh,
        )
        profile.validate()
        np.testing.assert_allclose(profile.z, self.zh)
        np.testing.assert_allclose(profile.umean, self.u_true(profile.z), atol=1e-8)


class TestCLIFromDriverWithProf(unittest.TestCase):
    def test_from_driver_with_prof_flag_uses_prof_grid_and_derives_zsize(self):
        from udprep.synthetic_inflow import main

        expnr = 777
        nprocy = 1
        jtot = 10
        ktot = 5
        jh, kh = 1, 1
        zsize = 5.0
        ylen = 10.0
        dt = 0.1
        nt = 8

        k = np.arange(ktot + 1)
        zh = zsize * (k / ktot) ** 1.5
        zf = 0.5 * (zh[:-1] + zh[1:])

        # A small travelling-wave fluctuation (period = jtot in y, varying
        # with t) gives length/time_scales_from_driver_files something
        # non-degenerate to estimate (a flat, zero-variance field trips its
        # "t_scale must be > 0" validation) -- its mean over one full y
        # period is exactly 0 to float round-off, so it does not perturb the
        # mean-profile check below.
        j = np.arange(jtot)
        tt = np.arange(nt)
        fluct = 0.05 * np.cos(2.0 * np.pi * (j[None, :, None] + tt[:, None, None]) / jtot)

        u = np.broadcast_to((0.2 * zf)[None, None, :], (nt, jtot, ktot)).copy() + fluct
        v = np.zeros((nt, jtot, ktot)) + fluct
        w = np.zeros((nt, jtot, ktot)) + fluct

        with tempfile.TemporaryDirectory() as d:
            directory = Path(d)
            t = np.arange(nt) * dt
            _write_driver_binary(directory, expnr, nprocy, jtot, ktot, jh, kh, t, u, v, w)

            # Minimal fake prof.inp.<expnr>: two header lines (content is
            # irrelevant, only the count matters, matching the Fortran and
            # centres_from_prof), then one row per level with column 1 = zf.
            prof_path = directory / f"prof.inp.{expnr}"
            with prof_path.open("w") as f:
                f.write("# SDBL flow\n")
                f.write("# z thl qt u v tke\n")
                for zfk in zf:
                    f.write(f"{zfk:.10f} 288.000000 0.000000 0.000000 0.000000 0.000000\n")

            out_dir = directory / "syntheticInflow_inputs"
            rc = main(
                [
                    "from-driver",
                    str(directory),
                    "--expnr", str(expnr),
                    "--nprocy", str(nprocy),
                    "--jtot", str(jtot),
                    "--ktot", str(ktot),
                    "--ylen", str(ylen),
                    "--dt", str(dt),
                    "--jh", str(jh),
                    "--kh", str(kh),
                    "--prof", str(prof_path),
                    "--output", str(out_dir),
                ]
            )
            self.assertEqual(rc, 0)

            rs = np.loadtxt(out_dir / "Reynolds_stress_profiles_velocity.txt", skiprows=1)
            self.assertEqual(rs.shape, (ktot + 1, 8))
            np.testing.assert_allclose(rs[:, 0], zh, atol=1e-8)
            np.testing.assert_allclose(rs[1:, 1], 0.2 * zh[1:], atol=1e-6)


class TestAutocorrelationHelpers(unittest.TestCase):
    def test_autocorr_of_white_noise_decays_fast(self):
        rng = np.random.default_rng(1)
        x = rng.standard_normal(2000)
        rho = _autocorr_1d(x)
        self.assertAlmostEqual(rho[0], 1.0, places=6)
        self.assertLess(abs(rho[1]), 0.2)

    def test_integral_scale_of_exponential_matches_tau(self):
        tau = 2.0
        dt = 0.1
        lag = np.arange(200) * dt
        rho = np.exp(-lag / tau)
        # never crosses zero analytically; the helper falls back to the full
        # window, so choose max_lag such that rho is already tiny there
        scale = _integral_scale_to_first_zero(rho, dt, max_lag=len(rho) - 1)
        self.assertAlmostEqual(scale, tau, delta=0.1)

    def test_periodic_length_scale_of_pure_cosine(self):
        n = 40
        dy = 0.5
        wavelength = 8.0
        y = np.arange(n) * dy
        g = np.cos(2 * np.pi * y / wavelength)
        field = np.tile(g, (10, 1))  # (t, y), identical rows
        L = integral_length_scale_periodic(field, dy)
        self.assertGreater(L, 0.0)
        self.assertLess(L, wavelength / 2.0)


if __name__ == "__main__":
    unittest.main()
