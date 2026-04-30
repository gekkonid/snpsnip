#!/usr/bin/env python3
"""
Unit tests for SNPSnip parsing functions.

Tests cover:
- _count_alt_alleles(): genotype string parsing
- _parse_gt_matrix(): genotype matrix parsing from strings and files
- _parse_dp_matrix(): depth matrix parsing from strings and files
- HDF5 roundtrip: write matrices to HDF5, read back, verify values
"""

import os
import sys
import tempfile
import unittest

import h5py
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snpsnip import _count_alt_alleles, _parse_gt_matrix, _parse_dp_matrix


class TestCountAltAlleles(unittest.TestCase):
    """Tests for _count_alt_alleles function."""

    def test_hom_ref_slash(self):
        self.assertEqual(_count_alt_alleles("0/0"), 0)

    def test_hom_ref_pipe(self):
        self.assertEqual(_count_alt_alleles("0|0"), 0)

    def test_het_slash(self):
        self.assertEqual(_count_alt_alleles("0/1"), 1)

    def test_het_pipe(self):
        self.assertEqual(_count_alt_alleles("0|1"), 1)

    def test_het_reversed_slash(self):
        self.assertEqual(_count_alt_alleles("1/0"), 1)

    def test_het_reversed_pipe(self):
        self.assertEqual(_count_alt_alleles("1|0"), 1)

    def test_hom_alt_slash(self):
        self.assertEqual(_count_alt_alleles("1/1"), 2)

    def test_hom_alt_pipe(self):
        self.assertEqual(_count_alt_alleles("1|1"), 2)

    def test_missing_dot(self):
        self.assertEqual(_count_alt_alleles("."), -1)

    def test_missing_slash(self):
        self.assertEqual(_count_alt_alleles("./."), -1)

    def test_missing_pipe(self):
        self.assertEqual(_count_alt_alleles(".|."), -1)

    def test_empty_string(self):
        self.assertEqual(_count_alt_alleles(""), -1)

    def test_none(self):
        self.assertEqual(_count_alt_alleles(None), -1)

    def test_multi_allelic_0_2(self):
        """0/2 should be het (1 alt allele)."""
        self.assertEqual(_count_alt_alleles("0/2"), 1)

    def test_multi_allelic_1_2(self):
        """1/2 should be het (2 alt alleles, but capped at 2)."""
        self.assertEqual(_count_alt_alleles("1/2"), 2)

    def test_multi_allelic_2_2(self):
        """2/2 should be hom-alt (2 alt alleles)."""
        self.assertEqual(_count_alt_alleles("2/2"), 2)

    def test_multi_allelic_0_3(self):
        """0/3 should be het (1 alt allele)."""
        self.assertEqual(_count_alt_alleles("0/3"), 1)

    def test_multi_allelic_3_3(self):
        """3/3 should be hom-alt (2 alt alleles, capped)."""
        self.assertEqual(_count_alt_alleles("3/3"), 2)

    def test_multi_allelic_triploid_0_0_1(self):
        """0/0/1 should be het (1 alt allele)."""
        self.assertEqual(_count_alt_alleles("0/0/1"), 1)

    def test_multi_allelic_triploid_0_1_1(self):
        """0/1/1 should be het (2 alt alleles)."""
        self.assertEqual(_count_alt_alleles("0/1/1"), 2)

    def test_multi_allelic_triploid_1_1_1(self):
        """1/1/1 should be hom-alt (3 alt alleles, capped at 2)."""
        self.assertEqual(_count_alt_alleles("1/1/1"), 2)

    def test_partial_missing_0_dot(self):
        """0/. should be missing."""
        self.assertEqual(_count_alt_alleles("0/."), -1)

    def test_partial_missing_dot_1(self):
        """./.1 should be missing."""
        self.assertEqual(_count_alt_alleles("./1"), -1)

    def test_invalid_string(self):
        """Non-numeric genotype should return -1."""
        self.assertEqual(_count_alt_alleles("a/b"), -1)


class TestParseGtMatrix(unittest.TestCase):
    """Tests for _parse_gt_matrix function."""

    def test_basic_biallelic_string(self):
        """Test basic biallelic genotypes from string."""
        data = "0/0\t0/1\t1/1\n0/0\t./.\t1/1\n"
        result = _parse_gt_matrix(data, n_snps=2, n_samples=3, source='string', dtype=np.int8)
        expected = np.array([
            [0, 1, 2],
            [0, -1, 2],
        ], dtype=np.int8)
        np.testing.assert_array_equal(result, expected)

    def test_pipe_separated(self):
        """Test pipe-separated (phased) genotypes."""
        data = "0|0\t0|1\t1|0\n1|1\t.|.\t0|0\n"
        result = _parse_gt_matrix(data, n_snps=2, n_samples=3, source='string', dtype=np.int8)
        expected = np.array([
            [0, 1, 1],
            [2, -1, 0],
        ], dtype=np.int8)
        np.testing.assert_array_equal(result, expected)

    def test_multi_allelic(self):
        """Test multi-allelic genotypes are correctly parsed."""
        data = "0/0\t0/2\t2/2\n1/2\t0/3\t3/3\n"
        result = _parse_gt_matrix(data, n_snps=2, n_samples=3, source='string', dtype=np.int8)
        expected = np.array([
            [0, 1, 2],
            [2, 1, 2],
        ], dtype=np.int8)
        np.testing.assert_array_equal(result, expected)

    def test_all_missing(self):
        """Test all missing genotypes."""
        data = "./.\t.|.\t.\n.\t./.\t.|.\n"
        result = _parse_gt_matrix(data, n_snps=2, n_samples=3, source='string', dtype=np.int8)
        expected = np.full((2, 3), -1, dtype=np.int8)
        np.testing.assert_array_equal(result, expected)

    def test_float_dtype_nan(self):
        """Test float dtype produces NaN for missing values."""
        data = "0/0\t./.\t1/1\n"
        result = _parse_gt_matrix(data, n_snps=1, n_samples=3, source='string', dtype=float)
        self.assertEqual(result[0, 0], 0.0)
        self.assertTrue(np.isnan(result[0, 1]))
        self.assertEqual(result[0, 2], 2.0)

    def test_custom_missing_sentinel(self):
        """Test custom missing sentinel value."""
        data = "0/0\t./.\t1/1\n"
        result = _parse_gt_matrix(data, n_snps=1, n_samples=3, source='string', dtype=np.int8, missing_sentinel=-9)
        expected = np.array([[0, -9, 2]], dtype=np.int8)
        np.testing.assert_array_equal(result, expected)

    def test_empty_lines_ignored(self):
        """Test empty lines are skipped."""
        data = "0/0\t0/1\n\n1/1\t./.\n\n"
        result = _parse_gt_matrix(data, n_snps=2, n_samples=2, source='string', dtype=np.int8)
        expected = np.array([
            [0, 1],
            [2, -1],
        ], dtype=np.int8)
        np.testing.assert_array_equal(result, expected)

    def test_trailing_tabs(self):
        """Test trailing tabs don't cause issues."""
        data = "0/0\t0/1\t\n1/1\t./.\t\n"
        result = _parse_gt_matrix(data, n_snps=2, n_samples=2, source='string', dtype=np.int8)
        expected = np.array([
            [0, 1],
            [2, -1],
        ], dtype=np.int8)
        np.testing.assert_array_equal(result, expected)

    def test_excess_columns_ignored(self):
        """Test columns beyond n_samples are ignored."""
        data = "0/0\t0/1\t1/1\t0/0\n"
        result = _parse_gt_matrix(data, n_snps=1, n_samples=3, source='string', dtype=np.int8)
        expected = np.array([[0, 1, 2]], dtype=np.int8)
        np.testing.assert_array_equal(result, expected)

    def test_from_file(self):
        """Test parsing from file."""
        data = "0/0\t0/1\t1/1\n./.\t1/0\t0/0\n"
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write(data)
            f.flush()
            result = _parse_gt_matrix(f.name, n_snps=2, n_samples=3, source='file', dtype=np.int8)
            expected = np.array([
                [0, 1, 2],
                [-1, 1, 0],
            ], dtype=np.int8)
            np.testing.assert_array_equal(result, expected)
            os.unlink(f.name)

    def test_shape(self):
        """Test output shape matches n_snps x n_samples."""
        data = "0/0\t0/1\n1/1\t./.\n0/0\t1/1\n"
        result = _parse_gt_matrix(data, n_snps=3, n_samples=2, source='string', dtype=np.int8)
        self.assertEqual(result.shape, (3, 2))


class TestParseDpMatrix(unittest.TestCase):
    """Tests for _parse_dp_matrix function."""

    def test_basic_string(self):
        """Test basic depth values from string."""
        data = "10\t20\t30\n5\t.\t15\n"
        result = _parse_dp_matrix(data, n_snps=2, n_samples=3, source='string')
        expected = np.array([
            [10, 20, 30],
            [5, -1, 15],
        ], dtype=np.int32)
        np.testing.assert_array_equal(result, expected)

    def test_missing_dot(self):
        """Test missing depth (.) becomes -1."""
        data = ".\t10\t.\n"
        result = _parse_dp_matrix(data, n_snps=1, n_samples=3, source='string')
        expected = np.array([[-1, 10, -1]], dtype=np.int32)
        np.testing.assert_array_equal(result, expected)

    def test_empty_values(self):
        """Test empty values become -1."""
        data = "10\t\t30\n"
        result = _parse_dp_matrix(data, n_snps=1, n_samples=3, source='string')
        expected = np.array([[10, -1, 30]], dtype=np.int32)
        np.testing.assert_array_equal(result, expected)

    def test_non_integer(self):
        """Test non-integer values become -1."""
        data = "10\tabc\t30\n"
        result = _parse_dp_matrix(data, n_snps=1, n_samples=3, source='string')
        expected = np.array([[10, -1, 30]], dtype=np.int32)
        np.testing.assert_array_equal(result, expected)

    def test_from_file(self):
        """Test parsing from file."""
        data = "100\t200\n.\t50\n"
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write(data)
            f.flush()
            result = _parse_dp_matrix(f.name, n_snps=2, n_samples=2, source='file')
            expected = np.array([
                [100, 200],
                [-1, 50],
            ], dtype=np.int32)
            np.testing.assert_array_equal(result, expected)
            os.unlink(f.name)

    def test_zero_depth(self):
        """Test zero depth is valid (not missing)."""
        data = "0\t10\t0\n"
        result = _parse_dp_matrix(data, n_snps=1, n_samples=3, source='string')
        expected = np.array([[0, 10, 0]], dtype=np.int32)
        np.testing.assert_array_equal(result, expected)

    def test_shape(self):
        """Test output shape matches n_snps x n_samples."""
        data = "10\t20\n30\t40\n50\t60\n"
        result = _parse_dp_matrix(data, n_snps=3, n_samples=2, source='string')
        self.assertEqual(result.shape, (3, 2))


class TestHDF5Roundtrip(unittest.TestCase):
    """Tests for HDF5 write/read roundtrip of GT and DP matrices."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.h5_path = os.path.join(self.tmpdir, "test.h5")

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmpdir)

    def test_gt_matrix_int8_roundtrip(self):
        """Test int8 GT matrix survives HDF5 roundtrip."""
        gt = np.array([
            [0, 1, 2],
            [-1, 0, 1],
            [2, -1, -1],
        ], dtype=np.int8)

        kw = dict(compression="gzip", compression_opts=4, shuffle=True, chunks=True)
        with h5py.File(self.h5_path, 'w') as h5:
            h5.create_dataset("gt_matrix", data=gt, **kw)

        with h5py.File(self.h5_path, 'r') as h5:
            read_gt = h5["gt_matrix"][()]

        np.testing.assert_array_equal(read_gt, gt)
        self.assertEqual(read_gt.dtype, np.int8)

    def test_dp_matrix_int32_roundtrip(self):
        """Test int32 DP matrix survives HDF5 roundtrip."""
        dp = np.array([
            [10, 20, 30],
            [-1, 5, 15],
            [100, -1, -1],
        ], dtype=np.int32)

        kw = dict(compression="gzip", compression_opts=4, shuffle=True, chunks=True)
        with h5py.File(self.h5_path, 'w') as h5:
            h5.create_dataset("dp_matrix", data=dp, **kw)

        with h5py.File(self.h5_path, 'r') as h5:
            read_dp = h5["dp_matrix"][()]

        np.testing.assert_array_equal(read_dp, dp)
        self.assertEqual(read_dp.dtype, np.int32)

    def test_gt_matrix_missing_sentinel_preserved(self):
        """Test -1 sentinel for missing genotypes is preserved in HDF5."""
        gt = np.full((100, 50), -1, dtype=np.int8)
        gt[0, 0] = 0
        gt[50, 25] = 1
        gt[99, 49] = 2

        kw = dict(compression="gzip", compression_opts=4, shuffle=True, chunks=True)
        with h5py.File(self.h5_path, 'w') as h5:
            h5.create_dataset("gt_matrix", data=gt, **kw)

        with h5py.File(self.h5_path, 'r') as h5:
            read_gt = h5["gt_matrix"][()]

        self.assertEqual(read_gt[0, 0], 0)
        self.assertEqual(read_gt[50, 25], 1)
        self.assertEqual(read_gt[99, 49], 2)
        self.assertEqual(read_gt[0, 1], -1)
        self.assertEqual(read_gt[1, 0], -1)
        self.assertEqual(np.sum(read_gt == -1), 100 * 50 - 3)

    def test_dp_matrix_missing_sentinel_preserved(self):
        """Test -1 sentinel for missing depth is preserved in HDF5."""
        dp = np.full((100, 50), -1, dtype=np.int32)
        dp[0, 0] = 10
        dp[50, 25] = 20
        dp[99, 49] = 30

        kw = dict(compression="gzip", compression_opts=4, shuffle=True, chunks=True)
        with h5py.File(self.h5_path, 'w') as h5:
            h5.create_dataset("dp_matrix", data=dp, **kw)

        with h5py.File(self.h5_path, 'r') as h5:
            read_dp = h5["dp_matrix"][()]

        self.assertEqual(read_dp[0, 0], 10)
        self.assertEqual(read_dp[50, 25], 20)
        self.assertEqual(read_dp[99, 49], 30)
        self.assertEqual(read_dp[0, 1], -1)
        self.assertEqual(read_dp[1, 0], -1)
        self.assertEqual(np.sum(read_dp == -1), 100 * 50 - 3)

    def test_gt_matrix_large_roundtrip(self):
        """Test large GT matrix roundtrip with mixed values."""
        np.random.seed(42)
        gt = np.random.choice([-1, 0, 1, 2], size=(500, 100), p=[0.05, 0.45, 0.3, 0.2]).astype(np.int8)

        kw = dict(compression="gzip", compression_opts=4, shuffle=True, chunks=True)
        with h5py.File(self.h5_path, 'w') as h5:
            h5.create_dataset("gt_matrix", data=gt, **kw)

        with h5py.File(self.h5_path, 'r') as h5:
            read_gt = h5["gt_matrix"][()]

        np.testing.assert_array_equal(read_gt, gt)

    def test_dp_matrix_large_roundtrip(self):
        """Test large DP matrix roundtrip with mixed values."""
        np.random.seed(42)
        dp = np.random.randint(-1, 100, size=(500, 100)).astype(np.int32)
        dp[dp == -1] = -1  # ensure some -1 values

        kw = dict(compression="gzip", compression_opts=4, shuffle=True, chunks=True)
        with h5py.File(self.h5_path, 'w') as h5:
            h5.create_dataset("dp_matrix", data=dp, **kw)

        with h5py.File(self.h5_path, 'r') as h5:
            read_dp = h5["dp_matrix"][()]

        np.testing.assert_array_equal(read_dp, dp)

    def test_r_style_missing_conversion(self):
        """Test R-style -1 to NA conversion simulation."""
        gt = np.array([
            [0, 1, -1],
            [-1, 2, 0],
            [1, -1, 2],
        ], dtype=np.int8)

        kw = dict(compression="gzip", compression_opts=4, shuffle=True, chunks=True)
        with h5py.File(self.h5_path, 'w') as h5:
            h5.create_dataset("gt_matrix", data=gt, **kw)

        with h5py.File(self.h5_path, 'r') as h5:
            read_gt = h5["gt_matrix"][()]

        gt_with_na = read_gt.astype(np.int16)
        gt_with_na[gt_with_na == -1] = -999

        self.assertEqual(gt_with_na[0, 0], 0)
        self.assertEqual(gt_with_na[0, 1], 1)
        self.assertEqual(gt_with_na[0, 2], -999)
        self.assertEqual(gt_with_na[1, 0], -999)
        self.assertEqual(gt_with_na[1, 1], 2)
        self.assertEqual(gt_with_na[1, 2], 0)

    def test_hdf5_full_structure(self):
        """Test full HDF5 file structure matches expected layout."""
        gt = np.array([[0, 1, 2], [-1, 0, 1]], dtype=np.int8)
        dp = np.array([[10, 20, 30], [-1, 5, 15]], dtype=np.int32)
        samples = np.array(["sample1", "sample2", "sample3"], dtype=object)
        str_dt = h5py.string_dtype()

        kw = dict(compression="gzip", compression_opts=4, shuffle=True, chunks=True)
        with h5py.File(self.h5_path, 'w') as h5:
            h5.attrs["snpsnip_version"] = "1.0.0"
            h5.create_dataset("has_dp_matrix", data=1)
            h5.create_dataset("samples", data=samples, dtype=str_dt)
            h5.create_dataset("gt_matrix", data=gt, **kw)
            h5.create_dataset("dp_matrix", data=dp, **kw)

            sg = h5.create_group("sample_stats")
            sg.create_dataset("id", data=samples, dtype=str_dt)
            sg.create_dataset("missing_rate", data=np.array([0.1, 0.2, 0.3]))

        with h5py.File(self.h5_path, 'r') as h5:
            self.assertEqual(h5.attrs["snpsnip_version"], "1.0.0")
            self.assertEqual(h5["has_dp_matrix"][()], 1)
            read_samples = h5["samples"][()]
            decoded = [s.decode() if isinstance(s, bytes) else s for s in read_samples]
            self.assertEqual(decoded, ["sample1", "sample2", "sample3"])
            np.testing.assert_array_equal(h5["gt_matrix"][()], gt)
            np.testing.assert_array_equal(h5["dp_matrix"][()], dp)
            self.assertIn("sample_stats/id", h5)
            self.assertIn("sample_stats/missing_rate", h5)


if __name__ == "__main__":
    unittest.main()
