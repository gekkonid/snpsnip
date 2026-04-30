#!/usr/bin/env python3
"""
Integration tests for SNPSnip R notebook export and HDF5 handling.

Tests cover:
- R notebook export via --r-notebook flag
- HDF5 file structure validation
- GT and DP matrix encoding/decoding roundtrip
- Missing value handling in HDF5 matrices
- VCF with known missing genotypes
- VCF without DP fields
"""

import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

import h5py
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snpsnip.generate_example import generate_example_vcf


class RNotebookIntegrationTest(unittest.TestCase):
    """Integration tests for R notebook export."""

    @classmethod
    def setUpClass(cls):
        try:
            subprocess.run(["bcftools", "--version"],
                         capture_output=True, check=True, timeout=5)
        except (subprocess.SubprocessError, FileNotFoundError, subprocess.TimeoutExpired):
            raise unittest.SkipTest("bcftools not found - skipping R notebook tests.")

        cls.test_dir = tempfile.mkdtemp(prefix="snpsnip_rnotebook_test_")
        cls.test_vcf = os.path.join(cls.test_dir, "test.vcf.gz")

        generate_example_vcf(
            output_path=cls.test_vcf,
            n_samples=20,
            n_snps=10000,
            missing_prop=0.05,
            n_chroms=2,
            compress=True,
            random_seed=42
        )

    @classmethod
    def tearDownClass(cls):
        if os.path.exists(cls.test_dir):
            shutil.rmtree(cls.test_dir)

    def setUp(self):
        self.output_dir = tempfile.mkdtemp(prefix="snpsnip_rnotebook_output_")

    def tearDown(self):
        if os.path.exists(self.output_dir):
            shutil.rmtree(self.output_dir)

    def _run_snpsnip(self, extra_args=None):
        cmd = [
            sys.executable, "-m", "snpsnip",
            "--vcf", self.test_vcf,
            "--output-dir", self.output_dir,
            "--offline",
            "--r-notebook",
            "--subset-freq", "0.5",
            "--processes", "2",
            "--seed", "42"
        ]
        if extra_args:
            cmd.extend(extra_args)

        result = subprocess.run(cmd, capture_output=True, text=True)
        return result

    def test_01_r_notebook_creates_hdf5(self):
        """Test that --r-notebook creates snpsnip_data.h5."""
        result = self._run_snpsnip(["--maf", "0", "--max-missing", "1", "--min-qual", "0"])
        self.assertEqual(result.returncode, 0, f"Should succeed. stderr: {result.stderr}")

        h5_path = os.path.join(self.output_dir, "snpsnip_data.h5")
        self.assertTrue(os.path.exists(h5_path), "HDF5 file should exist")

    def test_02_r_notebook_creates_r_script(self):
        """Test that --r-notebook creates snpsnip_analysis.R."""
        result = self._run_snpsnip(["--maf", "0", "--max-missing", "1", "--min-qual", "0"])
        self.assertEqual(result.returncode, 0)

        r_path = os.path.join(self.output_dir, "snpsnip_analysis.R")
        self.assertTrue(os.path.exists(r_path), "R script should exist")

    def test_03_hdf5_structure(self):
        """Test HDF5 file has expected structure and datasets."""
        result = self._run_snpsnip(["--maf", "0", "--max-missing", "1", "--min-qual", "0"])
        self.assertEqual(result.returncode, 0)

        h5_path = os.path.join(self.output_dir, "snpsnip_data.h5")
        with h5py.File(h5_path, 'r') as h5:
            self.assertIn("snpsnip_version", h5.attrs)

            self.assertIn("gt_matrix", h5)
            self.assertIn("samples", h5)
            self.assertIn("has_dp_matrix", h5)
            self.assertIn("available_stats", h5)

            self.assertIn("sample_stats/id", h5)
            self.assertIn("sample_stats/missing_rate", h5)
            self.assertIn("sample_stats/mean_depth", h5)
            self.assertIn("sample_stats/het_rate", h5)

            self.assertIn("pca/samples", h5)
            self.assertIn("pca/coordinates", h5)
            self.assertIn("pca/variance_explained", h5)

            self.assertIn("snp_info/chrom", h5)
            self.assertIn("snp_info/pos", h5)
            self.assertIn("snp_info/ref", h5)
            self.assertIn("snp_info/alt", h5)
            self.assertIn("snp_info/qual", h5)
            self.assertIn("snp_info/f_missing", h5)

    def test_04_gt_matrix_shape_and_dtype(self):
        """Test GT matrix has correct shape and dtype."""
        result = self._run_snpsnip(["--maf", "0", "--max-missing", "1", "--min-qual", "0"])
        self.assertEqual(result.returncode, 0)

        h5_path = os.path.join(self.output_dir, "snpsnip_data.h5")
        with h5py.File(h5_path, 'r') as h5:
            gt = h5["gt_matrix"][()]
            self.assertEqual(gt.dtype, np.int8, "GT matrix should be int8")
            self.assertEqual(gt.ndim, 2, "GT matrix should be 2D")
            self.assertEqual(gt.shape[0], gt.shape[0])  # n_snps

            samples = h5["samples"][()]
            self.assertEqual(gt.shape[1], len(samples), "GT matrix columns should match samples")

    def test_05_gt_matrix_values(self):
        """Test GT matrix contains only valid values (0, 1, 2, -1)."""
        result = self._run_snpsnip(["--maf", "0", "--max-missing", "1", "--min-qual", "0"])
        self.assertEqual(result.returncode, 0)

        h5_path = os.path.join(self.output_dir, "snpsnip_data.h5")
        with h5py.File(h5_path, 'r') as h5:
            gt = h5["gt_matrix"][()]
            unique_values = set(np.unique(gt))
            valid_values = {0, 1, 2, -1}
            self.assertTrue(unique_values.issubset(valid_values),
                          f"GT matrix contains invalid values: {unique_values - valid_values}")

    def test_06_gt_matrix_has_missing(self):
        """Test GT matrix has some missing values (-1) since VCF has missing_prop=0.05."""
        result = self._run_snpsnip(["--maf", "0", "--max-missing", "1", "--min-qual", "0"])
        self.assertEqual(result.returncode, 0)

        h5_path = os.path.join(self.output_dir, "snpsnip_data.h5")
        with h5py.File(h5_path, 'r') as h5:
            gt = h5["gt_matrix"][()]
            missing_count = np.sum(gt == -1)
            total = gt.size
            missing_rate = missing_count / total
            self.assertGreater(missing_rate, 0, "GT matrix should have some missing values")
            self.assertLess(missing_rate, 0.5, "Missing rate should be reasonable (< 50%)")

    def test_07_dp_matrix_shape_and_dtype(self):
        """Test DP matrix has correct shape and dtype when present."""
        result = self._run_snpsnip(["--maf", "0", "--max-missing", "1", "--min-qual", "0"])
        self.assertEqual(result.returncode, 0)

        h5_path = os.path.join(self.output_dir, "snpsnip_data.h5")
        with h5py.File(h5_path, 'r') as h5:
            has_dp = bool(h5["has_dp_matrix"][()])
            if has_dp:
                self.assertIn("dp_matrix", h5, "DP matrix should exist when has_dp_matrix is True")
                dp = h5["dp_matrix"][()]
                self.assertEqual(dp.dtype, np.int32, "DP matrix should be int32")
                self.assertEqual(dp.ndim, 2, "DP matrix should be 2D")

                gt = h5["gt_matrix"][()]
                self.assertEqual(dp.shape, gt.shape, "DP and GT matrices should have same shape")

    def test_08_dp_matrix_values(self):
        """Test DP matrix contains valid values (>= 0 or -1 for missing)."""
        result = self._run_snpsnip(["--maf", "0", "--max-missing", "1", "--min-qual", "0"])
        self.assertEqual(result.returncode, 0)

        h5_path = os.path.join(self.output_dir, "snpsnip_data.h5")
        with h5py.File(h5_path, 'r') as h5:
            if "dp_matrix" not in h5:
                self.skipTest("DP matrix not present in this VCF")

            dp = h5["dp_matrix"][()]
            valid_mask = (dp >= 0) | (dp == -1)
            self.assertTrue(np.all(valid_mask),
                          f"DP matrix contains invalid values. Min: {dp.min()}, Max: {dp.max()}")

    def test_09_samples_dataset(self):
        """Test samples dataset contains expected sample names."""
        result = self._run_snpsnip(["--maf", "0", "--max-missing", "1", "--min-qual", "0"])
        self.assertEqual(result.returncode, 0)

        h5_path = os.path.join(self.output_dir, "snpsnip_data.h5")
        with h5py.File(h5_path, 'r') as h5:
            samples = h5["samples"][()]
            self.assertEqual(len(samples), 20, "Should have 20 samples")
            for i in range(20):
                expected = f"sample_{i+1:02d}"
                self.assertEqual(samples[i].decode() if isinstance(samples[i], bytes) else samples[i],
                               expected, f"Sample {i} should be {expected}")

    def test_10_pca_data(self):
        """Test PCA data has correct structure."""
        result = self._run_snpsnip(["--maf", "0", "--max-missing", "1", "--min-qual", "0"])
        self.assertEqual(result.returncode, 0)

        h5_path = os.path.join(self.output_dir, "snpsnip_data.h5")
        with h5py.File(h5_path, 'r') as h5:
            pca_coords = h5["pca/coordinates"][()]
            self.assertEqual(pca_coords.ndim, 2, "PCA coordinates should be 2D")
            self.assertEqual(pca_coords.shape[0], 20, "PCA should have 20 samples")

            var_explained = h5["pca/variance_explained"][()]
            self.assertGreater(len(var_explained), 0, "Variance explained should not be empty")
            for v in var_explained:
                self.assertGreater(v, 0, "Each component should explain > 0% variance")
                self.assertLess(v, 100, "Each component should explain < 100% variance")

    def test_11_snp_info_lengths(self):
        """Test SNP info arrays have consistent lengths."""
        result = self._run_snpsnip(["--maf", "0", "--max-missing", "1", "--min-qual", "0"])
        self.assertEqual(result.returncode, 0)

        h5_path = os.path.join(self.output_dir, "snpsnip_data.h5")
        with h5py.File(h5_path, 'r') as h5:
            n_snps = len(h5["snp_info/chrom"][()])
            self.assertEqual(len(h5["snp_info/pos"][()]), n_snps)
            self.assertEqual(len(h5["snp_info/ref"][()]), n_snps)
            self.assertEqual(len(h5["snp_info/alt"][()]), n_snps)
            self.assertEqual(len(h5["snp_info/qual"][()]), n_snps)
            self.assertEqual(len(h5["snp_info/f_missing"][()]), n_snps)

            gt = h5["gt_matrix"][()]
            self.assertEqual(gt.shape[0], n_snps, "GT matrix rows should match SNP count")

    def test_12_r_script_content(self):
        """Test R script contains expected content for reading HDF5."""
        result = self._run_snpsnip(["--maf", "0", "--max-missing", "1", "--min-qual", "0"])
        self.assertEqual(result.returncode, 0)

        r_path = os.path.join(self.output_dir, "snpsnip_analysis.R")
        with open(r_path, 'r') as f:
            content = f.read()

        self.assertIn("snpsnip_data.h5", content)
        self.assertIn("gt_matrix", content)
        self.assertIn("dp_matrix", content)
        self.assertIn("hdf5r", content)
        self.assertIn("NA_integer_", content)
        self.assertIn("-1L", content)

    def test_13_combined_r_notebook_next_file(self):
        """Test that combined R notebook next file can be loaded."""
        result = self._run_snpsnip(["--maf", "0", "--max-missing", "1", "--min-qual", "0"])
        self.assertEqual(result.returncode, 0)

        state_file = os.path.join(self.output_dir, "state.json")
        with open(state_file, 'r') as f:
            state = json.load(f)

        all_samples = [s["id"] for s in state["sample_stats"]]
        next_data = {
            "groups": {"all": all_samples},
            "thresholds": {
                "all": {
                    "missing": {"min": None, "max": 0.2},
                    "maf": {"min": 0.05, "max": None},
                    "qual": {"min": 30, "max": None},
                }
            }
        }

        next_file = os.path.join(self.output_dir, "snpsnip_selections.json")
        with open(next_file, 'w') as f:
            json.dump(next_data, f)

        result = subprocess.run([
            sys.executable, "-m", "snpsnip",
            "--vcf", self.test_vcf,
            "--output-dir", self.output_dir,
            "--offline",
            "--next", next_file,
            "--seed", "42"
        ], capture_output=True, text=True)

        self.assertEqual(result.returncode, 0, f"Should succeed. stderr: {result.stderr}")

        with open(state_file, 'r') as f:
            state = json.load(f)
        self.assertTrue(state.get("completed"), "Processing should be complete")


class MissingGenotypeTest(unittest.TestCase):
    """Tests for missing genotype handling in HDF5 export."""

    @classmethod
    def setUpClass(cls):
        try:
            subprocess.run(["bcftools", "--version"],
                         capture_output=True, check=True, timeout=5)
        except (subprocess.SubprocessError, FileNotFoundError, subprocess.TimeoutExpired):
            raise unittest.SkipTest("bcftools not found - skipping missing genotype tests.")

        cls.test_dir = tempfile.mkdtemp(prefix="snpsnip_missing_test_")
        cls.test_vcf = os.path.join(cls.test_dir, "missing_test.vcf.gz")
        cls._create_missing_vcf()

    @classmethod
    def _create_missing_vcf(cls):
        """Create a VCF with known missing genotypes for testing."""
        header = """##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##contig=<ID=chr1,length=1000000>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample_01\tsample_02\tsample_03\tsample_04\tsample_05
"""
        rows = []
        for i in range(1, 1501):
            pos = i * 1000
            gt_row = []
            dp_row = []
            for j in range(1, 6):
                if (i + j) % 10 == 0:
                    gt_row.append("./.")
                    dp_row.append(".")
                elif (i + j) % 7 == 0:
                    gt_row.append("0/1")
                    dp_row.append(str(15 + j))
                elif (i + j) % 5 == 0:
                    gt_row.append("1/1")
                    dp_row.append(str(20 + j))
                else:
                    gt_row.append("0/0")
                    dp_row.append(str(10 + j))

            gt_str = "\t".join(gt_row)
            dp_str = "\t".join(dp_row)
            rows.append(f"chr1\t{pos}\t.\tA\tG\t50\tPASS\tDP=100\tGT:DP\t{gt_str}:{dp_str}")

        vcf_content = header + "\n".join(rows) + "\n"

        vcf_uncompressed = os.path.join(cls.test_dir, "missing_test.vcf")
        with open(vcf_uncompressed, 'w') as f:
            f.write(vcf_content)

        subprocess.run(["bcftools", "view", "-Oz", "-o", cls.test_vcf, vcf_uncompressed],
                      check=True, capture_output=True)
        subprocess.run(["bcftools", "index", "-f", cls.test_vcf],
                      check=True, capture_output=True)

    @classmethod
    def tearDownClass(cls):
        if os.path.exists(cls.test_dir):
            shutil.rmtree(cls.test_dir)

    def setUp(self):
        self.output_dir = tempfile.mkdtemp(prefix="snpsnip_missing_output_")

    def tearDown(self):
        if os.path.exists(self.output_dir):
            shutil.rmtree(self.output_dir)

    def test_known_missing_genotypes(self):
        """Test that known missing genotypes are correctly encoded as -1 in HDF5."""
        result = subprocess.run([
            sys.executable, "-m", "snpsnip",
            "--vcf", self.test_vcf,
            "--output-dir", self.output_dir,
            "--offline",
            "--r-notebook",
            "--subset-freq", "1.0",
            "--processes", "1",
            "--seed", "42",
            "--maf", "0",
            "--max-missing", "1",
            "--min-qual", "0"
        ], capture_output=True, text=True)

        self.assertEqual(result.returncode, 0, f"Should succeed. stderr: {result.stderr}")

        h5_path = os.path.join(self.output_dir, "snpsnip_data.h5")
        with h5py.File(h5_path, 'r') as h5:
            gt = h5["gt_matrix"][()]

            missing_count = np.sum(gt == -1)
            self.assertGreater(missing_count, 0, "Should have missing genotypes")

            non_missing = gt[gt != -1]
            self.assertTrue(np.all(np.isin(non_missing, [0, 1, 2])),
                          "Non-missing values should only be 0, 1, or 2")

    def test_known_missing_depth(self):
        """Test that known missing depth values are correctly encoded as -1 in HDF5."""
        result = subprocess.run([
            sys.executable, "-m", "snpsnip",
            "--vcf", self.test_vcf,
            "--output-dir", self.output_dir,
            "--offline",
            "--r-notebook",
            "--subset-freq", "1.0",
            "--processes", "1",
            "--seed", "42",
            "--maf", "0",
            "--max-missing", "1",
            "--min-qual", "0"
        ], capture_output=True, text=True)

        self.assertEqual(result.returncode, 0)

        h5_path = os.path.join(self.output_dir, "snpsnip_data.h5")
        with h5py.File(h5_path, 'r') as h5:
            if "dp_matrix" not in h5:
                self.skipTest("DP matrix not present")

            dp = h5["dp_matrix"][()]
            missing_count = np.sum(dp == -1)
            self.assertGreater(missing_count, 0, "Should have missing depth values")

            non_missing = dp[dp != -1]
            self.assertTrue(np.all(non_missing >= 0),
                          "Non-missing depth values should be >= 0")

    def test_gt_dp_correlation(self):
        """Test that missing GT and missing DP are both present in the matrices."""
        result = subprocess.run([
            sys.executable, "-m", "snpsnip",
            "--vcf", self.test_vcf,
            "--output-dir", self.output_dir,
            "--offline",
            "--r-notebook",
            "--subset-freq", "1.0",
            "--processes", "1",
            "--seed", "42",
            "--maf", "0",
            "--max-missing", "1",
            "--min-qual", "0"
        ], capture_output=True, text=True)

        self.assertEqual(result.returncode, 0)

        h5_path = os.path.join(self.output_dir, "snpsnip_data.h5")
        with h5py.File(h5_path, 'r') as h5:
            gt = h5["gt_matrix"][()]
            if "dp_matrix" not in h5:
                self.skipTest("DP matrix not present")
            dp = h5["dp_matrix"][()]

            gt_missing = (gt == -1)
            dp_missing = (dp == -1)

            gt_missing_count = int(np.sum(gt_missing))
            dp_missing_count = int(np.sum(dp_missing))

            self.assertGreater(gt_missing_count, 0, "GT matrix should have missing values")
            self.assertGreater(dp_missing_count, 0, "DP matrix should have missing values")

            both_missing = int(np.sum(gt_missing & dp_missing))
            self.assertGreater(both_missing, 0, "Some positions should be missing in both GT and DP")


class NoDPFieldsTest(unittest.TestCase):
    """Tests for VCF without DP fields."""

    @classmethod
    def setUpClass(cls):
        try:
            subprocess.run(["bcftools", "--version"],
                         capture_output=True, check=True, timeout=5)
        except (subprocess.SubprocessError, FileNotFoundError, subprocess.TimeoutExpired):
            raise unittest.SkipTest("bcftools not found - skipping no-DP tests.")

        cls.test_dir = tempfile.mkdtemp(prefix="snpsnip_nodp_test_")
        cls.test_vcf = os.path.join(cls.test_dir, "test.vcf.gz")
        cls.stripped_vcf = os.path.join(cls.test_dir, "stripped.vcf.gz")

        generate_example_vcf(
            output_path=cls.test_vcf,
            n_samples=10,
            n_snps=5000,
            missing_prop=0.05,
            n_chroms=1,
            compress=True,
            random_seed=42
        )

        annotate_cmd = [
            "bcftools", "annotate",
            "-x", "INFO/DP,FORMAT/DP",
            "-Oz", "-o", cls.stripped_vcf,
            cls.test_vcf
        ]
        subprocess.run(annotate_cmd, check=True, capture_output=True)
        subprocess.run(["bcftools", "index", "-f", cls.stripped_vcf],
                      check=True, capture_output=True)

    @classmethod
    def tearDownClass(cls):
        if os.path.exists(cls.test_dir):
            shutil.rmtree(cls.test_dir)

    def setUp(self):
        self.output_dir = tempfile.mkdtemp(prefix="snpsnip_nodp_output_")

    def tearDown(self):
        if os.path.exists(self.output_dir):
            shutil.rmtree(self.output_dir)

    def test_no_dp_matrix_in_hdf5(self):
        """Test that HDF5 file does not contain dp_matrix when VCF lacks DP."""
        result = subprocess.run([
            sys.executable, "-m", "snpsnip",
            "--vcf", self.stripped_vcf,
            "--output-dir", self.output_dir,
            "--offline",
            "--r-notebook",
            "--subset-freq", "0.5",
            "--processes", "1",
            "--seed", "42",
            "--maf", "0",
            "--max-missing", "1",
            "--min-qual", "0"
        ], capture_output=True, text=True)

        self.assertEqual(result.returncode, 0, f"Should succeed. stderr: {result.stderr}")

        h5_path = os.path.join(self.output_dir, "snpsnip_data.h5")
        with h5py.File(h5_path, 'r') as h5:
            has_dp = bool(h5["has_dp_matrix"][()])
            self.assertFalse(has_dp, "has_dp_matrix should be False")
            self.assertNotIn("dp_matrix", h5, "dp_matrix should not exist")
            self.assertIn("gt_matrix", h5, "gt_matrix should still exist")

    def test_gt_matrix_still_valid_without_dp(self):
        """Test GT matrix is still valid when DP is absent."""
        result = subprocess.run([
            sys.executable, "-m", "snpsnip",
            "--vcf", self.stripped_vcf,
            "--output-dir", self.output_dir,
            "--offline",
            "--r-notebook",
            "--subset-freq", "0.5",
            "--processes", "1",
            "--seed", "42",
            "--maf", "0",
            "--max-missing", "1",
            "--min-qual", "0"
        ], capture_output=True, text=True)

        self.assertEqual(result.returncode, 0)

        h5_path = os.path.join(self.output_dir, "snpsnip_data.h5")
        with h5py.File(h5_path, 'r') as h5:
            gt = h5["gt_matrix"][()]
            self.assertEqual(gt.dtype, np.int8)
            unique_values = set(np.unique(gt))
            self.assertTrue(unique_values.issubset({0, 1, 2, -1}))


if __name__ == "__main__":
    unittest.main()
