#!/usr/bin/env python3
"""
Integration tests for SNPSnip

These tests cover the full workflow including:
- VCF subset creation
- Sample filtering with and without sample lists
- Variant filtering
- Multiple groups
- Offline mode operation
"""

import json
import os
import shutil
import subprocess
import tempfile
import sys
import unittest
from pathlib import Path


class SNPSnipIntegrationTest(unittest.TestCase):
    """Integration tests for SNPSnip in offline mode"""

    @classmethod
    def setUpClass(cls):
        """Create test VCF file once for all tests"""
        # Check if bcftools is available
        try:
            subprocess.run(["bcftools", "--version"],
                         capture_output=True, check=True, timeout=5)
        except (subprocess.SubprocessError, FileNotFoundError, subprocess.TimeoutExpired):
            raise unittest.SkipTest("bcftools not found - skipping integration tests. "
                                   "Please install bcftools to run integration tests.")

        cls.test_dir = tempfile.mkdtemp(prefix="snpsnip_test_")
        cls.test_vcf = os.path.join(cls.test_dir, "test.vcf.gz")

        # Create a realistic test VCF with sufficient SNPs (10,000 variants)
        cls._create_test_vcf()

    @classmethod
    def tearDownClass(cls):
        """Clean up test directory"""
        if os.path.exists(cls.test_dir):
            shutil.rmtree(cls.test_dir)

    @classmethod
    def _create_test_vcf(cls):
        """Create a test VCF file with 10,000 variants and 20 samples"""
        vcf_content = """##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=chr1,length=10000000>
##contig=<ID=chr2,length=10000000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"""

        # Add 20 sample names
        samples = [f"sample_{i:02d}" for i in range(1, 21)]
        vcf_content += "\t" + "\t".join(samples) + "\n"

        # Generate 10,000 variants (5,000 per chromosome)
        variant_lines = []
        for chrom in ["chr1", "chr2"]:
            for i in range(5000):
                pos = (i + 1) * 1000
                variant_id = f"{chrom}_{pos}"
                qual = 30 + (i % 50)  # Quality varies from 30-80
                dp = 100 + (i % 100)  # Depth varies from 100-200
                af = 0.1 + (i % 40) / 100  # AF varies from 0.1 to 0.5

                # Generate genotypes with varying patterns
                genotypes = []
                for j in range(20):
                    # Create some variation in genotypes
                    if (i + j) % 7 == 0:
                        gt = "0/0"  # Homozygous reference
                    elif (i + j) % 5 == 0:
                        gt = "1/1"  # Homozygous alternate
                    elif (i + j) % 11 == 0:
                        gt = "./."  # Missing
                    else:
                        gt = "0/1"  # Heterozygous

                    sample_dp = 10 + ((i + j) % 30)
                    genotypes.append(f"{gt}:{sample_dp}")

                line = f"{chrom}\t{pos}\t{variant_id}\tA\tT\t{qual}\tPASS\tDP={dp};AF={af:.3f}\tGT:DP\t"
                line += "\t".join(genotypes)
                variant_lines.append(line)

        vcf_content += "\n".join(variant_lines)

        # Write uncompressed VCF first
        uncompressed_vcf = cls.test_vcf.replace(".gz", "")
        with open(uncompressed_vcf, 'w') as f:
            f.write(vcf_content)

        # Compress and index with bcftools
        subprocess.run(["bcftools", "view", "-Oz", "-o", cls.test_vcf, uncompressed_vcf], check=True)
        subprocess.run(["bcftools", "index", "-f", cls.test_vcf], check=True)

        # Clean up uncompressed file
        os.remove(uncompressed_vcf)

    def setUp(self):
        """Set up test output directory for each test"""
        self.output_dir = tempfile.mkdtemp(prefix="snpsnip_output_")

    def tearDown(self):
        """Clean up test output directory"""
        if os.path.exists(self.output_dir):
            shutil.rmtree(self.output_dir)

    def _run_snpsnip(self, extra_args=None):
        """Helper to run snpsnip with common arguments"""
        cmd = [
            sys.executable, "-m", "snpsnip",
            "--vcf", self.test_vcf,
            "--output-dir", self.output_dir,
            "--offline",
            "--subset-freq", "0.5",  # Use 50% to ensure we get enough SNPs
            "--processes", "2"
        ]
        if extra_args:
            cmd.extend(extra_args)

        result = subprocess.run(cmd, capture_output=True, text=True)
        return result

    def _check_state_file(self, expected_stage):
        """Check that state file exists and has expected stage"""
        state_file = os.path.join(self.output_dir, "state.json")
        self.assertTrue(os.path.exists(state_file), "State file should exist")

        with open(state_file, 'r') as f:
            state = json.load(f)

        self.assertEqual(state["stage"], expected_stage, f"Stage should be {expected_stage}")
        return state

    def test_01_initial_processing_no_sample_list(self):
        """Test initial processing without sample list"""
        result = self._run_snpsnip()

        # Should succeed
        self.assertEqual(result.returncode, 0, f"Should succeed. stderr: {result.stderr}")

        # Check state
        state = self._check_state_file("ready_for_sample_filtering")

        # Verify subset VCF was created
        self.assertIsNotNone(state.get("subset_vcf"), "Subset VCF should be created")
        self.assertTrue(os.path.exists(state["subset_vcf"]), "Subset VCF file should exist")

        # Verify sample stats were computed
        self.assertIsNotNone(state.get("sample_stats"), "Sample stats should exist")
        self.assertEqual(len(state["sample_stats"]), 20, "Should have stats for 20 samples")

        # Verify PCA was computed
        self.assertIsNotNone(state.get("pca"), "PCA should exist")
        self.assertIn("samples", state["pca"], "PCA should have samples")

        # Verify HTML was generated
        html_file = os.path.join(self.output_dir, "sample_filtering.html")
        self.assertTrue(os.path.exists(html_file), "Sample filtering HTML should exist")

        # Check that subset VCF has enough SNPs (should have ~5000 with 0.5 subset freq)
        subset_vcf = state["subset_vcf"]
        count_result = subprocess.run(
            ["bcftools", "index", "--stats", subset_vcf],
            capture_output=True, text=True, check=True
        )
        variant_count = sum(
            int(line.split('\t')[2]) for line in count_result.stdout.strip().split('\n')
            if line and not line.startswith('#') and len(line.split('\t')) >= 3
        )
        self.assertGreaterEqual(variant_count, 5000, "Subset should have at least 5000 SNPs")
        print(f"Test 1: Subset VCF created with {variant_count} SNPs")

    def test_02_initial_processing_with_sample_list(self):
        """Test initial processing with sample list"""
        # Create sample list with 10 valid samples and 2 invalid ones
        sample_list_file = os.path.join(self.output_dir, "samples.txt")
        with open(sample_list_file, 'w') as f:
            for i in range(1, 11):
                f.write(f"sample_{i:02d}\n")
            f.write("invalid_sample_1\n")
            f.write("invalid_sample_2\n")

        result = self._run_snpsnip(["--sample-list", sample_list_file])

        # Should succeed
        self.assertEqual(result.returncode, 0, f"Should succeed. stderr: {result.stderr}")

        # Check warnings about invalid samples
        self.assertIn("invalid_sample", result.stderr, "Should warn about invalid samples")

        # Check state
        state = self._check_state_file("ready_for_sample_filtering")

        # Verify subset VCF was created with only selected samples
        subset_vcf = state["subset_vcf"]
        samples_result = subprocess.run(
            ["bcftools", "query", "-l", subset_vcf],
            capture_output=True, text=True, check=True
        )
        subset_samples = samples_result.stdout.strip().split('\n')
        self.assertEqual(len(subset_samples), 10, "Should have 10 samples in subset")

        print(f"Test 2: Sample list filtering works correctly ({len(subset_samples)} samples)")

    def test_03_full_workflow_no_groups(self):
        """Test full workflow without sample groups"""
        # Stage 1: Initial processing
        result = self._run_snpsnip()
        self.assertEqual(result.returncode, 0)

        # Create a next file simulating user selection (all samples in one group)
        state_file = os.path.join(self.output_dir, "state.json")
        with open(state_file, 'r') as f:
            state = json.load(f)

        all_samples = [s["id"] for s in state["sample_stats"]]
        next_data = {
            "groups": {
                "all_samples": all_samples
            }
        }

        next_file = os.path.join(self.output_dir, "next_sample_groups.json")
        with open(next_file, 'w') as f:
            json.dump(next_data, f)

        # Stage 2: Process variant stats
        result = self._run_snpsnip(["--next", next_file])
        self.assertEqual(result.returncode, 0, f"Stage 2 should succeed. stderr: {result.stderr}")

        state = self._check_state_file("ready_for_variant_filtering")

        # Verify variant stats were computed
        self.assertIn("all_samples", state["variant_stats"], "Should have variant stats for group")
        group_stats = state["variant_stats"]["all_samples"]
        self.assertIn("qual", group_stats, "Should have quality stats")
        self.assertIn("depth", group_stats, "Should have depth stats")
        self.assertIn("af", group_stats, "Should have AF stats")

        # Verify HTML was generated
        html_file = os.path.join(self.output_dir, "variant_filtering.html")
        self.assertTrue(os.path.exists(html_file), "Variant filtering HTML should exist")

        # Create next file with variant thresholds
        next_data = {
            "thresholds": {
                "all_samples": {
                    "qual": {"min": 30, "max": None},
                    "depth": {"min": 50, "max": None},
                    "af": {"min": 0.1, "max": None}
                }
            }
        }

        next_file = os.path.join(self.output_dir, "next_variant_filters.json")
        with open(next_file, 'w') as f:
            json.dump(next_data, f)

        # Stage 3: Apply final filters
        result = self._run_snpsnip(["--next", next_file])
        self.assertEqual(result.returncode, 0, f"Stage 3 should succeed. stderr: {result.stderr}")

        state = self._check_state_file("ready_for_final")
        self.assertTrue(state.get("completed"), "Processing should be complete")

        # Verify output VCF was created
        output_files = state.get("output_files", [])
        self.assertEqual(len(output_files), 1, "Should have one output file")

        output_vcf = output_files[0]
        self.assertTrue(os.path.exists(output_vcf), "Output VCF should exist")
        self.assertTrue(os.path.exists(output_vcf + ".csi"), "Output VCF index should exist")

        # Verify output VCF has correct samples and variants
        samples_result = subprocess.run(
            ["bcftools", "query", "-l", output_vcf],
            capture_output=True, text=True, check=True
        )
        output_samples = samples_result.stdout.strip().split('\n')
        self.assertEqual(len(output_samples), 20, "Should have all 20 samples")

        print(f"Test 3: Full workflow completed successfully")

    def test_04_full_workflow_with_multiple_groups(self):
        """Test full workflow with multiple sample groups"""
        # Stage 1: Initial processing
        result = self._run_snpsnip()
        self.assertEqual(result.returncode, 0)

        # Create next file with multiple groups
        state_file = os.path.join(self.output_dir, "state.json")
        with open(state_file, 'r') as f:
            state = json.load(f)

        all_samples = [s["id"] for s in state["sample_stats"]]
        next_data = {
            "groups": {
                "group_A": all_samples[:10],  # First 10 samples
                "group_B": all_samples[10:],  # Last 10 samples
            }
        }

        next_file = os.path.join(self.output_dir, "next_sample_groups.json")
        with open(next_file, 'w') as f:
            json.dump(next_data, f)

        # Stage 2: Process variant stats
        result = self._run_snpsnip(["--next", next_file])
        self.assertEqual(result.returncode, 0)

        state = self._check_state_file("ready_for_variant_filtering")

        # Verify variant stats for both groups
        self.assertIn("group_A", state["variant_stats"], "Should have stats for group_A")
        self.assertIn("group_B", state["variant_stats"], "Should have stats for group_B")

        # Create next file with different thresholds for each group
        next_data = {
            "thresholds": {
                "group_A": {
                    "qual": {"min": 35, "max": None},
                    "depth": {"min": 60, "max": None},
                    "af": {"min": 0.15, "max": None}
                },
                "group_B": {
                    "qual": {"min": 30, "max": None},
                    "depth": {"min": 50, "max": None},
                    "af": {"min": 0.1, "max": 0.4}
                }
            }
        }

        next_file = os.path.join(self.output_dir, "next_variant_filters.json")
        with open(next_file, 'w') as f:
            json.dump(next_data, f)

        # Stage 3: Apply final filters
        result = self._run_snpsnip(["--next", next_file])
        self.assertEqual(result.returncode, 0)

        state = self._check_state_file("ready_for_final")
        self.assertTrue(state.get("completed"), "Processing should be complete")

        # Verify two output VCFs were created
        output_files = state.get("output_files", [])
        self.assertEqual(len(output_files), 2, "Should have two output files")

        # Verify both output files exist and have correct samples
        for output_vcf in output_files:
            self.assertTrue(os.path.exists(output_vcf), f"Output VCF {output_vcf} should exist")

            # Check which group this is
            if "group_A" in output_vcf:
                expected_samples = 10
                group_name = "group_A"
            else:
                expected_samples = 10
                group_name = "group_B"

            samples_result = subprocess.run(
                ["bcftools", "query", "-l", output_vcf],
                capture_output=True, text=True, check=True
            )
            output_samples = samples_result.stdout.strip().split('\n')
            self.assertEqual(len(output_samples), expected_samples,
                           f"{group_name} should have {expected_samples} samples")

        print(f"Test 4: Multiple groups workflow completed successfully")

    def test_05_insufficient_snps_error(self):
        """Test that insufficient SNPs causes appropriate error"""
        # Use very low subset frequency to trigger error
        cmd = [
            sys.executable, "-m", "snpsnip",
            "--vcf", self.test_vcf,
            "--output-dir", self.output_dir,
            "--offline",
            "--subset-freq", "0.001",  # Very low to get < 5000 SNPs
            "--processes", "2"
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)

        # Should fail with error
        self.assertNotEqual(result.returncode, 0, "Should fail with insufficient SNPs")

        # Check error message
        self.assertIn("minimum required: 5000", result.stderr,
                     "Should mention minimum requirement")
        self.assertIn("Increase --subset-freq", result.stderr,
                     "Should suggest increasing subset-freq")

        print(f"Test 5: Insufficient SNPs error handled correctly")

    def test_06_with_initial_filters(self):
        """Test with initial MAF and quality filters"""
        result = self._run_snpsnip([
            "--maf", "0.2",
            "--max-missing", "0.2",
            "--min-qual", "40"
        ])

        # Should succeed
        self.assertEqual(result.returncode, 0, f"Should succeed. stderr: {result.stderr}")

        # Check state
        state = self._check_state_file("ready_for_sample_filtering")

        # Verify subset was created (may have fewer SNPs due to filters)
        self.assertIsNotNone(state.get("subset_vcf"))

        print(f"Test 6: Initial filters applied successfully")

    def test_07_with_sample_list_and_groups(self):
        """Test workflow with both sample list and groups"""
        # Create sample list
        sample_list_file = os.path.join(self.output_dir, "samples.txt")
        with open(sample_list_file, 'w') as f:
            for i in range(1, 16):  # Select 15 samples
                f.write(f"sample_{i:02d}\n")

        # Stage 1: Initial processing with sample list
        result = self._run_snpsnip(["--sample-list", sample_list_file])
        self.assertEqual(result.returncode, 0)

        state_file = os.path.join(self.output_dir, "state.json")
        with open(state_file, 'r') as f:
            state = json.load(f)

        # Create groups (some samples may not be in sample list)
        next_data = {
            "groups": {
                "group_1": [f"sample_{i:02d}" for i in range(1, 8)],
                "group_2": [f"sample_{i:02d}" for i in range(8, 16)]
            }
        }

        next_file = os.path.join(self.output_dir, "next_sample_groups.json")
        with open(next_file, 'w') as f:
            json.dump(next_data, f)

        # Stage 2
        result = self._run_snpsnip(["--next", next_file, "--sample-list", sample_list_file])
        self.assertEqual(result.returncode, 0)

        # Stage 3 with filters
        next_data = {
            "thresholds": {
                "group_1": {"qual": {"min": 30, "max": None}},
                "group_2": {"qual": {"min": 30, "max": None}}
            }
        }

        next_file = os.path.join(self.output_dir, "next_variant_filters.json")
        with open(next_file, 'w') as f:
            json.dump(next_data, f)

        result = self._run_snpsnip([
            "--next", next_file,
            "--sample-list", sample_list_file
        ])
        self.assertEqual(result.returncode, 0)

        state_file = os.path.join(self.output_dir, "state.json")
        with open(state_file, 'r') as f:
            state = json.load(f)

        # Verify output files
        output_files = state.get("output_files", [])
        self.assertEqual(len(output_files), 2, "Should have two output files")

        # Verify samples in output match sample list
        for output_vcf in output_files:
            samples_result = subprocess.run(
                ["bcftools", "query", "-l", output_vcf],
                capture_output=True, text=True, check=True
            )
            output_samples = samples_result.stdout.strip().split('\n')
            # All samples should be from the sample list
            for sample in output_samples:
                self.assertIn(sample, [f"sample_{i:02d}" for i in range(1, 16)])

        print(f"Test 7: Sample list with groups completed successfully")


def run_tests():
    """Run all integration tests"""
    # Create test suite
    loader = unittest.TestLoader()
    suite = loader.loadTestsFromTestCase(SNPSnipIntegrationTest)

    # Run tests with verbose output
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)

    # Return exit code
    return 0 if result.wasSuccessful() else 1


if __name__ == "__main__":
    exit(run_tests())
