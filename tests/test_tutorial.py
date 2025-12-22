#!/usr/bin/env python3
"""
Integration test for tutorial.md

This test validates that the tutorial documentation works by running all bash
code blocks in sequence, simulating a user following the tutorial.
"""

import json
import os
import re
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


class TutorialTest(unittest.TestCase):
    """Test that tutorial.md commands work when run in sequence"""

    test_passed = False

    @classmethod
    def setUpClass(cls):
        """Set up test environment"""
        # Check if bcftools is available
        try:
            subprocess.run(["bcftools", "--version"],
                         capture_output=True, check=True, timeout=5)
        except (subprocess.SubprocessError, FileNotFoundError, subprocess.TimeoutExpired):
            raise unittest.SkipTest("bcftools not found - skipping tutorial tests.")

        # Find tutorial.md and repo root
        cls.repo_root = Path(__file__).parent.parent
        cls.tutorial_path = cls.repo_root / "tutorial.md"
        if not cls.tutorial_path.exists():
            raise unittest.SkipTest(f"tutorial.md not found at {cls.tutorial_path}")

        # Create test directory (separate from repo)
        cls.test_dir = tempfile.mkdtemp(prefix="snpsnip_tutorial_test_")
        cls.original_cwd = os.getcwd()

    @classmethod
    def tearDownClass(cls):
        """Clean up test directory only on success"""
        os.chdir(cls.original_cwd)
        if cls.test_passed and os.path.exists(cls.test_dir):
            shutil.rmtree(cls.test_dir)
            print(f"Test passed - cleaned up {cls.test_dir}")
        elif os.path.exists(cls.test_dir):
            print(f"\n{'='*60}")
            print(f"Test failed - keeping test directory for inspection:")
            print(f"  {cls.test_dir}")
            print(f"{'='*60}")

    def test_tutorial_commands(self):
        """Run all tutorial bash commands in sequence"""
        # Extract all bash code blocks
        with open(self.tutorial_path, 'r') as f:
            content = f.read()

        # Find all bash code blocks
        pattern = r'```bash\n(.*?)\n```'
        matches = re.findall(pattern, content, re.DOTALL)

        self.assertGreater(len(matches), 0, "No bash code blocks found in tutorial.md")
        print(f"\nFound {len(matches)} bash code blocks in tutorial.md")

        # Create a shell script that runs all commands
        script_lines = [
            "#!/bin/bash",
            "set -e  # Exit on error",
            ""
        ]

        for i, code_block in enumerate(matches):
            # Replace 'pip install -e .' with path to repo
            if 'pip install -e .' in code_block:
                code_block = code_block.replace('pip install -e .', f'pip install -e {self.repo_root}')
                print(f"Modified pip install to use repo path: {self.repo_root}")

            script_lines.append(f"# Block {i+1} from tutorial")
            script_lines.append(code_block)
            script_lines.append("")

        # Write the script
        script_path = os.path.join(self.test_dir, "run_tutorial.sh")
        with open(script_path, 'w') as f:
            f.write('\n'.join(script_lines))

        # Make it executable
        os.chmod(script_path, 0o755)

        print(f"\nRunning tutorial commands...")
        print(f"Test directory: {self.test_dir}")

        # Run the script
        result = subprocess.run(
            ["bash", script_path],
            capture_output=True,
            text=True,
            timeout=600,
            cwd=self.test_dir
        )

        # Print output
        print("\n" + "="*60)
        print("STDOUT:")
        print("="*60)
        print(result.stdout)

        if result.stderr:
            print("\n" + "="*60)
            print("STDERR:")
            print("="*60)
            print(result.stderr)

        # Check result
        self.assertEqual(result.returncode, 0,
                        f"Tutorial commands failed with exit code {result.returncode}")

        # Verify expected files were created
        expected_files = [
            "example.vcf.gz",
            "example.vcf.gz.csi",
            "snpsnip_sample_filters.json",
            "snpsnip_variant_filters.json",
        ]

        for filename in expected_files:
            filepath = os.path.join(self.test_dir, filename)
            self.assertTrue(os.path.exists(filepath),
                          f"Expected file not created: {filename}")

        # Verify offline mode outputs
        offline_dir = os.path.join(self.test_dir, "example_offline")
        if os.path.exists(offline_dir):
            # Check final VCF was created
            final_vcf = os.path.join(offline_dir, "all_samples_filtered.vcf.gz")
            self.assertTrue(os.path.exists(final_vcf),
                          f"Final filtered VCF not created: {final_vcf}")
            self.assertTrue(os.path.exists(final_vcf + ".csi"),
                          f"Final filtered VCF index not created: {final_vcf}.csi")

            # Verify final VCF has variants
            count_result = subprocess.run(
                ["bcftools", "view", "-H", final_vcf],
                capture_output=True,
                text=True,
                timeout=30
            )
            variant_count = len([l for l in count_result.stdout.strip().split('\n') if l])
            self.assertGreater(variant_count, 0, "Final VCF has no variants")
            print(f"Final VCF created with {variant_count} variants")

            # Check state.json exists and has expected structure
            state_file = os.path.join(offline_dir, "state.json")
            if os.path.exists(state_file):
                with open(state_file, 'r') as f:
                    state = json.load(f)

                # Verify state has expected keys
                self.assertTrue(state.get("completed", False),
                               f"Expected state completed=True, got '{state.get('completed')}'")
                self.assertIn("stage", state, "state.json missing 'stage' key")
                self.assertEqual(state["stage"], "completed",
                               f"Expected stage='completed', got '{state.get('stage')}'")

                print(f"State validated: stage={state['stage']}")

        print("\n" + "="*60)
        print("TUTORIAL VALIDATION SUCCESSFUL")
        print("="*60)

        # Mark test as passed for cleanup
        TutorialTest.test_passed = True


if __name__ == '__main__':
    unittest.main()
