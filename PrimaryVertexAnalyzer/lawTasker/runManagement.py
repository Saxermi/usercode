#!/usr/bin/env python3
import os
import re
import time
import shutil
import subprocess
import luigi
import law
 
#####################################################################
# Task to set up the experimental run directory.
#####################################################################
class SetupRunDirectoryTask(law.Task):
    # Base directory where experimental_run_X folders will be created.
    base_dir = luigi.Parameter(default="/work/frejalom/ba")
    
    def output(self):
        # The output is a file that stores the chosen run directory.
        return law.LocalFileTarget("run_directory.txt")
 
    def run(self):
        # Loop until a free experimental_run_<n> folder is found.
        n = 1
        while True:
            run_dir = os.path.join(self.base_dir, f"experimental_run_{n}")
            if not os.path.exists(run_dir):
                os.makedirs(run_dir)
                break
            n += 1
        with self.output().open("w") as outf:
            outf.write(run_dir)
        self.logger.info(f"Created run directory: {run_dir}")
 
 
#####################################################################
# Task to extract hardcoded parameters from a given C++ file.
# It searches for the two methods and extracts numeric literals found in their bodies.
#####################################################################
class ExtractParametersTask(law.Task):
    # The C++ file to scan (adjust the default path as needed).
    cpp_file = law.Parameter(default="/t3home/frejalom/cmssw/CMSSW_14_2_0_pre4/src/RecoVertex/PrimaryVertexProducer/src/DAClusterizerInZ_vect.cc")
    
    def requires(self):
        # Depends on the run directory being set up.
        return SetupRunDirectoryTask()
 
    def output(self):
        # Write the parameter log into the run directory (name includes the run directory name).
        run_dir = open(self.input().path).read().strip()
        output_path = os.path.join(run_dir, f"parameter_log_{os.path.basename(run_dir)}.txt")
        return law.LocalFileTarget(output_path)
 
    def run(self):
        with open(self.cpp_file, "r") as f:
            content = f.read()
 
        methods = [
            "DAClusterizerInZ_vect::vertices_in_blocks",
            "DAClusterizerInZ_vect::vertices_no_blocks",
        ]
        results = []
        for func in methods:
            # A simple regex that finds the function body between { and }.
            # (This assumes the body is on one “level”; adjust if needed.)
            pattern = re.compile(rf"{re.escape(func)}\s*\(.*?\)\s*\{{(.*?)\}}", re.DOTALL)
            match = pattern.search(content)
            if match:
                body = match.group(1)
                # Extract numeric literals (integers and decimals).
                literals = re.findall(r"(?<![\w.])\d+(?:\.\d+)?(?![\w.])", body)
                results.append(f"{func}: " + ", ".join(literals))
            else:
                results.append(f"{func}: <not found>")
        with self.output().open("w") as outf:
            outf.write("\n".join(results))
        self.logger.info("Extracted parameters written.")
 
 
#####################################################################
# Dynamic job task that submits a job using your bash script via lawslurm.
# This replaces your original submit_job() function.
#####################################################################
class RunJobTask(law.Task):
    sample    = law.Parameter()      # e.g. "ttbar_new_1_3_split"
    overlap   = law.FloatParameter()   # e.g. 0 or 0.4
    blocksize = law.IntParameter()     # e.g. 256 or 512
    notify    = law.BoolParameter(default=False)
 
    def requires(self):
        # Need the run directory to build the job-specific output path.
        return SetupRunDirectoryTask()
 
    def output(self):
        run_dir = open(self.input().path).read().strip()
        # Create a subdirectory structure as in your original code.
        overlap_dir = f"n{int(abs(self.overlap)*10)*10}" if self.overlap < 0 else f"{int(self.overlap*10)*10}"
        job_dir = os.path.join(run_dir, self.sample, overlap_dir, str(self.blocksize))
        if not os.path.exists(job_dir):
            os.makedirs(job_dir)
        # Dummy output file indicating job completion.
        return law.LocalFileTarget(os.path.join(job_dir, "job_done.txt"))
 
    def run(self):
        # Choose the appropriate bash script.
        bash_script = "pvslurmmaster_notifyme_frejalom.bsh" if self.notify else "pvslurmmaster_frejalom.bsh"
        run_dir = open(self.input().path).read().strip()
        overlap_dir = f"n{int(abs(self.overlap)*10)*10}" if self.overlap < 0 else f"{int(self.overlap*10)*10}"
        job_dir = os.path.join(run_dir, self.sample, overlap_dir, str(self.blocksize))
        if not os.path.exists(job_dir):
            os.makedirs(job_dir)
        # Build the command line.
        cmd = [
            "bash", bash_script,
            "-d", self.sample,
            "-o", str(self.overlap),
            "-p", job_dir,
            "-b", str(self.blocksize),
            "-l", "pv",
        ]
        self.logger.info("Submitting job with command: " + " ".join(cmd))
        # Instead of calling sbatch, lawslurm (if configured) will submit this command as a job.
        ret = subprocess.run(cmd, capture_output=True, text=True)
        if ret.returncode != 0:
            self.logger.error("Job failed: " + ret.stderr)
            raise RuntimeError("Job submission failed")
        # Write dummy output to mark completion.
        with self.output().open("w") as outf:
            outf.write("Job completed successfully\n")
        self.logger.info(f"Job for {self.sample} (overlap {self.overlap}, blocksize {self.blocksize}) submitted.")
 
 
#####################################################################
# Task to run the ROOT plotting code.
# (Here a temporary Python file is written and executed.
# Replace the dummy plotting code with your full ROOT code as needed.)
#####################################################################
class PlottingTask(law.Task):
    def requires(self):
        return SetupRunDirectoryTask()
 
    def output(self):
        run_dir = open(self.input().path).read().strip()
        plot_dir = os.path.join(run_dir, "plot_output")
        if not os.path.exists(plot_dir):
            os.makedirs(plot_dir)
        return law.LocalFileTarget(os.path.join(plot_dir, "plot_done.txt"))
 
    def run(self):
        run_dir = open(self.input().path).read().strip()
        self.logger.info(f"Using run directory: {run_dir}")

        dest_folders = {}
        for ov in overlaps:
            for bs in blocksizes:
                dest_folder = os.path.join(run_dir, f"o{ov*100}bs{bs}")
                os.makedirs(dest_folder, exist_ok=True)
                dest_folders[(ov, bs)] = dest_folder
                self.logger.info(f"Created destination folder {dest_folder} for overlap {ov} and blocksize {bs}")

        # Walk through the run directory tree and copy CSV files from job directories.
        # Expected job directories are structured as: run_dir/sample/overlap_folder/blocksize/
        for root, dirs, files in os.walk(run_dir):
            for file in files:
                if file.endswith(".csv"):
                    # Compute the relative path from the run directory.
                    rel = os.path.relpath(root, run_dir)
                    parts = rel.split(os.sep)
                    # Expect at least three levels: sample, overlap folder, blocksize folder.
                    if len(parts) >= 3:
                        # The overlap folder is the third-last element and blocksize is second-last.
                        overlap_dir = parts[-3]
                        block_dir = parts[-2]
                        if (overlap_dir, block_dir) in dest_folders:
                            source_file = os.path.join(root, file)
                            dest_file = os.path.join(dest_folders[(overlap_dir, block_dir)], file)
                            shutil.copy(source_file, dest_file)
                            self.logger.debug(f"Copied {source_file} to {dest_file}")
        self.logger.info("CSV files copied for all overlap and blocksize combinations.")

        plot_file = "/t3home/frejalom/cmssw/CMSSW_14_2_0_pre4/src/usercode/download/plotter.py"
        self.logger.info(f"Running existing plotter script at {plot_file}")
        ret = subprocess.run(["python3", plot_file], capture_output=True, text=True)
        if ret.returncode != 0:
            self.logger.error("Plotting failed: " + ret.stderr)
            raise RuntimeError("Plotting failed")
        self.logger.info("Plotting output:\n" + ret.stdout)

        with self.output().open("w") as outf:
            outf.write("Plotting completed successfully\n")
        self.logger.info("Plotting completed.")
 
 
#####################################################################
# Task to archive .out files from a given source directory into an 'outarchive' folder.
#####################################################################
class ArchiveOutTask(law.Task):
    def requires(self):
        return SetupRunDirectoryTask()
 
    def output(self):
        run_dir = open(self.input().path).read().strip()
        archive_dir = os.path.join(run_dir, "outarchive")
        if not os.path.exists(archive_dir):
            os.makedirs(archive_dir)
        return law.LocalFileTarget(os.path.join(archive_dir, "archive_done.txt"))
 
    def run(self):
        source_dir = "/t3home/frejalom/cmssw/CMSSW_14_2_0_pre4/src/usercode/PrimaryVertexAnalyzer/test/jobs"
        run_dir = open(self.input().path).read().strip()
        archive_dir = os.path.join(run_dir, "outarchive")
        for fname in os.listdir(source_dir):
            if fname.endswith(".out"):
                shutil.copy(os.path.join(source_dir, fname), archive_dir)
        with self.output().open("w") as outf:
            outf.write("Archive completed successfully\n")
        self.logger.info("Archived .out files.")
 
 
#####################################################################
# The master (wrapper) task which ties everything together.
# It requires the setup, extraction, dynamic job submissions, plotting, and archiving tasks.
#####################################################################
class MasterWorkflow(law.WrapperTask):
    def workflow_requires(self):
        reqs = dict()
        # First set up the run directory and extract the parameters.
        reqs["setup"] = SetupRunDirectoryTask()
        reqs["extract"] = ExtractParametersTask()
        # Dynamically create a list of job tasks.
        samples = [
            "ttbar_new_1_3_split",
            # "ttbar_new_2_3_split",
            # "ttbar_new_3_3_split",
            # "ttbar_new_4_3_split",
            # "ttbar_new_5_3_split",
            # "ttbar_new_6_3_split",
            # "ttbar_new_7_3_split",
            # "ttbar_new_8_3_split",
            # "ttbar_new_9_3_split",
            # "ttbar_new_10_3_split",
            # "ttbar_new_11_3_split",
            # "ttbar_new_12_3_split",
            # "ttbar_new_13_3_split",
            # "ttbar_new_14_3_split",
            # "ttbar_new_15_3_split",
            # "ttbar_new_16_3_split",
            # "ttbar_new_17_3_split",
            # "ttbar_new_18_3_split",
            # "ttbar_new_19_3_split",
            # "ttbar_new_20_3_split",
            # "ttbar_new_21_3_split",
            # "ttbar_new_22_3_split",
            # "ttbar_new_23_3_split",
            # "ttbar_new_24_3_split",
            # "ttbar_new_25_3_split",
            # "ttbar_new_26_3_split",
            # "ttbar_new_27_3_split",
            # "ttbar_new_28_3_split",
            # "ttbar_new_29_3_split",
            # "ttbar_new_30_3_split",
        ]
        overlaps = [0, 0.4]
        blocksizes = [256]
        jobs = []
        for sample in samples:
            for overlap in overlaps:
                for blocksize in blocksizes:
                    jobs.append(RunJobTask(sample=sample, overlap=overlap, blocksize=blocksize, notify=False))
        reqs["jobs"] = jobs
        # Then add plotting and archive tasks.
        reqs["plotting"] = PlottingTask()
        reqs["archive"] = ArchiveOutTask()
        return reqs
 
if __name__ == "__main__":
    law.run()