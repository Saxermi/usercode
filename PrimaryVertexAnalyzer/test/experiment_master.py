import os
import subprocess
import itertools
import time
import logging
from datetime import datetime

# Global variable to control test mode
TEST_MODE = False  # Set to False to submit all jobs

# Get the current timestamp for the log filename
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
base_path = "experimental_run_18"   stop
log_filename = f"{base_path}_{timestamp}.log"
log_path = os.path.join("/work/msaxer", log_filename)

# Set up logging to file and console
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(log_path),
        logging.StreamHandler()
    ]
)

def create_directory(path):
    """
    Function to create a directory if it doesn't exist.
    
    Args:
        path (str): The path of the directory to create.
    """
    path = os.path.join("/work/msaxer", path)
    if not os.path.exists(path):
        os.makedirs(path)
        logging.info(f"Directory created: {path}")
        print(f"Directory created: {path}")
    else:
        logging.info(f"Directory already exists: {path}")
        print(f"Directory already exists: {path}")

def submit_job(sample, overlap, blocksize, iterating_blocksize=False, notify=False, queue="standard"):
    """
    Function to submit the slurm job using the corresponding bash script.

    Args:
        sample (str): The sample name, e.g., Subset_SToMuMu_01.txt.
        overlap (float): The overlap factor, e.g., 0.5.
        blocksize (int): The block size, e.g., 256.
        notify (bool): Whether to use the notification bash script.
        queue (str): The queue to submit the job to, e.g., "short", "standard", "long".
    """
    # Strip the '.txt' extension from the sample to create the directory name
    # sample_name = os.path.splitext(sample)[0]

    # Determine which script to use
    bash_script = "pvslurmmaster_notifyme.bsh" if notify else "pvslurmmaster.bsh"

    # Calculate the directory name based on overlap (convert to an integer percentage)
    #overlap_dir = int(abs(overlap) * 10) * 10  # Converts 0.1 to 10, 0.2 to 20, ..., 0.9 to 90
    overlap_dir = int(abs(overlap) * 100)  # for higher precision

    # Add "n" if overlap is negative
    overlap_dir_name = f"n{overlap_dir}" if overlap < 0 else f"{overlap_dir}"

    # Hardcoded base path for experimental runs
    base_path = "experimental_run_18"
    dir_create_path = "/work/msaxer/"

    # Create the full path based on the sample, overlap, and blocksize
    if iterating_blocksize:
        path = os.path.join(base_path, sample, overlap_dir_name, str(blocksize))
    else:
        path = os.path.join(base_path, sample, overlap_dir_name)

    # Create the directory if it doesn't exist
    create_directory(path)

    # Define the command to run
    cmd = [
        "sbatch"
    ]

    # Add the queue option if provided
    if queue:
        cmd.extend(["-p", queue])

    cmd.extend([
        bash_script,
        "-d",
        sample,
      #  "-n",
       # str(1000),
        "-o",
        str(overlap),
        "-p",
        path,  # Pass the path to store output in the correct directory
        "-b",
        str(blocksize),
        "-l",
        "pv",
    ])

    # Submit the job using subprocess and capture the output
    try:
        result = subprocess.run(
            cmd,
            check=True,
            universal_newlines=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        # Print and log the output from sbatch (which typically includes the job ID)
        logging.info(result.stdout)
        logging.info(f"Job for {sample} with overlap {overlap} and blocksize {blocksize} submitted successfully.")
        print(result.stdout)
        print(f"Job for {sample} with overlap {overlap} and blocksize {blocksize} submitted successfully.")
    except subprocess.CalledProcessError as e:
        logging.error(f"Failed to submit job for {sample} with overlap {overlap} and blocksize {blocksize}. Error: {e}")
        print(f"Failed to submit job for {sample} with overlap {overlap} and blocksize {blocksize}. Error: {e}")

def main():
    # Samples to be used (from TTbar_01 to TTbar_15)
    # samples = [f"TTbar_{str(i).zfill(2)}" for i in range(1, 16)]
    
    #subsets = [
     #   "Subset_SToMuMu_01", "Subset_ZMM_03", "Subset_SToMuMu_02", 
      #  "Subset_TTbar_01", "Subset_TTbar_02", "Subset_HiggsGluonFusion_01", 
      #  "Subset_HiggsGluonFusion_02", "Subset_ZMM_01", "Subset_ZMM_02"
    #]
    subsets = [
        
        "Subset_TTbar_01", "Subset_TTbar_02"
        
    ]
    #subsets = [
     #   "Subset_ZMM_02"
    #]
    #subsets = ["Subset_SToMuMu_01"]
    # Overlap values from 0.0 to 0.9 in 0.1 increments (include negative values if needed)
    overlaps = [0,0.3,0.35,0.4]  # Add negative overlaps here

    # Block sizes to iterate over
    blocksizes = [512]
    #blocksizes = [ 512]

    # If in test mode, only submit two jobs, one with notify and one without
    if TEST_MODE:
        logging.info("Running in test mode...")
        print("Running in test mode...")
        submit_job(subsets[0], overlaps[0], blocksize=blocksizes[0], iterating_blocksize=True, notify=False)
        time.sleep(10)
        submit_job(subsets[0], overlaps[0], blocksize=blocksizes[0], iterating_blocksize=True, notify=True)
    else:
        # Iterate over all combinations of samples, overlaps, and block sizes
        for idx, (sample, overlap, blocksize) in enumerate(itertools.product(subsets, overlaps, blocksizes)):
            # Use notify script every 10th job
            # notify = (idx + 1) % 10 == 0
            # notify = True
            # Wait for 1 second before submitting each job
            time.sleep(1)
            submit_job(sample, overlap, blocksize, iterating_blocksize=True, notify=True, queue="long")

if __name__ == "__main__":
    main()
