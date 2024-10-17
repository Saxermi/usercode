import os
import subprocess
import itertools
import time

# Global variable to control test mode
TEST_MODE = True  # Set to False to submit all jobs


def submit_job(sample, overlap, notify=False):
    """
    Function to submit the slurm job using the corresponding bash script.

    Args:
        sample (str): The sample name, e.g., TTbar_01.
        overlap (float): The overlap factor, e.g., 0.5.
        notify (bool): Whether to use the notification bash script.
    """
    # Determine which script to use
    bash_script = "pvslurmmaster_notifyme.bsh" if notify else "pvslurmmaster.bsh"

    # Calculate the directory name based on overlap (convert to an integer percentage)
    overlap_dir = (
        int(overlap * 10) * 10
    )  # Converts 0.1 to 10, 0.2 to 20, ..., 0.9 to 90
    path = f"experimental_run_1/{overlap_dir}"

    # Define the command to run
    cmd = [
        "sbatch",
        bash_script,
        "-n",
        "10",
        "-d",
        sample,
        "-o",
        str(overlap),
        "-p",
        path,  # Pass the path to store output in the correct directory
        "-l",
        "pvBlock",
    ]

    # Submit the job using subprocess and capture the output
    try:
        result = subprocess.run(
            cmd,
            check=True,
            universal_newlines=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        # Print the output from sbatch (which typically includes the job ID)
        print(result.stdout)
        print(f"Job for {sample} with overlap {overlap} submitted successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Failed to submit job for {sample} with overlap {overlap}. Error: {e}")


def main():
    # Samples to be used (from TTbar_01 to TTbar_15)
    samples = [f"TTbar_{str(i).zfill(2)}" for i in range(1, 16)]

    # Overlap values from 0.0 to 0.9 in 0.1 increments
    overlaps = [round(x * 0.1, 1) for x in range(10)]

    # If in test mode, only submit two jobs, one with notify and one without
    if TEST_MODE:
        print("Running in test mode...")
        submit_job(samples[0], overlaps[0], notify=False)
        submit_job(samples[0], overlaps[0], notify=True)
    else:
        # Iterate over all combinations of samples and overlap factors
        for idx, (sample, overlap) in enumerate(itertools.product(samples, overlaps)):
            # Use notify script every 10th job
            notify = (idx + 1) % 10 == 0

            # Wait for 1 second before submitting each job
            time.sleep(1)

            submit_job(sample, overlap, notify=notify)


if __name__ == "__main__":
    main()
