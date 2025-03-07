import os
import paramiko
from scp import SCPClient

#test

def create_ssh_client(server, user, password):
    ssh = paramiko.SSHClient()
    ssh.load_system_host_keys()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect(server, username=user, password=password)
    return ssh

def download_files(ssh, remote_path, local_path, os_type):
    with SCPClient(ssh.get_transport()) as scp:
        stdin, stdout, stderr = ssh.exec_command(f"find {remote_path} -type f -name '*.root'")
        files = stdout.readlines()
        
        for file in files:
            file = file.strip()
            # Preserving directory structure
            relative_path = file.replace(remote_path, "").lstrip('/')
            
            if os_type == "Windows":
                relative_path = relative_path.replace('/', '\\')
            
            local_file_path = os.path.join(local_path, relative_path)
            local_dir = os.path.dirname(local_file_path)
            
            if not os.path.exists(local_dir):
                os.makedirs(local_dir)
                
            print(f"Downloading {file} to {local_file_path}")
            scp.get(file, local_file_path)

def main():
    # Hardcoded server and remote path
    server = "t3ui01.psi.ch"  # Replace with the actual server address or hostname
    remote_path = "cmssw/CMSSW_14_1_0_pre7/src/usercode/PrimaryVertexAnalyzer/test"  # Replace with the actual remote path where .root files are stored

    # Ask user for OS, username, and password
    os_type = input("Which OS are you using? (Linux/Windows): ").strip()
    
    if os_type not in ['Linux', 'Windows']:
        print("Invalid OS type.")
        return

    # Prompt for username and password
    username = input("Enter your SSH username: ").strip()
    password = input("Enter your SSH password: ").strip()

    # Local directory input
    local_path = input("Enter the local directory where files should be saved: ").strip()
    
    if not os.path.exists(local_path):
        print(f"Local directory {local_path} does not exist.")
        return
    
    # Establish SSH connection and download files
    try:
        ssh = create_ssh_client(server, username, password)
        download_files(ssh, remote_path, local_path, os_type)
        print("Download completed successfully.")
    except Exception as e:
        print(f"An error occurred: {e}")
    finally:
        ssh.close()

if __name__ == "__main__":
    main()
