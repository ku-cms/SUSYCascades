import os
import ROOT

# Directory where subdirectories and root files are located
base_dir = "/ospool/cms-user/zflowers/NTUPLES/Processing/Summer23BPix_130X_Essentials/"

# List of branch names to check
branches_to_check = ["weight", "XSec", "Nweight", "Nevent"]

# Loop through subdirectories in the base directory
for subdir in os.listdir(base_dir):
    subdir_path = os.path.join(base_dir, subdir)
    
    if os.path.isdir(subdir_path):  # Check if it's a subdirectory
        # Assuming the root file is in the subdirectory (replace with correct filename if needed)
        root_file_path = os.path.join(subdir_path, subdir_path)
        root_file_path += "/"+subdir+"_0_0.root"
        
        if os.path.exists(root_file_path):
            # Open the root file
            root_file = ROOT.TFile.Open(root_file_path)
            
            if root_file and not root_file.IsZombie():
                # Get the tree from the file
                tree = root_file.Get("KUAnalysis")
                
                if tree:
                    # Load the branches in the tree
                    tree.SetBranchStatus("*", 0)  # Disable all branches to start
                    for branch in branches_to_check:
                        tree.SetBranchStatus(branch, 1)  # Enable the branches we want to check
                    if tree.GetEntries() == 0: continue
                    # Get the first entry (event) and check the values
                    tree.GetEntry(0)
                    
                    # Check the value of each branch for the first event
                    for branch in branches_to_check:
                        branch_value = getattr(tree, branch)
                        
                        if branch_value == 0:
                            print(f"Warning: In file {root_file_path}, branch '{branch}' has value 0 for the first event.")
                        
                    # Close the root file after processing
                    root_file.Close()
                else:
                    print(f"Warning: Tree 'KUAnalysis' not found in {root_file_path}.")
            else:
                print(f"Warning: Could not open root file {root_file_path}.")

