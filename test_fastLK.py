# Create script files for comparing log-likelihood values
import subprocess as sub
import sys
import os

if os.path.isdir('/project/exaptation/'):
    # print ("exaptation")
    test_data_path = "/project/exaptation/fast_approx_pruning/test_data/"
    tool_path = "/project/exaptation/Tools/"
    fastLK_command_path = "/project/exaptation/fast_approx_pruning/phylo_likelihood_mut_diff/fastLK"
elif os.path.isdir('/home/kalaghat/'):
    # print ("voyager")
    test_data_path = "/home/kalaghat/Projects/gisaid_nicolaDM/test_data/"
    tool_path = "/home/kalaghat/Tools/"
    fastLK_command_path = "/home/kalaghat/Projects/phylo_likelihood_mut_diff/fastLK"

# paths and file names
test_prefix = "test_4a"
# test_prefix = sys.argv[1]

phyml_command = tool_path +  "phyml-mpi -i " + test_data_path + test_prefix + ".phy -o n -u " +  test_data_path + test_prefix + ".fastree -f 0.25,0.25,0.25,0.25 -m 0.04,0.3,0.1,0.02,1.0"
sub.call (phyml_command, shell=True)
iqtree_command = tool_path + "iqtree2 -s " + test_data_path + test_prefix + ".fa -st DNA -te " + test_data_path + test_prefix + ".fastree "  + " -pre "+ test_prefix + "-blfix -keep-ident -m GTR\{0.04,0.3,0.1,0.02,1.0,1.0\}+FQ -redo -nt 1"
sub.call (iqtree_command, shell=True)
fastLK_cpp_command = fastLK_command_path + " -r " + test_data_path + "EPI_ISL_402124_lowercase.fasta -s + " + test_data_path + test_prefix + ".fa -m " + test_data_path + test_prefix + ".mut_diff -w mut_diff -t " + test_data_path + test_prefix + ".fastree"
sub.call (fastLK_cpp_command, shell=True)
fastLK_python_command = "python3 " + tool_path + "calculateFastLK.py --path " + test_data_path + " --reference EPI_ISL_402124_lowercase.fasta  --diffFile " + test_prefix + ".mut_diff --treeFile " + test_prefix + ".fastree"
sub.call (fastLK_python_command, shell=True)

# Parse output files of various tools and compile loglikelihood values into a csv file