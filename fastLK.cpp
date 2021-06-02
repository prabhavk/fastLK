#include "fastLK.hpp"
#include <string>

// using namespace std;

int main(int argc, char * argv[])
{
    string path_to_reference_file = "";
    string path_to_mut_diff_file = "";
    string path_to_tree_file = "";
    string path_to_log_file = "";
    string workflow_type = "";
    float clv_threshold = 0;
    string path_to_sequence_alignment_file = "";

    if (argc < 2) {        
        cerr << "Usage: " << argv[0] << " -r reference_file -s multiple_sequence_alignment -e clv_threshold -m mutation_diff_file -w workflow_type -t tree_file -l log_file" << endl;
        return (-1);
    } else {        
        // parse arguments            
        for (int i = 1; i < argc ; i++) {
        // path to reference fasta file
            if (strcmp(argv[i], "-r") == 0) {
                if (i < argc -1) {
                    path_to_reference_file = argv[++i];                    
                }
        // path to multiple sequence alignment file
            } else if (strcmp(argv[i], "-s") == 0) {
                if (i < argc -1) {
                    path_to_sequence_alignment_file = argv[++i];
                }
        // path to mut diff file
            } else if (strcmp(argv[i], "-m") == 0) {
                if (i < argc -1) {
                    path_to_mut_diff_file = argv[++i];
                }
        // path to input tree (unrooted tree with branch lengths)
            } else if (strcmp(argv[i], "-t") == 0) {              
                if (i < argc -1) {
                    path_to_tree_file = argv[++i];
                }
        // path to conditional likelihood vector threshold        
            } else if (strcmp(argv[i], "-e") == 0) {              
                if (i < argc -1) {
                    clv_threshold = stof(argv[++i]);
                }
            } else if (strcmp(argv[i], "-l") == 0) {
                if (i < argc - 1) {
                    path_to_log_file = argv[++i];
                }                
            } else if (strcmp(argv[i], "-w") == 0) {
                if (i < argc - 1) {
                    workflow_type = argv[++i];
                }                
            }
        }
    }

    fastLK_overview * fastLK_manager = new fastLK_overview(path_to_reference_file, path_to_mut_diff_file, path_to_sequence_alignment_file, path_to_tree_file, path_to_log_file, workflow_type, clv_threshold);

    delete fastLK_manager;

    return (0); 
}