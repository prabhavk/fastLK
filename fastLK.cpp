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
    string path = "";
    bool verbose_flag_bool = true;
    string verbose_flag_string;
    float clv_threshold = 0;
    string path_to_sequence_alignment_file = "";
    
    cout << "Command to be executed is ";
    for (int i = 1; i < argc ; i++) {
        cout << argv[i] << "\t";
    }
    cout << endl;    

    if (argc < 2) {        
        cerr << "Usage: " << argv[0] << " -r reference_file -p path -s multiple_sequence_alignment -e clv_threshold -m mutation_diff_file -w workflow_type -t tree_file -l log_file -v verbose_flag" << endl;
        return (-1);
    } else {
        if (verbose_flag_bool) {        
            cout << "[";
        for (int i = 1; i < argc -1 ; i++) {
            cout << "\"" << argv[i] << "\"" << ",\t";
        }
         cout << "\"" << argv[argc-1] << "\"";
        cout << endl;        
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
        // verbose flag
            } else if (strcmp(argv[i], "-v") == 0) {
                if (i < argc -1) {
                    verbose_flag_string = argv[++i];
                    if (strcmp(verbose_flag_string.c_str(), "1") == 0) {                        
                        verbose_flag_bool = true;
                    } else {
                        verbose_flag_bool = false;
                    }
                }
        // path to directory containing input files
            } else if (strcmp(argv[i], "-p") == 0) {
                if (i < argc -1) {
                    path = argv[++i];
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
        // path to conditional likelihood vector threshold        ss
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
    if (verbose_flag_bool) {
        cout << "verbose_flag_bool is True" << endl;
    } else {
        cout << "verbose_flag_bool is False" << endl;
    }
    
    // add "path" to input files
    cout << "path to input files is " << path << endl;
    path_to_reference_file  = path + path_to_reference_file;
    path_to_mut_diff_file  = path + path_to_mut_diff_file;
    path_to_sequence_alignment_file  = path + path_to_sequence_alignment_file;
    path_to_tree_file  = path + path_to_tree_file;
    path_to_log_file = path + path_to_log_file;

    fastLK_overview * fastLK_manager = new fastLK_overview(path_to_reference_file, path_to_mut_diff_file, path_to_sequence_alignment_file, path_to_tree_file, path_to_log_file, workflow_type, clv_threshold,verbose_flag_bool);

    delete fastLK_manager;    
    }

    return (0); 
}