#ifndef FASTLK_H
#define FASTLH_H
#include <vector>
#include <memory>
#include <map>
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <iomanip>
#include <regex>
#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/MatrixFunctions>
#include <boost/algorithm/string.hpp>
#include <algorithm>
#include <numeric>
#include <chrono>
// #include <boost/algorithm/string/replace.hpp>
using namespace Eigen;
// using namespace std;

// global variables

const unsigned int RESERVE_genome_list_elem = 200;
const unsigned int GENOME_LENGTH = 29891;


struct genome_list_elem {
    unsigned char dna;
    unsigned int start_pos;
    unsigned int n_pos;
    float length = 0;
    bool flag_for_removal = false;
    void Set_dna_pos_and_length(unsigned char dna, unsigned int start_pos, unsigned int n_pos, float length);
    // array <double, 4> clv;
    // double log_scaling_factor;
    genome_list_elem() {

    }
    ~genome_list_elem() {
        // delete this;

    }

};

void genome_list_elem::Set_dna_pos_and_length(unsigned char dna, unsigned int start_pos, unsigned int n_pos, float length) {
    this->dna = dna;
    this->start_pos = start_pos;
    this->n_pos = n_pos;
    this->length = length;
}

class node {
public:
    std::vector < node * > neighbors;
    std::vector < node * > children;
    bool verbose = true;     
    int degree = 0;
    int in_degree = 0;
    int out_degree = 0;
    int times_visited = 0;  
    std::string concatenated_descendant_names = "";
    bool leaf = false;
    bool empty_genome_list = true;
    node * parent = this;
    float initial_edge_length;
    std::array <double, 4> clv;
    double log_scaling_factor = 0;
    std::string name;
    std::vector <unsigned char> complete_sequence;
    std::vector <unsigned char> compressed_sequence;

    std::vector <genome_list_elem *> genome_list;        
    std::map <genome_list_elem *, std::array <double, 4>> clv_for_list_elem;
    
    void Add_parent(node * p);
    void Add_child(node * c);
    void Set_pos_for_reference_characters();
    node (std::string node_name) {        
        this->name = node_name;        
        genome_list.reserve(RESERVE_genome_list_elem);
    }

    ~node() {        
        for (genome_list_elem * list_elem : this->genome_list) {
            delete list_elem;
        }

    }
};

void node::Set_pos_for_reference_characters() {
    // if length of genome list is empty
    if (this->empty_genome_list) {
        genome_list_elem * list_elem_ref = new genome_list_elem();
        list_elem_ref->dna = 4; // 4 is the code for reference sequence
        list_elem_ref->length = this->initial_edge_length;
        list_elem_ref->start_pos = 1;
        list_elem_ref->n_pos = GENOME_LENGTH;
        this->genome_list.push_back(list_elem_ref);
        this->empty_genome_list = false;
    } else {
        int num_genome_list_elems = this->genome_list.size();    
        genome_list_elem * list_elem;    
        for (int i = 0; i < num_genome_list_elems; i++) {
            list_elem = this->genome_list[i];
            if (list_elem->dna == 4) { // 4 is code for refrence character
                if (i == 0) {
                    // reference list element is the first genome list element
                    list_elem->start_pos = 1;
                    list_elem->n_pos = this->genome_list[i+1]->start_pos - list_elem->start_pos;
                } else if (i == num_genome_list_elems -1) {
                    // reference list element is the last genome list element
                    list_elem->start_pos = this->genome_list[i-1]->start_pos + this->genome_list[i-1]->n_pos;
                    list_elem->n_pos = GENOME_LENGTH - list_elem->start_pos + 1;
                } else {
                    // reference list element is an intermediate genome list element
                    list_elem->start_pos = this->genome_list[i-1]->start_pos + this->genome_list[i-1]->n_pos;
                    list_elem->n_pos = this->genome_list[i+1]->start_pos - list_elem->start_pos;
                }
            }
            if (verbose) {assert (list_elem->n_pos > 0);}
        }
    }    
}

void node::Add_parent(node * p) {
    this->in_degree += 1;
    this->parent = p;
}

void node::Add_child(node * c) {
    this->out_degree += 1;
    this->children.push_back(c);
}

class tree {
public:
    node * root;
    void Set_descendant_names();
    double log_likelihood = 0.0;
    double clv_threshold;
    int h_ind = 0;
    bool verbose = false;
    bool debug = false;
    bool logging = true;    
    bool rooted;
    std::vector <unsigned char> reference_sequence;
    int num_char_patterns;
    std::string tree_file_name;
    std::string alignment_file_name;
    std::string workflow_type;
    std::string mut_diff_file_name;
    std::string reference_file_name;
    std::string log_file_name;
    std::ofstream log_file;
    void Read_newick_file();
    // void WriteNewickFile(string tree_file_name, bool rooted);
    void Root_unrooted_tree();
    std::map <std::string, unsigned char> DNA_to_char;
    std::map <std::string, node*> node_list;
    std::vector <node *> leaves;
    void Set_leaves();
    void Set_root();
    void Add_fasta_sequences();
    void Add_mut_diff_sequences();
    void Populate_DNA_to_char_map();
    void Compress_sequences();
    void Add_ref_nuc_counts_based_on_genome_coordinates();
    std::map <unsigned int, std::array <unsigned int, 4>> cum_ref_nuc_counts_map; // store interval counts of a, c, g, t
    std::vector <int> character_pattern_weights;
    void Set_model_parameters();
    void Set_pi_rho_using_ref_seq();
    void Optimize_UNREST();
    void Optimize_branch_lengths();
    float Compute_scaling_factor(Matrix4f Q);    
    void Read_reference_sequence();
    void Compute_loglikelihood_using_standard_pruning_algorithm();    
    void Compute_loglikelihood_using_fast_pruning_algorithm();
    std::array <double, 4> Get_clv_for_dna(unsigned char dna_obs);
    void Root_tree_along_edge(node * u, node * v, float dist_from_u); 
    void Set_nodes_for_postorder_traversal();
    void Add_node(std::string node_name);
    void Reset_times_visited();
    node * Get_node(std::string node_name);
    bool Contains_node(std::string node_name);
    void Add_undirected_edge(node * u, node * v, float length);
    void Add_directed_edge(node * u, node * v, float length);
    void OptimizeBranchLengths();
    std::vector <node *> nodes_for_postorder_traversal;
    float Get_undirected_edge_length(node * u, node * v);
    float Get_directed_edge_length(node * p, node * c);
    void Reset_log_scaling_factors_and_clvs();
    std::map <std::pair <node *, node *>, float> undirected_edge_length_map;
    std::map <std::pair <node *, node *>, float> directed_edge_length_map;
    Matrix4f Q;
    void Update_list_elements(node* parent, node * left_child, node * right_child, genome_list_elem * list_elem_parent, genome_list_elem * list_elem_left_child, genome_list_elem * list_elem_right_child);
    void Combine_ref_type_list_elements(node * n);
    void Compute_log_likelihood_score_for_tree_using_genome_list();
    std::array <float, 4> pi_rho;
    // TESTING
    void verbose_combining_genome_ref_elements();

    // Matrix4f Q_UNREST;

    tree (){
        // Only hardcoded value should be set here
        this->Populate_DNA_to_char_map();
        // root = new node("h_root");
        // this->node_list["h_root"] = root;         
    }

    ~tree(){
        // this->node_list.clear();
        for (std::pair <std::string, node *> elem : this->node_list) {
            delete elem.second;
        }
    }
};


void tree::OptimizeBranchLengths() {
    float t_curr;
    for (std::pair<std::pair <node*, node*>, float> nodePair_Length: this->directed_edge_length_map) {
        std::pair <node*, node*> edge = nodePair_Length.first;
        t_curr = nodePair_Length.second;
        // estimate 1 + (q_xx)*(t_curr)
        float max_val = 0;
        for (int x = 0; x < 4; x++) {
            if (max_val < abs(Q(x,x)*t_curr)) {
                max_val = abs(Q(x,x)*t_curr);
            }
        }
        std::cout << "Maximum value of (q_xx)*(t_curr) is" << max_val << std::endl;
    }
}

void tree::Root_unrooted_tree() {
    node * l;
    if (this->verbose) {
        if (this->Contains_node("EPI_ISL_1208981")) {
                l = this->Get_node("EPI_ISL_1208981");
            } else {
                l = this->leaves[0];
            }
    }
    // std::cout << "node for rooting is " << l->name << std::endl;
    this->Root_tree_along_edge(l, l->neighbors[0], 0.0);
    this->Set_nodes_for_postorder_traversal();
}

void tree::Set_descendant_names() {
    node * c_l; node * c_r;
    for (node * n : this->nodes_for_postorder_traversal) {
        if (n->leaf) {
            n->concatenated_descendant_names = n->name;
        } else {
            c_l = n->children[0];
            c_r = n->children[1];
            n->concatenated_descendant_names = c_l->concatenated_descendant_names + "," + c_r->concatenated_descendant_names;
            // std::cout << n->concatenated_descendant_names << std::endl;
        }
    }
}

void tree::Reset_log_scaling_factors_and_clvs() {
    for (std::pair <std::string, node *> elem: this->node_list) {
        elem.second->clv[0] = 0; elem.second->clv[1] = 0; elem.second->clv[2] = 0; elem.second->clv[3] = 0;
        elem.second->log_scaling_factor = 0;
    }
}

std::array <double, 4> tree::Get_clv_for_dna(unsigned char dna_obs) {
    std::array <double, 4> clv;    
    // if (verbose) {assert (dna_obs < 4); }    
    clv[0] = 0.0; clv[1] = 0.0; clv[2] = 0.0; clv[3] = 0.0;
    if (dna_obs < 4) { // a, c, g, t
        clv[dna_obs] = 1.0;
    } else if (dna_obs == 5) { // -, n
        clv[0] = 1.0; clv[1] = 1.0; clv[2] = 1.0; clv[3] = 1.0;
    } else if (dna_obs == 6) { // w = a, t
        clv[0] = 1.0; clv[3] = 1.0;
    } else if (dna_obs == 7) { // s = c, g
        clv[1] = 1.0; clv[2] = 1.0;
    } else if (dna_obs == 8) { // m = a, c
        clv[0] = 1.0; clv[1] = 1.0;
    } else if (dna_obs == 9) { // k = g, t
        clv[2] = 1.0; clv[3] = 1.0;
    } else if (dna_obs == 10) { // r = a, g
        clv[0] = 1.0; clv[2] = 1.0;
    } else if (dna_obs == 11) { // y = c, t
        clv[1] = 1.0; clv[3] = 1.0;
    } else if (dna_obs == 12) { // b = c, g, t
        clv[1] = 1.0; clv[2] = 1.0; clv[3] = 1.0;
    } else if (dna_obs == 13) { // d = a, g, t
        clv[0] = 1.0; clv[2] = 1.0; clv[3] = 1.0;
    } else if (dna_obs == 14) { // h = a, c, t
        clv[0] = 1.0; clv[1] = 1.0; clv[3] = 1.0;
    } else if (dna_obs == 15) { // v = a, c, g
        clv[0] = 1.0; clv[1] = 1.0; clv[2] = 1.0;
    }        
    if (this->verbose) {assert(*std::max_element(clv.begin(),clv.end()) > 0.0);}   
    return (clv);
}

void tree::Reset_times_visited() {
    for (std::pair <std::string, node *> node_elem: this->node_list) {
        node_elem.second->times_visited = 0;
    }
}

void tree::Set_nodes_for_postorder_traversal() {
    this->nodes_for_postorder_traversal.clear();
    for (std::pair <std::string, node *> elem : this->node_list) {
        elem.second->times_visited = 0;
    }
    std::vector <node *> nodes_to_visit;
    for (node * n: this->leaves) {
        nodes_to_visit.push_back(n);
        this->nodes_for_postorder_traversal.push_back(n);
    }
    int num_nodes_to_visit = nodes_to_visit.size();
    node * v;
    while (num_nodes_to_visit > 0) {
        v = nodes_to_visit[num_nodes_to_visit -1];
        num_nodes_to_visit --;
        nodes_to_visit.pop_back();
        v->parent->times_visited ++;
        if (v->parent->times_visited == 2) {
            nodes_to_visit.push_back(v->parent);
            this->nodes_for_postorder_traversal.push_back(v->parent);
            num_nodes_to_visit ++;
        }
    }
    // std::cout << "Nodes for postorder traversal are" << std::endl;
    // for (node * n: this->nodes_for_postorder_traversal) {        
    //     std::cout << n->name << std::endl;
    // }
    // std::cout << "---------------------------------" << std::endl;    
}

void tree::Add_directed_edge(node * p, node * c, float edge_length) {
    c->Add_parent(p);
    p->Add_child(c);
    this->directed_edge_length_map[std::make_pair(p,c)] = edge_length;
}


void tree::Root_tree_along_edge(node * u, node * v, float dist_from_u) {
    this->root = new node("h_root");
    this->node_list["h_root"] = this->root;
    std::pair <node *, node *> edge;
    if (u < v) {
        edge = std::make_pair(u,v);
    } else {
        edge = std::make_pair(v,u);
    }
    if (this->verbose) {
        std::cout << u->name << " is a leaf with degree " << u->degree << std::endl;
        std::cout << v->name << " is not a leaf with degree " << v->degree << std::endl;
    }

    if (verbose) {assert(this->undirected_edge_length_map.find(edge) != this->undirected_edge_length_map.end());}    
    std::vector <node *> nodes_visited;
    std::vector <node *> nodes_to_visit; 
    nodes_to_visit.push_back(u);
    nodes_to_visit.push_back(v);
    nodes_visited.push_back(u);
    nodes_visited.push_back(v);
    float edge_length;
    edge_length = this->Get_undirected_edge_length(u,v);
    if (verbose) {assert(edge_length >= dist_from_u);}

    this->Add_directed_edge(this->root,u,dist_from_u);
    this->Add_directed_edge(this->root,v,edge_length - dist_from_u);
    int num_nodes_to_visit = nodes_to_visit.size();
    node * p;
    while (num_nodes_to_visit > 0) {
        p = nodes_to_visit[num_nodes_to_visit - 1];        
        nodes_to_visit.pop_back();
        nodes_visited.push_back(p);
        num_nodes_to_visit --;
        for (node * c: p->neighbors) {
            if (find(nodes_visited.begin(),nodes_visited.end(),c) == nodes_visited.end()) {
                nodes_to_visit.push_back(c);
                num_nodes_to_visit ++;                
                edge_length = this->Get_undirected_edge_length(p,c);
                this->Add_directed_edge(p,c,edge_length);
            }       
        }
    }    
    // Preorder traversal operations
}

void tree::Compute_loglikelihood_using_standard_pruning_algorithm() {
    Matrix4f Q = this->Q;
    Matrix4f Q_scaled_l; Matrix4f Q_scaled_r;
    Matrix4f P_l; Matrix4f P_r;
    float t_l; float t_r;
	float scaling_factor;
    node * c_l; node * c_r;
    unsigned char ch;    
	double partial_likelihood_l;
    double partial_likelihood_r;
    double log_likelihood_contri_from_root_seq = 0;
    double largest_elem;
	this->log_likelihood = 0;
	double site_likelihood;
    // bool verbose = true;
    unsigned int num_char_patterns = this->character_pattern_weights.size();
    // std::cout << "Iterating over character patterns" << std::endl;
	for (unsigned int site = 0; site < num_char_patterns; site++) {
    // for (unsigned int site = 0; site < 1; site++) {
        this->Reset_log_scaling_factors_and_clvs();
        if (this->verbose) {
            std::cout << "log scaling factors and clvs resetted" << std::endl;
        }        
        for (node * n : this->nodes_for_postorder_traversal) {
            if (n->leaf) {
                // set conditional likelihood vector using observed character  
                if (this->verbose) {
                    std::cout << "leaf node name is\t" << n->name << std::endl;
                    std::cout << " dna for " << " site " << site << " is " << (int) n->compressed_sequence[site] << std::endl;
                }                                              
                n->clv = this->Get_clv_for_dna(n->compressed_sequence[site]);
            } else {
                if (this->verbose) {
                    std::cout << "Non-leaf node name is\t" << n->name << std::endl;                    
                }
                // compute conditional likelihood vector by multiplying partial likelihoods
                c_l = n->children[0];
                c_r = n->children[1];
                t_l = this->Get_undirected_edge_length(n,c_l);
                t_r = this->Get_undirected_edge_length(n,c_r);
                Q_scaled_l = Q*t_l; Q_scaled_r = Q*t_r;
                P_l = Q_scaled_l.exp(); P_r = Q_scaled_r.exp();
                largest_elem = 0;
                for (unsigned dna_p = 0; dna_p < 4; dna_p ++) {
                    partial_likelihood_l = 0; partial_likelihood_r = 0;
                    for (unsigned dna_c = 0; dna_c < 4; dna_c ++) {
                        partial_likelihood_l += P_l(dna_p,dna_c) * c_l->clv[dna_c];
                        partial_likelihood_r += P_r(dna_p,dna_c) * c_r->clv[dna_c];
                    }
                    n->clv[dna_p] = partial_likelihood_l * partial_likelihood_r;
                    if (largest_elem < n->clv[dna_p]) {
                        largest_elem = n->clv[dna_p];
                    }
                }
                // rescale clvs
                for (unsigned dna_p = 0; dna_p < 4; dna_p ++) {
                    n->clv[dna_p] /= largest_elem;
                }                
                n->log_scaling_factor += log(largest_elem);
                n->log_scaling_factor += c_l->log_scaling_factor + c_r->log_scaling_factor;
                c_l->log_scaling_factor = 0; c_r->log_scaling_factor = 0;
            }
        }
        
		site_likelihood = 0;
		for (unsigned char dna = 0; dna < 4; dna++) {			
			site_likelihood += this->pi_rho[dna] * this->root->clv[dna];
		}
        log_likelihood_contri_from_root_seq += log(site_likelihood) * this->character_pattern_weights[site];
        if (this->verbose) {
            std::cout << "log_likelihood_contri_from_root_seq\t" << std::setprecision(8) << log_likelihood_contri_from_root_seq << std::endl;
        }        
		this->log_likelihood += (this->root->log_scaling_factor + log(site_likelihood)) * this->character_pattern_weights[site];
        this->root->log_scaling_factor = 0;
	}
}


void tree::Set_leaves() {
    this->leaves.clear();
    for (std::pair <std::string, node *> elem : this->node_list) {
        if (elem.second->degree == 1) {
            this->leaves.push_back(elem.second);
            elem.second->leaf = true;
        } else {
            elem.second->leaf = false;
        }
    }
    if (this->verbose) {
        std::cout << "Number of leaves is " << this->leaves.size() << std::endl;
    }    
}

void tree::Set_root() {
    for (std::pair <std::string, node *> elem : this->node_list) {
        if (elem.second->in_degree == 0) {
            this->root = elem.second;
            // std::cout << "Root has been set" << std::endl;
        }
    }    
}

void tree::verbose_combining_genome_ref_elements() {
    node * n = new node("test");
    unsigned char dna_ref = 4;
    unsigned char dna_diff = 2;
    genome_list_elem * list_elem_ref_1a = new genome_list_elem();    
    genome_list_elem * list_elem_ref_1b = new genome_list_elem();
    genome_list_elem * list_elem_ref_1c = new genome_list_elem();
    genome_list_elem * list_elem_diff_1 = new genome_list_elem();
    genome_list_elem * list_elem_ref_2a = new genome_list_elem();
    genome_list_elem * list_elem_ref_2b = new genome_list_elem();        
    genome_list_elem * list_elem_ref_2c = new genome_list_elem(); 
    list_elem_ref_1a->Set_dna_pos_and_length(dna_ref, 1, 10, 0.01);
    list_elem_ref_1b->Set_dna_pos_and_length(dna_ref, 11, 10, 0.01);
    list_elem_ref_1c->Set_dna_pos_and_length(dna_ref, 21, 10, 0.01);
    list_elem_diff_1->Set_dna_pos_and_length(dna_diff, 31, 10, 0.0);
    list_elem_ref_2a->Set_dna_pos_and_length(dna_ref, 41, 10, 0.01);
    list_elem_ref_2b->Set_dna_pos_and_length(dna_ref, 51, 10, 0.01);
    list_elem_ref_2c->Set_dna_pos_and_length(dna_ref, 61, 10, 0.00);

    n->genome_list.push_back(list_elem_ref_1a), n->genome_list.push_back(list_elem_ref_1b),n->genome_list.push_back(list_elem_ref_1c);
    n->genome_list.push_back(list_elem_diff_1), n->genome_list.push_back(list_elem_ref_2a),n->genome_list.push_back(list_elem_ref_2b),n->genome_list.push_back(list_elem_ref_2c);

    this->Combine_ref_type_list_elements(n);
    delete n;
}

void tree::Combine_ref_type_list_elements(node * n) {
    // collapse contiguous list elements of type reference
    int num_list_elem;
    bool continue_search_for_elements_to_combine = false;
    std::vector <genome_list_elem *>::iterator start_pos_iterator = n->genome_list.begin();
    int ind_of_elem_extended;
    if (verbose) {
        for (int i = 0; i < n->genome_list.size() - 1; i++) {
            assert(n->genome_list[i]->start_pos + n->genome_list[i]->n_pos == n->genome_list[i+1]->start_pos);
        }
    }
    if (false) {
        std::cout << "Num elements before combining " << n->genome_list.size() << std::endl;
        for (genome_list_elem * list_elem : n->genome_list) {
            std::cout << "dna: " << (int) list_elem->dna;
            std::cout << ", start_pos: " << list_elem->start_pos;
            std::cout << ", n_pos: " << list_elem->n_pos;
            std::cout << ", length: " << list_elem->length << std::endl;        
        }
        // std::cout << "Combining list elements for " << n->name << std::endl;
    }    
    
    std::vector <genome_list_elem *> list_elems_to_remove;
    std::vector <genome_list_elem *> list_elems_to_add;
    // vector <int> ind_of_list_elem_to_extend;
    // float thresh = pow(10,-10);
    // map <genome_list_elem *, vector <int>> list_elems_to_join;        
    std::vector <genome_list_elem *>::iterator genome_list_end = n->genome_list.end();
    std::vector <genome_list_elem *>::iterator genome_list_begin = n->genome_list.begin();
    std::vector <genome_list_elem *>::iterator it;
    genome_list_elem * list_elem_to_extend;
    genome_list_elem * list_elem_i;
    genome_list_elem * list_elem_i_plus_one;
    bool list_elem_to_extend_found = false;
    bool combine_elements = false;
    float length_thresh =  pow(10,-10);
    for (it = genome_list_begin;  it < genome_list_end - 1; it++) {
        if ((*it)->dna == 4 && (*(it+1))->dna == 4 && abs((*it)->length - (*(it+1))->length) < length_thresh) {            
            if (!list_elem_to_extend_found) {
                list_elem_to_extend = *it;
                list_elems_to_add.clear();
                list_elem_to_extend_found = true;
                // list_elems_to_add.push_back(list_elem_to_extend);
                list_elems_to_add.push_back(*(it+1));
                list_elems_to_remove.push_back(*(it+1));
            } else {
                list_elems_to_add.push_back(*(it+1));
                list_elems_to_remove.push_back(*(it+1));
            }
            // std::cout << "dna: " << (int) (*it)->dna;
            // std::cout << ", start_pos: " << (*it)->start_pos;
            // std::cout << ", n_pos: " << (*it)->n_pos;
            // std::cout << ", length: " << (*it)->length << std::endl;
        } else {
            if (list_elem_to_extend_found) {
                for (genome_list_elem * list_elem : list_elems_to_add) {
                    list_elem_to_extend->n_pos += list_elem->n_pos;
                }
                list_elem_to_extend_found = false;
            }
        }
    }
    if (list_elem_to_extend_found) {
        // std::cout << "list elem to extend found" << std::endl;
        for (genome_list_elem * list_elem : list_elems_to_add) {
            list_elem_to_extend->n_pos += list_elem->n_pos;
            // std::cout << "dna: " << (int) list_elem->dna;
            // std::cout << ", start_pos: " << list_elem->start_pos;
            // std::cout << ", n_pos: " << list_elem->n_pos;
            // std::cout << ", length: " << list_elem->length << std::endl;
        }
        list_elem_to_extend_found = false;
    }

    for (genome_list_elem * list_elem: list_elems_to_remove) {
        n->genome_list.erase(remove(n->genome_list.begin(),n->genome_list.end(),list_elem),n->genome_list.end());
    }

    
    if (false) {
        std::cout << "Num elements after combining " << n->genome_list.size() << std::endl;
        for (genome_list_elem * list_elem : n->genome_list) {
            std::cout << "dna: " << (int) list_elem->dna;
            std::cout << ", start_pos: " << list_elem->start_pos;
            std::cout << ", n_pos: " << list_elem->n_pos;
            std::cout << ", length: " << list_elem->length << std::endl;        
        }
    }    
    if (verbose) {
        for (int i = 0; i < n->genome_list.size() - 1; i++) {
            assert(n->genome_list[i]->start_pos + n->genome_list[i]->n_pos == n->genome_list[i+1]->start_pos);
        }
    }
}

void tree::Update_list_elements(node * parent, node * left_child, node * right_child, genome_list_elem * list_elem_parent, genome_list_elem * list_elem_left_child, genome_list_elem * list_elem_right_child) {
    // dna_l and dna_r are characters observed for list_elem_left and list_elem_right        
    int n_pos;
    assert (list_elem_left_child->start_pos == list_elem_right_child->start_pos);
    // Set start pos for parent element    
    list_elem_parent->start_pos = list_elem_left_child->start_pos;
    if (list_elem_left_child->n_pos < list_elem_right_child->n_pos) {
        n_pos = list_elem_left_child->n_pos;
    } else {
        n_pos = list_elem_right_child->n_pos;
    }    
    // Set n pos for parent element
    list_elem_parent->n_pos = n_pos;    

    // Adjust the starting positions and n_pos of left and right child
    list_elem_left_child->n_pos -= n_pos;
    list_elem_left_child->start_pos += n_pos;
    list_elem_right_child->n_pos -= n_pos;
    list_elem_right_child->start_pos += n_pos;
    
    unsigned char dna_l = list_elem_left_child->dna;
    unsigned char dna_r = list_elem_right_child->dna;
    std::array <unsigned int, 4> ref_nuc_counts;
    std::array <double, 4> clv_parent;
    std::array <double, 4> clv_left_child;
    std::array <double, 4> clv_right_child;
    double largest_elem; double partial_likelihood_l; double partial_likelihood_r;
    unsigned int index_largest_elem;
    int num_clv_elems_larger_than_threshold;
    // if (list_elem_parent->start_pos == 6730) {
    //     std::cout << " ************************** ";
    //     std::cout << "Updating list element for " << parent->name << std::endl;
    // }
    // Iterate over combinations of dna in left child and right child and compute the clv for parent node    
    if (dna_l == 5) {        
        if (dna_r == 5) {            
            // 1) dna_l is N and dna_r is N
            list_elem_parent->dna = 5;            
        } else {
            // 2) dna_l is N and dna_r is not N
            // set dna of parent as dna of right child            
            list_elem_parent->dna = dna_r;
            // set branch length for parent as branch length of right child            
            list_elem_parent->length = list_elem_right_child->length;
            // copy clv if necessary
            if (right_child->clv_for_list_elem.find(list_elem_right_child) != right_child->clv_for_list_elem.end()) {
                parent->clv_for_list_elem[list_elem_parent] = right_child->clv_for_list_elem[list_elem_right_child];
                if (verbose) {assert(*std::max_element(parent->clv_for_list_elem[list_elem_parent].begin(),parent->clv_for_list_elem[list_elem_parent].end()) > 0.0);}
            }
        }
    } else {
        if (dna_r == 5) {
            // 3) dna_r is N and dna_l is not N
            // set dna of parent as dna of left child            
            list_elem_parent->dna = dna_l;
            // set branch length for parent as branch length of right child            
            list_elem_parent->length = list_elem_left_child->length;
            // copy clv if necessary
            if (left_child->clv_for_list_elem.find(list_elem_left_child) != left_child->clv_for_list_elem.end()) {
                parent->clv_for_list_elem[list_elem_parent] = left_child->clv_for_list_elem[list_elem_left_child];
                if (verbose) {assert(*std::max_element(parent->clv_for_list_elem[list_elem_parent].begin(),parent->clv_for_list_elem[list_elem_parent].end()) > 0.0);}
            }
        } else {
            // dna_l is not N and dna_r is not N            
            if (dna_l == dna_r && dna_l < 5) {
                // dna_l = dna_r and dna_l is in {a,c,g,t,x} (x is ref)            
                if (dna_l < 4) {
                    // 4) dna_l and dna_r are in {a,c,t,g}. n_pos is one
                    // if (n_pos != 1) {
                    //     std::cout << "ref pos is " << list_elem_parent->start_pos << std::endl;
                    //     std::cout << "n pos is " << list_elem_parent->n_pos << std::endl;                        
                    //     std::cout << "dna_l is " << (int) dna_l << "\t" << "dna_r is " << (int) dna_r << std::endl;
                    //     std::cout << "left child name is " << left_child->name << std::endl;
                    //     std::cout << "right child name is " << right_child->name << std::endl;
                    //     // std::cout << "clv left child " << dna_l << " is " << clv_left_child[dna] << "\t";
                    //     // std::cout << "clv right child " << dna_r << " is " << clv_right_child[dna] << "\t";
                    //     // std::cout << "clv parent " << dna << " is " << clv_parent[dna] << std::endl;
                    // }
                    if (verbose) {assert(n_pos == 1);}
                    // set dna of parent as dna of left child
                    list_elem_parent->dna = dna_l;
                    // set branch length of parent as 0 (default)
                    // compute likelihood as P(i,i,l) (probability that i remains i after time l)
                    // as P(i,i,l) = 1 + Q(i,i)l (because l is very small)
                    // because we are storing log-likelihood
                    // we can compute log (P(i,i,l)) ~ log (1 + Q(i,i)l) ~ Q(i,i)l
                    // thus contribution towards log_likelihood is Q(i,i)(l_c + l_r)
                    parent->log_scaling_factor += this->Q(dna_l,dna_l) * (list_elem_left_child->length + list_elem_right_child->length);
                } else if (dna_l == 4) {
                    // 5) dna_l and dna_r are reference nucleobases. n_pos is probably large
                    // std::cout << "log_likelihood added here 1" << std::endl;
                    // std::cout << "n_pos is " << list_elem_parent->n_pos << std::endl;
                    // set dna of parent as dna of left child
                    list_elem_parent->dna = dna_l;
                    // set branch length of parent as 0 (default)
                    // we can compute log (P(i,i,l)^{n_i}) ~ n_i log (1 + Q(i,i)l) ~ n_iQ(i,i)l 
                    // where n_i is the number of ref nucs that is i
                    // thus contribution towards log_likelihood is sum_i n_iQ(i,i)(l_c + l_r)
                    this->cum_ref_nuc_counts_map[list_elem_parent->start_pos - 1];
                    for (unsigned int dna = 0; dna < 4; dna ++) {
                        ref_nuc_counts[dna] = this->cum_ref_nuc_counts_map[list_elem_parent->start_pos + list_elem_parent->n_pos -1][dna];
                        ref_nuc_counts[dna] -= this->cum_ref_nuc_counts_map[list_elem_parent->start_pos -1][dna];
                        // std::cout << " ref_nuc_counts[" << (int) dna << "] is " << (int) ref_nuc_counts[dna] << std::endl;
                        // std::cout << " Q(" << (int) dna << "," << (int) dna << ") is " << this->Q(dna,dna) << std::endl;
                        // std::cout << " sum of edge lengths is " << list_elem_left_child->length + list_elem_right_child->length << std::endl;
                        parent->log_scaling_factor += ref_nuc_counts[dna] * this->Q(dna,dna) * (list_elem_left_child->length + list_elem_right_child->length);                        
                    }
                    // std::cout << "log likelihood added is " << parent->log_scaling_factor << std::endl;
                }
            } else {
                // 6) dna_l and dna_r are (i) ambiguous and identical, (ii) ambiguous and distinct, or (iii) non-ambiguous and distinct
                // set branch length for parent as 0 (default)
                // if (list_elem_parent->start_pos == 6730) {
                // std::cout << "dna_l is " << (int) dna_l << "\t" << "dna_r is " << (int) dna_r << std::endl;
                // }
                // compute clv of parent using clv of children
                // std::cout << "********************************************************************************" << std::endl;
                // std::cout << "case 6: dna_l and dna_r are (i) ambiguous and identical, (ii) ambiguous and distinct, or (iii) non-ambiguous and distinct" << std::endl;
                // std::cout << "********************************************************************************" << std::endl;
                if (left_child->clv_for_list_elem.find(list_elem_left_child) != left_child->clv_for_list_elem.end()) {
                    if (dna_l < 6) {
                        if (this->verbose) {
                            std::cout << "dna_l is " << int(dna_l) << std::endl;
                            std::cout << "clv entry for dna_l is " << left_child->clv_for_list_elem[list_elem_left_child][dna_l] << std::endl;
                        }                        
                    }
                    if (this->debug) {assert(dna_l > 5);}
                    clv_left_child = left_child->clv_for_list_elem[list_elem_left_child];
                    if (this->debug) {assert(*std::max_element(clv_left_child.begin(),clv_left_child.end()) > 0.0);}
                } else {
                    // std::cout << "I'm here now" << std::endl;
                    if (this->debug) {assert(dna_l < 5);}
                    if (dna_l < 4) {
                        clv_left_child = this->Get_clv_for_dna(dna_l);
                    } else {
                        // std::cout << "Hello! I'm here now" << std::endl;
                        // dna = 4 (reference)                        
                        // std::cout << "dna for reference is " << (int) this->reference_sequence[list_elem_parent->start_pos -1] << std::endl;
                        clv_left_child = this->Get_clv_for_dna(this->reference_sequence[list_elem_parent->start_pos -1]);
                        // std::cout << "max clv element is " << *std::max_element(clv_left_child.begin(),clv_left_child.end()) << std::endl;
                        if (this->debug) {assert(*std::max_element(clv_left_child.begin(),clv_left_child.end()) > 0.0);}
                    }
                }

                if (right_child->clv_for_list_elem.find(list_elem_right_child) != right_child->clv_for_list_elem.end()) {
                    // std::cout << "****\t" << (int) dna_r << std::endl; 
                    if (this->debug) {assert(dna_r > 5);}
                    clv_right_child = right_child->clv_for_list_elem[list_elem_right_child];
                    if (this->debug) {assert(*std::max_element(clv_right_child.begin(),clv_right_child.end()) > 0.0);}
                } else {
                    if (this->debug) {assert(dna_r < 5);}
                    if (dna_r < 4) {
                        clv_right_child = this->Get_clv_for_dna(dna_r);
                    } else {                        
                        // std::cout << "dna for reference is " << (int) this->reference_sequence[list_elem_parent->start_pos -1] << std::endl;
                        clv_right_child = this->Get_clv_for_dna(this->reference_sequence[list_elem_parent->start_pos -1]);
                        if (this->debug) {assert(*std::max_element(clv_right_child.begin(),clv_right_child.end()) > 0.0);}
                    }
                }
                if (this->debug) {assert(*std::max_element(clv_left_child.begin(),clv_left_child.end()) > 0.0);}
                if (this->debug) {assert(*std::max_element(clv_right_child.begin(),clv_right_child.end()) > 0.0);}
                // assert(list_elem_left_child->length > 0);
                // assert(list_elem_right_child->length > 0);
                clv_parent[0] = 1.0; clv_parent[1] = 1.0; clv_parent[2] = 1.0; clv_parent[3] = 1.0;
                largest_elem = 0;
                for (unsigned char dna_p = 0; dna_p < 4; dna_p ++) {
                    partial_likelihood_l = 0; partial_likelihood_r = 0;
                    for (unsigned char dna_c = 0; dna_c < 4; dna_c ++) {
                        if (dna_p == dna_c) {
                            // std::cout << (int) dna_p << "\t" << clv_left_child[dna_p] << "\t" << clv_right_child[dna_p] << std::endl;
                            partial_likelihood_l += (1.0 + (this->Q(dna_p,dna_p) * list_elem_left_child->length)) * clv_left_child[dna_c];
                            partial_likelihood_r += (1.0 + (this->Q(dna_p,dna_p) * list_elem_right_child->length)) * clv_right_child[dna_c];
                        } else {
                            partial_likelihood_l += (this->Q(dna_p,dna_c) * list_elem_left_child->length) * clv_left_child[dna_c];
                            partial_likelihood_r += (this->Q(dna_p,dna_c) * list_elem_right_child->length) * clv_right_child[dna_c];
                        }                        
                    }
                    clv_parent[dna_p] = partial_likelihood_l * partial_likelihood_r;
                    // std::cout << "updated clv is " << clv_parent[dna_p] << std::endl;
                    if (largest_elem < clv_parent[dna_p]) {
                        index_largest_elem = dna_p;
                        largest_elem = clv_parent[dna_p];
                    }
                    // std::cout << "largest element is " << largest_elem << std::endl;
                    // std::cout << "index of largest element is " << (int) index_largest_elem << std::endl;
                }                
                // rescale clvs
                for (unsigned char dna_p = 0; dna_p < 4; dna_p++) {
                    clv_parent[dna_p] /= largest_elem;                    
                }
                // std::cout << "log_likelihood added here 2" << std::endl;                
                parent->log_scaling_factor += log(largest_elem);
                parent->clv_for_list_elem[list_elem_parent] = clv_parent;
                if (this->debug) {assert(*std::max_element(clv_parent.begin(),clv_parent.end()) > 0.0);}
                // set dna of parent as type "O" if second largest element is larger than 10^-6 (user-defined threshold)
                num_clv_elems_larger_than_threshold = 0;
                for (unsigned char dna_p = 0; dna_p < 4; dna_p ++) {
                    if (clv_parent[dna_p] > this->clv_threshold) {
                        num_clv_elems_larger_than_threshold += 1;
                    }
                }
                if (num_clv_elems_larger_than_threshold > 1) {
                    list_elem_parent->dna = 16; // type "O"
                } else {
                    if (this->reference_sequence[list_elem_parent->start_pos - 1] == index_largest_elem) {
                        list_elem_parent->dna = 4; // reference
                        parent->clv_for_list_elem.erase(list_elem_parent);
                    } else {
                        // std::cout << "dna assigned to parent " << parent->name << " is " << (int) index_largest_elem << std::endl;
                        list_elem_parent->dna = index_largest_elem;
                        parent->clv_for_list_elem.erase(list_elem_parent);
                    }
                }
            }
        }
    }    
}


void tree::Compute_log_likelihood_score_for_tree_using_genome_list() {
    if (this->verbose) {
        std::cout << "-------------------------------------------------------" << std::endl;
        std::cout << "Computing log-likehood using genome list for the root" << std::endl;
    }    
    this->log_likelihood = this->root->log_scaling_factor;
    double log_likelihood_contri_from_root_seq = this->log_likelihood;
    int n_pos; unsigned int dna_r;
    std::array <double, 4> clv_root;
    double site_likelihood;
    std::array <unsigned int, 4> ref_nuc_counts;
    int tot_pos = 0;
    if (this->verbose) {
        std::cout << "Number of list elements in genome list for the root " << this->root->genome_list.size() << std::endl;
    }    
    for (genome_list_elem * list_elem_root : this->root->genome_list) {
        if (this->verbose) {
            std::cout << "########################################################" << std::endl;
        }
        n_pos = list_elem_root->n_pos;
        dna_r = list_elem_root->dna;
        if (this->verbose) {
            std::cout << "list_elem_root dna is " << (int) list_elem_root->dna << std::endl;        
            std::cout << "list_elem_root start pos is " << (int) list_elem_root->start_pos << std::endl;
            std::cout << "list_elem_root n_pos is " << (int) list_elem_root->n_pos << std::endl;
        }

        if (dna_r == 4) { // reference
            if (this->verbose) {
                std::cout << "***********************************************************" << std::endl;
            }        
            for (unsigned char dna = 0; dna < 4; dna ++) {                
                ref_nuc_counts[dna] = this->cum_ref_nuc_counts_map[list_elem_root->start_pos + list_elem_root->n_pos -1][dna];
                ref_nuc_counts[dna] -= this->cum_ref_nuc_counts_map[list_elem_root->start_pos -1][dna];
                // if (ref_nuc_counts[dna] == 0) {
                //     std::cout << "list elem start pos is " << list_elem_root->start_pos << std::endl;
                //     std::cout << "list elem n_pos is " << list_elem_root->n_pos << std::endl;
                //     std::cout << "I should be printed" << std::endl;
                //     std::cout << "cumulative counts for current list elem is ";
                //     std::cout << this->cum_ref_nuc_counts_map[list_elem_root->start_pos + list_elem_root->n_pos -1][dna] << std::endl;
                //     std::cout << "cumulative counts for previous list elem is " << this->cum_ref_nuc_counts_map[list_elem_root->start_pos -1][dna] << std::endl;
                // }
                // assert(ref_nuc_counts[dna] > 0);
                if (this->verbose) {
                    std::cout << "cum counts is " << (int) ref_nuc_counts[dna] << std::endl;
                    std::cout << "pi_rho is " << this->pi_rho[dna] << std::endl;
                }
                this->log_likelihood += ref_nuc_counts[dna] * log(this->pi_rho[dna]);
            }
            if (this->verbose) {            
                std::cout << "***********************************************************" << std::endl;
            }
        } else if (dna_r != 5) { // skipping N
            if (list_elem_root->dna < 4) {
                clv_root = this->Get_clv_for_dna(list_elem_root->dna);
            } else if (list_elem_root->dna > 5) {
                // std::cout << (int) list_elem_root->dna << std::endl;
                if (this->verbose) {assert (this->root->clv_for_list_elem.find(list_elem_root) != this->root->clv_for_list_elem.end()); }
                clv_root = this->root->clv_for_list_elem[list_elem_root];
            }
            if (this->verbose) { assert (list_elem_root->n_pos == 1); }
            site_likelihood = 0;
            for (unsigned char dna = 0; dna < 4; dna ++) {
                site_likelihood += this->pi_rho[dna] * clv_root[dna];
            }
            this->log_likelihood += log(site_likelihood);
        }        
    }
    log_likelihood_contri_from_root_seq = this->log_likelihood - log_likelihood_contri_from_root_seq;
    if (this->verbose) {
        std::cout << "log_likelihood_contri_from_root_seq\t" << std::setprecision(8) << log_likelihood_contri_from_root_seq << std::endl;
    }
}

void tree::Compute_loglikelihood_using_fast_pruning_algorithm() {
    if (this->verbose) {
        std::cout << "------------------------------------------------------------------" << std::endl;
    }

    // Apply pruning algorithm to genome lists
    Matrix4f Q = this->Q;
    Matrix4f Q_scaled_l; Matrix4f Q_scaled_r;
    Matrix4f P_l; Matrix4f P_r;
    float t_l; float t_r;
	float scaling_factor;
    node * left_child; node * right_child;
    unsigned int ch;
	double partial_likelihood_l;
    double partial_likelihood_r;
    double largest_elem;
	this->log_likelihood = 0;
	double site_likelihood;    
    int list_elem_ind;
    genome_list_elem * list_elem_left_child;
    genome_list_elem * list_elem_right_child;    
    genome_list_elem * list_elem_parent;
    int ind_left;
    int ind_right;
    int num_list_elem_left_child;
    int num_list_elem_right_child;
    int end_pos_left; int end_pos_right;
    bool verbose_continue = false;
    bool verbose = false;    
    // if (this->logging) {
        // std::cout << "logging is switched on here" << std::endl;
    // }
    for (node * parent : this->nodes_for_postorder_traversal) {
        if (this->verbose) {
            std::cout << parent->name << std::endl;
            std::cout << "*****" << std::endl;
        }
        if (!parent->leaf) { // iterate over ancestral nodes
            left_child = parent->children[0];
            right_child = parent->children[1];
            ind_left = 0; ind_right = 0;
            num_list_elem_left_child = left_child->genome_list.size();
            num_list_elem_right_child = right_child->genome_list.size();            
            // std::cout << "Updating list elements for " << left_child->name << " and " << right_child->name << std::endl;
            // Add genome wide log scaling factors of children to parent
            if (this->verbose) {
                std::cout << "num_list_elem_left_child is " << num_list_elem_left_child << "\t" << "num_list_elem_right_child " << num_list_elem_right_child << std::endl;
            }
            parent->log_scaling_factor = left_child->log_scaling_factor + right_child->log_scaling_factor;
            if (this->verbose) {
                std::cout << parent->name << "\t" << parent->concatenated_descendant_names << "\n";
                this->log_file << parent->name << "\t" << parent->concatenated_descendant_names << "\n";
            }
            while (ind_left < num_list_elem_left_child && ind_right < num_list_elem_right_child) {
                if (this->verbose) {
                    std::cout << "in while loop" << std::endl;
                }
                // verbose_continue = false;
                // std::cout << "Check 1" << std::endl;
                list_elem_parent = new genome_list_elem();
                // std::cout << "Check 2" << std::endl;
                parent->genome_list.push_back(list_elem_parent);
                // std::cout << "Check 3" << std::endl;
                list_elem_left_child = left_child->genome_list[ind_left];
                // std::cout << "Check 4" << std::endl;
                list_elem_right_child = right_child->genome_list[ind_right];        
                // std::cout << "Check 5" << std::endl;
                if (this->verbose) {
                    std::cout << std::endl << "Before update" << std::endl ;
                    std::cout << "Left child:" ;
                    std::cout << "\tname:\t" << left_child->name;
                    std::cout << "\tdna:\t" << (int) list_elem_left_child->dna;
                    std::cout << "\tstart pos:\t" << (int) list_elem_left_child->start_pos;
                    std::cout << "\tn pos:\t" << (int) list_elem_left_child->n_pos << std::endl;
                    
                    std::cout << "Right child:" ;
                    std::cout << "\tname:\t" << right_child->name;
                    std::cout << "\tdna:\t" << (int) list_elem_right_child->dna;
                    std::cout << "\tstart pos:\t" << (int) list_elem_right_child->start_pos;
                    std::cout << "\tn pos:\t" << (int) list_elem_right_child->n_pos << std::endl;
                }
                // std::cout << "log scaling factor for " << parent->name << "\tbefore update is\t" << parent->log_scaling_factor << std::endl;                    
                this->Update_list_elements(parent, left_child, right_child, list_elem_parent, list_elem_left_child, list_elem_right_child);        
                // std::cout << "log scaling factor for " << parent->name << "\tafter update is\t" << parent->log_scaling_factor << std::endl;    
                if (parent->parent != parent) {
                    // add length of branch from grandparent to parent to list_elem_parent->length
                    list_elem_parent->length += this->Get_directed_edge_length(parent->parent, parent);
                }
                if (this->verbose || list_elem_parent->dna > 16) {
                    std::cout << "After update" << std::endl ;
                    std::cout << "left child:" ;
                    std::cout << "\tname:\t" << left_child->name;
                    std::cout << "\tdna:\t" << (int) list_elem_left_child->dna;
                    std::cout << "\tstart pos:\t" << (int) list_elem_left_child->start_pos;
                    std::cout << "\tn pos:\t" << (int) list_elem_left_child->n_pos << std::endl;
                    
                    
                    std::cout << "right child:" ;
                    std::cout << "\tname:\t" << right_child->name;
                    std::cout << "\tdna:\t" << (int) list_elem_right_child->dna;
                    std::cout << "\tstart pos:\t" << (int) list_elem_right_child->start_pos;
                    std::cout << "\tn pos:\t" << (int) list_elem_right_child->n_pos << std::endl;
                    

                    std::cout << "parent:" ;
                    std::cout << "\tname:\t" << parent->name;
                    std::cout << "\tdna:\t" << (int) list_elem_parent->dna;
                    std::cout << "\tstart pos:\t" << (int) list_elem_parent->start_pos;
                    std::cout << "\tn pos:\t" << (int) list_elem_parent->n_pos << std::endl;
                    assert (list_elem_parent->n_pos > 0);
                }
                // if (verbose) {assert (list_elem_parent->dna < 17);}

                if (list_elem_left_child->n_pos == 0) {
                    if (this->verbose) {
                        std::cout << "========= shift left index =========" << std::endl;
                    }                    
                    ind_left += 1;
                }
                if (list_elem_right_child->n_pos == 0) {
                    if (this->verbose) {
                        std::cout << "========= shift right index ========" << std::endl;
                    }                    
                    ind_right += 1;
                }
            }
            if (logging) {
                if (this->verbose) {
                    std::cout << "logging " << std::endl;                
                }
                this->log_file << "\t" << parent->name << "\t" << parent->log_scaling_factor << std::endl;
                if (this->verbose) {
                    std::cout << "Check 6" << std::endl;
                }
            }
            if (this->verbose) {            
                std::cout << "parent name is " << parent->name << std::endl;
                std::cout << "parent genome list size is " << parent->genome_list.size() << std::endl;
            }
            assert((*(parent->genome_list.end()-1))->start_pos + (*(parent->genome_list.end()-1))->n_pos -1 == GENOME_LENGTH);
            if (this->verbose) {
                std::cout << "ind_left is\t" << ind_left << "\tnum list elem in left child is\t" << num_list_elem_left_child << std::endl;
                std::cout << "ind_right is\t" << ind_right << "\tnum list elem in right child is\t" << num_list_elem_right_child << std::endl;
            }            
            if (this->verbose) {assert (ind_left == num_list_elem_left_child && ind_right == num_list_elem_right_child);}  
            // std::cout << "total log likelihood added to " << parent->name << "\tafter update is\t" << parent->log_scaling_factor << std::endl;              
            // break;
            // collapse contiguous reference type list elements into a single list element
            this->Combine_ref_type_list_elements(parent);
            if (this->verbose) {
                assert((*(parent->genome_list.end()-1))->start_pos + (*(parent->genome_list.end()-1))->n_pos -1 == GENOME_LENGTH);
            }
            // assert(*(parent->genome_list.end()-1)->start_pos + *(parent->genome_list.end()-1)->n_pos -1 = GENOME_LENGTH);
        }        
    }
    // compute total log likelihood score
    if (this->verbose) {
        std::cout << "Computing log-likelihood score for tree using genome list" << std::endl;
    }
    this->Compute_log_likelihood_score_for_tree_using_genome_list();
    if (this->verbose) {
        std::cout << "Completed computing log-likelihood score for tree using genome list" << std::endl;
    }        
}

void tree::Optimize_branch_lengths() {

// P = I+Qt

// try 
// t = -(sum_x (sum_y (n_xy))/(sum_x(n_xx q_xx))

}

void tree::Set_model_parameters() {
    bool empirical = false;
    if (empirical) {
        // pi_r has been set using reference sequence        
        // traverse tree and set UNREST rate matrix as F - I
        // where F is the frequentist estimate of the transition matrix that sums over all branches 
        // count changes  
    }
    else {
        // for (int i = 0; i < 4; i++) {
        //     this->pi_rho[i] = 0.25;
        // }    
        std::vector <float> exchangeability_params;
        exchangeability_params.push_back(0.04); // A<->C beta
        exchangeability_params.push_back(0.3);  // A<->G alpha
        exchangeability_params.push_back(0.1);  // A<->T gamma
        exchangeability_params.push_back(0.02); // C<->G delta
        exchangeability_params.push_back(1.0);  // C<->T eta
        exchangeability_params.push_back(1.0);  // G<->T epsilon

        this->Q(0,1) = this->pi_rho[1] * exchangeability_params[0]; // pi[C] * R[A<->C]
        this->Q(0,2) = this->pi_rho[2] * exchangeability_params[1]; // pi[G] * R[A<->G]
        this->Q(0,3) = this->pi_rho[3] * exchangeability_params[2]; // pi[T] * R[A<->T]
        this->Q(1,0) = this->pi_rho[0] * exchangeability_params[0]; // pi[A] * R[A<->C]
        this->Q(1,2) = this->pi_rho[2] * exchangeability_params[3]; // pi[G] * R[C<->G]
        this->Q(1,3) = this->pi_rho[3] * exchangeability_params[4]; // pi[T] * R[C<->T]
        this->Q(2,0) = this->pi_rho[0] * exchangeability_params[1]; // pi[G] * R[A<->G]
        this->Q(2,1) = this->pi_rho[1] * exchangeability_params[3]; // pi[C] * R[C<->G]
        this->Q(2,3) = this->pi_rho[3] * exchangeability_params[5]; // pi[T] * R[G<->T]
        this->Q(3,0) = this->pi_rho[0] * exchangeability_params[2]; // pi[A] * R[A<->T]
        this->Q(3,1) = this->pi_rho[1] * exchangeability_params[4]; // pi[C] * R[C<->T]
        this->Q(3,2) = this->pi_rho[2] * exchangeability_params[5]; // pi[G] * R[G<->T]

        this->Q(0,0) = - (this->Q(0,1) + this->Q(0,2) + this->Q(0,3));
        this->Q(1,1) = - (this->Q(1,0) + this->Q(1,2) + this->Q(1,3));
        this->Q(2,2) = - (this->Q(2,0) + this->Q(2,1) + this->Q(2,3));
        this->Q(3,3) = - (this->Q(3,0) + this->Q(3,1) + this->Q(3,2));
    }        
    Q /= this->Compute_scaling_factor(Q);        

    // std::cout << "GTR rate matrix is " << std::endl << this->Q << std::endl;
    // Set root and compute vertex list for post order traversal
    // node * l;
    // if (this->verbose) {
    //     if (this->Contains_node("EPI_ISL_1208981")) {
    //             l = this->Get_node("EPI_ISL_1208981");
    //         } else {
    //             l = this->leaves[0];
    //         }
    // }
    // std::cout << "node for rooting is " << l->name << std::endl;
    // this->Root_tree_along_edge(l, l->neighbors[0], 0.0);
    // this->Set_nodes_for_postorder_traversal();
}

float tree::Compute_scaling_factor(Matrix4f Q) {
    MatrixXf Q_aug = ArrayXXf::Zero(4,5);
	for (int row = 0; row < 4; row++){
		for (int col = 0; col < 4; col++){
			Q_aug(row, col) = Q(row, col);
		}
	}
	for (int row = 0; row < 4; row++){
		Q_aug(row, 4) = 1;
	}	
	MatrixXf b = ArrayXXf::Zero(5,1);
	for (int row = 0; row < 4; row++){
		b(row,0) = 0;
	}
	b(4,0) = 1;	
	MatrixXf pi = ArrayXXf::Zero(1,4);
	pi = Q_aug.transpose().colPivHouseholderQr().solve(b).transpose();
	float scalingFactor = 0;
	for (int i = 0; i < 4; i++){
		scalingFactor -= pi(0,i) * Q(i,i);
	}
	return scalingFactor;
}

void tree::Populate_DNA_to_char_map() {
    // IUPAC code
    // std::cout << "populating dna to char map" << std::endl;
    // std::cout << this->DNA_to_char.size() << std::endl;
    this->DNA_to_char["a"] = 0;    
    this->DNA_to_char["c"] = 1;
    this->DNA_to_char["g"] = 2;
    this->DNA_to_char["t"] = 3;
    this->DNA_to_char["x"] = 4; // x stands for reference nucleobase
    this->DNA_to_char["-"] = 5;
    this->DNA_to_char["n"] = 5;
    this->DNA_to_char["w"] = 6;
    this->DNA_to_char["s"] = 7;
    this->DNA_to_char["m"] = 8;
    this->DNA_to_char["k"] = 9;
    this->DNA_to_char["r"] = 10;
    this->DNA_to_char["y"] = 11;
    this->DNA_to_char["b"] = 12;
    this->DNA_to_char["d"] = 13;
    this->DNA_to_char["h"] = 14;
    this->DNA_to_char["v"] = 15;
    this->DNA_to_char["o"] = 16; // o stands for ambiguous type    
}

void tree::Add_node(std::string u_name) {
    node * u = new node(u_name);
    this->node_list.insert(std::pair<std::string,node *>(u_name,u));
    if (this->verbose) {
        std::cout << "Added node " << u->name << std::endl;
    }
}

node * tree::Get_node(std::string u_name) {
    bool node_present = this->Contains_node(u_name);
    if (!node_present) {
        if (this->verbose) {
            std::cout << "node " << u_name << " is missing " << std::endl;
        }
    }    
    if (verbose) {assert(node_present);}
    return(this->node_list[u_name]);
}

bool tree::Contains_node(std::string u_name) {
    bool node_present;
    std::
    map <std::string, node *>::iterator it = this->node_list.find(u_name);
    if (it != this->node_list.end()) {
        node_present = true;
    } else {
        node_present = false;
    }
    return (node_present);
}

void tree::Add_undirected_edge(node * u, node * v, float length) {    
    u->neighbors.push_back(v);
    v->neighbors.push_back(u);
    u->degree += 1;
    v->degree += 1;
    if (u < v) {
        this->undirected_edge_length_map.insert(std::make_pair(std::make_pair(u,v),length));
    } else {
        this->undirected_edge_length_map.insert(std::make_pair(std::make_pair(v,u),length));
    }
}

float tree::Get_undirected_edge_length(node * u, node * v) {
    float length;
    if (u < v) {
        length = this->undirected_edge_length_map[std::make_pair(u,v)];
    } else {
        length = this->undirected_edge_length_map[std::make_pair(v,u)];
    }
    return (length);
}

float tree::Get_directed_edge_length(node * p, node * c) {
    return (this->directed_edge_length_map[std::make_pair(p,c)]);
}

void tree::Read_newick_file() {
    // default value for rooted is true
    this->leaves.clear();
    if (this->verbose) {
        std::cout << "Tree file name is " << tree_file_name << std::endl;    
    }
    std::string node_name; std::string h_name;
    node * h; node * n; node * l;
    float length; float v_length;
    std::string newick_string;
    std::string sibling_string;
    std::string sibling_string_without_parenthesis;
    std::string newick_string_split;
    std::vector <std::string> split_sibling_string;
    std::vector <std::string> node_name_and_length;
    std::ifstream inputFile(this->tree_file_name.c_str());
    getline(inputFile, newick_string);    
    if (this->verbose) {
        std::cout << "Newick string is " << newick_string << std::endl;
    }
    std::regex sibling_pattern ("\\([^\\(\\)]+\\)");
    std::smatch cherry_matches;
    bool continue_search = true;
    float min_length = pow(10,-8);
    while (continue_search) {
        if (regex_search(newick_string, cherry_matches, sibling_pattern)) {
            h_name = "h_" + std::to_string(this->h_ind);
            this->Add_node(h_name);
            h = this->Get_node(h_name);
            this->h_ind++;
            sibling_string = cherry_matches[0];
            sibling_string_without_parenthesis = "";
            for (int i = 1; i < sibling_string.size() -1; i++) {
                sibling_string_without_parenthesis += sibling_string[i];
            }
            boost::split(split_sibling_string, sibling_string_without_parenthesis, [](char c){return c == ',';});
            if (split_sibling_string.size() == 2 || split_sibling_string.size() == 3) {
                for (std::string sibling_string : split_sibling_string) {
                    boost::split(node_name_and_length, sibling_string, [](char c){return c == ':';});
                    node_name = node_name_and_length[0];
                    length = std::stof(node_name_and_length[1]);
                    if (length < min_length) {
                        length = min_length;
                    }
                    if (!this->Contains_node(node_name)) {
                        this->Add_node(node_name);                        
                        this->leaves.push_back(this->Get_node(node_name));
                    }
                    n = this->Get_node(node_name);
                    this->Add_directed_edge(h, n, length);
                    this->Add_undirected_edge(n, h, length);
                }
            } else {
                std::cerr << "sibling string has not been properly parsed" << std::endl;
                // exit(-1);
            }
            newick_string = std::regex_replace(newick_string, sibling_pattern, h_name, std::regex_constants::format_first_only);
        } else {
            continue_search = false;
        }        
    }
    unsigned int num_leaves = this->leaves.size();
    unsigned int num_nodes = this->node_list.size();
    if (this->verbose) {
        std::cout << "Number of edges is " << this->undirected_edge_length_map.size() << std::endl;
        std::cout << "Number of leaves is " << num_leaves << std::endl;
        std::cout << "Number of nodes is " << num_nodes << std::endl;    
    }    
    if (num_nodes == 2*(num_leaves) -2) {
        this->rooted = false;
        // std::cout << "Input tree is unrooted" << std::endl;
    } else {
        assert(num_nodes == 2*(num_leaves) -1);
        this->rooted = true;
        // std::cout << "Input tree is rooted" << std::endl;
    }
}

// void tree::WriteNewickFile(string tree_file_name, bool rooted) {
//     // Hardcoded rooted to false
//     rooted = false;
//     string newick_string = "";    

// }

void tree::Read_reference_sequence() {
    unsigned char c;
    this->reference_sequence.clear();	  
	std::ifstream inputFile(this->reference_file_name.c_str());
	std::string seq_name;
	std::string seq = "";    
	for (std::string line; getline(inputFile, line);) {
		if (line[0]=='>') {
			seq_name = line.substr(1,line.length());
		}
		else {
		    seq += line;
		}
	}
    inputFile.close();
    // set pi_rho    
    std::array <int, 4> nuc_counts;
    for (int i = 0; i < 4; i++) {
        nuc_counts[i] = 0;
        this->pi_rho[i] = 0.0;
    }

	for (char const dna: seq) {
        c = this->DNA_to_char[std::string(1,tolower(dna))];
		this->reference_sequence.push_back(c);
        if (c < 4) {
            nuc_counts[c] += 1;
        }
	}
    float total_counts = 0.0;
    for (int i = 0; i < 4; i++) {
        total_counts += (float) nuc_counts[i];
    }
    for (int i = 0; i < 4; i ++) {
        this->pi_rho[i] = (float) nuc_counts[i] / total_counts;
        // std::cout << "pi[" << i << "] is " << this->pi_rho[i] << "\t";
    }
    // std::cout << std::endl;
}

void tree::Add_fasta_sequences() {
    std::vector <unsigned char> recoded_sequence;
    node * v;
	recoded_sequence.clear();
	unsigned int site = 0;
	std::ifstream inputFile(this->alignment_file_name.c_str());
	std::string seq_name;
	std::string seq = "";	
	for(std::string line; getline(inputFile, line);) {        
		if (line[0] == '>') {
			if (seq != "") {
				for (char const dna: seq) {                    
					recoded_sequence.push_back(this->DNA_to_char[std::string(1,tolower(dna))]);
					site += 1;
					}                
                v = this->Get_node(seq_name);
                v->complete_sequence = recoded_sequence;                
				recoded_sequence.clear();
			}
			seq_name = line.substr(1,line.length());
			seq = "";
			site = 0;
		}
		else {
		    seq += line ;
		}
	}
	for (char const dna: seq) {
		recoded_sequence.push_back(this->DNA_to_char[std::string(1,tolower(dna))]);
		site += 1;
	}
    std::cout << "All but last sequences have been added" << std::endl;
    std::cout << "Adding node " << seq_name << std::endl;
	v = this->Get_node(seq_name);
    v->complete_sequence = recoded_sequence;
	recoded_sequence.clear();
	inputFile.close();
}

void tree::Add_ref_nuc_counts_based_on_genome_coordinates() {    
    unsigned int end_pos;
    this->cum_ref_nuc_counts_map.clear();
    std::array <unsigned int, 4> cum_nuc_counts;
    cum_nuc_counts[0] = 0; cum_nuc_counts[1] = 0; cum_nuc_counts[2] = 0; cum_nuc_counts[3] = 0;
    for (node * leaf : this->leaves) {
        for (genome_list_elem * list_elem : leaf->genome_list) {            
            end_pos = list_elem->start_pos + list_elem->n_pos - 1;            
            if (this->cum_ref_nuc_counts_map.find(end_pos) == this->cum_ref_nuc_counts_map.end()) {
                // std::cout << "position added in cum ref count map is " << end_pos << std::endl;
                this->cum_ref_nuc_counts_map[end_pos] = cum_nuc_counts;
            }
        }
    }
    // std::cout << "initial sum of cum counts is " << cum_nuc_counts[0] + cum_nuc_counts[1] + cum_nuc_counts[2] + cum_nuc_counts[3] << std::endl;
    // std::cout << "###########################################" << std::endl;
    unsigned int c;
    int counter = 0;    
    for (unsigned int pos = 1; pos < GENOME_LENGTH + 1; pos++) {
        counter ++;
        c = this->reference_sequence[pos -1];        
        if (c > 4 || c < 0) {
            std::cout << (int) c << std::endl;
            std::cout << pos << std::endl;
        }
        assert (c < 4 && c >= 0);
        cum_nuc_counts[c] ++;
        // cum_nuc_counts[c] += 1;        
        // std::cout << "cum counts for pos " << pos << " is " << cum_nuc_counts[0] + cum_nuc_counts[1] + cum_nuc_counts[2] + cum_nuc_counts[3] << " and character is " << (int) c << std::endl;
        // for (int i = 0; i < 4; i++) {
        //     std::cout << "cum counts for " << i << " is " << (int) cum_nuc_counts[i] << std::endl;
        // }
        if (this->cum_ref_nuc_counts_map.find(pos) != this->cum_ref_nuc_counts_map.end()) {
            // std::cout << "cum counts for pos " << pos << " is " << cum_nuc_counts[0] + cum_nuc_counts[1] + cum_nuc_counts[2] + cum_nuc_counts[3] << std::endl;
            for (unsigned int dna = 0; dna < 4; dna++) {
                this->cum_ref_nuc_counts_map[pos][dna] = cum_nuc_counts[dna];
            }            
        }        
    }
    for (unsigned char dna = 0; dna < 4; dna++) {
        this->pi_rho[dna] = cum_nuc_counts[dna]/float(GENOME_LENGTH);
    }    
    // std::cout << "###########################################" << std::endl;
    cum_nuc_counts[0] = 0; cum_nuc_counts[1] = 0; cum_nuc_counts[2] = 0; cum_nuc_counts[3] = 0;
    // for the case where a stretch of ref nucs starts at 1
    this->cum_ref_nuc_counts_map[0] = cum_nuc_counts;
    // std::cout << "size of cum ref nuc counts map is " << this->cum_ref_nuc_counts_map.size() << std::endl;
    std::array <unsigned int, 4> previous_count;
    previous_count[0] = 0; previous_count[1] = 0; previous_count[2] = 0; previous_count[3] = 0;
    int previous_total_count = 0;
    int current_total_count = 0;
    int pos = 0;
    for (std::pair<unsigned int, std::array<unsigned int, 4UL>> ref_count_map_elem : this->cum_ref_nuc_counts_map) {        
        pos = ref_count_map_elem.first;
        // std::cout << "Position is " << pos << std::endl;        
        // for (unsigned int dna = 0; dna < 4; dna ++) {
            // std::cout << "Number of counts for dna " << (int) dna << " is " << (int) ref_count_map_elem.second[dna] << std::endl;
            // current_total_count += ref_count_map_elem.second[dna];
        // }
        // if (pos > 0) {
        //     std::cout << "total previous count is " << previous_total_count << std::endl;
        //     std::cout << "total current count is " << current_total_count << std::endl;
        //     std::cout << "diff is " << current_total_count - previous_total_count << std::endl;
        //     assert(current_total_count - previous_total_count > 0);
        // }
        // previous_total_count = current_total_count;
    }
// assert that consecutive
}

void tree::Add_mut_diff_sequences() {
    // vector <unsigned char> recoded_sequence;
    node * n;
	// recoded_sequence.clear();
	// unsigned int site = 0;
	std::ifstream inputFile(this->mut_diff_file_name.c_str());
	std::string seq_name = "";
	// string seq = "";
    std::vector <std::string> split_line;
    split_line.reserve(3);
    std::string dna;
    std::array <double, 4> clv;
    genome_list_elem * list_elem_diff; // Read from mut diff file
    genome_list_elem * list_elem_ref; // List element for positions identical to reference sequence
    genome_list_elem * list_elem_prev_diff;
    unsigned int genome_list_size;    
    // std::cout << "Reading mutdiff file" << std::endl;
    // this->verbose = true;
	for (std::string line; getline(inputFile, line );) {
		if (line[0]=='>') {
            if (seq_name != "") {                              
                n->Set_pos_for_reference_characters();
            }
            seq_name = line.substr(1,line.length());
            if (this->verbose) {
                std::cout << "seq_name\t" << seq_name << std::endl;
            }
            n = this->Get_node(seq_name);
            n->initial_edge_length = this->Get_directed_edge_length(n->parent,n);                        
		}
		else {
            boost::split(split_line, line, [](char c){return c == '\t';});
            list_elem_diff = new genome_list_elem();
            dna = split_line[0];
            list_elem_diff->dna = this->DNA_to_char[dna];
            list_elem_diff->length = n->initial_edge_length;
            // assert(list_elem_diff->length > 0);
		    if (split_line.size() == 3) {
                if (verbose) {assert (split_line[0] == "n" || split_line[0] == "-" );}
                list_elem_diff->start_pos = stoi(split_line[1]);
                list_elem_diff->n_pos = stoi(split_line[2]);
            } else {
                if (verbose) {assert (split_line.size() == 2);}
                list_elem_diff->start_pos = stoi(split_line[1]);
                list_elem_diff->n_pos = 1;
            }
            split_line.clear();
            if (list_elem_diff->dna > 5) { // ambiguous IUPAC code and not N or -
                clv = Get_clv_for_dna(list_elem_diff->dna);
                n->clv_for_list_elem[list_elem_diff] = clv;
            }
            // if (n->name == "EPI_ISL_1589845") {
            //     std::cout << "dna is " << (int) list_elem_diff->dna << std::endl;
            //     std::cout << "start pos is " << list_elem_diff->start_pos << std::endl;
            //     std::cout << "n pos is " << list_elem_diff->n_pos << std::endl;
            // }
            if (n->empty_genome_list && list_elem_diff->start_pos > 1) {
                list_elem_ref = new genome_list_elem();
                list_elem_ref->dna = 4; // 4 is the code for reference sequence
                list_elem_ref->length = n->initial_edge_length;
                n->genome_list.push_back(list_elem_ref);
                n->empty_genome_list = false;
            }
            // erase ref element that is sandwiched between contiguous diff elements            
            genome_list_size = n->genome_list.size();
            if (genome_list_size > 1) {
                list_elem_prev_diff =  n->genome_list[genome_list_size - 2];                
                // std::cout << (int) list_elem_prev_diff->dna << "\t" << list_elem_prev_diff->n_pos << std::endl;
                if (verbose) {assert (list_elem_prev_diff->dna != 4);}
                if (list_elem_prev_diff->start_pos + list_elem_prev_diff->n_pos == list_elem_diff->start_pos) {
                    list_elem_ref = n->genome_list[genome_list_size-1];
                    if (verbose) {assert(list_elem_ref->dna == 4);}
                    n->genome_list.pop_back();
                    delete list_elem_ref;                    
                }
            }

            n->genome_list.push_back(list_elem_diff);
            n->empty_genome_list = false;
            if (list_elem_diff->start_pos + list_elem_diff->n_pos - 1 < GENOME_LENGTH) {
                list_elem_ref = new genome_list_elem();
                list_elem_ref->dna = 4; // 4 is the code for reference sequence
                list_elem_ref->length = this->Get_directed_edge_length(n->parent, n);
                n->genome_list.push_back(list_elem_ref);            
            }
            // Add a list element comprising positions identical to the reference sequence
            // add start pos and end pos based on pos in neighboring elems (do this after populating the entire list of elems)
		}
	}    
    n->Set_pos_for_reference_characters();
	inputFile.close();
    // this->verbose = true;    
}
void tree::Compress_sequences() {
    
    this->character_pattern_weights.clear();
    std::map <std::vector<unsigned char>,int> patterns_to_weights;
	int number_of_sites = this->leaves[0]->complete_sequence.size();
    int number_of_leaves = this->leaves.size();
    std::cout << "number of sites is " << number_of_sites << std::endl;
	std::vector <unsigned char> character_pattern;
    for (node * l : this->leaves) {
        l->compressed_sequence.clear();
    }
	for (int site = 0; site < number_of_sites; site++) {
		character_pattern.clear();
		for (int leaf_ind = 0; leaf_ind < number_of_leaves; leaf_ind ++) {
			character_pattern.push_back(this->leaves[leaf_ind]->complete_sequence[site]);
            }
		if (patterns_to_weights.find(character_pattern) != patterns_to_weights.end()) {
			patterns_to_weights[character_pattern] += 1;
		} else {
            patterns_to_weights[character_pattern] = 1;                        
			for (int leaf_ind = 0; leaf_ind < number_of_leaves; leaf_ind ++) {
                this->leaves[leaf_ind]->compressed_sequence.push_back(character_pattern[leaf_ind]);
			}
		}
	}    
    this->num_char_patterns = this->leaves[0]->compressed_sequence.size();    
    std::cout << "Number of character patterns is " << this->num_char_patterns << std::endl;
    for (int site = 0; site < this->num_char_patterns; site++) {
		character_pattern.clear();
		for (int leaf_ind = 0; leaf_ind < number_of_leaves; leaf_ind ++) {
            // std::cout << (int)this->leaves[leaf_ind]->compressed_sequence[site];
			character_pattern.push_back(this->leaves[leaf_ind]->compressed_sequence[site]);
            }
        // std::cout << std::endl;
        this->character_pattern_weights.push_back(patterns_to_weights[character_pattern]);
    }    
}

class fastLK_overview {
public:    
    tree * T;
    void Run_workflow(std::string workflow_type);
    fastLK_overview(std::string path_to_reference_file, std::string path_to_mut_diff_file, std::string path_to_sequence_alignment_file, std::string path_to_tree_file, std::string path_to_log_file, std::string workflow_type, float clv_threshold, bool verbose_flag_bool){        
        std::chrono::system_clock::time_point start_time = std::chrono::high_resolution_clock::now();		        
        this->T = new tree;
        this->T->tree_file_name = path_to_tree_file;
        this->T->alignment_file_name = path_to_sequence_alignment_file;
        this->T->reference_file_name = path_to_reference_file;
        this->T->mut_diff_file_name = path_to_mut_diff_file;
        if (path_to_log_file != "") {
            this->T->log_file_name = path_to_log_file;
        } else {
            this->T->log_file_name = path_to_mut_diff_file + ".log";
        }
        this->T->log_file.open(this->T->log_file_name);
        this->T->clv_threshold = clv_threshold;
        if (workflow_type == "") {
            this->T->workflow_type = "mut_diff";
        } else {
            this->T->workflow_type = workflow_type;
        }
        if (verbose_flag_bool) {
            this->T->verbose = true;
        } else {
            this->T->verbose = false;
        }
        // workflow_type = "verbose";
        //workflow_type = "standard"; 
        this->Run_workflow(this->T->workflow_type);
        // std::cout << "workflow type is " << this->T->workflow_type << std::endl;
        std::chrono::system_clock::time_point end_time = std::chrono::high_resolution_clock::now();
        std::cout << "Total CPU time used is " << std::chrono::duration_cast<std::chrono::seconds>(end_time-start_time).count() << " second(s)\n";		
        this->T->log_file << "Total CPU time used is " << std::chrono::duration_cast<std::chrono::seconds>(end_time-start_time).count() << " second(s)\n";
        this->T->log_file.close();
    }

    ~fastLK_overview(){
        delete this->T;
    }
};

void fastLK_overview::Run_workflow(std::string workflow_type){
    // Read tree file
    // Add sequences (perform global site pattern compression)        
    this->T->Read_newick_file();
    // std::cout << "Setting root" << std::endl;
    this->T->Set_root();
    // std::cout << "Setting leaves" << std::endl;
    this->T->Set_leaves();
    // std::cout << "Setting nodes for postorder traversal" << std::endl;
    this->T->Set_nodes_for_postorder_traversal();    
    // std::cout << "Reading reference sequence" << std::endl;
    this->T->Read_reference_sequence();
    if (!this->T->rooted) {
        if (this->T->verbose) {
            std::cout << "Input tree is unrooted" << std::endl;
        }
        this->T->Root_unrooted_tree();
    }    
    // std::cout << "Setting model parameters" << std::endl;         
    this->T->Set_model_parameters();

    if (workflow_type == "standard") {
        if (this->T->verbose) {
            std::cout << "Adding fasta sequences" << std::endl;
        }
        this->T->Add_fasta_sequences();
        if (this->T->verbose) {       
            std::cout << "Compressing sequences" << std::endl; 
        }
        this->T->Compress_sequences();  
        if (this->T->verbose) {         
            std::cout << "Computing log-likelihood using standard approach" << std::endl;      
        }
        this->T->Compute_loglikelihood_using_standard_pruning_algorithm();        
        std::cout << "log likelihood score computed using standard pruning algorithm is " << std::setprecision(8) << this->T->log_likelihood << std::endl;
    } else if (workflow_type == "mut_diff") {
        // std::cout << "Setting descendant names" << std::endl;
        // this->T->Set_descendant_names();
        // std::cout << "Adding mutdiff sequences" << std::endl;        
        this->T->Add_mut_diff_sequences();
        // for (node * n : this->T->leaves) {
        //     if (n->clv_for_list_elem.size() > 0) {
        //         std::cout << n->name << std::endl;
        //     }
        // }
        // std::cout << "Adding ref nuc counts" << std::endl;
        this->T->Add_ref_nuc_counts_based_on_genome_coordinates();       
        // std::cout << "Computing log-likelihood using fast pruning algorithm" << std::endl;
        this->T->Compute_loglikelihood_using_fast_pruning_algorithm();
        // cout << "Completed computing log-likelihood score" << endl;
        cout << "log likelihood score computed using fast pruning algorithm is " << setprecision(8) << this->T->log_likelihood << endl;
        cout << "Optimizing branch lengths of fully labeled tree" << endl;
        this->T->OptimizeBranchLengths();
    
    } else if (workflow_type == "verbose") {        
        this->T->verbose_combining_genome_ref_elements();
    }    
    }

    // ************************************************ //
    // Compute log likelihood for GTR rate matrix using //
    // standard Felsenstein's pruning algorithm         //
    // ************************************************ //

    // std::cout << "Log likelihood score computing using standard pruning algorithm is " << std::setprecisionecision(8) << this->T->log_likelihood << std::endl;


    // ************************************************ //
    // Compute log likelihood for GTR rate matrix using //
    // approximate Felsenstein's pruning algorithm      //
    // ************************************************ //

#endif