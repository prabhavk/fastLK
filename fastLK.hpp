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
// #include <boost/algorithm/string/replace.hpp>
using namespace Eigen;
using namespace std;

// global variables

const unsigned short RESERVE_genome_list_elem = 200;
const unsigned short GENOME_LENGTH = 29891;


struct genome_list_elem {
    unsigned char dna;
    unsigned int start_pos;
    unsigned int n_pos;
    float length = 0;
    // array <double, 4> clv;
    // double log_scaling_factor;
    genome_list_elem() {

    }
    ~genome_list_elem() {

    }

};


class node {
public:
    vector < node * > neighbors;
    vector < node * > children;
    bool debug = true;     
    int degree = 0;
    int in_degree = 0;
    int out_degree = 0;
    int times_visited = 0;
    string concated_descendant_tip_names = "";
    bool leaf = false;    
    node * parent = this;
    array <double, 4> clv;
    double log_scaling_factor = 0;
    string name;
    vector <unsigned char> complete_sequence;
    vector <unsigned char> compressed_sequence;

    vector <genome_list_elem *> genome_list;        
    map <genome_list_elem *, array <double, 4>> clv_for_list_elem;
    
    void Add_parent(node * p);
    void Add_child(node * c);
    void Set_pos_for_reference_characters();
    node (string node_name) {        
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
        if (debug) {assert (list_elem->n_pos > 0);}
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
    double log_likelihood = 0.0;
    double clv_threshold;
    int h_ind = 0;
    bool debug = true;    
    vector <unsigned char> reference_sequence;
    int num_char_patterns;
    string tree_file_name;
    string alignment_file_name;
    string mut_diff_file_name;
    string reference_file_name;
    void Read_newick_file(bool rooted);
    // void WriteNewickFile(string tree_file_name, bool rooted);
    map <string, unsigned char> DNA_to_char;
    map <string, node*> node_list;
    vector <node *> leaves;
    void Set_leaves();
    void Add_fasta_sequences();
    void Add_mut_diff_sequences();
    void Populate_DNA_to_char_map();
    void Compress_sequences();
    void Add_ref_nuc_counts_based_on_genome_coordinates();
    map <unsigned short, array <unsigned char, 4>> cum_ref_nuc_counts_map; // store interval counts of a, c, g, t
    vector <int> character_pattern_weights;
    void Set_model_parameters();
    float Compute_scaling_factor(Matrix4f Q);    
    void Read_reference_sequence();
    void Compute_loglikelihood_using_standard_pruning_algorithm();    
    void Compute_loglikelihood_using_fast_pruning_algorithm();
    array <double, 4> Get_clv_for_dna(unsigned char dna_obs);
    void Root_tree_along_edge(node * u, node * v, float dist_from_u); 
    void Set_nodes_for_postorder_traversal();
    void Add_node(string node_name);
    void Reset_times_visited();
    node * Get_node(string node_name);
    bool Contains_node(string node_name);
    void Add_undirected_edge(node * u, node * v, float length);
    void Add_directed_edge(node * u, node * v, float length);
    vector <node *> nodes_for_postorder_traversal;
    float Get_undirected_edge_length(node * u, node * v);
    float Get_directed_edge_length(node * p, node * c);
    void Reset_log_scaling_factors_and_clvs();
    map <pair <node *, node *>, float> undirected_edge_length_map;
    map <pair <node *, node *>, float> directed_edge_length_map;
    Matrix4f Q_GTR;
    void Update_list_elements(node* parent, node * left_child, node * right_child, genome_list_elem * list_elem_parent, genome_list_elem * list_elem_left_child, genome_list_elem * list_elem_right_child);
    void Combine_ref_type_list_elements(node * n);
    void Compute_log_likelihood_score_for_tree_using_genome_list();
    array <float, 4> pi_rho;    

    // Matrix4f Q_UNREST;

    tree (){
        // Only hardcoded value should be set here
        this->Populate_DNA_to_char_map();
        root = new node("h_root");
        this->node_list["h_root"] = root;         
    }

    ~tree(){
        for (pair <string, node *> elem : this->node_list) {
            delete elem.second;
        }

    }
};

void tree::Reset_log_scaling_factors_and_clvs() {
    for (pair <string, node *> elem: this->node_list) {
        elem.second->clv[0] = 0; elem.second->clv[1] = 0; elem.second->clv[2] = 0; elem.second->clv[3] = 0;
        elem.second->log_scaling_factor = 0;
    }
}

array <double, 4> tree::Get_clv_for_dna(unsigned char dna_obs) {
    array <double, 4> clv;    
    // if (debug) {assert (dna_obs < 4); }    
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
    if (debug) {assert(*max_element(clv.begin(),clv.end()) > 0.0);}    
    return (clv);
}

void tree::Reset_times_visited() {
    for (pair <string, node *> node_elem: this->node_list) {
        node_elem.second->times_visited = 0;
    }
}

void tree::Set_nodes_for_postorder_traversal() { 
    this->nodes_for_postorder_traversal.clear();
    for (pair <string, node *> elem : this->node_list) {
        elem.second->times_visited = 0;
    }
    vector <node *> nodes_to_visit;
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
}

void tree::Add_directed_edge(node * p, node * c, float edge_length) {
    c->Add_parent(p);
    p->Add_child(c);
    this->directed_edge_length_map[make_pair(p,c)] = edge_length;
}


void tree::Root_tree_along_edge(node * u, node * v, float dist_from_u) {
    pair <node *, node *> edge;
    if (u < v) {
        edge = make_pair(u,v);
    } else {
        edge = make_pair(v,u);
    }
    // cout << "u is a leaf with degree " << u->degree << endl;
    // cout << "v is not a leaf with degree " << v->degree << endl;
    if (debug) {assert(this->undirected_edge_length_map.find(edge) != this->undirected_edge_length_map.end());}    
    vector <node *> nodes_visited;
    vector <node *> nodes_to_visit; 
    nodes_to_visit.push_back(u);
    nodes_to_visit.push_back(v);
    nodes_visited.push_back(u);
    nodes_visited.push_back(v);
    float edge_length;
    edge_length = this->Get_undirected_edge_length(u,v);
    if (debug) {assert(edge_length >= dist_from_u);}

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
    Matrix4f Q = this->Q_GTR;
    Matrix4f Q_scaled_l; Matrix4f Q_scaled_r;
    Matrix4f P_l; Matrix4f P_r;
    float t_l; float t_r;
	float scaling_factor;
    node * c_l; node * c_r;
    unsigned char ch;    
	double partial_likelihood_l;
    double partial_likelihood_r;
    double largest_elem;
	this->log_likelihood = 0;
	double site_likelihood;
    unsigned int num_char_patterns = this->character_pattern_weights.size();
	for (unsigned int site = 0; site < num_char_patterns; site++) {
    // for (unsigned int site = 0; site < 1; site++) {
        this->Reset_log_scaling_factors_and_clvs();
        for (node * n : this->nodes_for_postorder_traversal) {
            if (n->leaf) {
                // set conditional likelihood vector using observed character                
                // cout << "node name is\t" << n->name << endl;
                // cout << " dna for " << " site " << site << " is " << (int) n->compressed_sequence[site] << endl;
                n->clv = this->Get_clv_for_dna(n->compressed_sequence[site]);
            } else {
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
		this->log_likelihood += (this->root->log_scaling_factor + log(site_likelihood)) * this->character_pattern_weights[site];
        this->root->log_scaling_factor = 0;
	}
}

void tree::Set_leaves() {
    this->leaves.clear();
    for (pair <string, node *> elem : this->node_list) {
        if (elem.second->degree == 1) {
            this->leaves.push_back(elem.second);
            elem.second->leaf = true;
        }
    }
    cout << "Number of leaves is " << this->leaves.size() << endl;
}

void tree::Combine_ref_type_list_elements(node * n) {
    // collapse contiguous list elements of type reference
    int num_list_elem;
    bool continue_search_for_elements_to_combine = true;
    genome_list_elem * list_elem_i;
    genome_list_elem * list_elem_i_plus_one;
    vector <genome_list_elem *>::iterator start_pos_iterator = n->genome_list.begin();
    int ind_of_elem_extended;
    if (debug) {
        for (int i = 0; i < n->genome_list.size() - 1; i++) {
            assert(n->genome_list[i]->start_pos + n->genome_list[i]->n_pos == n->genome_list[i+1]->start_pos);
        }
    }
    // cout << "Combining list elements for " << n->name << endl;
    // cout << "num elements before combining: " << n->genome_list.size() << endl;
    while (continue_search_for_elements_to_combine) {
        continue_search_for_elements_to_combine = false;
        num_list_elem = n->genome_list.size();
        for (int i = 0; i < num_list_elem - 1; i++) {
            list_elem_i = n->genome_list[i];
            list_elem_i_plus_one = n->genome_list[i+1];
            if (list_elem_i->dna == 4 && list_elem_i_plus_one->dna == 4) {
                continue_search_for_elements_to_combine = true;
                ind_of_elem_extended = i;
            }
        }
        if (continue_search_for_elements_to_combine) {
            // cout << "elements to combine found" << endl;
            list_elem_i = n->genome_list[ind_of_elem_extended];
            list_elem_i_plus_one = n->genome_list[ind_of_elem_extended + 1];
            if (debug) {assert(list_elem_i->dna == 4 && list_elem_i_plus_one->dna == 4);}
            list_elem_i->n_pos += list_elem_i_plus_one->n_pos;
            n->genome_list.erase(start_pos_iterator + ind_of_elem_extended + 1);
            delete list_elem_i_plus_one;
        }
    }
    // cout << "num elements after combining:  " << n->genome_list.size() << endl;;
    // cout << "===========================================================" << endl;
    if (debug) {
        for (int i = 0; i < n->genome_list.size() - 1; i++) {
            assert(n->genome_list[i]->start_pos + n->genome_list[i]->n_pos == n->genome_list[i+1]->start_pos);
        }
    }
}

void tree::Update_list_elements(node * parent, node * left_child, node * right_child, genome_list_elem * list_elem_parent, genome_list_elem * list_elem_left_child, genome_list_elem * list_elem_right_child) {
    // dna_l and dna_r are characters observed for list_elem_left and list_elem_right        
    int n_pos;
    if (debug) {assert (list_elem_left_child->start_pos == list_elem_right_child->start_pos);}
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
    array <unsigned char, 4> ref_nuc_counts;
    array <double, 4> clv_parent;
    array <double, 4> clv_left_child;
    array <double, 4> clv_right_child;
    double largest_elem; double partial_likelihood_l; double partial_likelihood_r;
    unsigned char index_largest_elem;
    int num_clv_elems_larger_than_threshold;
    // if (list_elem_parent->start_pos == 6730) {
    //     cout << " ************************** ";
    //     cout << "Updating list element for " << parent->name << endl;
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
                if (debug) {assert(*max_element(parent->clv_for_list_elem[list_elem_parent].begin(),parent->clv_for_list_elem[list_elem_parent].end()) > 0.0);}
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
                if (debug) {assert(*max_element(parent->clv_for_list_elem[list_elem_parent].begin(),parent->clv_for_list_elem[list_elem_parent].end()) > 0.0);}
            }
        } else {
            // dna_l is not N and dna_r is not N            
            if (dna_l == dna_r && dna_l < 5) {
                // dna_l = dna_r and dna_l is in {a,c,g,t,x} (x is ref)            
                if (dna_l < 4) {
                    // 4) dna_l and dna_r are in {a,c,t,g}. n_pos is one
                    // if (n_pos != 1) {
                    //     cout << "ref pos is " << list_elem_parent->start_pos << endl;
                    //     cout << "n pos is " << list_elem_parent->n_pos << endl;                        
                    //     cout << "dna_l is " << (int) dna_l << "\t" << "dna_r is " << (int) dna_r << endl;
                    //     cout << "left child name is " << left_child->name << endl;
                    //     cout << "right child name is " << right_child->name << endl;
                    //     // cout << "clv left child " << dna_l << " is " << clv_left_child[dna] << "\t";
                    //     // cout << "clv right child " << dna_r << " is " << clv_right_child[dna] << "\t";
                    //     // cout << "clv parent " << dna << " is " << clv_parent[dna] << endl;
                    // }
                    if (debug) {assert(n_pos == 1);}
                    // set dna of parent as dna of left child
                    list_elem_parent->dna = dna_l;
                    // set branch length of parent as 0 (default)
                    // compute likelihood as P(i,i,l) (probability that i remains i after time l)
                    // as P(i,i,l) = 1 + Q(i,i)l (because l is very small)
                    // because we are storing log-likelihood
                    // we can compute log (P(i,i,l)) ~ log (1 + Q(i,i)l) ~ Q(i,i)l
                    // thus contribution towards log_likelihood is Q(i,i)(l_c + l_r)
                    parent->log_scaling_factor += this->Q_GTR(dna_l,dna_l) * (list_elem_left_child->length + list_elem_right_child->length);
                } else if (dna_l == 4) {
                    // 5) dna_l and dna_r are reference nucleobases. n_pos is probably large
                    // cout << "log_likelihood added here 1" << endl;
                    // cout << "n_pos is " << list_elem_parent->n_pos << endl;
                    // set dna of parent as dna of left child
                    list_elem_parent->dna = dna_l;
                    // set branch length of parent as 0 (default)
                    // we can compute log (P(i,i,l)^{n_i}) ~ n_i log (1 + Q(i,i)l) ~ n_iQ(i,i)l 
                    // where n_i is the number of ref nucs that is i
                    // thus contribution towards log_likelihood is sum_i n_iQ(i,i)(l_c + l_r)
                    this->cum_ref_nuc_counts_map[list_elem_parent->start_pos - 1];
                    for (unsigned char dna = 0; dna < 4; dna ++) {
                        ref_nuc_counts[dna] = this->cum_ref_nuc_counts_map[list_elem_parent->start_pos + list_elem_parent->n_pos -1][dna];
                        ref_nuc_counts[dna] -= this->cum_ref_nuc_counts_map[list_elem_parent->start_pos -1][dna];
                        // cout << " ref_nuc_counts[" << (int) dna << "] is " << (int) ref_nuc_counts[dna] << endl;
                        // cout << " Q(" << (int) dna << "," << (int) dna << ") is " << this->Q_GTR(dna,dna) << endl;
                        // cout << " sum of edge lengths is " << list_elem_left_child->length + list_elem_right_child->length << endl;
                        parent->log_scaling_factor += ref_nuc_counts[dna] * this->Q_GTR(dna,dna) * (list_elem_left_child->length + list_elem_right_child->length);                        
                    }
                    // cout << "log likelihood added is " << parent->log_scaling_factor << endl;
                }
            } else {
                // 6) dna_l and dna_r are (i) ambiguous and identical, (ii) ambiguous and distinct, or (iii) non-ambiguous and distinct
                // set branch length for parent as 0 (default)
                // if (list_elem_parent->start_pos == 6730) {
                // cout << "dna_l is " << (int) dna_l << "\t" << "dna_r is " << (int) dna_r << endl;
                // }
                // compute clv of parent using clv of children
                // cout << "********************************************************************************" << endl;
                // cout << "case 6: dna_l and dna_r are (i) ambiguous and identical, (ii) ambiguous and distinct, or (iii) non-ambiguous and distinct" << endl;
                // cout << "********************************************************************************" << endl;
                if (left_child->clv_for_list_elem.find(list_elem_left_child) != left_child->clv_for_list_elem.end()) {
                    if (dna_l < 6) {
                        cout << "dna_l is " << int(dna_l) << endl;
                        cout << "clv entry for dna_l is " << left_child->clv_for_list_elem[list_elem_left_child][dna_l] << endl;
                    }
                    if (debug) {assert(dna_l > 5);}
                    clv_left_child = left_child->clv_for_list_elem[list_elem_left_child];
                    if (debug) {assert(*max_element(clv_left_child.begin(),clv_left_child.end()) > 0.0);}
                } else {
                    // cout << "I'm here now" << endl;
                    if (debug) {assert(dna_l < 5);}
                    if (dna_l < 4) {
                        clv_left_child = this->Get_clv_for_dna(dna_l);
                    } else {
                        // cout << "Hello! I'm here now" << endl;
                        // dna = 4 (reference)                        
                        // cout << "dna for reference is " << (int) this->reference_sequence[list_elem_parent->start_pos -1] << endl;
                        clv_left_child = this->Get_clv_for_dna(this->reference_sequence[list_elem_parent->start_pos -1]);
                        // cout << "max clv element is " << *max_element(clv_left_child.begin(),clv_left_child.end()) << endl;
                        if (debug) {assert(*max_element(clv_left_child.begin(),clv_left_child.end()) > 0.0);}
                    }
                }

                if (right_child->clv_for_list_elem.find(list_elem_right_child) != right_child->clv_for_list_elem.end()) {
                    // cout << "****\t" << (int) dna_r << endl; 
                    if (debug) {assert(dna_r > 5);}
                    clv_right_child = right_child->clv_for_list_elem[list_elem_right_child];
                    if (debug) {assert(*max_element(clv_right_child.begin(),clv_right_child.end()) > 0.0);}
                } else {
                    if (debug) {assert(dna_r < 5);}
                    if (dna_r < 4) {
                        clv_right_child = this->Get_clv_for_dna(dna_r);
                    } else {                        
                        // cout << "dna for reference is " << (int) this->reference_sequence[list_elem_parent->start_pos -1] << endl;
                        clv_right_child = this->Get_clv_for_dna(this->reference_sequence[list_elem_parent->start_pos -1]);
                        if (debug) {assert(*max_element(clv_right_child.begin(),clv_right_child.end()) > 0.0);}
                    }
                }
                if (debug) {assert(*max_element(clv_left_child.begin(),clv_left_child.end()) > 0.0);}
                if (debug) {assert(*max_element(clv_right_child.begin(),clv_right_child.end()) > 0.0);}
                // assert(list_elem_left_child->length > 0);
                // assert(list_elem_right_child->length > 0);
                clv_parent[0] = 1.0; clv_parent[1] = 1.0; clv_parent[2] = 1.0; clv_parent[3] = 1.0;
                largest_elem = 0;
                for (unsigned char dna_p = 0; dna_p < 4; dna_p ++) {
                    partial_likelihood_l = 0; partial_likelihood_r = 0;
                    for (unsigned char dna_c = 0; dna_c < 4; dna_c ++) {
                        if (dna_p == dna_c) {
                            // cout << (int) dna_p << "\t" << clv_left_child[dna_p] << "\t" << clv_right_child[dna_p] << endl;
                            partial_likelihood_l += (1.0 + (this->Q_GTR(dna_p,dna_p) * list_elem_left_child->length)) * clv_left_child[dna_c];
                            partial_likelihood_r += (1.0 + (this->Q_GTR(dna_p,dna_p) * list_elem_right_child->length)) * clv_right_child[dna_c];
                        } else {
                            partial_likelihood_l += (this->Q_GTR(dna_p,dna_c) * list_elem_left_child->length) * clv_left_child[dna_c];
                            partial_likelihood_r += (this->Q_GTR(dna_p,dna_c) * list_elem_right_child->length) * clv_right_child[dna_c];
                        }                        
                    }
                    clv_parent[dna_p] = partial_likelihood_l * partial_likelihood_r;
                    // cout << "updated clv is " << clv_parent[dna_p] << endl;
                    if (largest_elem < clv_parent[dna_p]) {
                        index_largest_elem = dna_p;
                        largest_elem = clv_parent[dna_p];
                    }
                    // cout << "largest element is " << largest_elem << endl;
                    // cout << "index of largest element is " << (int) index_largest_elem << endl;
                }                
                // rescale clvs
                for (unsigned char dna_p = 0; dna_p < 4; dna_p++) {
                    clv_parent[dna_p] /= largest_elem;                    
                }
                // cout << "log_likelihood added here 2" << endl;                
                parent->log_scaling_factor += log(largest_elem);
                parent->clv_for_list_elem[list_elem_parent] = clv_parent;
                if (debug) {assert(*max_element(clv_parent.begin(),clv_parent.end()) > 0.0);}
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
                        // cout << "**************************" << endl;
                        // cout << "Unlikely case encountered" << endl;
                        // cout << "ref pos is " << list_elem_parent->start_pos << endl;
                        // cout << "n pos is " << list_elem_parent->n_pos << endl;                        
                        // cout << "dna_l is " << (int) dna_l << "\t" << "dna_r is " << (int) dna_r << endl;
                        // for (int dna = 0; dna < 4; dna ++) {
                        //     cout << "clv left child " << dna << " is " << clv_left_child[dna] << "\t";
                        //     cout << "clv right child " << dna << " is " << clv_right_child[dna] << "\t";
                        //     cout << "clv parent " << dna << " is " << clv_parent[dna] << endl;
                        // }
                        // cout << "dna assigned to parent " << parent->name << " is " << (int) index_largest_elem << endl;
                        list_elem_parent->dna = index_largest_elem;
                        parent->clv_for_list_elem.erase(list_elem_parent);
                    }
                }       
            }
        }
    }    
}


void tree::Compute_log_likelihood_score_for_tree_using_genome_list() {
    this->log_likelihood = this->root->log_scaling_factor;
    int n_pos; unsigned char dna_r;
    array <double, 4> clv_root;
    double site_likelihood;
    array <unsigned char, 4> ref_nuc_counts;
    int tot_pos = 0;
    for (genome_list_elem * list_elem_root : this->root->genome_list) {
        n_pos = list_elem_root->n_pos;        
        dna_r = list_elem_root->dna;
        if (dna_r == 4) {
            for (unsigned char dna = 0; dna < 4; dna ++) {
                ref_nuc_counts[dna] = this->cum_ref_nuc_counts_map[list_elem_root->start_pos + list_elem_root->n_pos -1][dna];
                ref_nuc_counts[dna] -= this->cum_ref_nuc_counts_map[list_elem_root->start_pos -1][dna];
                this->log_likelihood += ref_nuc_counts[dna] * log(pi_rho[dna]);
            }
        } else {
            if (list_elem_root->dna > 4){
                assert (this->root->clv_for_list_elem.find(list_elem_root) != this->root->clv_for_list_elem.end());
                clv_root = this->root->clv_for_list_elem[list_elem_root];
            } else {
                clv_root = this->Get_clv_for_dna(list_elem_root->dna);
            }                        
            assert (list_elem_root->n_pos == 1);
            site_likelihood = 0;
            for (unsigned char dna = 0; dna < 4; dna ++) {
                site_likelihood += pi_rho[dna] * clv_root[dna];
            }
            this->log_likelihood += log(site_likelihood);
        }        
    }
}

void tree::Compute_loglikelihood_using_fast_pruning_algorithm() {
    // Apply pruning algorithm to genome lists
    Matrix4f Q = this->Q_GTR;
    Matrix4f Q_scaled_l; Matrix4f Q_scaled_r;
    Matrix4f P_l; Matrix4f P_r;
    float t_l; float t_r;
	float scaling_factor;
    node * left_child; node * right_child;
    unsigned char ch;
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
    bool debug_continue = false;
    bool verbose = false;
    for (node * parent : this->nodes_for_postorder_traversal) {
        if (!parent->leaf) { // iterate over ancestral nodes
            left_child = parent->children[0];
            right_child = parent->children[1];
            ind_left = 0; ind_right = 0;
            num_list_elem_left_child = left_child->genome_list.size();
            num_list_elem_right_child = right_child->genome_list.size();
            // cout << "Updating list elements for " << left_child->name << " and " << right_child->name << endl;
            // Add genome wide log scaling factors of children to parent
            parent->log_scaling_factor = left_child->log_scaling_factor + right_child->log_scaling_factor;
            while (ind_left < num_list_elem_left_child && ind_right < num_list_elem_right_child) {
                // debug_continue = false;
                list_elem_parent = new genome_list_elem;
                parent->genome_list.push_back(list_elem_parent);
                list_elem_left_child = left_child->genome_list[ind_left];
                list_elem_right_child = right_child->genome_list[ind_right];         
                if (verbose) {
                    cout << endl << "Before update" << endl ;
                    cout << "Left child:" ;
                    cout << "\tname:\t" << left_child->name;
                    cout << "\tdna:\t" << (int) list_elem_left_child->dna;
                    cout << "\tstart pos:\t" << (int) list_elem_left_child->start_pos;
                    cout << "\tn pos:\t" << (int) list_elem_left_child->n_pos << endl;
                    
                    cout << "Right child:" ;
                    cout << "\tname:\t" << right_child->name;
                    cout << "\tdna:\t" << (int) list_elem_right_child->dna;
                    cout << "\tstart pos:\t" << (int) list_elem_right_child->start_pos;
                    cout << "\tn pos:\t" << (int) list_elem_right_child->n_pos << endl;
                }
                // cout << "log scaling factor for " << parent->name << "\tbefore update is\t" << parent->log_scaling_factor << endl;    
                this->Update_list_elements(parent, left_child, right_child, list_elem_parent, list_elem_left_child, list_elem_right_child);        
                // cout << "log scaling factor for " << parent->name << "\tafter update is\t" << parent->log_scaling_factor << endl;    
                if (parent->parent != parent) {
                    // add length of branch from grandparent to parent to list_elem_parent->length
                    list_elem_parent->length += this->Get_directed_edge_length(parent->parent, parent);
                }
                if (verbose || list_elem_parent->dna>16) {
                    cout << "After update" << endl ;
                    cout << "left child:" ;
                    cout << "\tname:\t" << left_child->name;
                    cout << "\tdna:\t" << (int) list_elem_left_child->dna;
                    cout << "\tstart pos:\t" << (int) list_elem_left_child->start_pos;
                    cout << "\tn pos:\t" << (int) list_elem_left_child->n_pos << endl;
                    
                    cout << "right child:" ;
                    cout << "\tname:\t" << right_child->name;
                    cout << "\tdna:\t" << (int) list_elem_right_child->dna;
                    cout << "\tstart pos:\t" << (int) list_elem_right_child->start_pos;
                    cout << "\tn pos:\t" << (int) list_elem_right_child->n_pos << endl;

                    cout << "parent:" ;
                    cout << "\tname:\t" << parent->name;
                    cout << "\tdna:\t" << (int) list_elem_parent->dna;
                    cout << "\tstart pos:\t" << (int) list_elem_parent->start_pos;
                    cout << "\tn pos:\t" << (int) list_elem_parent->n_pos << endl;
                }
                if (debug) {assert (list_elem_parent->dna < 17);}

                if (list_elem_left_child->n_pos == 0) {
                    if (verbose) {
                        cout << "========= shift left index =========" << endl;
                    }                    
                    ind_left += 1;
                }
                if (list_elem_right_child->n_pos == 0) {
                    if (verbose) {
                        cout << "========= shift right index ========" << endl;
                    }                    
                    ind_right += 1;
                }
            }
            // cout << "ind_left is\t" << ind_left << "\tnum list elem in left child is\t" << num_list_elem_left_child << endl;
            // cout << "ind_right is\t" << ind_right << "\tnum list elem in right child is\t" << num_list_elem_right_child << endl;
            if (debug) {assert (ind_left == num_list_elem_left_child && ind_right == num_list_elem_right_child);}  
            // cout << "total log likelihood added to " << parent->name << "\tafter update is\t" << parent->log_scaling_factor << endl;              
            // break;
            // collapse contiguous reference type list elements into a single list element
            this->Combine_ref_type_list_elements(parent);            
        }
        // compute total log likelihood score
        this->Compute_log_likelihood_score_for_tree_using_genome_list();        
    }        
}

void tree::Set_model_parameters() {
    bool iqtree_params = true;
    if (iqtree_params) {
        for (int i = 0; i < 4; i ++) {
            this->pi_rho[i] = 0.25;
        }
    }
    vector <float> exchangeability_params;
    exchangeability_params.push_back(0.04); //A<->C beta
    exchangeability_params.push_back(0.3); //A<->G alpha
    exchangeability_params.push_back(0.1); //A<->T gamma
    exchangeability_params.push_back(0.02); //C<->G delta
    exchangeability_params.push_back(1.0); //C<->T eta
    exchangeability_params.push_back(1.0); //G<->T epsilon

    this->Q_GTR(0,1) = this->pi_rho[1] * exchangeability_params[0]; // pi[C] * R[A<->C]
    this->Q_GTR(0,2) = this->pi_rho[2] * exchangeability_params[1]; // pi[G] * R[A<->G]
    this->Q_GTR(0,3) = this->pi_rho[3] * exchangeability_params[2]; // pi[T] * R[A<->T]
    this->Q_GTR(1,0) = this->pi_rho[0] * exchangeability_params[0]; // pi[A] * R[A<->C]
    this->Q_GTR(1,2) = this->pi_rho[2] * exchangeability_params[3]; // pi[G] * R[C<->G]
    this->Q_GTR(1,3) = this->pi_rho[3] * exchangeability_params[4]; // pi[T] * R[C<->T]
    this->Q_GTR(2,0) = this->pi_rho[0] * exchangeability_params[1]; // pi[G] * R[A<->G]
    this->Q_GTR(2,1) = this->pi_rho[1] * exchangeability_params[3]; // pi[C] * R[C<->G]
    this->Q_GTR(2,3) = this->pi_rho[3] * exchangeability_params[5]; // pi[T] * R[G<->T]
    this->Q_GTR(3,0) = this->pi_rho[0] * exchangeability_params[2]; // pi[A] * R[A<->T]
    this->Q_GTR(3,1) = this->pi_rho[1] * exchangeability_params[4]; // pi[C] * R[C<->T]
    this->Q_GTR(3,2) = this->pi_rho[2] * exchangeability_params[5]; // pi[G] * R[G<->T]

    this->Q_GTR(0,0) = - (this->Q_GTR(0,1) + this->Q_GTR(0,2) + this->Q_GTR(0,3));
    this->Q_GTR(1,1) = - (this->Q_GTR(1,0) + this->Q_GTR(1,2) + this->Q_GTR(1,3));
    this->Q_GTR(2,2) = - (this->Q_GTR(2,0) + this->Q_GTR(2,1) + this->Q_GTR(2,3));
    this->Q_GTR(3,3) = - (this->Q_GTR(3,0) + this->Q_GTR(3,1) + this->Q_GTR(3,2));
    
    Q_GTR /= this->Compute_scaling_factor(Q_GTR);

    // cout << "GTR rate matrix is " << endl << this->Q_GTR << endl;
    // Set root and compute vertex list for post order traversal
    node * l = this->leaves[0];
    this->Root_tree_along_edge(l, l->neighbors[0], 0.0);
    this->Set_nodes_for_postorder_traversal();
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
    // cout << "populating dna to char map" << endl;
    // cout << this->DNA_to_char.size() << endl;
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

void tree::Add_node(string u_name) {
    node * u = new node(u_name);
    this->node_list.insert(pair<string,node *>(u_name,u));
}

node * tree::Get_node(string u_name) {
    bool node_present = this->Contains_node(u_name);
    if (!node_present) {
        cout << "node " << u_name << " is missing " << endl;
    }    
    if (debug) {assert(node_present);}
    return(this->node_list[u_name]);
}

bool tree::Contains_node(string u_name) {
    bool node_present;
    map <string, node *>::iterator it = this->node_list.find(u_name);
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
        this->undirected_edge_length_map.insert(make_pair(make_pair(u,v),length));
    } else {
        this->undirected_edge_length_map.insert(make_pair(make_pair(v,u),length));
    }
}

float tree::Get_undirected_edge_length(node * u, node * v) {
    float length;
    if (u < v) {
        length = this->undirected_edge_length_map[make_pair(u,v)];
    } else {
        length = this->undirected_edge_length_map[make_pair(v,u)];
    }
    return (length);
}

float tree::Get_directed_edge_length(node * p, node * c) {
    return (this->directed_edge_length_map[make_pair(p,c)]);
}

void tree::Read_newick_file(bool rooted) {
    // Hardcoded rooted to false
    cout << "Tree file name is " << tree_file_name << endl;
    rooted = false;
    string node_name; string h_name;
    node * h; node * n;
    float length; float v_length;
    string newick_string;
    string sibling_string;
    string sibling_string_without_parenthesis;
    string newick_string_split;
    vector <string> split_sibling_string;
    vector <string> node_name_and_length;
    ifstream inputFile(this->tree_file_name.c_str());
    getline(inputFile, newick_string);
    regex sibling_pattern ("\\([^\\(\\)]+\\)");
    smatch cherry_matches;    
    bool continue_search = true;
    while (continue_search) {
        if (regex_search(newick_string, cherry_matches, sibling_pattern)) {
            h_name = "h_" + to_string(this->h_ind);
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
                for (string sibling_string : split_sibling_string) {
                    boost::split(node_name_and_length, sibling_string, [](char c){return c == ':';});
                    node_name = node_name_and_length[0];
                    length = stof(node_name_and_length[1]);                    
                    if (!this->Contains_node(node_name)) {
                        this->Add_node(node_name);
                    }
                    n = this->Get_node(node_name);
                    this->Add_undirected_edge(n, h, length);
                }
            } else {
                cerr << "sibling string has not been properly parsed" << endl;
            }            
            newick_string = regex_replace(newick_string, sibling_pattern, h_name, regex_constants::format_first_only);                
        } else {
            continue_search = false;
        }
        // cout << "Number of edges is " << this->undirected_edge_length_map.size() << endl;
        // cout << "Number of vertices is " << this->node_list.size() << endl;
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
	ifstream inputFile(this->reference_file_name.c_str());
	string seq_name;
	string seq = "";    
	for (string line; getline(inputFile, line);) {
		if (line[0]=='>') {
			seq_name = line.substr(1,line.length());
		}
		else {
		    seq += line;
		}
	}
    inputFile.close();
    // set pi_rho    
    array <int, 4> nuc_counts;
    for (int i = 0; i < 4; i++) {
        nuc_counts[i] = 0;
        this->pi_rho[i] = 0.0;
    }

	for (char const dna: seq) {
        c = this->DNA_to_char[string(1,tolower(dna))];
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
        // cout << "pi[" << i << "] is " << this->pi_rho[i] << "\t";
    }
    // cout << endl;
}

void tree::Add_fasta_sequences() {
    vector <unsigned char> recoded_sequence;
    node * v;
	recoded_sequence.clear();
	unsigned int site = 0;
	ifstream inputFile(this->alignment_file_name.c_str());
	string seq_name;
	string seq = "";	
	for(string line; getline(inputFile, line );) {
		if (line[0]=='>') {
			if (seq != "") {
				for (char const dna: seq) {                    
					recoded_sequence.push_back(this->DNA_to_char[string(1,tolower(dna))]);
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
		recoded_sequence.push_back(this->DNA_to_char[string(1,tolower(dna))]);
		site += 1;
	}
	v = this->Get_node(seq_name);
    v->complete_sequence = recoded_sequence;
	recoded_sequence.clear();
	inputFile.close();
}

void tree::Add_ref_nuc_counts_based_on_genome_coordinates() {
    // map <unsigned short, array <unsigned char, 4>> cum_ref_nuc_counts_map;
    unsigned short end_pos;
    this->cum_ref_nuc_counts_map.clear();
    array <unsigned char, 4> cum_nuc_counts;
    cum_nuc_counts[0] = 0; cum_nuc_counts[1] = 0; cum_nuc_counts[2] = 0; cum_nuc_counts[3] = 0;
    for (node * leaf : this->leaves) {
        for (genome_list_elem * list_elem : leaf->genome_list) {            
            end_pos = list_elem->start_pos + list_elem->n_pos - 1;
            if (this->cum_ref_nuc_counts_map.find(end_pos) == this->cum_ref_nuc_counts_map.end()) {
                this->cum_ref_nuc_counts_map[end_pos] = cum_nuc_counts;
            }
        }
    }
    unsigned char c;
    for (unsigned short pos = 1; pos < GENOME_LENGTH; pos++) {
        c = this->reference_sequence[pos];        
        // assert (c < 4);
        cum_nuc_counts[c] += 1;
        if (this->cum_ref_nuc_counts_map.find(pos) != this->cum_ref_nuc_counts_map.end()) {
            this->cum_ref_nuc_counts_map[pos] = cum_nuc_counts;
        }
    }
    cum_nuc_counts[0] = 0; cum_nuc_counts[1] = 0; cum_nuc_counts[2] = 0; cum_nuc_counts[3] = 0;
    // for the case where a stretch of ref nucs starts at 1
    this->cum_ref_nuc_counts_map[0] = cum_nuc_counts;
}

void tree::Add_mut_diff_sequences() {
    // vector <unsigned char> recoded_sequence;
    node * n;
	// recoded_sequence.clear();
	// unsigned int site = 0;
	ifstream inputFile(this->mut_diff_file_name.c_str());
	string seq_name = "";
	// string seq = "";
    vector <string> split_line;
    split_line.reserve(3);
    string dna;
    array <double, 4> clv;
    genome_list_elem * list_elem_diff; // Read from mut diff file
    genome_list_elem * list_elem_ref; // List element for positions identical to reference sequence
    genome_list_elem * list_elem_prev_diff;
    unsigned int genome_list_size;
    bool empty_genome_list = true;
	for (string line; getline(inputFile, line );) {
		if (line[0]=='>') {
            if (seq_name != "") {                              
                n->Set_pos_for_reference_characters();
            }
            seq_name = line.substr(1,line.length());
            n = this->Get_node(seq_name);            
            empty_genome_list = true;
		}
		else {
            boost::split(split_line, line, [](char c){return c == '\t';});
            list_elem_diff = new genome_list_elem();
            dna = split_line[0];
            list_elem_diff->dna = this->DNA_to_char[dna];
            list_elem_diff->length = this->Get_directed_edge_length(n->parent, n);
            // assert(list_elem_diff->length > 0);
		    if (split_line.size() == 3) {
                if (debug) {assert (split_line[0] == "n" || split_line[0] == "-" );}
                list_elem_diff->start_pos = stoi(split_line[1]);
                list_elem_diff->n_pos = stoi(split_line[2]);
            } else {
                if (debug) {assert (split_line.size() == 2);}
                list_elem_diff->start_pos = stoi(split_line[1]);
                list_elem_diff->n_pos = 1;
            }
            split_line.clear();
            if (list_elem_diff->dna > 5) { // ambiguous IUPAC code and not N or -
                clv = Get_clv_for_dna(list_elem_diff->dna);
                n->clv_for_list_elem[list_elem_diff] = clv;
            }
            // if (n->name == "EPI_ISL_1589845") {
            //     cout << "dna is " << (int) list_elem_diff->dna << endl;
            //     cout << "start pos is " << list_elem_diff->start_pos << endl;
            //     cout << "n pos is " << list_elem_diff->n_pos << endl;
            // }
            if (empty_genome_list && list_elem_diff->start_pos > 1) {
                list_elem_ref = new genome_list_elem();
                list_elem_ref->dna = 4; // 4 is the code for reference sequence
                list_elem_ref->length = this->Get_directed_edge_length(n->parent, n);
                n->genome_list.push_back(list_elem_ref);
                empty_genome_list = false;
            }
            // erase ref element that is sandwiched between contiguous diff elements            
            genome_list_size = n->genome_list.size();
            if (genome_list_size > 1) {
                list_elem_prev_diff =  n->genome_list[genome_list_size - 2];
                // cout << (int) list_elem_prev_diff->dna << "\t" << list_elem_prev_diff->n_pos << endl;
                if (debug) {assert (list_elem_prev_diff->dna != 4);}
                if (list_elem_prev_diff->start_pos + list_elem_prev_diff->n_pos == list_elem_diff->start_pos) {
                    list_elem_ref = n->genome_list[genome_list_size-1];
                    if (debug) {assert(list_elem_ref->dna == 4);}
                    n->genome_list.pop_back();
                    delete list_elem_ref;                    
                }
            }

            n->genome_list.push_back(list_elem_diff);
            empty_genome_list = false;
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
}
void tree::Compress_sequences() {
    
    this->character_pattern_weights.clear();
    map <vector<unsigned char>,int> patterns_to_weights;
	int number_of_sites = this->leaves[0]->complete_sequence.size();
    int number_of_leaves = this->leaves.size();
    cout << "number of sites is " << number_of_sites << endl;
	vector <unsigned char> character_pattern;
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
    cout << "Number of character patterns is " << this->num_char_patterns << endl;
    for (int site = 0; site < this->num_char_patterns; site++) {
		character_pattern.clear();
		for (int leaf_ind = 0; leaf_ind < number_of_leaves; leaf_ind ++) {
            // cout << (int)this->leaves[leaf_ind]->compressed_sequence[site];
			character_pattern.push_back(this->leaves[leaf_ind]->compressed_sequence[site]);
            }
        // cout << endl;
        this->character_pattern_weights.push_back(patterns_to_weights[character_pattern]);
    }    
}

class fastLK_overview {
public:    
    tree * T;
    void Run_workflow(string workflow_type);
    fastLK_overview(string path_to_reference_file, string path_to_mut_diff_file, string path_to_sequence_alignment_file, string path_to_tree_file, float clv_threshold){        
        this->T = new tree;
        this->T->tree_file_name = path_to_tree_file;
        this->T->alignment_file_name = path_to_sequence_alignment_file;
        this->T->reference_file_name = path_to_reference_file;
        this->T->mut_diff_file_name = path_to_mut_diff_file;
        this->T->clv_threshold = clv_threshold;
        string workflow_type = "mut_diff";
        // string workflow_type = "standard"; 
        this->Run_workflow(workflow_type);
    }

    ~fastLK_overview(){
        delete this->T;
    }
};

void fastLK_overview::Run_workflow(string workflow_type){

    // Read tree file
    // Add sequences (perform global site pattern compression)    
    this->T->Read_newick_file(false);    
    this->T->Set_leaves();
    this->T->Read_reference_sequence();
    this->T->Set_model_parameters();
    if (workflow_type == "standard") {
        this->T->Add_fasta_sequences();
        this->T->Compress_sequences();
        this->T->Compute_loglikelihood_using_standard_pruning_algorithm();        
        cout << "log likelihood score computed using standard pruning algorithm is " << setprecision(8) << this->T->log_likelihood << endl;
    } else if (workflow_type == "mut_diff") {
        this->T->Add_mut_diff_sequences();
        for (node * n : this->T->leaves) {
            if (n->clv_for_list_elem.size() > 0) {
                // cout << n->name << endl;
            }
        }        
        this->T->Add_ref_nuc_counts_based_on_genome_coordinates();        
        this->T->Compute_loglikelihood_using_fast_pruning_algorithm();
        cout << "log likelihood score computed using fast pruning algorithm is " << setprecision(8) << this->T->log_likelihood << endl;
    }    

    // ************************************************ //
    // Compute log likelihood for GTR rate matrix using //
    // standard Felsenstein's pruning algorithm         //
    // ************************************************ //

    // cout << "Log likelihood score computing using standard pruning algorithm is " << setprecision(8) << this->T->log_likelihood << endl;


    // ************************************************ //
    // Compute log likelihood for GTR rate matrix using //
    // approximate Felsenstein's pruning algorithm      //
    // ************************************************ //
}

#endif