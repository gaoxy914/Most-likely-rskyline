#ifndef __MRSKY_H__
#define __MRSKY_H__

#include "Tuple.h"
#include "Region.h"
#include "Heap.h"

class MRSKY {
private:
    vector<int> layer0;
    vector<vector<int> > groups_at_layer0;

    unordered_map<int, unordered_set<int> > graph;
    vector<int> in_degree;
    unordered_map<int, unordered_set<int> > reverse_graph;

    BigFloat prob_res_at_layer0;
    vector<int> res_at_layer0;

    long int num_candidates;

    int dim;
    int n;
    double center;
    string data_path;
    vector<Tuple> tuples;

    bool fdominates(const vector<int>& rskyline, const int& t, const Region& F);
    bool fdominates(const unordered_set<int>& rskyline, const int& t, const Region& F);
    bool fdominates(const vector<vector<double> >& rskyline, const int& t, const Region& F);
    bool fdominates(const vector<double>& s, const int& t, const Region& F);
    bool dominates(const vector<double>& s, const int& t, const Region& F);
    void print_graph();
    pair<int, double> compute_ratio(const unordered_set<int>& rskyline, const unordered_map<int, int>& count, const Region& F);
    void check_count(const unordered_set<int>& partial_result, const unordered_map<int, int>& count, const Region& F);
    bool valid_swap(const int& t, const int& s, const unordered_set<int>& rskyline, const unordered_map<int, int>& count, const Region& F);

public:

    MRSKY(const int& dim, const int& n, const double& center);
    virtual ~MRSKY();

    void gen_ind_data();
    void gen_anti_data();
    void load_data(const char* data_path);
    void load_car_data();
    void load_nba_data();
    void write_data(const char* data_path);
    void print_data();

    void construct_layers(const Region& F);

    int point_reduction(const Region& F);
    int path_reduction(const Region& F);
    int region_reduction(const Region& F);
    int region_reduction_2d(const Region& F);
    int virtual_region_reduction(const Region& F);
    int data_reduction_USKY(const Region& F);
    
    void group_layer0();

    BigFloat DSA(vector<int>& result, const Region& F);
    BigFloat IBBA(vector<int>& result, const Region& F);
    void IBBA_rec(const vector<int>& rskyline, const BigFloat& prob, vector<pair<int, double> > candidates, vector<int>& opt_result, BigFloat& opt_prob, const int& depth);
    BigFloat EBBA(vector<int>& result, const Region& F);
    void EBBA_rec(const vector<int>& rskyline, const BigFloat& prob, vector<pair<int, double> > candidates, vector<int>& opt_result, BigFloat& opt_prob, const int& depth);
    BigFloat greedy(vector<int>& result, const Region& F);
    BigFloat half_rounding(vector<int>& result, const Region& F);
    int init_S(const Region& F, unordered_set<int>& S);
    BigFloat LSA(vector<int>& result, const Region& F, const double& tau);
    BigFloat LSA_PLUS(vector<int>& result, const Region& F, const double& tau);
    
    void ARSP(vector<pair<int, double> >& result, const Region& F);

    void check_prob(const vector<int>& rskyline, const BigFloat& prob, const Region& F);
};

#endif