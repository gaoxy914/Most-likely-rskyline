#ifndef __TUPLE_H__
#define __TUPLE_H__

#include "BigFloat.h"

class Tuple {
public:
    int id;
    int dim;
    double *coord;
    double prob;

    Tuple();
    Tuple(const int& id, const int& dim, const double* coord, const double& prob);
    Tuple(const Tuple& other);
    virtual ~Tuple();
    Tuple& operator= (const Tuple& other);
    double operator[] (const int& index);
    double score(const vector<double>& weight) const;
    double min_score(const vector<vector<double> > weights) const;
    bool fdominate(const vector<vector<double> >& weights, const Tuple& other);
    friend ostream& operator<< (ostream& out, const Tuple& t);

    struct Comparator {
        vector<double> weight;
        Comparator(const vector<double>& weight) : weight(weight) {}
        
        bool operator ()(const Tuple& t, const Tuple& s) const {
            return t.score(weight) < s.score(weight);
        }
    };
};

#endif