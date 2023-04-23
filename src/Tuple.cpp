#include "Tuple.h"

Tuple::Tuple() {
    id = -1;
    dim = 0;
    coord = nullptr;
    prob = 0;
}

Tuple::Tuple(const int& id, const int& dim, const double* coord, const double& prob) {
    this->id = id;
    this->dim = dim;
    this->prob = prob;
    this->coord = new double[this->dim];
    if (coord != nullptr) {
        memcpy(this->coord, coord, this->dim*sizeof(double));
    } else {
        memset(this->coord, 0, this->dim*sizeof(double));
    }
}

Tuple::Tuple(const Tuple& other) {
    id = other.id;
    dim = other.dim;
    prob = other.prob;
    coord = new double[dim];
    memcpy(coord, other.coord, dim*sizeof(double));
}

Tuple::~Tuple() {
    if (coord != nullptr) {
        delete[] coord;
        coord = nullptr;
    }
    dim = 0;
    prob = 0;
}

Tuple& Tuple::operator= (const Tuple& other) {
    if (&other != this) {
        id = other.id;
        prob = other.prob;
        if (dim != other.dim) {
            dim = other.dim;
            if (coord != nullptr) {
                delete[] coord;
            }
            coord = new double[dim];
        }
        memcpy(coord, other.coord, dim*sizeof(double));
    }
    return *this;
}

double Tuple::operator[] (const int& index) {
    return coord[index];
}

double Tuple::score(const vector<double>& weight) const {
    double s = 0;
    for (int i = 0; i < dim; ++ i) {
        s += weight[i]*coord[i];
    }
    return s;
}

double Tuple::min_score(const vector<vector<double> > weights) const {
    double min_score = dim;
    for (auto weight : weights) {
        min_score = min(min_score, this->score(weight));
    }
    return min_score;
}

bool Tuple::fdominate(const vector<vector<double> >& weights, const Tuple& other) {
    double score1 = 0, score2 = 0;
    bool unique = false;
    for (int i = 0; i < weights.size(); ++ i) {
        score1 = this->score(weights[i]);
        score2 = other.score(weights[i]);
        if (score1 > score2) {
            return false;
        } else if (score1 < score2) {
            unique = true;
        }
    }
    return unique;
}

ostream& operator<< (ostream& out, const Tuple& t) {
    out << "[" << t.id;
    out << " (";
    for (int i = 0; i < t.dim - 1; ++ i) {
        out << t.coord[i] << ", ";
    }
    out << t.coord[t.dim - 1] << "), " << t.prob << "]";
    return out;
}