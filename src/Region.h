#ifndef __REGION_H__
#define __REGION_H__

#include "Tuple.h"

#include "libqhullcpp/Qhull.h"
#include "libqhullcpp/QhullError.h"
#include "libqhullcpp/RboxPoints.h"
#include "libqhullcpp/QhullFacetList.h"
#include "libqhullcpp/QhullFacetSet.h"
#include "libqhullcpp/QhullVertexSet.h"
#include "libqhullcpp/QhullUser.h"
#include "libqhullcpp/QhullVertex.h"

using orgQhull::Coordinates;
using orgQhull::Qhull;
using orgQhull::QhullFacetListIterator;
using orgQhull::QhullHyperplane;

class HalfSpace {
public:
    int dim;
    vector<double> coef;
    bool side;

    HalfSpace() {}
    virtual ~HalfSpace() {}
};

class Region {
public:
    int dim;
    int m;
    vector<HalfSpace> spaces;
    vector<vector<double> > vertices;
    vector<double> inner;

    Region();
    Region(const int& dim, const int& m);
    Region(const Region& other);
    virtual ~Region();
    Region& operator= (const Region& other);

    void gen_weak_query();
    void gen_inter_query();
    void load_query(const char* query_path);
    void write_query(const char* query_path);
    void print_query();
    void compute_vertex();

};

#endif