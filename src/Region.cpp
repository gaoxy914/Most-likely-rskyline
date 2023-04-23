#include "Region.h"

Region::Region() {
    dim = 0;
    m = 0;
}

Region::Region(const int& dim, const int& m) {
    this->dim = dim;
    this->m = m;
    inner.resize(dim);
    spaces.resize(m);
    for (int i = 0; i < m; ++ i) {
        spaces[i].dim = dim - 1;
        spaces[i].coef.resize(dim);
    }
}

Region::~Region() {}

Region::Region(const Region& other) {
    dim = other.dim;
    m = other.m;
    inner = other.inner;
    spaces = other.spaces;
    vertices = other.vertices;
}

Region& Region::operator= (const Region& other) {
    if (&other != this) {
        dim = other.dim;
        m = other.m;
        inner = other.inner;
        spaces = other.spaces;
        vertices = other.vertices;       
    }
    return *this;
}

void Region::gen_weak_query() {
    assert(m = dim - 1);
    double sum = (1.0 + dim)*dim/2.0;
    for (int i = 0; i < dim; ++ i) {
        inner[i] = (dim - i)/sum;
    }
    for (int i = 0; i < m - 1; ++ i) {
        for (int j = 0; j < dim; ++ j) {
            if (j == i) spaces[i].coef[j] = 1;
            if (j == i + 1) spaces[i].coef[j] = -1;
        }
        spaces[i].side = true;
    }
    for (int i = 0; i < dim; ++ i) {
        if (i == dim - 2) spaces[m - 1].coef[i] = 2;
        else spaces[m - 1].coef[i] = 1;
    }
    spaces[m - 1].side = true;
}

void Region::gen_inter_query() {
    double sum = 0;
    for (int i = 0; i < dim; ++ i) {
        inner[i] = rand_uniform(0, 1);
        sum += inner[i];
    }
    for (int i = 0; i < dim; ++ i) {
        inner[i] /= sum;
    }
    double t[dim], s[dim], side;
    for (int i = 0; i < m; ++ i) {
        side = 0;
        for (int j = 0; j < dim; ++ j) {
            t[j] = drand48();
            s[j] = drand48();
            spaces[i].coef[j] = t[j] - s[j];
            side += spaces[i].coef[j]*inner[j];
        }
        for (int j = 0; j < dim - 1; ++ j) {
            spaces[i].coef[j] -= spaces[i].coef[dim - 1];
        }
        spaces[i].coef[dim - 1] *= -1;
        spaces[i].side = (side >= 0);
    }
}
    

void Region::load_query(const char* query_path) {
    ifstream file((string(query_path) + to_string(dim) + string(".qry")).c_str(), ios::in);
    int mMax = 0;
    file >> mMax;
    if (m > mMax) {
        printf("Number of constraints exceeds max.\n");
        exit(1);
    }
    for (int i = 0; i < dim; ++ i) {
        file >> inner[i];
    }
    for (int i = 0; i < m; ++ i) {
        for (int j = 0; j < dim; ++ j) {
            file >> spaces[i].coef[j];
        }
        file >> spaces[i].side;
    }
    file.close();
}

void Region::write_query(const char* query_path) {
    ofstream file ((string(query_path) + to_string(dim) + string(".qry")).c_str(), ios::out);
    file << m << " ";
    for (int i = 0; i < dim; ++ i) {
        file << inner[i] << " ";
    }
    for (int i = 0; i < m; ++ i) {
        for (int j = 0; j < dim; ++ j) {
            file << spaces[i].coef[j] << " ";
        }
        file << spaces[i].side << " ";
    }
    file.close();
}

void Region::print_query() {
    /* printf("Query Parameters: dim = %d, m = %d.\n", dim, m);
    cout << "Inner point: (";
    for (int i = 0; i < dim - 1; ++ i) {
        cout << inner[i] << ",";
    }
    cout << inner[dim - 1] << ")" << endl;
    for (int i = 0; i < m; ++ i) {
        for (int j = 0; j < dim - 2; ++ j) {
            cout << spaces[i].coef[j] << " x w[" << j << "] + ";
        }
        cout << spaces[i].coef[dim - 2] << " x w[" << dim - 2 << "]";
        cout << (spaces[i].side ? " >= " : " <= ") << spaces[i].coef[dim - 1] << endl;
    } */
    cout << "Vertices : " << vertices.size() << endl;
    /* for (int i = 0; i < vertices.size(); ++ i) {
        cout << "(";
        for (int j = 0; j < dim - 1; ++ j) {
            cout << vertices[i][j] << ",";
        }
        cout << vertices[i][dim - 1] << ")" << endl;
    } */
}

void Region::compute_vertex() {
    if (dim > 2) {
        int n = m + 2*(dim - 1) + 1;
        int size = dim*n;
        double normals[size];
        int index = 0;
        for (int i = 0; i < dim; ++ i) {
            if (i == dim - 1) normals[index ++] = -1;
            else if (i < dim - 1) normals[index ++] = 1;
        }
        for (int i = 0; i < dim - 1; ++ i) {
            for (int j = 0; j < dim; ++ j) {
                if (i == j) normals[index ++] = 1;
                else normals[index ++] = 0;
            }
            normals[index - 1] = -1;
            for (int j = 0; j < dim; ++ j) {
                if (i == j) normals[index ++] = -1;
                else normals[index ++] = 0;
            }
        }
        for (int i = 0; i < m; ++ i) {
            if (spaces[i].side) {
                for (int j = 0; j < dim; ++ j) {
                    if (j == dim - 1) normals[index ++] = spaces[i].coef[j];
                    else if (j < dim - 1) normals[index ++] = -spaces[i].coef[j];
                }
            } else {
                for (int j = 0; j < dim; ++ j) {
                    if (j == dim - 1) normals[index ++] = -spaces[i].coef[j];
                    else if (j < dim - 1) normals[index ++] = spaces[i].coef[j];
                }
            }
        }
        Coordinates feasible;
        for (int i = 0; i < dim - 1; ++ i) feasible << inner[i];
        Qhull q;
        q.setFeasiblePoint(feasible);
        q.runQhull("", dim, n, normals, "H");
        vertices.resize(q.facetCount(), vector<double>(dim, 0));
        QhullFacetListIterator it(q.facetList());
        index = 0;
        while (it.hasNext()) {
            QhullHyperplane plane = it.next().hyperplane();
            vertices[index][dim - 1] = 1;
            for (int i = 0; i < dim - 1; ++ i) {
                vertices[index][i] = -plane[i]/plane.offset() + inner[i];
                if (vertices[index][i] < 1e-5) vertices[index][i] = 0;
                vertices[index][dim - 1] -= vertices[index][i];
            }
            if (vertices[index][dim - 1] < 1e-5) vertices[index][dim - 1] = 0;
            ++ index;
        }
    } else if (dim == 2) {
        double left = 0, right = 1;
        for (int i = 0; i < m; ++ i) {
            if (spaces[i].coef[0] == 0) continue;
            double beta = spaces[i].coef[1]/spaces[i].coef[0];
            if (spaces[i].side) {
                if (spaces[i].coef[0] > 0) {
                    left = max(left, beta);
                } else if (spaces[i].coef[0] < 0) {
                    right = min(right, beta);
                }
            } else {
                if (spaces[i].coef[0] > 0) {
                    right = min(right, beta);
                } else if (spaces[i].coef[0] < 0) {
                    left = max(left, beta);
                }
            }
        }
        vertices.resize(2, vector<double>(2, 0));
        vertices[0][0] = left;
        vertices[0][1] = 1 - left;
        vertices[1][0] = right;
        vertices[1][1] = 1 - right;
    }
}