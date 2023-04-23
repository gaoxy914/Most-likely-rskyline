#include "MRSKY.h"

MRSKY::MRSKY(const int& dim, const int& n, const double& center) {
    srand((unsigned)time(nullptr));
    srand48((unsigned)time(nullptr));
    this->dim = dim;
    this->n = n;
    this->center = center;
    this->tuples.resize(n, Tuple(-1, dim, nullptr, 0));
    stringstream stream;
    stream << fixed << setprecision(2) << center;
    data_path = to_string(dim) + "_" + stream.str();
    this->num_candidates = 0;
    this->prob_res_at_layer0 = BigFloat(1);
    this->in_degree.resize(n, 0);
}

MRSKY::~MRSKY() {

}

void MRSKY::gen_ind_data() {
    double delta = min(center, 1 - center);
    double lower = center - delta, upper = center + delta;
    for (int i = 0; i < n; ++ i) {
        tuples[i].id = i;
        for (int j = 0; j < dim; ++ j) {
            tuples[i].coord[j] = rand_uniform(0, 1);
        }
        tuples[i].prob = rand_uniform(lower, upper);
    }
}

void MRSKY::gen_anti_data() {
    double delta = min(center, 1 - center);
    double lower = center - delta, upper = center + delta;
    double x[dim];
    for (int i = 0; i < n; ++ i) {
        tuples[i].id = i;
        double range = 0.5*dim + rand_normal(0, 0.05);
        tuples[i].coord[0] = rand_uniform(0, 1)*min(1.0, range);
        for (int j = 1; j < dim; ++ j) {
            range -= tuples[i].coord[j - 1];
            tuples[i].coord[j] = rand_uniform(0, 1)*min(1.0, range);
            if (j == dim - 1) tuples[i].coord[j] = min(1.0, range);
        }
        tuples[i].prob = rand_uniform(lower, upper);
    }
}

void MRSKY::load_data(const char* data_path) {
    ifstream file((string(data_path) + this->data_path + ".dat").c_str(), ios::in);
    if (!file.is_open()) {
        printf("Fail in opening data file.\n");
        exit(1);
    }
    int nMax = 0;
    file >> nMax;
    if (n > nMax) {
        printf("n exceeds nMax.\n");
        exit(1);
    }
    for (int i = 0; i < n; ++ i) {
        tuples[i].id = i;
        for (int j = 0; j < dim; ++ j) {
            file >> tuples[i].coord[j];
        }
        file >> tuples[i].prob;
    }
    file.close();
}

void MRSKY::load_car_data() {
    ifstream file("data/car.dat", ios::in);
    if (!file.is_open()) {
        printf("Fail in opening data file.\n");
        exit(1);
    }
    int nMax = 0;
    file >> nMax;
    if (n > nMax) {
        printf("n exceeds nMax.\n");
        exit(1);
    }
    for (int i = 0; i < n; ++ i) {
        tuples[i].id = i;
        for (int j = 0; j < dim; ++ j) {
            file >> tuples[i].coord[j];
        }
        file >> tuples[i].prob;
    }
    file.close();
}

void MRSKY::load_nba_data() {
    ifstream file("data/nba.dat", ios::in);
    if (!file.is_open()) {
        printf("Fail in opening data file.\n");
        exit(1);
    }
    int nMax = 0;
    file >> nMax;
    if (n > nMax) {
        printf("n exceeds nMax.\n");
        exit(1);
    }
    for (int i = 0; i < n; ++ i) {
        tuples[i].id = i;
        for (int j = 0; j < dim; ++ j) {
            file >> tuples[i].coord[j];
        }
        file >> tuples[i].prob;
    }
    file.close();
}

void MRSKY::write_data(const char* data_path) {
    ofstream file((string(data_path) + this->data_path + ".dat").c_str(), ios::out);
    file << n << " ";
    for (int i = 0; i < n; ++ i) {
        for (int j = 0; j < dim; ++ j) {
            file << tuples[i].coord[j] << " ";
        }
        file << tuples[i].prob << " ";
    }
    file.close();
}

void MRSKY::print_data() {
    for (int i = 0; i < tuples.size(); ++ i) {
        cout << i << "," << tuples[i] << endl;
    }
}

void MRSKY::construct_layers(const Region& F) {
    sort(tuples.begin(), tuples.end(), Tuple::Comparator(F.inner));
    vector<vector<int> > layers;
    vector<int> in_degree(tuples.size(), 0);
    unordered_map<int, vector<int> > edges;
    layers.push_back(vector<int>{0});
    for (int i = 1; i < tuples.size(); ++ i) {
        int l = 0, r = layers.size() - 1, j = 0;
        while (l < r) {
            j = (l + r)/2;
            if (fdominates(layers[j], i, F)) {
                l = j + 1;
            } else {
                r = j;
            }
        }
        if (fdominates(layers[l], i, F)) {
            ++ l;
            if (l == layers.size()) {
                layers.push_back(vector<int>{});
            }
        }
        layers[l].push_back(i);
        if (l > 0) {
            for (auto s : layers[l - 1]) {
                if (tuples[s].fdominate(F.vertices, tuples[i])) {
                    edges[s].push_back(i);
                    ++ in_degree[i];
                }
            }
            // check layer
            if (in_degree[i] == 0) {
                cout << "wrong layer " << tuples[i].id << endl;
                for (int i = 0; i < layers.size(); ++ i) {
                    cout << "L" << i << ":\n";
                    for (auto t : layers[i]) {
                        cout << tuples[t].id << ",";
                    }
                    cout << endl;
                }
                exit(1);
            }
        }
    }
    // output G
    for (int i = 0; i < layers.size(); ++ i) {
        cout << "L" << i << ":\n";
        for (auto t : layers[i]) {
            cout <<"[";
            cout << "(" << tuples[t].id << "," << tuples[t].prob << "):";
            for (auto s : edges[t]) {
                cout << tuples[s].id << ", ";
            }
            cout <<"]";
        }
        cout << endl;
    }
}

int MRSKY::point_reduction(const Region& F) {
    sort(tuples.begin(), tuples.end(), Tuple::Comparator(F.inner));
    vector<int> PS;
    for (int i = 0; i < tuples.size(); ++ i) {
        if (!fdominates(layer0, i, F)) {
            layer0.push_back(i);
            graph[i] = unordered_set<int>{};
            reverse_graph[i] = unordered_set<int>{};
            in_degree[i] = 0;
            if (tuples[i].prob > 0.5) {
                PS.push_back(i);
            }
        } else if (fdominates(PS, i, F)) {
            -- n;
            in_degree[i] = -1;
        } else {
            graph[i] = unordered_set<int>{};
            reverse_graph[i] = unordered_set<int>{};
            in_degree[i] = 0;
            for (int j = 0; j < i; ++ j) {
                if (in_degree[j] != -1 && tuples[j].fdominate(F.vertices, tuples[i])) {
                    graph[j].insert(i);
                    reverse_graph[i].insert(j);
                    ++ in_degree[i];
                }
            }
            if (tuples[i].prob > 0.5) {
                PS.push_back(i);
            }
        }
    }
    return n;
}

int MRSKY::path_reduction(const Region& F) {
    sort(tuples.begin(), tuples.end(), Tuple::Comparator(F.inner));
    vector<int> PS;
    vector<double> beta(tuples.size(), 0);
    for (int i = 0; i < tuples.size(); ++ i) {
        if (!fdominates(layer0, i, F)) {
            layer0.push_back(i);
            graph[i] = unordered_set<int>{};
            reverse_graph[i] = unordered_set<int>{};
            in_degree[i] = 0;
            beta[i] = tuples[i].prob/(1 - tuples[i].prob);
            if (tuples[i].prob > 0.5) {
                PS.push_back(i);
            }
        } else if (fdominates(PS, i, F)) {
            -- n;
            in_degree[i] = -1;
        } else {
            graph[i] = unordered_set<int>{};
            reverse_graph[i] = unordered_set<int>{};
            in_degree[i] = 0;
            if (tuples[i].prob > 0.5) {
                PS.push_back(i);
                for (int j = 0; j < i; ++ j) {
                    if (in_degree[j] != -1 && tuples[j].fdominate(F.vertices, tuples[i])) {
                        graph[j].insert(i);
                        reverse_graph[i].insert(j);
                        ++ in_degree[i];
                    }
                }
            } else {
                for (int j = 0; j < i; ++ j) {
                    if (tuples[j].fdominate(F.vertices, tuples[i])) {
                        if (in_degree[j] != -1) {
                            graph[j].insert(i);
                            reverse_graph[i].insert(j);
                            ++ in_degree[i];    
                        }
                        if (beta[i] < 1)
                            beta[i] = max(beta[i], beta[j]/(1 - tuples[i].prob));
                    }
                }
                if (beta[i] >= 1) {
                    PS.push_back(i);
                }
            }
        }
    }
    return n;
}

int MRSKY::region_reduction(const Region& F) {
    sort(tuples.begin(), tuples.end(), Tuple::Comparator(F.inner));
    vector<int> PS;
    for (int i = 0; i < tuples.size(); ++ i) {
        if (!fdominates(layer0, i, F)) {
            layer0.push_back(i);
            graph[i] = unordered_set<int>{};
            reverse_graph[i] = unordered_set<int>{};
            in_degree[i] = 0;
            if (tuples[i].prob > 0.5) {
                PS.push_back(i);
            }
        } else if (fdominates(PS, i, F)) {
            -- n;
            in_degree[i] = -1;
        } else {
            graph[i] = unordered_set<int>{};
            reverse_graph[i] = unordered_set<int>{};
            in_degree[i] = 0;
            if (tuples[i].prob > 0.5) {
                PS.push_back(i);
                for (int j = 0; j < i; ++ j) {
                    if (in_degree[j] != -1 && tuples[j].fdominate(F.vertices, tuples[i])) {
                        graph[j].insert(i);
                        reverse_graph[i].insert(j);
                        ++ in_degree[i];
                    }
                }
            } else {
                double r = 0;
                for (int j = 0; j < i; ++ j) {
                    if (tuples[j].fdominate(F.vertices, tuples[i])) {
                        if (in_degree[j] != -1) {
                            graph[j].insert(i);
                            reverse_graph[i].insert(j);
                            ++ in_degree[i];
                        }
                        if (r < 1) {
                            r = tuples[j].prob/(1 - tuples[j].prob);
                            r /= (1 - tuples[i].prob);
                            for (int k = j + 1; k < i; ++ k) {
                                if (tuples[j].fdominate(F.vertices, tuples[k])\
                                    && tuples[k].fdominate(F.vertices, tuples[i])) {
                                    r /= (1 - tuples[k].prob);
                                }
                            }
                        }
                    }
                }
                if (r >= 1) {
                    PS.push_back(i);
                }
            }
        }
    }
    return n;
}

int MRSKY::region_reduction_2d(const Region& F) {
    map<pair<double, double>, vector<double> > grid; // vector<double> [0] prob, [1] h, [2] v, [3] rho
    vector<double> X, Y;
    pair<double, double> p = make_pair(0, 0);
    for (int i = 0; i < tuples.size(); ++ i) {
        X.push_back(tuples[i][0]);
        Y.push_back(tuples[i][1]);
        p.first = tuples[i][0];
        for (int j = 0; j < tuples.size(); ++ j) {
            p.second = tuples[j][1];
            grid[p] = vector<double>(4, 0);
            if (i == j) grid[p][0] = tuples[i].prob;
            else grid[p][0] = 0;
        }
    }
    sort(X.begin(), X.end());
    sort(Y.begin(), Y.end());
    pair<double, double> l = make_pair(0, 0); // left to p
    for (auto y : Y) {
        p.second = y;
        l.second = y;
        for (int i = 0; i < X.size(); ++ i) {
            p.first = X[i];
            if (i == 0) grid[p][1] = 1;
            else {
                l.first = X[i - 1];
                grid[p][1] = grid[l][1]*(1 - grid[l][0]);
            }
        }
    }
    pair<double, double> b = make_pair(0, 0); // bottom to p
    for (int i = 0; i < X.size(); ++ i) {
        p.first = X[i];
        b.first = X[i];
        for (int j = 0; j < Y.size(); ++ j) {
            p.second = Y[j];
            if (j == 0) grid[p][2] = 1;
            else {
                b.second = Y[j - 1];
                grid[p][2] = grid[b][2]*(1 - grid[b][0]);
            }
            if (i == 0) {
                grid[p][3] = grid[p][2];
            } else if (j == 0) {
                grid[p][3] = grid[p][1];
            } else {
                l.first = X[i - 1];
                l.second = Y[j];
                grid[p][3] = grid[l][3]*(1 - grid[l][0])*grid[b][2]*(1 - grid[b][0]);
            }
        }
    }
    vector<int> PS;
    sort(tuples.begin(), tuples.end(), Tuple::Comparator(F.inner));
    for (int i = 0; i < tuples.size(); ++ i) {
        if (fdominates(PS, i, F)) {
            -- n;
        } else {
            double r = tuples[i].prob/(1 - tuples[i].prob);
            if (r >= 1) {
                PS.push_back(i);
            } else {
                pair<double, double> t = make_pair(tuples[i][0], tuples[i][1]);
                for (int j = 0; j < i; ++ j) {
                    if (tuples[j].fdominate(F.vertices, tuples[i])) {
                        pair<double, double> s = make_pair(tuples[j][0], tuples[j][1]);
                        pair<double, double> ts = make_pair(tuples[i][0], tuples[j][1]);
                        pair<double, double> st = make_pair(tuples[j][0], tuples[i][1]);
                        double r = grid[t][3]*grid[s][3]/(grid[ts][3]*grid[st][3]);
                        r *= grid[ts][1]/grid[s][1];
                        r *= grid[st][2]/grid[s][2];
                        if (r >= 1) {
                            cout << i << '\t' << j << '\t' << r << endl;
                            PS.push_back(i);
                            break;
                        }
                    }
                }
            }
        }
    }
    return n;
}

int MRSKY::virtual_region_reduction(const Region& F) {
    sort(tuples.begin(), tuples.end(), Tuple::Comparator(F.inner));
    vector<int> PS;
    vector<vector<double> > virtual_PS;
    vector<double> beta(tuples.size(), 0);
    for (int i = 0; i < tuples.size(); ++ i) {
        if (!fdominates(layer0, i, F)) {
            layer0.push_back(i);
            graph[i] = unordered_set<int>{};
            reverse_graph[i] = unordered_set<int>{};
            in_degree[i] = 0;
            beta[i] = tuples[i].prob/(1 - tuples[i].prob);
            if (tuples[i].prob > 0.5) {
                PS.push_back(i);
            } else { // look for a virtual point
                vector<double> virtual_point(F.vertices.size(), 0);
                double r = tuples[i].prob/(1 - tuples[i].prob);
                for (int j = i + 1; j < tuples.size(); ++ j) {
                    if (tuples[i].fdominate(F.vertices, tuples[j])) {
                        r /= (1 - tuples[j].prob);
                        for (int k = 0; k < F.vertices.size(); ++ k) {
                            virtual_point[k] = max(virtual_point[k], tuples[j].score(F.vertices[k]));
                        }
                        if (r >= 1) {
                            virtual_PS.push_back(virtual_point);
                            break;
                        }
                    }
                }
            }
        } else if (fdominates(PS, i, F)) {
            -- n;
            in_degree[i] = -1;
        } else if (fdominates(virtual_PS, i, F)) {
            -- n;
            in_degree[i] = -1;
        } else {
            graph[i] = unordered_set<int>{};
            reverse_graph[i] = unordered_set<int>{};
            in_degree[i] = 0;
            if (tuples[i].prob > 0.5) {
                PS.push_back(i);
                for (int j = 0; j < i; ++ j) {
                    if (in_degree[j] != -1 && tuples[j].fdominate(F.vertices, tuples[i])) {
                        graph[j].insert(i);
                        reverse_graph[i].insert(j);
                        ++ in_degree[i];
                    }
                }
            } else {
                for (int j = 0; j < i; ++ j) {
                    if (tuples[j].fdominate(F.vertices, tuples[i])) {
                        if (in_degree[j] != -1) {
                            graph[j].insert(i);
                            reverse_graph[i].insert(j);
                            ++ in_degree[i];    
                        }
                        if (beta[i] < 1)
                            beta[i] = max(beta[i], beta[j]/(1 - tuples[i].prob));
                    }
                }
                if (beta[i] >= 1) {
                    PS.push_back(i);
                } else {
                    double r = 0;
                    for (auto s : reverse_graph[i]) {
                        r = tuples[s].prob/(1 - tuples[s].prob);
                        r /= (1 - tuples[i].prob);
                        for (auto x : graph[s]) {
                            if (reverse_graph[i].find(x) != reverse_graph[i].end()) {
                                r /= (1 - tuples[x].prob);
                                if (r >= 1) {
                                    break;
                                }
                            }
                        }
                        if (r >= 1) {
                            PS.push_back(i);
                            break;
                        }
                    }
                }
            }
        }
    }
    /* for (int i = 0; i < tuples.size(); ++ i) {
        if (in_degree[i] != -1) {
            cout << i << '\t' << tuples[i] << endl;
            cout << "graph: ";
            for (auto t : graph[i]) {
                cout << t << '\t';
            }
            cout << endl;
            cout << "reverse_graph: ";
            for (auto t : reverse_graph[i]) {
                cout << t << '\t';
            }
            cout << endl;
        }
    } */
    return n;
}

int MRSKY::data_reduction_USKY(const Region& F) {
    vector<int> PS;
    sort(tuples.begin(), tuples.end(), Tuple::Comparator(F.inner));
    for (int i = 0; i < tuples.size(); ++ i) {
        if (fdominates(PS, i, F)) {
            --n;
        } else {
            double r = tuples[i].prob/(1 - tuples[i].prob);
            if (r >= 1) {
                PS.push_back(i);
            } else {
                for (int j = 0; j < i; ++ j) {
                    if (tuples[j].fdominate(F.vertices, tuples[i])) {
                        r = tuples[i].prob/(1 - tuples[i].prob);
                        r /= (1 - tuples[j].prob);
                        for (int k = j + 1; k < i; ++ k) {
                            if (tuples[j].fdominate(F.vertices, tuples[k])\
                             && tuples[k].fdominate(F.vertices, tuples[i])) {
                                r /= (1 - tuples[k].prob);
                            }
                        }
                        if (r >= 1) {
                            cout << i << '\t' << j << '\t' << r << endl;
                            PS.push_back(i);
                            break;
                        }
                    }
                }
            }
        }
    }
    return n;
}

void MRSKY::group_layer0() {
    vector<int> nodes;
    unordered_map<int, unordered_set<int> > edges;
    unordered_set<int> visited;
    for (auto t : layer0) {
        if (graph[t].size() == 0) {
            if (tuples[t].prob > 0.5) {
                res_at_layer0.push_back(t);
                prob_res_at_layer0 = prob_res_at_layer0*tuples[t].prob;
            } else {
                prob_res_at_layer0 = prob_res_at_layer0*(1 - tuples[t].prob);
            }
        } else {
            nodes.push_back(t);
            edges[t] = unordered_set<int>{};
        }
    }
    for (int i = 0; i < nodes.size(); ++ i) {
        for (int j = i + 1; j < nodes.size(); ++ j) {
            bool exists_edge = false;
            for (auto t : graph[nodes[i]]) {
                if (graph[nodes[j]].find(t) != graph[nodes[j]].end()) {
                    exists_edge = true;
                    break;
                }
            }
            if (exists_edge) {
                edges[nodes[i]].insert(nodes[j]);
                edges[nodes[j]].insert(nodes[i]);
            }
        }
    }
    for (auto u : nodes) {
        if (visited.find(u) == visited.end()) {
            groups_at_layer0.push_back(vector<int>{});
            stack<int> st;
            st.push(u);
            while (!st.empty()) {
                int v = st.top();
                st.pop();
                if (visited.find(v) == visited.end()) {
                    visited.insert(v);
                    groups_at_layer0.back().push_back(v);
                    for (auto w : edges[v]) {
                        if (visited.find(w) == visited.end()) {
                            st.push(w);
                        }
                    }
                }
            }
        }
    }
}

BigFloat MRSKY::DSA(vector<int>& result, const Region& F) {
    struct Candidate {
        vector<int> rskyline;
        BigFloat prob;
        double score;
    };

    BigFloat prob = prob_res_at_layer0;
    result.insert(result.end(), res_at_layer0.begin(), res_at_layer0.end());
    for (auto group : groups_at_layer0) {
        vector<int> partial_result;
        BigFloat partial_prob = BigFloat(0);
        priority_queue<pair<double, int>, vector<pair<double, int> >, greater<pair<double, int> > > sorted_tuples;
        stack<int> st;
        unordered_set<int> visited;
        for (auto t : group) {
            st.push(t);
        }
        while (!st.empty()) {
            int u = st.top();
            st.pop();
            if (visited.find(u) == visited.end()) {
                visited.insert(u);
                double min_score = dim;
                for (auto f : F.vertices) {
                    min_score = min(min_score, tuples[u].score(f));
                }
                sorted_tuples.push(make_pair(min_score, u));
                for (auto v : graph[u]) {
                    if (visited.find(v) == visited.end()) {
                        st.push(v);
                    }
                }
            }
        }
        vector<Candidate> candidates;
        candidates.push_back(Candidate{vector<int>{}, BigFloat(1), 1.0*dim});
        while (!sorted_tuples.empty()) {
            double min_score = sorted_tuples.top().first;
            int t = sorted_tuples.top().second;
            sorted_tuples.pop();
            double max_score = 0;
            for (auto f : F.vertices) {
                max_score = max(max_score, tuples[t].score(f));
            }
            int candidates_size = candidates.size();
            for (int i = 0; i < candidates_size; ++ i) {
                if (candidates[i].score <= min_score) {
                    if (candidates[i].prob > partial_prob) {
                        partial_prob = candidates[i].prob;
                        partial_result = candidates[i].rskyline;
                    }
                    swap(candidates[i], candidates[candidates_size - 1]);
                    swap(candidates[candidates_size - 1], candidates.back());
                    -- candidates_size;
                    -- i;
                    candidates.pop_back();
                } else if (!fdominates(candidates[i].rskyline, t, F)) {
                    Candidate new_candidate(candidates[i]);
                    new_candidate.rskyline.push_back(t);
                    new_candidate.prob = new_candidate.prob*tuples[t].prob;
                    if (new_candidate.prob > partial_prob) {
                        new_candidate.score = min(new_candidate.score, max_score);
                        candidates.push_back(new_candidate);
                        // this->num_candidates ++;
                    }
                    candidates[i].prob = candidates[i].prob*(1 - tuples[t].prob);
                    if (candidates[i].prob < partial_prob) {
                        swap(candidates[i], candidates[candidates_size - 1]);
                        swap(candidates[candidates_size - 1], candidates.back());
                        -- candidates_size;
                        -- i;
                        candidates.pop_back();
                    }
                }
            }
        }
        for (auto candidate : candidates) {
            if (candidate.prob > partial_prob) {
                partial_prob = candidate.prob;
                partial_result = candidate.rskyline;
            }
        }
        vector<Candidate>().swap(candidates);
        result.insert(result.end(), partial_result.begin(), partial_result.end());
        prob = prob*partial_prob;
    }
    // cout << this->num_candidates << endl;
    // check_prob(result, prob, F);
    return prob;
}

BigFloat MRSKY::IBBA(vector<int>& result, const Region& F) {
    BigFloat prob = prob_res_at_layer0;
    result.insert(result.end(), res_at_layer0.begin(), res_at_layer0.end());
    int i = 0;
    /* sort(groups_at_layer0.begin(), groups_at_layer0.end(), [](const vector<int>& x, const vector<int>& y){
        return x.size() < y.size();
    }); */
    for (auto group : groups_at_layer0) {
        // cout << group.size() << endl;
        vector<pair<int, double> > candidates;
        for (auto t : group) {
            candidates.push_back(make_pair(t, tuples[t].prob));
        }
        vector<int> partial_result;
        BigFloat partial_prob = BigFloat(0);
        IBBA_rec(vector<int>{}, 1, candidates, partial_result, partial_prob, 0);
        prob = prob*partial_prob;
        result.insert(result.end(), partial_result.begin(), partial_result.end());
    }
    // cout << "candidates num: " << this->num_candidates << endl;
    for (auto t : result) {
        cout << tuples[t] << endl;
    }
    return prob;
}

void MRSKY::IBBA_rec(const vector<int>& rskyline, const BigFloat& prob, vector<pair<int, double> > candidates, vector<int>& opt_result, BigFloat& opt_prob, const int& depth) {
    // this->num_candidates ++;
    if (candidates.empty()) {
        if (prob > opt_prob) {
            opt_prob = prob;
            opt_result = rskyline;
        }
        return;
    }
    sort(candidates.begin(), candidates.end(), [](const pair<int, double>& t, const pair<int, double>& s){
        return t.second > s.second;
    });
    vector<int> new_rskyline(rskyline);
    BigFloat new_prob = prob;
    while (!candidates.empty() && candidates.begin()->second > 0.5) {
        new_rskyline.push_back(candidates.begin()->first);
        new_prob = new_prob*candidates.begin()->second;
        candidates.erase(candidates.begin());
    }
    vector<pair<int, double> > new_candidates(candidates);
    for (int i = 0; i < candidates.size() + 1; ++ i) {
        if (i != 0) {
            int t = new_rskyline.back();
            new_rskyline.pop_back();
            new_prob = new_prob/tuples[t].prob;
            new_prob = new_prob*(1 - tuples[t].prob);
            for (auto s : graph[t]) {
                -- in_degree[s];
                if (in_degree[s] == 0) {
                    new_candidates.push_back(make_pair(s, tuples[s].prob));
                }
            }
        }
        if (i != candidates.size()) {
            new_candidates.erase(new_candidates.begin());
            new_rskyline.push_back(candidates[i].first);
            new_prob = new_prob*candidates[i].second;
        }
        BigFloat upper_prob = new_prob;
        for (auto t : new_candidates) {
            upper_prob = upper_prob*max(t.second, 1 - t.second);
        }
        if (upper_prob > opt_prob) {
            IBBA_rec(new_rskyline, new_prob, new_candidates, opt_result, opt_prob, depth + 1);
        }
    }
    for (auto t : candidates) {
        for (auto s : graph[t.first]) {
            ++ in_degree[s];
        }
    }
}

BigFloat MRSKY::EBBA(vector<int>& result, const Region& F) {
    BigFloat prob = prob_res_at_layer0;
    result.insert(result.end(), res_at_layer0.begin(), res_at_layer0.end());
    for (auto group : groups_at_layer0) {
        vector<pair<int, double> > candidates;
        for (auto t : group) {
            candidates.push_back(make_pair(t, tuples[t].prob));
        }
        vector<int> partial_result;
        BigFloat partial_prob = BigFloat(0);
        EBBA_rec(vector<int>{}, 1, candidates, partial_result, partial_prob, 0);
        prob = prob*partial_prob;
        result.insert(result.end(), partial_result.begin(), partial_result.end());
    }
    // cout << this->num_candidates << endl;
    return prob;
}

void MRSKY::EBBA_rec(const vector<int>& rskyline, const BigFloat& prob, vector<pair<int, double> > candidates, vector<int>& opt_result, BigFloat& opt_prob, const int& depth) {
    // this->num_candidates ++;
    if (candidates.empty()) {
        if (prob > opt_prob) {
            opt_prob = prob;
            opt_result = rskyline;
        }
        return;
    }
    sort(candidates.begin(), candidates.end(), [](const pair<int, double>& t, const pair<int, double>& s){
        return t.second > s.second;
    });
    vector<int> new_rskyline(rskyline);
    BigFloat new_prob = prob;
    for (auto t : candidates) {
        new_rskyline.push_back(t.first);
        new_prob = new_prob*t.second;
    }
    if (new_prob > opt_prob) {
        opt_prob = new_prob;
        opt_result = new_rskyline;
    }
    int i = candidates.size() - 1;
    vector<pair<int, double> > back_candidates;
    while (i >= 0 && candidates[i].second <= 0.5) {
        auto t = new_rskyline.back();
        new_rskyline.pop_back();
        new_prob = new_prob/tuples[t].prob;
        new_prob = new_prob*(1 - tuples[t].prob);
        vector<pair<int, double> > new_candidates(back_candidates);
        for (auto s : graph[t]) {
            -- in_degree[s];
            if (in_degree[s] == 0) {
                new_candidates.push_back(make_pair(s, tuples[s].prob));
            }
        }
        BigFloat upper_prob = new_prob;
        for (auto t : new_candidates) {
            upper_prob = upper_prob*max(t.second, 1 - t.second);
            // cout << depth << '\t' << upper_prob << '\t' << opt_prob << endl;
        }
        if (upper_prob > opt_prob) {
            EBBA_rec(new_rskyline, new_prob, new_candidates, opt_result, opt_prob, depth + 1);
        }
        new_prob = new_prob/(1 - tuples[t].prob);
        back_candidates.push_back(make_pair(t, tuples[t].prob));
        for (auto s : graph[t]) {
            ++ in_degree[s];
        }
        -- i;
    }
}

BigFloat MRSKY::half_rounding(vector<int>& result, const Region& F) {
    BigFloat prob = BigFloat(1);
    for (int i = 0; i < tuples.size(); ++ i) {
        if (tuples[i].prob > 0.5) {
            if (!fdominates(result, i, F)) {
                result.push_back(i);
                prob = prob*tuples[i].prob;
            }
        }
    }
    for (int i = 0; i < tuples.size(); ++ i) {
        if (!fdominates(result, i, F)) {
            prob = prob*(1 - tuples[i].prob);
        }
    }
    return prob;
}

int MRSKY::init_S(const Region& F, unordered_set<int>& S) {
    sort(tuples.begin(), tuples.end(), Tuple::Comparator(F.inner));
    vector<int> PS;
    for (int i = 0; i < tuples.size(); ++ i) {
        if (!fdominates(layer0, i, F)) {
            layer0.push_back(i);
            graph[i] = unordered_set<int>{};
            reverse_graph[i] = unordered_set<int>{};
            in_degree[i] = 0;
            if (tuples[i].prob > 0.5) {
                PS.push_back(i);
                S.insert(i);
            }
        } else if (fdominates(PS, i, F)) {
            -- n;
            in_degree[i] = -1;
        } else {
            graph[i] = unordered_set<int>{};
            reverse_graph[i] = unordered_set<int>{};
            in_degree[i] = 0;
            if (tuples[i].prob > 0.5) {
                PS.push_back(i);
                S.insert(i);
                for (int j = 0; j < i; ++ j) {
                    if (in_degree[j] != -1 && tuples[j].fdominate(F.vertices, tuples[i])) {
                        graph[j].insert(i);
                        reverse_graph[i].insert(j);
                        ++ in_degree[i];
                    }
                }
            } else {
                double r = 0;
                int t = -1;
                for (int j = 0; j < i; ++ j) {
                    if (tuples[j].fdominate(F.vertices, tuples[i])) {
                        if (in_degree[j] != -1) {
                            graph[j].insert(i);
                            reverse_graph[i].insert(j);
                            ++ in_degree[i];
                        }
                        if (r < 1) {
                            r = tuples[j].prob/(1 - tuples[j].prob);
                            r /= (1 - tuples[i].prob);
                            for (int k = j + 1; k < i; ++ k) {
                                if (tuples[j].fdominate(F.vertices, tuples[k])\
                                    && tuples[k].fdominate(F.vertices, tuples[i])) {
                                    r /= (1 - tuples[k].prob);
                                }
                            }
                            if (r >= 1) {
                                t = j;
                            }
                        }
                    }
                }
                if (r >= 1) {
                    PS.push_back(i);
                    S.insert(t);
                }
            }
        }
    }
    return n;
}

BigFloat MRSKY::LSA(vector<int>& result, const Region& F, const double& tau) {
    unordered_set<int> S;
    init_S(F, S);
    group_layer0();
    BigFloat prob = prob_res_at_layer0;
    result.insert(result.end(), res_at_layer0.begin(), res_at_layer0.end());
    // int i = 0;
    for (auto group : groups_at_layer0) {
        // cout << "group " << i ++ << endl;
        stack<int> st;
        unordered_set<int> nodes;
        for (auto t : group) {
            st.push(t);
        }
        while (!st.empty()) {
            int u = st.top();
            st.pop();
            if (nodes.find(u) == nodes.end()) {
                nodes.insert(u);
                for (auto v : graph[u]) {
                    if (nodes.find(v) == nodes.end()) {
                        st.push(v);
                    }
                }
            }
        }
        vector<pair<int, double> > sort_nodes;
        unordered_set<int> partial_result;
        BigFloat partial_prob = BigFloat(1);
        for (auto t : nodes) {
            if (S.find(t) != S.end()) {
                sort_nodes.push_back(make_pair(t, tuples[t].score(F.inner)));
            }
        }
        sort(sort_nodes.begin(), sort_nodes.end(), [](const pair<int, double>& t, const pair<int, double>& s){
            return t.second < s.second;
        });
        for (auto p : sort_nodes) {
            auto t = p.first;
            if (!fdominates(partial_result, t, F)) {
                partial_result.insert(t);
                partial_prob = partial_prob*tuples[t].prob;
            }
        }
        for (auto t : nodes) {
            if (!fdominates(partial_result, t, F)) {
                partial_prob = partial_prob*(1 - tuples[t].prob);
            }
        }
        unordered_map<int, int> count; // nodes - partial_result
        for (auto t : nodes) {
            count[t] = 0;
            if (partial_result.find(t) == partial_result.end()) {
                for (auto s : partial_result) {
                    if (graph[t].find(s) != graph[t].end()) { // t < s
                        -- count[t];
                    }
                    if (reverse_graph[t].find(s) != reverse_graph[t].end()) { // s < t
                        ++ count[t];
                    }
                }
            }
        }
        auto t = compute_ratio(partial_result, count, F);
        while (t.second > tau) {
            // cout << t.second << '\t' << count[t.first] << endl;
            if (count[t.first] > 0) {
                for (auto s : reverse_graph[t.first]) {
                    if (partial_result.erase(s)) {
                        for (auto x : graph[s]) {
                            -- count[x];
                        }
                        for (auto x : reverse_graph[s]) {
                            ++ count[x];
                        }
                    }
                }
            } else if (count[t.first] < 0) {
                for (auto s : graph[t.first]) {
                    if (partial_result.erase(s)) {
                        for (auto x : graph[s]) {
                            -- count[x];
                        }
                        for (auto x : reverse_graph[s]) {
                            ++ count[x];
                        }
                    }
                }
            }
            partial_prob = partial_prob*t.second;
            partial_result.insert(t.first);
            for (auto s : graph[t.first]) {
                ++ count[s];
            }
            for (auto s : reverse_graph[t.first]) {
                -- count[s];
            }
            count[t.first] = 0;
            t = compute_ratio(partial_result, count, F);
        }
        // cout << t.second << '\t' << t.first << endl;

        // post-deletion
        vector<int> temp;
        for (auto x : partial_result) {
            temp.push_back(x);
        }
        bool removable = true;
        for (auto t : temp) {
            removable = true;
            for (auto s : graph[t]) {
                if (count[s] == 1) {
                    removable = false;
                    break;
                }
            }
            if (removable && tuples[t].prob < 0.5) {
                partial_result.erase(t);
                partial_prob = partial_prob/tuples[t].prob;
                partial_prob = partial_prob*(1 - tuples[t].prob);
                for (auto s : graph[t]) {
                    -- count[s];
                }
                for (auto s : reverse_graph[t]) {
                    ++ count[s];
                }
            }
        }

        prob = prob*partial_prob;
        result.insert(result.end(), partial_result.begin(), partial_result.end());
    }
    /* for (auto t : result) {
        cout << tuples[t] << endl;
    } */
    return prob;
}

BigFloat MRSKY::LSA_PLUS(vector<int>& result, const Region& F, const double& tau) {
    unordered_set<int> S;
    init_S(F, S);
    group_layer0();
    BigFloat prob = prob_res_at_layer0;
    result.insert(result.end(), res_at_layer0.begin(), res_at_layer0.end());
    int i = 0;
    for (auto group : groups_at_layer0) {
        stack<int> st;
        unordered_set<int> nodes;
        for (auto t : group) {
            st.push(t);
        }
        while (!st.empty()) {
            int u = st.top();
            st.pop();
            if (nodes.find(u) == nodes.end()) {
                nodes.insert(u);
                for (auto v : graph[u]) {
                    if (nodes.find(v) == nodes.end()) {
                        st.push(v);
                    }
                }
            }
        }
        vector<pair<int, double> > sort_nodes;
        unordered_set<int> partial_result;
        // vector<int> temp;
        BigFloat partial_prob = BigFloat(1);
        for (auto t : nodes) {
            if (S.find(t) != S.end()) {
                // partial_result.insert(t);
                // partial_prob = partial_prob*tuples[t].prob;
                sort_nodes.push_back(make_pair(t, tuples[t].score(F.inner)));
            }
        }
        sort(sort_nodes.begin(), sort_nodes.end(), [](const pair<int, double>& t, const pair<int, double>& s){
            return t.second < s.second;
        });
        for (auto p : sort_nodes) {
            auto t = p.first;
            if (!fdominates(partial_result, t, F)) {
                partial_result.insert(t);
                // temp.push_back(t);
                partial_prob = partial_prob*tuples[t].prob;
            }
        }
        for (auto t : nodes) {
            if (!fdominates(partial_result, t, F)) {
                partial_prob = partial_prob*(1 - tuples[t].prob);
            }
        }
        /* unordered_map<int, int> count; // nodes - partial_result
        for (auto t : nodes) {
            count[t] = 0;
            if (partial_result.find(t) == partial_result.end()) {
                for (auto s : partial_result) {
                    if (graph[t].find(s) != graph[t].end()) { // t < s
                        -- count[t];
                    }
                    if (reverse_graph[t].find(s) != reverse_graph[t].end()) { // s < t
                        ++ count[t];
                    }
                }
            }
        }
        auto t = compute_ratio(partial_result, count, F);
        while (t.second > tau) {
            if (count[t.first] > 0) {
                for (auto s : reverse_graph[t.first]) {
                    if (partial_result.erase(s)) {
                        for (auto x : graph[s]) {
                            -- count[x];
                        }
                        for (auto x : reverse_graph[s]) {
                            ++ count[x];
                        }
                    }
                }
            } else if (count[t.first] < 0) {
                for (auto s : graph[t.first]) {
                    if (partial_result.erase(s)) {
                        for (auto x : graph[s]) {
                            -- count[x];
                        }
                        for (auto x : reverse_graph[s]) {
                            ++ count[x];
                        }
                    }
                }
            }
            partial_prob = partial_prob*t.second;
            partial_result.insert(t.first);
            for (auto s : graph[t.first]) {
                ++ count[s];
            }
            for (auto s : reverse_graph[t.first]) {
                -- count[s];
            }
            count[t.first] = 0;
            // check_count(partial_result, count, F);
            t = compute_ratio(partial_result, count, F);
        }

        // post-deletion
        vector<int> temp;
        for (auto x : partial_result) {
            temp.push_back(x);
        }
        bool removable = true;
        for (auto t : temp) {
            removable = true;
            for (auto s : graph[t]) {
                if (count[s] == 1) {
                    removable = false;
                    break;
                }
            }
            if (removable && tuples[t].prob < 0.5) {
                partial_result.erase(t);
                partial_prob = partial_prob/tuples[t].prob;
                partial_prob = partial_prob*(1 - tuples[t].prob);
                for (auto s : graph[t]) {
                    -- count[s];
                }
                for (auto s : reverse_graph[t]) {
                    ++ count[s];
                }
            }
        } */

        // find swap
        /* bool swap = true;
        while (swap) {
            swap = false;
            for (auto t : nodes) {
                if (partial_result.find(t) == partial_result.end()) {
                    for (auto s : partial_result) {
                        if (valid_swap(t, s, partial_result, count, F)) {
                            // swap s and t
                            partial_result.erase(s);
                            partial_prob = partial_prob/tuples[s].prob;
                            partial_prob = partial_prob*(1 - tuples[s].prob);
                            partial_prob = partial_prob*tuples[t].prob;
                            for (auto x : graph[s]) {
                                -- count[x];
                            }
                            for (auto x : reverse_graph[s]) {
                                ++ count[x];
                            }
                            for (auto x : graph[t]) {
                                ++ count[x];
                            }
                            for (auto x : reverse_graph[t]) {
                                -- count[x];
                            }
                            swap = true;
                        }
                    }
                }
            }
        } */
        prob = prob*partial_prob;
        result.insert(result.end(), partial_result.begin(), partial_result.end());
    }
    return prob;
}

void MRSKY::ARSP(vector<pair<int, double> >& result, const Region& F) {
    result.resize(tuples.size());
    sort(tuples.begin(), tuples.end(), Tuple::Comparator(F.inner));
    for (int i = 0; i < tuples.size(); ++ i) {
        result[i].first = i;
        result[i].second = tuples[i].prob;
        for (int j = 0; j < i; ++ j) {
            if (tuples[j].fdominate(F.vertices, tuples[i])) {
                // if (tuples[i].id == 568) cout << tuples[j].id << endl;
                result[i].second *= (1 - tuples[j].prob);
            }
        }
    }
    sort(result.begin(), result.end(), [](const pair<int, double>& t, const pair<int, double>& s){
        return t.second > s.second;
    });
    for (int i = 0; i <= 20; ++ i) {
        cout << tuples[result[i].first] << '\t' << result[i].first << '\t' << result[i].second << endl;
    }
}

bool MRSKY::fdominates(const vector<int>& rskyline, const int& t, const Region& F) {
    for (auto i : rskyline) {
        if (i == t) return true;
        if (tuples[i].fdominate(F.vertices, tuples[t])) {
            return true;
        }
    }
    return false;
}

bool MRSKY::fdominates(const unordered_set<int>& rskyline, const int& t, const Region& F) {
    for (auto i : rskyline) {
        if (i == t) return true;
        if (tuples[i].fdominate(F.vertices, tuples[t])) {
            return true;
        }
    }
    return false;
}

bool MRSKY::fdominates(const vector<vector<double> >& rskyline, const int& t, const Region& F) {
    for (auto s : rskyline) {
        if (dominates(s, t, F)) {
            return true;
        }
    }
    return false;
}

bool MRSKY::fdominates(const vector<double>& s, const int& t, const Region& F) {
    double score1 = 0, score2 = 0;
    bool unique = true;
    for (int i = 0; i < F.vertices.size(); ++ i) {
        for (int j = 0; j < dim; ++ j) {
            score1 += s[j]*F.vertices[i][j];
        }
        score2 = tuples[t].score(F.vertices[i]);
        if (score1 > score2) { 
            return false;
        }
        unique = score1 < score2;
    }
    return unique;
}

bool MRSKY::dominates(const vector<double>& s, const int& t, const Region& F) {
    bool unique = false;
    double score;
    for (int i = 0; i < F.vertices.size(); ++ i) {
        score = tuples[t].score(F.vertices[i]);
        if (s[i] > score) { 
            return false;
        } else if (s[i] < score) {
            unique = true;
        }
    }
    return unique;
}

void MRSKY::print_graph() {
    stack<int> st;
    unordered_set<int> visited;
    for (auto t : layer0) {
        st.push(t);
    }
    while (!st.empty()) {
        int u = st.top();
        st.pop();
        if (visited.find(u) == visited.end()) {
            cout << tuples[u].id << ',' << in_degree[u] << endl;
            visited.insert(u);
            for (auto v : graph[u]) {
                if (visited.find(v) == visited.end()) {
                    st.push(v);
                }
            }
        }
    }
}

pair<int, double> MRSKY::compute_ratio(const unordered_set<int>& rskyline, const unordered_map<int, int>& count, const Region& F) {
    double ratio_max = 0, ratio = 0;
    int t_max = 0;
    for (auto t : count) {
        if(rskyline.find(t.first) != rskyline.end()) continue;
        if (t.second == 0) {
            ratio = tuples[t.first].prob/(1 - tuples[t.first].prob);
            for (auto s : graph[t.first]) {
                if (count.at(s) == 0) {
                    ratio /= (1 - tuples[s].prob);
                }
            }
        } else if (t.second > 0) {
            unordered_set<int> temp(rskyline);
            ratio = tuples[t.first].prob;
            for (auto s : reverse_graph[t.first]) {
                if (rskyline.find(s) != rskyline.end()) {
                    ratio *= (1 - tuples[s].prob)/tuples[s].prob;
                    temp.erase(s);
                }
            }
            temp.insert(t.first);
            for (auto p : count) {
                if (p.second > 0) {
                    auto x = p.first;
                    if (!fdominates(temp, x, F)) {
                        ratio *= (1 - tuples[x].prob);
                    }
                }
            }
        } else {
            ratio = tuples[t.first].prob/(1 - tuples[t.first].prob);
            for (auto s : graph[t.first]) {
                if (rskyline.find(s) != rskyline.end()) {
                    ratio /= tuples[s].prob;
                } else if (count.at(s) <= 0) {
                    ratio /= (1 - tuples[s].prob);
                }
            }
        }
        if (ratio > ratio_max) {
            ratio_max = ratio;
            t_max = t.first;
        }
    }
    // cout << ratio_max << endl;
    return make_pair(t_max, ratio_max);
}

void MRSKY::check_count(const unordered_set<int>& partial_result, const unordered_map<int, int>& count, const Region& F) {
    for (auto iter : count) {
        int num = 0;
        if (partial_result.find(iter.first) == partial_result.end()) {
            for (auto s : partial_result) {
                if (tuples[s].fdominate(F.vertices, tuples[iter.first])) {
                    ++ num;
                } else if (tuples[iter.first].fdominate(F.vertices, tuples[s])) {
                    -- num;
                }
            }
        }
        if (num != iter.second) {
            cout << "wrong count: " << iter.first << '\t' << num << '\t' << iter.second << endl;
        }
    }
}

void MRSKY::check_prob(const vector<int>& rskyline, const BigFloat& prob, const Region& F) {
    BigFloat check_prob = BigFloat(1);
    unordered_map<int, bool> in_rskyline;
    for (auto t : rskyline) {
        check_prob = check_prob*tuples[t].prob;
        in_rskyline[t] = true;
    }
    for (int i = 0; i < tuples.size(); ++ i) {
        if (in_rskyline.find(i) == in_rskyline.end()) {
            if (!fdominates(rskyline, i, F)) {
                if (in_degree[i] == -1) {
                    cout << tuples[i].id << " wrong\n";
                }
                check_prob = check_prob*(1 - tuples[i].prob);
            }
        }
    }
    cout << "compute prob: " << prob << ", check prob: " << check_prob << endl;
}

bool MRSKY::valid_swap(const int& t, const int& s, const unordered_set<int>& rskyline, const unordered_map<int, int>& count, const Region& F) {
    if (count.at(t) != 0) return false;

    
}