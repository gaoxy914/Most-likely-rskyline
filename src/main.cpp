#include "MRSKY.h"

int main(int argc, char const *argv[]) {
    if (strcmp(argv[1], "genDQ") == 0) { // gene all datasets and queries for exps
        for (int d = 2; d <= 10; d += 1) {
            for (double alpha = 0.2; alpha <= 0.8; alpha += 0.15) {
                MRSKY mrsky1(d, 4000000, alpha);
                mrsky1.gen_ind_data();
                mrsky1.write_data("data/ind/");
                MRSKY mrsky2(d, 4000000, alpha);
                mrsky2.gen_anti_data();
                mrsky2.write_data("data/anti/");
            }
            Region region1(d, d - 1);
            region1.gen_weak_query();
            region1.write_query("query/weak/");
            Region region2(d, 7);
            region2.gen_inter_query();
            region2.write_query("query/inter/");
        }
    } else if (strcmp(argv[1], "cmpRR") == 0) {
        int n = atoi(argv[3]);
        int dim = atoi(argv[4]);
        double alpha = atof(argv[5]);
        int c = atoi(argv[7]);
        string data_path = "data/";
        data_path += string(argv[2]);
        data_path += "/";
        string query_path = "query/";
        query_path += string(argv[6]);
        query_path += "/";
        Region F = Region(dim, c);
        F.load_query(query_path.c_str());
        F.compute_vertex();
        F.print_query();
        long long algmtime, reducemtime, seconds, useconds;
        struct timeval start, end;   
        MRSKY mrsky1(dim, n, alpha);
        mrsky1.load_data(data_path.c_str());
        gettimeofday(&start, nullptr);
        int nReduce = mrsky1.path_reduction(F);
        gettimeofday(&end, nullptr);
        seconds = end.tv_sec - start.tv_sec;
        useconds = end.tv_usec - start.tv_usec;
        reducemtime = seconds*1000000 + useconds;
        cout << "-------------------------------" << endl;
        cout << "n: " << n << " d: " << dim << " alpha: " << alpha << endl;
        cout << "path reduction time: " << reducemtime << " data size: " << nReduce << " ratio: " << 1.0*nReduce/n << endl;
        
        MRSKY mrsky2(dim, n, alpha);
        mrsky2.load_data(data_path.c_str());
        gettimeofday(&start, nullptr);
        nReduce = mrsky2.region_reduction(F);
        gettimeofday(&end, nullptr);
        seconds = end.tv_sec - start.tv_sec;
        useconds = end.tv_usec - start.tv_usec;
        reducemtime = seconds*1000000 + useconds;
        cout << "-------------------------------" << endl;
        cout << "n: " << n << " d: " << dim << " alpha: " << alpha << endl;
        cout << "region reduction time: " << reducemtime << " data size: " << nReduce << " ratio: " << 1.0*nReduce/n << endl;
        cout << "-------------------------------" << endl;
    } else {
        int n = atoi(argv[3]);
        int dim = atoi(argv[4]);
        double alpha = atof(argv[5]);
        int c = atoi(argv[7]);
        MRSKY mrsky(dim, n, alpha);
        string data_path = "data/";
        string query_path = "query/";
        data_path += string(argv[2]);
        data_path += "/";
        if (strcmp(argv[2], "car") == 0) {
            mrsky.load_car_data();
        } else if (strcmp(argv[2], "nba") == 0) {
            mrsky.load_nba_data();
        } else {
            mrsky.load_data(data_path.c_str());
        }
        Region F = Region(dim, c);
        query_path += string(argv[6]);
        query_path += "/";
        F.load_query(query_path.c_str());
        F.compute_vertex();
        F.print_query();
        BigFloat prob;
        vector<int> result;
        long long algmtime, reducemtime, seconds, useconds;
        struct timeval start, end;     
        if (strcmp(argv[1], "LSA") == 0) {
            double tau = atof(argv[8]);
            gettimeofday(&start, nullptr);
            prob = mrsky.LSA(result, F, tau);
            gettimeofday(&end, nullptr);
            seconds = end.tv_sec - start.tv_sec;
            useconds = end.tv_usec - start.tv_usec;    
            algmtime = seconds*1000000 + useconds;
            cout << "-------------------------------" << endl;
            cout << "n: " << n << " d: " << dim << " alpha: " << alpha << " c: " << c << endl;
            cout << argv[1] << " prob: " << prob << " total time: " << algmtime << endl;
            mrsky.check_prob(result, prob, F);
            cout << "-------------------------------" << endl;
        } else if (strcmp(argv[1], "LSA+") == 0) {
            double tau = atof(argv[8]);
            gettimeofday(&start, nullptr);
            prob = mrsky.LSA_PLUS(result, F, tau);
            gettimeofday(&end, nullptr);
            seconds = end.tv_sec - start.tv_sec;
            useconds = end.tv_usec - start.tv_usec;    
            algmtime = seconds*1000000 + useconds;
            cout << "-------------------------------" << endl;
            cout << "n: " << n << " d: " << dim << " alpha: " << alpha << " c: " << c << endl;
            cout << argv[1] << " prob: " << prob << " total time: " << algmtime << endl;
            mrsky.check_prob(result, prob, F);
            cout << "-------------------------------" << endl;
        } else {
            gettimeofday(&start, nullptr);
            int nReduce = mrsky.virtual_region_reduction(F);
            mrsky.group_layer0();
            gettimeofday(&end, nullptr);
            seconds = end.tv_sec - start.tv_sec;
            useconds = end.tv_usec - start.tv_usec;
            reducemtime = seconds*1000000 + useconds;
            cout << "-------------------------------" << endl;
            cout << "n: " << n << " d: " << dim << " alpha: " << alpha << " c: " << c << endl;
            cout << "reduction time: " << reducemtime << " data size: " << nReduce << " ratio: " << 1.0*nReduce/n << endl;
            if (strcmp(argv[1], "DSA") == 0) {
                gettimeofday(&start, nullptr);
                prob = mrsky.DSA(result, F);
                gettimeofday(&end, nullptr);
                seconds = end.tv_sec - start.tv_sec;
                useconds = end.tv_usec - start.tv_usec;
                algmtime = seconds*1000000 + useconds;
            } else if (strcmp(argv[1], "IBBA") == 0) {
                gettimeofday(&start, nullptr);
                prob = mrsky.IBBA(result, F);
                gettimeofday(&end, nullptr);
                seconds = end.tv_sec - start.tv_sec;
                useconds = end.tv_usec - start.tv_usec;
                algmtime = seconds*1000000 + useconds;
            } else if (strcmp(argv[1], "EBBA") == 0) {
                gettimeofday(&start, nullptr);
                prob = mrsky.EBBA(result, F);
                gettimeofday(&end, nullptr);
                seconds = end.tv_sec - start.tv_sec;
                useconds = end.tv_usec - start.tv_usec;
                algmtime = seconds*1000000 + useconds;
            } else {
                cout << "wrong algorithm.\n";
                exit(1);
            }
            cout << argv[1] << " time: " << algmtime << " prob: " << prob << " total time: " << algmtime + reducemtime << endl;
            mrsky.check_prob(result, prob, F);
            cout << "-------------------------------" << endl;
        }
    }
    return 0;
}