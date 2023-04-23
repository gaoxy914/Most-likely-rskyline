#ifndef __BIGFLOAT_H__
#define __BIGFLOAT_H__

#include <iomanip>
#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <deque>
#include <random>
#include <fstream>
#include <map>
#include <unordered_map>
#include <cstring>
#include <sys/time.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sstream>
#include <queue>
#include <iostream>
#include <set>
#include <stack>
#include <unordered_set>

#define PREC 6

using namespace std;

double roundoff(const double& value, unsigned char prec);
double rand_uniform(const double& a, const double& b);
double rand_normal(const double& mean, const double& var);

class BigFloat {
private:
    int exponent;
    double decimal;

public:

    BigFloat();
    BigFloat(const double& num);
    BigFloat(const BigFloat& other);
    BigFloat& operator= (const BigFloat& other);
    BigFloat& operator= (const double& num);
    bool operator> (const BigFloat& other) const;
    bool operator< (const BigFloat& other) const;
    friend BigFloat operator* (const BigFloat& left, const BigFloat& right);
    friend BigFloat operator* (const BigFloat& left, const double& right);
    friend BigFloat operator/ (const BigFloat& left, const BigFloat& right);
    friend BigFloat operator/ (const BigFloat& left, const double& right);
    friend ostream& operator<< (ostream& out, const BigFloat& num);
};

#endif