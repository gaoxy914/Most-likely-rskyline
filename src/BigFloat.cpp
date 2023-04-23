#include "BigFloat.h"

double roundoff(const double& value, unsigned char prec) {
    double pow_10 = pow(10.0, (double)prec);
    return round(value*pow_10)/pow_10;
}
double rand_uniform(const double& a, const double& b) {
    return roundoff(drand48()*(b - a) + a, PREC);
}

double rand_normal(const double& mean, const double& var) {
    double x, y, s;
    do {
        x = drand48();
        if (rand()%2) x = -x;
        y = drand48();
        if (rand()%2) y = -y;
        s = x*x + y*y;
    } while (s >= 1);
    double num = x*sqrt(-2.*log(s)/s);
    num = num*var + mean;
    return roundoff(num, PREC);
}

BigFloat::BigFloat() {
    exponent = 0;
    decimal = 0;
}

BigFloat::BigFloat(const double& num) {
    // assert(num >= 0 && num <= 1);
    if (num == 0) {
        decimal = 0;
        exponent = -1;
    } else if (num == 1) {
        decimal = 1;
        exponent = -1;
    } else {
        decimal = num;
        exponent = 0;
        while (decimal < 0.1) {
            ++ exponent;
            decimal *= 10;
        }
    }
}
 
BigFloat::BigFloat(const BigFloat& other) {
    exponent = other.exponent;
    decimal = other.decimal;
}

BigFloat& BigFloat::operator= (const BigFloat& other) {
    if (&other != this) {
        exponent = other.exponent;
        decimal = other.decimal;
    }
    return *this;
}
 
BigFloat& BigFloat::operator= (const double& num) {
    *this = BigFloat(num);
    return *this;
}
 
bool BigFloat::operator> (const BigFloat& other) const {
    if (other.exponent == -1) {
        if (other.decimal == 0) return true;
        else if (other.decimal == 1) return false;
    }
    if (exponent == other.exponent) {
        return decimal > other.decimal;
    }
    return exponent < other.exponent;
}

bool BigFloat::operator< (const BigFloat& other) const {
    if (other.exponent == -1) {
        if (other.decimal == 0) return false;
        else if (other.decimal == 1) return true;
    }
    if (exponent == other.exponent) {
        return decimal < other.decimal;
    }
    return exponent > other.exponent;
}
 
BigFloat operator* (const BigFloat& left, const BigFloat& right) {
    if (right.exponent == -1 && right.decimal == 0) return BigFloat(0);
    if (left.exponent == -1 && left.decimal == 0) return BigFloat(0);
    if (right.exponent == -1 && right.decimal == 1) return left;
    if (left.exponent == -1 && left.decimal == 1) return right;
    BigFloat res(left);
    res.exponent += right.exponent;
    res.decimal *= right.decimal;
    while (res.decimal < 0.1) {
        ++ res.exponent;
        res.decimal *= 10;
    }
    while (res.decimal > 1) {
        -- res.exponent;
        res.decimal /= 10;
    }
    return res;
}
 
BigFloat operator* (const BigFloat& left, const double& right) {
    return left*BigFloat(right);
}

BigFloat operator/ (const BigFloat& left, const BigFloat& right) {
    assert(left.exponent >= right.exponent);
    BigFloat res(left);
    res.exponent -= right.exponent;
    res.decimal /= right.decimal;
    while (res.decimal < 0.1) {
        ++ res.exponent;
        res.decimal *= 10;
    }
    while (res.decimal > 1) {
        -- res.exponent;
        res.decimal /= 10;
    }
    assert(res.decimal >= 0);
    return res;
}

BigFloat operator/ (const BigFloat& left, const double& right) {
    return left/BigFloat(right);
}

ostream& operator<< (ostream& out, const BigFloat& num) {
    if (num.exponent == -1) {
        if (num.decimal == 1) out << "1";
        else if (num.decimal == 0) out << "0";
    } else {
        if (num.exponent > 0) out << num.decimal << "*10^" << num.exponent;
        else out << num.decimal;
    }
    return out;
}