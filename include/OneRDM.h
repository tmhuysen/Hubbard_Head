//
// Created by Wulfix on 19/11/2017.
//

#ifndef HUBBARDCLEAN_ONERDM_H
#define HUBBARDCLEAN_ONERDM_H

#include <boost/math/special_functions.hpp>
#include <vector>
#include <AddressingMatrix.h>
#include <armadillo>
#include "bitset.h"

class OneRDM {
public:
    OneRDM(std::vector<double> coefs, unsigned nup, unsigned ndown, unsigned sites);
    void print();




private:
    double precision;
    std::vector<double> coefs;
    unsigned nup;
    unsigned ndown;
    unsigned sites;
    unsigned long n_bf_alpha;
    unsigned long n_bf_beta;

    unsigned long oneRDMDim;
    AddressingMatrix ad_a;
    AddressingMatrix ad_b;

    arma::mat rdmAlpha;
public:
    const arma::mat &getRdmAlpha() const;

    const arma::mat &getRdmBeta() const;
private:
    arma::mat rdmBeta;


    void rdmCalc();
    double oneRDMAlpha(unsigned long i, unsigned long j);
    double oneRDMBeta(unsigned long i, unsigned long j);
};


#endif //HUBBARDCLEAN_ONERDM_H
