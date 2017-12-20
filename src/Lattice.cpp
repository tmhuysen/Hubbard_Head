//
// Created by Wulfix on 16/11/2017.
//

#include "Lattice.h"


Lattice::Lattice(const arma::mat &H) {
    this->H = H;
    this->L = static_cast<unsigned>(H.n_cols);
}

Lattice::Lattice(const arma::umat &A, double t, double U) {
    this->L = static_cast<unsigned>(A.n_cols);
    this->H = U * arma::eye(this->L, this->L) - t * A;
}

unsigned Lattice::getLength()const {
    return this->L;
}

double Lattice::getElement(unsigned site1, unsigned site2)const {
    return H(site1,site2);
}

const arma::mat &Lattice::getH() const {
    return H;
}
