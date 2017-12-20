//
// Created by Wulfix on 16/11/2017.
//

#ifndef HUBBARDCLEAN_LATTICE_H
#define HUBBARDCLEAN_LATTICE_H


#include <armadillo>

class Lattice {
private:
    unsigned L;     // The number of lattice sites
    arma::mat H;
public:
    const arma::mat &getH() const;
    // The hopping matrix representing the values of the on-site interaction and the hopping parameters

public:
    /** Constructor based on the adjacency matrix A, and uniform t and U parameters.
     *  What is meant by uniform parameters, is that the interaction on every site has the same U, and every hopping between two neighbouring sites is equal to t
     */
    Lattice(const arma::umat& A, double t, double U);


    /** Constructor based on the hopping matrix H
     */
    Lattice(const arma::mat& H);
    Lattice(unsigned dim):H (arma::zeros(dim,dim)) {this->L=dim;};

    /**
     * Ask for
     */
    double getElement(unsigned site1, unsigned site2)const;
    unsigned getLength()const;



};



#endif //HUBBARDCLEAN_LATTICE_H
