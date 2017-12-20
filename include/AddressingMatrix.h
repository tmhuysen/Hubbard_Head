//
// Created by Wulfix on 17/11/2017.
//

#ifndef HUBBARDMODEL_ADDRESSINGMATRIX_H
#define HUBBARDMODEL_ADDRESSINGMATRIX_H

#include "bitset.h"
#include <iostream>
#include <Eigen/Dense>
class AddressingMatrix {
private:
    int n_sites;
    int n_electrons;
    Eigen::MatrixXi addressMatrix;

    void generateMatrix();

public:
    /**
     * Constructors
     */
    AddressingMatrix()= default;
    AddressingMatrix(unsigned int n_sites, unsigned int n_electrons);

    /**
     *
     * @param bitVector enter a bitset
     * @return get the order of the bitset based on permutations, of a total of n_sites bits with fixed n_electrons 1 bits.
     *
     */
    unsigned long fetchAddress(boost::dynamic_bitset<> bitVector);

    /**
     *
     * @param address of desired bitset
     * @return bitset based on permutations of a bitset with n_sites bits and fixed n_electrons 1 bits.
     */
    boost::dynamic_bitset<> generateBinaryVector(unsigned long address);
    /**
     * prints the matrix (tests)
     */
    void print();




};


#endif //HUBBARDMODEL_ADDRESSINGMATRIX_H
