//
// Created by Wulfix on 19/11/2017.
//

#include <utility>

#include "OneRDM.h"


OneRDM::OneRDM(std::vector<double> coefs, unsigned nup, unsigned ndown, unsigned sites) {
    precision = 1e-15;

    this -> nup = nup;
    this -> ndown = ndown;
    this -> sites = sites;
    this -> coefs = std::move(coefs);
    this -> ad_a = AddressingMatrix(sites,nup);
    this -> ad_b = AddressingMatrix(sites,ndown);

    auto n_bf_alpha_ = boost::math::binomial_coefficient<double>(this->sites, nup);
    if (n_bf_alpha_ > 4294967295.0) {
        // before casting into unsigned long, we have to make sure that it can fit
        throw std::overflow_error("The number of basis functions for the separated sector is too high to be cast into unsigned long.");
    }
    this->n_bf_alpha = static_cast<unsigned long>(n_bf_alpha_);

    auto n_bf_beta_ = boost::math::binomial_coefficient<double>(sites, ndown);
    if (n_bf_beta_ > 4294967295.0) {
        // before casting into unsigned long, we have to make sure that it can fit
        throw std::overflow_error("The number of basis functions for the separated sector is too high to be cast into unsigned long.");
    }

    this->n_bf_beta = static_cast<unsigned long>(n_bf_beta_);
    rdmBeta = arma::zeros(sites,sites);
    rdmAlpha = arma::zeros(sites,sites);

    rdmCalc();



}

/**
 * alpha and beta evaluations for each site,
 * The calculations are separated into two functions,
 * because a single function would be far less elegant.
 *
 */

void OneRDM::rdmCalc(){
    for (unsigned long i = 0; i<sites; i ++){
        for(unsigned long j = i; j<sites; j++){

            rdmAlpha(i,j) = oneRDMAlpha(i,j);
            rdmBeta(i,j) = oneRDMBeta(i,j);


        }
    }
    rdmAlpha = arma::symmatu(rdmAlpha);
    rdmBeta = arma::symmatu(rdmBeta);

}


/**
 *
 * @param i annihilation operator
 * @param j creation operator
 * @return return the sum of coefficient pairs for the evaluation
 */

double OneRDM::oneRDMAlpha(unsigned long i, unsigned long j) {
    double coefsum = 0;
    //incase we want to parallelise on this level
    boost::dynamic_bitset<> alpha_set = ad_a.generateBinaryVector(0);

    for (unsigned long l = 0; l<n_bf_alpha ; l++ ){
        boost::dynamic_bitset<> alpha_target = boost::dynamic_bitset<>(alpha_set);
        //operator evaluations
        if (annihilation(alpha_target, i) && creation(alpha_target, j)) {
            double phase_factor = phaseCheck(alpha_set, alpha_target);
            unsigned long address = ad_a.fetchAddress(alpha_target);
            //evaluations doesn't change since operators only affect one spin
            //Hence the evaluations can be kept for while iterating over the total basis
            //and adjusting the addresses accordingly

            for(unsigned long k = 0; k<n_bf_beta ; k++){
                double element_set = coefs.at(l+k*n_bf_beta);
                double element_target = coefs.at(address + k*n_bf_beta);
                coefsum += element_target*element_set*phase_factor;
            }
        }

        alpha_set = next_bitset_permutation(alpha_set);

    }

    if (abs(coefsum) < precision){
        coefsum = 0;
    }

    return coefsum;
}

/**
 *
 * @param i annihilation operator
 * @param j creation operator
 * @return return the sum of coefficient pairs for the evaluation
 */


double OneRDM::oneRDMBeta(unsigned long i, unsigned long j){
    double coefsum = 0;
    boost::dynamic_bitset<> beta_set = ad_b.generateBinaryVector(0);
    for (unsigned long l = 0; l<n_bf_beta ; l++ ){
        boost::dynamic_bitset<> beta_target = boost::dynamic_bitset<>(beta_set);
        if (annihilation(beta_target, i) && creation(beta_target, j)) {
            double phase_factor = phaseCheck(beta_set, beta_target);
            unsigned long address = ad_b.fetchAddress(beta_target);
            for(unsigned long k = 0; k<n_bf_alpha ; k++){
                double element_target = coefs.at(address*n_bf_beta+k);
                double element_set = coefs.at(l*n_bf_beta+k);
                coefsum += element_target*element_set*phase_factor;
            }
        }
        beta_set = next_bitset_permutation(beta_set);

    }
    if (abs(coefsum) < precision){
        coefsum = 0;
    }
    return coefsum;
}

void OneRDM::print() {
    std::cout<<std::endl<<rdmAlpha<<std::endl<<rdmBeta;


}

const arma::mat &OneRDM::getRdmAlpha() const {
    return rdmAlpha;
}

const arma::mat &OneRDM::getRdmBeta() const {
    return rdmBeta;
}


