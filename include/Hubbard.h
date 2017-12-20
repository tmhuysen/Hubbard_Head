//
// Created by Wulfix on 16/11/2017.
//

#ifndef HUBBARDCLEAN_HUBBARD_H
#define HUBBARDCLEAN_HUBBARD_H


#include "Lattice.h"
#include "AddressingMatrix.h"
#include <boost/math/special_functions.hpp>
#include "io.h"

struct State {
    double eigenValue;          // The energy of the solution, a.k.a. the eigenvalue
    arma::vec eigenVector; // The coefficients of the solution with respect to the given basis, a.k.a. the eigenvector corresponding to the eigenvalue
    int S_z;            // Total projected spin !!!TIMES 2!!!
    //      For total projected spin = 1/2, S_z = 1
    //      For total projected spin = -1, S_z = -2
};



class Hubbard
{
    //private classes
private:
    /**
     * nesting to prevent circular inclusion. A hubbard instance has fixed set of parameters,
     * these can be used to calculate different spin sectors,
     * by distributing the total amount of electrons over up and down spins.
     *
     * A nested class can access the private members of the outer class if an instance of this class
     * is given as a member of our
     */
    class SpinSector {
        //private variable members
    private:
        Hubbard& p; //parent
        int S_z;
        arma::vec eigenvalues;
        arma::mat eigenvectors;

        // Total projected spin !!!TIMES 2!!!
        //      For total projected spin = 1/2, S_z = 1
        //      For total projected spin = -1, S_z = -2
        unsigned N_alpha;   // Number of alpha electrons
        unsigned N_beta;    // Number of beta electrons

        unsigned long n_bf; // number of total basis functions
        unsigned long n_bf_alpha;
        unsigned long n_bf_beta;

        arma::mat hamiltonian;
        //public methods
    public:
        /**
         * getters and setters.
         */
        const arma::mat & getHamiltonian() const;
        void setEigenvectors(const arma::mat &eigenvectors);
        void setEigenvalues(const arma::vec &eigenvalues);
        int getS_z() const;


        /** Default constructor
         *
         */
        //SpinSector();

        /** Create a spin sector with
         *      S_z:    integer representation of total projected spin (times 2)
         *                  total projected spin = 1/2 -> S_z = 1
         *                  total projected spin = -1 -> S_z = -2
         */
        SpinSector(Hubbard& p, int S_z);


        /** Helper function to print the Hamiltonian belonging to the spin sector
        */
        void print_hamiltonian() const;

        //private methods
    private:

        /**
         * Calculates the hamiltonian for given sector, it has a start and end to easily implement
         * parallelisation
         */
        void calculateSector(double start, double end);

        /**
         *
         * This function has been created to easily be able to implement different libraries
         * without having to change multiple istances of adding to matrix of a given library.
         *
         */
        void addToHamiltonian(double value, size_t index1, size_t index2);
        /*
         * idem
         */
        void symmetryFill();


    };

    //private members
private:
    unsigned N;
public:
    Lattice &getLattice() const;

private:
    // The number of electrons to place in the lattice
    Lattice& lattice;
public:
    void setLattice(Lattice &lattice);

private:
    // REFERENCE to a lattice. A reference is used to avoid copying.

    // Hubbard will serve as SpinSector factory.
    // this->spin_sectors is a std::vector of SpinSectors.
    //      They are ordered with increasing total projected spin
    //          spin_sectors[0]     is the spin sector with the lowest number of alpha electrons (0)
    //          spin_sectors[N]   is the spin sector with the highest number of alpha electrons (N)
    std::vector<SpinSector> spinSectors;


    // For each distribution of electrons we will need a different addressingsMatrix these may be same for several spin sectors
    // So we storing them allows us to only generate them once.
    // Currently not used because of error.
    std::vector<AddressingMatrix> addressing_list;
    std::vector<State> groundstates;
    //public methods
public:
    /**
     * getters & setters
     */
    std::vector<State> getGroundstates();
    arma::mat getSector(int index);

    /**
     * Constructor
     */
    Hubbard(unsigned N, Lattice& lattice);

    void print();

    //private methods
private:
    /**
     * creates all spinsectors
     */
    void calculateSectors();
    /**
     * creates all AddressingMatrixes for the given electrons and lattice dimension.
     */
    void generateAddressingMatrix();
    /**
     * This will simply diagonalize the hamiltonian
     * of a given sector and check if it has the lowest groundstates.
     */
    void solveSectors();
    /** Adds a state to the groundstate list
    *     if the state is lower than the current energy of all assumed "groundstates"
    *     the list is cleared and the new lowest state is added.
    *     if the state is equal in energy it is added as a degenerate state.
    */
    void groundStates(State state);




};

//Compare tools for the State struct.
bool compareState(const State &o1, const State &o2);
bool areSame(const State &o1, const State &o2);


#endif //HUBBARDCLEAN_HUBBARD_H
