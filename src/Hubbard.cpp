//
// Created by Wulfix on 16/11/2017.
//

#include "Hubbard.h"

Hubbard::Hubbard(unsigned N, Lattice &lattice): lattice(lattice) {
    if (N > (2 * this->lattice.getLength())) {

        throw std::invalid_argument("The lattice is not compatible with the given number of electrons.");
    } else {
        this->N = N;
    }
    generateAddressingMatrix();

    calculateSectors();
    solveSectors();
}

void Hubbard::generateAddressingMatrix() {
    for(size_t electrons = 0; electrons<N+1;electrons++){
        //emplace takes constructor parameters and constructs and instance of <Class> of the vector.
        addressing_list.emplace_back(lattice.getLength(),electrons);
    }

}

void Hubbard::calculateSectors() {
    spinSectors = std::vector<SpinSector>();
    for(unsigned i= 0; i<N+1;i++){
        spinSectors.emplace_back(*this,-static_cast<int>(N) + 2 * i);
    }


}


void Hubbard::print(){
    for(const auto &x: spinSectors){
        std::cout<<std::endl;
        x.print_hamiltonian();
        std::cout<<std::endl;
        std::cout<<x.getS_z();
    }





}

std::vector<State> Hubbard::getGroundstates(){
        return groundstates;
}
/**
 * changes the lattice of a given hubbard model, keeping all other parameters intact.
 * One must only set a new lattice if this lattice has the same dimension as the previous lattice
 * otherwise one should just create a new hubbard instance.
 */
void Hubbard::setLattice(Lattice &lattice) {
    if (this->lattice.getLength() != lattice.getLength()){
        throw std::logic_error("lattice with different dimensions warrants a complete new hubbard instance");
    }else{
        Hubbard::lattice = lattice;
        calculateSectors();
        solveSectors();
    }

}

void Hubbard::solveSectors() {
    groundstates = { State {std::numeric_limits<double>::max(),arma::vec()} };
    // Diagonalize the Hamiltonian for every spin sector, and assign the eigenvalues and eigenvectors to the corresponding SpinSector instance
    for (auto& spin_sector : this->spinSectors) {      // use reference, otherwise we would make a copy...
        auto H = spin_sector.getHamiltonian();
        //std::clock_t start = std::clock();
        arma::vec eigenvalues;
        arma::mat eigenvectors;
        //std::clock_t end = std::clock();
        //std::cout<<std::endl<<(end-start)<<" time for solver";
        arma::eig_sym(eigenvalues, eigenvectors, H);
        //We store all eigen solutions
        spin_sector.setEigenvalues(eigenvalues);
        spin_sector.setEigenvectors(eigenvectors);
        //We extract only the groundstate
        for (int i = 0; i<eigenvalues.size(); i++) {
            groundStates(State {eigenvalues[i], eigenvectors.col(i), spin_sector.getS_z()});
        }
    }


}
void Hubbard::groundStates(State state) {
    if(areSame(state,this->groundstates.at(0))){
        groundstates.push_back(state);
    }
    else{
        if(compareState(state,groundstates.at(0))){
            groundstates = std::vector<State> {state};
        }
    }

}

arma::mat Hubbard::getSector(int index) {
    return spinSectors.at(index).getHamiltonian();
}

Lattice &Hubbard::getLattice() const {
    return lattice;
};


bool compareState(const State &o1, const State &o2){
    return o1.eigenValue < o2.eigenValue;
}
bool areSame(const State &o1, const State &o2) {
    double precision = 10000000;

    double ELIPSON = (fabs(o1.eigenValue) > fabs(o2.eigenValue)) ?  o2.eigenValue/precision : o1.eigenValue/precision   ;

    return fabs(o1.eigenValue - o2.eigenValue) < fabs(ELIPSON);
}