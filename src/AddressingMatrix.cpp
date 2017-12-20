//
// Created by Wulfix on 17/11/2017.
//

#include "AddressingMatrix.h"


/*
 * Idea from
 * Trygve Helgaker, Poul Jorgensen, Jeppe Olsen Molecular Electronic-Structure Theory . August 2000.
 *
 */
AddressingMatrix::AddressingMatrix(unsigned int sites, unsigned int electrons) : n_sites(sites), n_electrons(electrons){
    generateMatrix();


}

void AddressingMatrix::generateMatrix(){
    //create zero filled matrix
    addressMatrix = Eigen::MatrixXi::Zero(n_sites+1, n_electrons+1);

    //fill first x elements of the first column with 1 (x being sites-electrons+1) (the +1 is because we want a root starting with value1)
    // ex: 2 electrons and 5 sites.
    //[ 1 0 0 ]
    //[ 1 0 0 ]
    //[ 1 0 0 ]
    //[ 1 0 0 ]
    //[ 0 0 0 ]
    //[ 0 0 0 ]

    for(int i = 0;i<n_sites-n_electrons+1; i++){
        addressMatrix(i,0) = 1;
    }


    //every element is the sum of the value of the element vertically above and the value of the element left diagonally above.
    // ex: 2 electrons and 5 sites.
    //[ 1 0 0 ]
    //[ 1 1 0 ]
    //[ 1 2 1 ]
    //[ 1 3 3 ]
    //[ 0 4 6 ]
    //[ 0 0 10]
    for(int i = 0;i<n_sites; i++){
        for(int j = i-(n_sites-n_electrons); j<n_electrons;j++){
            if(j >= 0) {
                addressMatrix(i + 1,j + 1) = addressMatrix(i,j) + addressMatrix(i,j + 1);
            }
        }
    }
    //It is clear we could technically scrap the first row, however this the idea is better presented in the current form.
}

unsigned long AddressingMatrix::fetchAddress(boost::dynamic_bitset<> bitVector) {
    unsigned long address = 0;
    unsigned int electronIndex = 0;

    //Every index of the bitset is a move in the diagram, 0 is vertically downards, 1 is right-diagonally downwards
    //When we move diagonally we add the difference of our starting element and our destination to the index address.
    // ex: 2 electrons and 5 sites.
    // 01010 vector
    // 0+1+0+3+0 = 4
    // 00011 0
    // 00101 1
    // 00110 2
    // 01001 3
    // 01010 4 correct!

    for(unsigned int i = 1; i<n_sites+1;i++){
        electronIndex += (bitVector[0]);
        address += (addressMatrix(i,electronIndex)
                    - addressMatrix((i-1)*(bitVector[0]),((electronIndex-1)*(bitVector[0]))))*((bitVector[0]));
        bitVector >>= 1;

    }
    return address;
}

boost::dynamic_bitset<> AddressingMatrix::generateBinaryVector(unsigned long address){
    if(n_electrons == 0){
        return boost::dynamic_bitset<>(n_sites);
    }
    boost::dynamic_bitset<> biVector = boost::dynamic_bitset<>(n_sites);
    int eIndex = n_electrons;
    // Enter an address, it will then look for the biggest diagonal difference between 2 set of elements in the diagram
    // that is smaller or equal to the address, when it can move diagonally it adds 1 to bitvector at the final index.
    // It will then remove the weight(aka the diffrence between diagonals) from the address.
    // if it cannot move diagonally it will move vertically up in the diagram and add 0 instead and shift to next highest index.
    // ex: 2 electrons and 5 sites.
    // address 5
    // 10-4 = 6 : 0 : 5
    // 6-3 = 3 : 01 : 2
    // 3-1 = 2 : 011 : 0
    // no more weight left in address => 01100
    // 00011 0
    // 00101 1
    // 00110 2
    // 01001 3
    // 01010 4
    // 01100 5 correct!
    for (int sIndex = n_sites; sIndex>0; sIndex--){
        int weight = (addressMatrix(sIndex,eIndex) - addressMatrix(sIndex-1,eIndex-1));
        if (weight <= address){
            address -= weight;
            biVector[sIndex-1] = true;
            eIndex--;
            if (eIndex == 0){
                break;
            }
        }

    }

    return biVector;

}

void AddressingMatrix::print() {
    std::cout<<addressMatrix;

};
