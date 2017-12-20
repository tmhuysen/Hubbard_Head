#define BOOST_TEST_MODULE "Hubbard"

#include "Hubbard.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain




BOOST_AUTO_TEST_CASE ( speed_test ) {
    std::clock_t start = std::clock();
    std::string path = "/Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/HubbardLibs/Hubbard_Head";
    std::string filename = path + "/tests/input_data/randomized_fourring_input.txt";
    std::string fileout = path +"/tests/test_output/randomized_fourring_output.txt";;
    std::ifstream is (filename);
    int counter = 0;
    std::ofstream outfile (fileout);
    outfile<<std::setprecision(12);
    arma::uword dim = 4;
    unsigned long elec = 4;
    unsigned long dimalt = 4;
    Lattice def = Lattice(dimalt);
    Hubbard instant = Hubbard(elec,def);
    char delimiter_char = ',';
    if (is.is_open()) {
        std::string line;
        while(std::getline(is, line)){
            try {
                counter++;

                arma::mat M = arma::zeros(dim, dim);

                // Read in the upper triangular part of the hopping matrix.
                // It is located on one line, so we can read in that line, and then break it down.
                std::stringstream linestream(
                        line); // convert the read line into a stringstream to perform std::getline() on it


                // The actual loop to read the upper triangular part of the hopping matrix
                for (arma::uword i = 0; i < dim; i++) {
                    for (arma::uword j = i; j < dim; j++) {
                        std::string value_as_string;
                        std::getline(linestream, value_as_string, delimiter_char);
                        M(i, j) = std::stod(value_as_string);
                    }
                }
                M = arma::symmatu(M);
                Lattice thisIn = Lattice(M);
                //std::clock_t start = std::clock();
                instant.setLattice(thisIn);
                //std::clock_t end = std::clock();
                //std::cout<<std::endl<<(end-start)<<" time for all";
                std::vector<State> ground = instant.getGroundstates();
                outfile << ground.at(0).eigenValue << std::endl;
            }catch(const std::invalid_argument& ia){
                std::cerr << "Invalid argument: " << ia.what() << '\n';
                std::cout<<std::endl<<counter<<" WTF";

            }








        }
        is.close();
        outfile.close();
    } else std::cout << "Unable to open file";


    std::clock_t end = std::clock();
    std::cout<<std::endl<<(end-start)<<" time for all";



}
