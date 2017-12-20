#define BOOST_TEST_MODULE "Hubbard"

#include "Hubbard.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain


/** Read a matrix from a given file name, and put the values in a given matrix.
 */
void read_matrix_from_file(const std::string& filename, arma::mat& M) {
    arma::uword dim = M.n_rows;

    std::ifstream is (filename);
    if (is.is_open()) {
        for (arma::uword i = 0; i < dim; i++) {
            for (arma::uword j = 0; j < dim; j++) {
                is >> M(i,j);
            }
        }
        is.close();
    }
}


BOOST_AUTO_TEST_CASE ( reference_to_lattice ) {

    // We don't want any copies to be made when creating a HubbardClass instance.
    arma::umat A = {{0, 1, 1},
                    {1, 0, 1},
                    {1, 1, 0}};
    double t = 1.0;
    double U = 4.0;
    Lattice lattice(A, t, U);

    size_t N = 2;

    Hubbard h(N, lattice);    // h1: even number of electrons

    // We can test that there is no copy by requiring that the addresses of the lattices are equal.
    // CONFIRMED: This test fails if in HubbardClass.hpp
    //      const Lattice& lattice;     is replaced with        const Lattice lattice;
    BOOST_CHECK_EQUAL(&lattice,&h.getLattice());
}
