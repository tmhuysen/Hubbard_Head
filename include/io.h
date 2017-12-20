//
// Created by Wulfix on 19/11/2017.
//

#ifndef HUBBARDCLEAN_IO_H
#define HUBBARDCLEAN_IO_H

#include <armadillo>
#include <Eigen/Dense>



/** Read and return the hopping matrix from a file that contains its upper triangular representation (since it is symmetric).
 *
 * The first line of such a file should be the number of sites.
 * The second line then is the upper triangle of the hopping matrix, delimited by <CHAR>
 *
 * Example input file:
 *      3
 *      1.0<CHAR>-0.5<CHAR>-0.2<CHAR>1.5<CHAR>-0.7<CHAR>2.3
 */
arma::mat read_hopping_matrix_from_file(const std::string& filename, char delimiter_char);


/** Read and process a file with (the upper triangular part of hopping matrices).
 *
 * The first line is the number of sites
 * Each further line is an upper triangle of the hopping matrix, delimited by <CHAR>
 *
 *
 * Creates a file named "<filename>_output.data",
 * with every line being the ground state energy of the Hubbard model with the specified hopping matrix on the corresponding line of <filename>
 *
 *
 * Example input file:
 *      4
 *      0.49504950495,1.82188218822,1.82188218822,1.82188218822,0.49504950495,1.82188218822,1.82188218822,0.49504950495,1.82188218822,0.49504950495
 *      0.49404940494,1.47344734473,1.47344734473,1.47344734473,0.49404940494,1.47344734473,1.47344734473,0.49404940494,1.47344734473,0.49404940494
 */
void process_file_with_hopping_matrices(const std::string& filename, char delimiter_char);

void writeOutTriMat(std::ofstream &out, arma::mat mat);

void Symmatu(Eigen::MatrixXd &mat);
#endif //HUBBARDCLEAN_IO_H
