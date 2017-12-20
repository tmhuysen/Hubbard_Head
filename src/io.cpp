//
// Created by Wulfix on 19/11/2017.
//

#include "io.h"


/** Read and return the hopping matrix from a file that contains its upper triangular representation (since it is symmetric).
 *
 * The first line of such a file should be the number of sites.
 * The second line then is the upper triangle of the hopping matrix, delimited by <CHAR>
 *
 * Example input file:
 *      3
 *      1.0<CHAR>-0.5<CHAR>-0.2<CHAR>1.5<CHAR>-0.7<CHAR>2.3
 */
arma::mat read_hopping_matrix_from_file(const std::string& filename, char delimiter_char) {


        // Make an uninitialized matrix to avoid going out-of-scope (M is only really initialized in the if)
    arma::mat M;
    std::ifstream is (filename);

    if (is.is_open()) {
            // The dimension of the resulting matrix is found as the first entry in the file name
        std::string dim_as_string;
        std::getline(is, dim_as_string);
        auto dim = static_cast<arma::uword>(std::stoi(dim_as_string));

            // Initialize a zero hopping matrix
        M = arma::zeros(dim, dim);

            // Read in the upper triangular part of the hopping matrix.
            // It is located on one line, so we can read in that line, and then break it down.
        std::string line;
        std::getline(is, line); // read in the actual line
        std::stringstream linestream (line); // convert the read line into a stringstream to perform std::getline() on it


            // The actual loop to read the upper triangular part of the hopping matrix
        for (arma::uword i = 0; i < dim; i++) {
            for (arma::uword j = i; j < dim; j++) {
                std::string value_as_string;
                std::getline(linestream, value_as_string, delimiter_char);
                M(i, j) = std::stod(value_as_string);
            }
        }
        is.close();
    }

        // After this loop, we only have the upper triangular form of the total hopping matrix
        // We can reflect the upper triangle to the lower triangle to obtain our result
    return arma::symmatu(M);;
}

void Symmatu(Eigen::MatrixXd &mat){
    for(int x =0; x<mat.innerSize();x++){
        for (int y = x+1; y<mat.innerSize();y++){
            mat(y,x) = mat(x,y);
        }
    }
}

void writeOutTriMat(std::ofstream &out, arma::mat mat){
    for(int i = 0; i<mat.n_cols;i++){
        for(int j = i; j<mat.n_cols;j++){
            if(i != mat.n_cols-1 || j != i){
                out<<mat(i,j)<<", ";
            }else{
                out<<mat(i,j)<<std::endl;
            }
        }
    }
}