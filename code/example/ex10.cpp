#include <eigen3/Eigen/Dense>
#include <vector>
#include <fstream>
#include "jlpg.hpp"
using namespace Eigen;
using namespace std;

// load_csv function refers from https://stackoverflow.com/questions/34247057/how-to-read-csv-file-and-assign-to-eigen-matrix
template<typename M>
M load_csv (const std::string & path) {
    std::ifstream indata;
    indata.open(path);
    std::string line;
    std::vector<double> values;
    uint rows = 0;
    while (std::getline(indata, line)) {
        std::stringstream lineStream(line);
        std::string cell;
        while (std::getline(lineStream, cell, ',')) {
            values.push_back(std::stod(cell));
        }
        ++rows;
    }
    return Map<const Matrix<typename M::Scalar, M::RowsAtCompileTime, M::ColsAtCompileTime, RowMajor>>(values.data(), rows, values.size()/rows);
}


Real masked_Frob(Mat X, Mat A){
    int Arow = A.rows(), Acol = A.cols();
    Real res = 0;
    for(int i = 0; i < Arow; i++){
        for(int j = 0; j< Acol; j++){
            if(A(i,j)>20) res+=(X(i,j)-A(i,j))*(X(i,j)-A(i,j));
        }
    }
    return res*5e-3;
}
Mat masked_Frob_grad(Mat X, Mat A){
    int Arow = A.rows(), Acol = A.cols();
    Mat U = Mat::Zero(Arow,Acol);
    for(int i = 0; i < Arow; i++){
        for(int j = 0; j< Acol; j++){
            if(A(i,j)>20) U(i,j) = X(i,j) - A(i,j);
        }
    }
    return U*1e-2;
}
int main(){
    clock_t  Tstart, Tend;
	Tstart = clock();
    Mat A = load_csv<Mat>("../../jester-data-1.csv");
    Tend = clock();
	Real during = (double)(Tend - Tstart)/CLOCKS_PER_SEC;
    cout << "time elasped: " << during << endl;
    cout << "Recommedation of JESTER dataset" << endl;
    grad_pair<Mat> mc(bind(masked_Frob, placeholders::_1, A), bind(masked_Frob_grad, placeholders::_1, A));
    Problem<Mat> p(mc, NUCLEAR_NORM);
    p.mu = 10; 
    Options opts(10000, 0.01, 0.01, 1e-0, 5e-1);
    opts.setClassical();
    ContOptions con_opts(opts, 10, 0.1, 100, 1e2, 1e2, 1e-1, 1e-1);

    Outputs out;
    
    Mat x = Mat::Zero(A.rows(),A.cols());
    Mat x1 = x, x2 = x, x3 = x;
    printf("type & cont. & iters & cputime & fval & optimality & flag \\\\\n");
    
    // Armijo + continuation
    opts.setArmijo(0.1);
    x1 = cont_pgm(p, x1, con_opts, out);
    printf("Armijo & 1 & %d & %g & %g & %g & %d \\\\\n", out.iter, out.cputime, out.F_cur, out.nrmG, out.Flag);

    // Nonmonotone + continuation
    opts.setNonmonotone(0.1, 10);
    x2 = cont_pgm(p, x2, con_opts, out);
    printf("Nonmonotone & 1 & %d & %g & %g & %g & %d \\\\\n", out.iter, out.cputime, out.F_cur, out.nrmG, out.Flag);

    // Classical + continuation
    opts.setClassical();
    x3 = cont_pgm(p, x3, con_opts, out);
    printf("Classical & 1 & %d & %g & %g & %g & %d \\\\\n", out.iter, out.cputime, out.F_cur, out.nrmG, out.Flag);
}