#include <eigen3/Eigen/Dense>
#include <vector>
#include <fstream>
#include "jlpg.hpp"
using namespace Eigen;
using namespace std;
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
    return res*0.5;
}
Mat masked_Frob_grad(Mat X, Mat A){
    int Arow = A.rows(), Acol = A.cols();
    Mat U = Mat::Zero(Arow,Acol);
    for(int i = 0; i < Arow; i++){
        for(int j = 0; j< Acol; j++){
            if(A(i,j)>20) U(i,j) = X(i,j) - A(i,j);
        }
    }
    return U;
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
    Options opts(10000, 1e-1, 1e-1, 1e-0, 5e-1);
    opts.setClassical();
    ContOptions con_opts(opts, 10, 0.1, 100, 1e2, 1e2, 1e-1, 1e-1);

    Outputs out;
    
    Mat X = Mat::Zero(A.rows(),A.cols());
    X = cont_pgm(p, X, con_opts, out);
    // x = pgm(p, x, opts, out);
    //cout << x << endl;
    cout << "cputime=" << out.cputime << endl;
    cout << "f=" << out.F_cur << "; nrmG=" << out.nrmG << "; Flag=" << out.Flag << endl;
}