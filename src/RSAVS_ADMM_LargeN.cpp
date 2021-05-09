#include <RcppEigen.h>
#include <cmath>
#include <iostream>
#include <omp.h>
#define NUM_THREADS 4
// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp; 
using namespace Eigen; 


/*
    This is the core script of applying ADMM for subgroup analysis problem when n(number of observations) is relatively large.
    The s vector in the computation would be length n * (n - 1) / 2.
    And we need to iterate over lam1_vec and lam2_vec.
    So the whole size for s_mat would be 
        (n * (n - 1) / 2) \times (len(lam1) * len(lam2))
    In R, for n = 1000, len(lam1) = 40, len(lam2) = 50, the matrix takes 7.4GB of RAM.
    For n = 2000 and others keeping unchanged, the matrix would take 29.8 GB
    In this script, the BIC value is computed directly during the iteration, not after the whold iteration finished.
    Hence there is no need to keep all the s_vec into that s_mat.
    Only the last one(provided as initial value for current iteration), current one, and the (currently) best one are kept.
*/

// pre-define functions
int Compute_Index(const int i, const int j, const int i_length, const int j_length){
/*
  The matrix has i_length rows and j_length columns.
  Current index is (i, j), both starting from 0.
  Return the current index(also starting from 0) if we consider the matrix as column-major, i.e. vec(Matrix).
  Note: Actually in this function, j_length is unnecessary(The ranges of i and j are not checked.)
*/
    int res = j * i_length + i;
    return(res);
}

Eigen::SparseMatrix<double> Generate_D_Matrix(int n){
/*
    This function generate the D matrix(pairwise difference matrix) in sparse matrix form.
    Find details of D in the paper.
    Find details about sparse matrix in Eigen's manual.
    If you want to return the result back to R, you will probably need package 'Matrix'.
    Args: n: number of observations.
    Returns: res: a (n * (n - 1) / 2) \times n matirx in sparse form.
*/
    Eigen::SparseMatrix<double> res(n * (n - 1) / 2, n);
    std::vector<Eigen::Triplet<double> > d_triplet;
    int row_index = 0;
    int pci = 0;    // column index of the positive entries
    int nci = 0;    // column index of the negative entries
    for(pci = 0; pci < (n - 1); pci ++){
        for(nci = pci + 1; nci < n; nci++){
            d_triplet.push_back(Eigen::Triplet<double>(row_index, pci, 1.0));
            d_triplet.push_back(Eigen::Triplet<double>(row_index, nci, -1.0));
            row_index = row_index + 1;
        }
    }
    res.setFromTriplets(d_triplet.begin(), d_triplet.end());
    res.makeCompressed();
    return res;
}

double SoftThresholding(const double z, const double r) {
/*
    This function is used to perform softthresholding on z at r
    res = z - r    if z > fabs(r)
          z + r    if z < -fabs(r)
          0        o.w.
    Args: z: input value.
          r: thresholding point, use its absolute value.
    Returns: res: softhresholding of z at r
*/        
    double holding_point;
    holding_point = fabs(r);
    if(z > holding_point) {
        return(z - holding_point);
    } 
    if(z < -holding_point) {
            return(z + holding_point);
    }
    return(0.0);     
}

/*
    These are the loss functions of the algorithm
*/
Eigen::VectorXd Loss_L1(const Eigen::VectorXd &invec, const Eigen::VectorXd& loss_param){
    // int n = invec.size();
    Eigen::VectorXd res = invec.array().abs();
    return(res);
}

Eigen::VectorXd Loss_L2(const Eigen::VectorXd &invec, const Eigen::VectorXd& loss_param){
    Eigen::VectorXd res = invec.array().pow(2);
    return(res);
}

Eigen::VectorXd Loss_Huber(const Eigen::VectorXd& invec, const Eigen::VectorXd& loss_param){
    const double huber_c = loss_param[0];
    const int n = invec.size();
    Eigen::VectorXd res = Eigen::MatrixXd::Zero(n, 1);
    for(int i = 0; i < n; i++){
        if(fabs(invec[i]) <= huber_c){
            res[i] = 0.5 * invec[i] * invec[i];
        }else{
            res[i] = huber_c * fabs(invec[i]) - 0.5 * huber_c * huber_c;
        }
    }
    return(res);
}

/*
    These are penalty functions of the algorithm 
*/
Eigen::VectorXd Penalty_Lasso(const Eigen::VectorXd &invec, 
                              const Eigen::VectorXd &param, 
                              const bool &derivative){
    const int n = invec.size();
    const double lam = fabs(param[0]);
    Eigen::VectorXd res = Eigen::MatrixXd::Zero(n, 1);
    if(derivative){
        for(int i = 0; i < n; i++){
            if(invec[i] >= 0){
                res[i] = lam;
            }else{
                res[i] = -lam;
            }
        }
    }else{
        for(int i = 0; i < n; i++){
            res[i] = lam * fabs(invec[i]);
        }
    }
    return(res);
}

Eigen::VectorXd Penalty_SCAD(const Eigen::VectorXd &invec, 
                             const Eigen::VectorXd &param, 
                             const bool &derivative){
    const int n = invec.size();
    const double lam = fabs(param[0]);
    const double gam = fabs(param[1]);
    Eigen::VectorXd res = Eigen::MatrixXd::Zero(n, 1);
    double tmp;
    if(derivative){
        for(int i = 0; i < n; i++){
            tmp= invec[i];
            if(tmp > lam * gam){
                res[i] = 0;
            }else if(tmp >= lam){
                res[i] = (lam * gam - tmp) / (gam - 1);
            }else if(tmp >= 0){
                res[i] = lam;
            }else if(tmp > -lam){
                res[i] = -lam;
            }else if(tmp > -lam * gam){
                res[i] = (-gam * lam - tmp) / (gam - 1);
            }else{
                res[i] = 0;
            }
        }
    }else{
        for(int i = 0; i < n; i++){
            tmp = invec[i];
            if(tmp> lam * gam){
                res[i] = (gam + 1) * lam * lam / 2;
            }else if(tmp >= lam){
                res[i] = (-tmp * tmp / 2 + gam * lam * tmp) / (gam - 1) - lam * lam / 2 / (gam - 1);
            }else if(tmp > -lam){
                res[i] = lam * fabs(tmp);
            }else if(tmp > -lam * gam){
                res[i] = (-tmp * tmp / 2 - gam * lam * tmp) / (gam - 1) - lam * lam / 2 / (gam - 1);
            }else{
                res[i] = (gam + 1) * lam * lam / 2;
            }
        }
    }
    return(res);
}

Eigen::VectorXd Penalty_MCP(const Eigen::VectorXd &invec, 
                            const Eigen::VectorXd &param, 
                            const bool &derivative){
    const int n = invec.size();
    const double lam = fabs(param[0]);
    const double gam = fabs(param[1]);
    Eigen::VectorXd res = Eigen::MatrixXd::Zero(n, 1);
    double tmp;
    if(derivative){
        for(int i = 0; i < n; i++){
            tmp= invec[i];
            if(tmp > gam * lam){
                res[i] = 0;
            }else if(tmp >= 0){
                res[i] = lam - tmp / gam;
            }else if(tmp > -lam * gam){
                res[i] = -lam - tmp / gam;
            }else{
                res[i] = 0;
            }
        }
    }else{
        for(int i = 0; i < n; i++){
            tmp = invec[i];
            if(tmp > gam * lam){
                res[i] = gam * lam * lam / 2;
            }else if(tmp >= 0){
                res[i] = lam * tmp - tmp * tmp / 2 / gam;
            }else if(tmp > -gam * lam){
                res[i] = -lam * tmp - tmp * tmp / 2 / gam;   
            }else{
                res[i] = gam * lam * lam / 2;
            }
        }
    }
    return(res);
}





/*  
    These are the updates of the algorithm of different types.
    The detail of types still needs modification.
    The update function for s vector and w vector are seperated because we use OpenMP during the update of s vector.
*/

void UpdateZ_L1(Eigen::VectorXd &z_invec, const Eigen::VectorXd &loss_param, const double &r){
/*
    Update the z vector for L1 loss
    Args: z_invec: the inputting z vector. 
                   During the algorithm its value is actually y - \mu - x^T\beta + q1 / r1
          loss_param: The parameter for the loss function. 
                      For L1, this is unused.
          r: quadratic parameter in the algorithm
    Returns: no return, but the z_invec is updated
*/
    int n = z_invec.size();
    for(int i = 0; i < n; i++){
        z_invec[i] = SoftThresholding(z_invec[i], 1.0 / n / r);
        // z_invec[i] = SoftThresholding(z_invec[i], 1.0 / 2.0 / r);    // use \sum\loss / 2, not \sum\loss / n
    }
}

void UpdateZ_L2(Eigen::VectorXd &z_invec, const Eigen::VectorXd &loss_param, const double &r){
/*
    Update the z vector for L2 loss.
    For L2 loss, another possibile solution is to omit z and directly update \mu and \beta.
    Args: z_invec: the inputting z vector. 
                   During the algorithm its value is actually y - \mu - x^T\beta + q1 / r1
          loss_param: The parameter for the loss function. 
                      For L2, this is unused.
          r: quadratic parameter in the algorithm
    Returns: no return, but the z_invec is updated
*/
    int n = z_invec.size();
    for(int i = 0; i < n; i++){
        z_invec[i] = z_invec[i] / (1.0 + 2.0 / n / r);
        // z_invec[i] = z_invec[i] / (1.0 + 2.0 / 2.0 / r);    // use \sum\loss / 2, not \sum\loss / n
    }
}

void UpdateZ_Huber(Eigen::VectorXd &z_invec, const Eigen::VectorXd &loss_param, const double &r){
/*
    Update the z vector for Huber loss
    Args: z_invec: the inputting z vector. 
                   During the algorithm its value is actually y - \mu - x^T\beta + q1 / r1
          loss_param: The parameter for the loss function. 
                      For Huber loss, this is the c.
          r: quadratic parameter in the algorithm
    Returns: no return, but the z_invec is updated
*/
    double c = loss_param[0];
    int n = z_invec.size();
    double nr = c / n / r;    // \sum\loss / n
    // nr = c / 2.0 / r;    // use \sum\loss / 2, not \sum\loss / n
    double z;
    for(int i = 0; i < n; i++){
        z = z_invec[i];
        if(z > (nr + c)){
            z_invec[i] = z - nr;
        }else{
            if(z < (-nr - c)){
                z_invec[i] = z + nr;
            }else{
                z_invec[i] = z / (1.0 + nr / c);
            }
        }
    }    
}

void UpdateS_Lasso(Eigen::VectorXd &s_invec, const Eigen::VectorXd &penalty_param, const double &r){
/*
    Update the s vector for the Lasso penalty
    Args:s_invec: the inputting s vector.
                  During the algorithm its value is d_mat * \mu + q2 / r2
         penalty_param: the parameters for the penalty
                        For Lasso, lambda = penlaty_param[0].
         r: the quadratic term in the algorithm
    Return: no return values, but the s_invec is updated
*/    
    double lambda = penalty_param[0];
    
    omp_set_num_threads(NUM_THREADS);
    #pragma omp parallel for
    for(int i = 0; i < s_invec.size(); i++){
        s_invec[i] = SoftThresholding(s_invec[i], lambda / r);
    }
}


void UpdateS_SCAD(Eigen::VectorXd &s_invec, const Eigen::VectorXd &penalty_param, const double &r){
/*
    Update the s vector for the SCAD penalty
    Args:s_invec: the inputting s vector.
                  During the algorithm its value is d_mat * \mu + q2 / r2
         penalty_param: the parameters for the penalty
                        For SCAD, lambda = penlaty_param[0] and gamma = penalty_param[1].
         r: the quadratic term in the algorithm
    Return: no return values, but the s_invec is updated
*/
    double lambda = penalty_param[0];
    double gamma = penalty_param[1];
    // double s_abs;
    
    omp_set_num_threads(NUM_THREADS);
    #pragma omp parallel for
    for(int i = 0; i < s_invec.size(); i++){
        // s_abs = fabs(s_invec[i]);
        if(fabs(s_invec[i]) <= ((1.0 + 1.0 / r) * lambda)){
            s_invec[i] = SoftThresholding(s_invec[i], lambda / r);
        } else{
            if(fabs(s_invec[i]) <= (lambda * gamma)){
                s_invec[i] = SoftThresholding(s_invec[i], gamma * lambda / r / (gamma - 1)) / (1.0 - 1.0 / r / (gamma - 1));
            }
        }
    }
}

void UpdateS_MCP(Eigen::VectorXd &s_invec, const Eigen::VectorXd &penalty_param, const double &r){
/*
    Update the s vector for the MCP penalty
    Args:s_invec: the inputting s vector.
                  During the algorithm its value is d_mat * \mu + q2 / r2
         penalty_param: the parameters for the penalty
                        For MCP, lambda = penlaty_param[0] and gamma = penalty_param[1].
         r: the quadratic term in the algorithm
    Return: no return values, but the s_invec is updated
*/    
    double lambda = penalty_param[0];
    double gamma = penalty_param[1];
    
    omp_set_num_threads(NUM_THREADS);
    #pragma omp parallel for
    for(int i = 0; i < s_invec.size(); i++){
        if(fabs(s_invec[i]) <= (lambda * gamma)){
            s_invec[i] = SoftThresholding(s_invec[i], lambda / r) / (1.0 - 1.0 / r / gamma);
        }
    }
}

void UpdateW_Lasso(Eigen::VectorXd &w_invec, const Eigen::VectorXd &penalty_param, const double &r){
/*
    Update the w vector for the Lasso penalty
    Args:w_invec: the inputting w vector.
                  During the algorithm its value is \beta + q3 / r3
         penalty_param: the parameters for the penalty
                        For Lasso, lambda = penlaty_param[0].
         r: the quadratic term in the algorithm
    Return: no return values, but the w_invec is updated
*/   
    double lambda = penalty_param[0];
    
    for(int i = 0; i < w_invec.size(); i++){
        w_invec[i] = SoftThresholding(w_invec[i], lambda / r);
    }  
}

void UpdateW_SCAD(Eigen::VectorXd &w_invec, const Eigen::VectorXd &penalty_param, const double &r){
/*
    Update the w vector for the SCAD penalty
    Args:w_invec: the inputting w vector.
                  During the algorithm its value is \beta + q3 / r3
         penalty_param: the parameters for the penalty
                        For SCAD, lambda = penlaty_param[0], gamma = penalty_param[1].
         r: the quadratic term in the algorithm
    Return: no return values, but the w_invec is updated
*/       
    double lambda = penalty_param[0];
    double gamma = penalty_param[1];
    double w_abs;
    
    for(int i = 0; i < w_invec.size(); i++){
        w_abs = fabs(w_invec[i]);
        if(w_abs <= ((1.0 + 1.0 / r) * lambda)){
            w_invec[i] = SoftThresholding(w_invec[i], lambda / r);
        } else{
            if(w_abs <= (lambda * gamma)){
                w_invec[i] = SoftThresholding(w_invec[i], gamma * lambda / r / (gamma - 1)) / (1.0 - 1.0 / r / (gamma - 1));
            }
        }
    }
}

void UpdateW_MCP(Eigen::VectorXd &w_invec, const Eigen::VectorXd &penalty_param, const double &r){
/*
    Update the w vector for the MCP penalty
    Args:w_invec: the inputting w vector.
                  During the algorithm its value is \beta + q3 / r3
         penalty_param: the parameters for the penalty
                        For MCP, lambda = penlaty_param[0] and gamma = penalty_param[1].
         r: the quadratic term in the algorithm
    Return: no return values, but the s_invec is updated
*/    
    double lambda = penalty_param[0];
    double gamma = penalty_param[1];

    for(int i = 0; i < w_invec.size(); i++){
        if(fabs(w_invec[i]) <= (lambda * gamma)){
            w_invec[i] = SoftThresholding(w_invec[i], lambda / r) / (1.0 - 1.0 / r / gamma);
        }
    }    
}



Eigen::VectorXd RSAVS_Compute_Loss_Value_Cpp(const Eigen::VectorXd& y_vec, const Eigen::MatrixXd& x_mat, const int& n, const int& p, 
                                             const std::string& l_type, const Eigen::VectorXd& l_param,
                                             const char& p1_type, const Eigen::VectorXd p1_param, 
                                             const char& p2_type, const Eigen::VectorXd p2_param, 
                                             const Eigen::VectorXd& const_r123, const Eigen::VectorXd& const_abc, 
                                             const Eigen::VectorXd& mu_vec, const Eigen::VectorXd& beta_vec, 
                                             const Eigen::VectorXd& z_vec, const Eigen::VectorXd& s_vec, const Eigen::VectorXd& w_vec, 
                                             const Eigen::VectorXd& q1_vec, const Eigen::VectorXd& q2_vec, const Eigen::VectorXd& q3_vec){
    
    Eigen::VectorXd loss_vec = Eigen::MatrixXd::Zero(7, 1);
    double loss, loss_p1, loss_p2, loss_p3, loss_aug1, loss_aug2, loss_aug3;
    Eigen::VectorXd tmp;
    Eigen::SparseMatrix<double> d_mat = Generate_D_Matrix(n);
    
    // ------ setup loss function ------
    Eigen::VectorXd (*Loss_Function)(const Eigen::VectorXd &, const Eigen::VectorXd &);
    Loss_Function = Loss_L1;
    if(l_type == "L2"){
        Loss_Function = Loss_L2;
    }
    if(l_type == "Huber"){
        Loss_Function = Loss_Huber;
    }
    
    // ------ setup p1 penalty function ------
    Eigen::VectorXd (*P1_Function)(const Eigen::VectorXd &, const Eigen::VectorXd &, const bool &);
    P1_Function = Penalty_Lasso;
    if(p1_type == "S"){
        P1_Function = Penalty_SCAD;
    }
    if(p1_type == "M"){
        P1_Function = Penalty_MCP;
    }
    
    // ------ setup p2 penalty function ------
    Eigen::VectorXd (*P2_Function)(const Eigen::VectorXd &, const Eigen::VectorXd &, const bool &);
    P2_Function = Penalty_Lasso;
    if(p2_type == "S"){
        P2_Function = Penalty_SCAD;
    }
    if(p2_type == "M"){
        P2_Function = Penalty_MCP;
    }
    
    // ------ compute loss value ------
    // --- main loss function ---
    loss_p1 = 1.0 / const_abc[0] * (Loss_Function(z_vec, loss_param).sum());
    loss_p2 = const_abc[1] * (P1_Function(s_vec, p1_param, false).sum());
    loss_p3 = const_abc[2] * (P2_Function(w_vec, p2_param, false).sum());
    
    // --- augmented lagrangian part ---
    tmp = y_vec - mu_vec - x_mat * beta_vec - z_vec;
    loss_aug1 = const_r123[0] / 2.0 * tmp.dot(tmp) + tmp.dot(q1_vec);
    
    tmp = d_mat * mu_vec - s_vec;
    loss_aug2 = const_r123[1] / 2.0 * tmp.dot(tmp) + tmp.dot(q2_vec);
    
    tmp = beta_vec - w_vec;
    loss_aug3 = const_r123[2] / 2.0 * tmp.dot(tmp) + tmp.dot(q3_vec);
    
    loss = loss_p1 + loss_p2 + loss_p3 + loss_aug1 + loss_aug2 + loss_aug3;
    
    loss_vec[0] = loss;
    loss_vec[1] = loss_p1;
    loss_vec[2] = loss_p2;
    loss_vec[3] = loss_p3;
    loss_vec[4] = loss_aug1;
    loss_vec[5] = loss_aug2;
    loss_vec[6] = loss_aug3;
    
    return(loss_vec);
}

// main body 
// [[Rcpp::export]]
Rcpp::List RSAVS_LargeN_Rcpp(const Eigen::MatrixXd x_mat, const Eigen::VectorXd y_vec, const int n, const int p, 
                             const std::string& loss_type, const Eigen::VectorXd loss_param, 
                             const char p1_type, Eigen::VectorXd p1_param, 
                             const char p2_type, Eigen::VectorXd p2_param, 
                             const Eigen::VectorXd lam1_vec, const Eigen::VectorXd lam2_vec,
                             const double r1, const double r2, const double r3, 
                             const double phi, const double tol, const double max_iter){
/*
    This is the main ADMM algorithm for solving RSAVA when n is large.
    In this algorithm, we have lam1(for penalty 1, hence \mu) and lam2(for penalty 2, hence \beta) to iterate.
    Outside the function, it's designed as a matrix:
            +--- lam_2 ---+
            ||    -->     |
            ||    -->     |
            |+----------> |
          lam_1           |
            |             |
            |             |
            |             |
            +-------------+
    But in this function, it's converted to a vector in a column-major pattern.
    Refer to the function Compute_Index for the details.
    Args:
    Returns:
*/          

// prepare some variables
    int i, j, j_start, i_old, j_old, ind, ind_old, k, best_ind, best_i, best_j;
    /*
        i, j: iterator over lam1_vec and lam2_vec
        j_start: starting index of j. When i = 0, j starts from 1, otherwise j starts from 0
        i_old, j_old: the index of initial values for current (i, j)
        ind, ind_old: the index in the result matrix corresponding to coordinate (i, j), (i_old, j_old) of the lambda matrix.
        k: inner loop counter for the algorithm
        best_ind: the row index of the resulting matrix for the best solution
        best_i, best_j: the row and column index of the lambda matrix for the best solution
    */
    double lam1, lam2;
    double diff = 1.0;
    double primal_residual, dual_residual, res_part1, res_part2, res_part3;
    const int lam1_len = lam1_vec.size();
    const int lam2_len = lam2_vec.size();
    Eigen::MatrixXd beta_mat = Eigen::MatrixXd::Zero(lam1_len * lam2_len, p);
    Rcpp::Function r_median("median");    // Access the R function median
    double y_median = Rcpp::as<double>(r_median(y_vec));    // For the initial value of mu vector
    Eigen::MatrixXd mu_mat = Eigen::MatrixXd::Constant(lam1_len * lam2_len, n, y_median);
    Eigen::MatrixXd z_mat = Eigen::MatrixXd::Zero(lam1_len * lam2_len, n);
    // Eigen::MatrixXd s_mat = Eigen::MatrixXd::Zero(lam1_len * lam2_len, n * (n - 1) / 2);    // danger for large n
    Eigen::MatrixXd w_mat = Eigen::MatrixXd::Zero(lam1_len * lam2_len, p);
    Eigen::MatrixXd q1_mat = Eigen::MatrixXd::Zero(lam1_len * lam2_len, n);
    // Eigen::MatrixXd q2_mat = Eigen::MatrixXd::Zero(lam1_len * lam2_len, n * (n - 1) / 2);    // danger for large n
    Eigen::MatrixXd q3_mat = Eigen::MatrixXd::Zero(lam1_len * lam2_len, p);
    Eigen::MatrixXd bic_mat = Eigen::MatrixXd::Zero(lam1_len, lam2_len);
    Eigen::MatrixXd k_mat = Eigen::MatrixXd::Zero(lam1_len, lam2_len);    // this matrix can tell us whether or not the algorithm converged at that lam1 and lam2
    void (*UpdateZ)(Eigen::VectorXd &, const Eigen::VectorXd &, const double &);
    void (*UpdateS)(Eigen::VectorXd &, const Eigen::VectorXd &, const double &);
    void (*UpdateW)(Eigen::VectorXd &, const Eigen::VectorXd &, const double &);
    
    Rcpp::Rcout << "Basic variables initialized." << std::endl;
    // variables needed during algorithm
    Eigen::VectorXd beta_vec = Eigen::MatrixXd::Zero(p, 1); 
    Eigen::VectorXd beta_old = Eigen::MatrixXd::Zero(p, 1);
    Rcpp::Rcout << "Beta finished." << std::endl;
    Eigen::VectorXd mu_vec = Eigen::MatrixXd::Constant(n, 1, y_median); 
    Eigen::VectorXd mu_old = Eigen::MatrixXd::Constant(n, 1, y_median);
    Rcpp::Rcout << "Mu finished." << std::endl;
    Eigen::VectorXd z_vec = Eigen::MatrixXd::Zero(n, 1);
    Eigen::VectorXd z_old = Eigen::MatrixXd::Zero(n, 1);
    Rcpp::Rcout << "Z finished." << std::endl;
    Eigen::VectorXd s_vec = Eigen::MatrixXd::Zero(n * (n - 1) / 2, 1);
    Eigen::VectorXd s_old = Eigen::MatrixXd::Zero(n * (n - 1) / 2, 1);
    Eigen::VectorXd s_anchor = Eigen::MatrixXd::Zero(n * (n - 1) / 2, 1);
    Eigen::VectorXd s_best = Eigen::MatrixXd::Zero(n * (n - 1) / 2, 1);
    Rcpp::Rcout << "S finished." << std::endl;
    Eigen::VectorXd w_vec = Eigen::MatrixXd::Zero(p, 1);
    Eigen::VectorXd w_old = Eigen::MatrixXd::Zero(p, 1);
    Rcpp::Rcout << "W finished." << std::endl;
    Eigen::VectorXd q1_vec = Eigen::MatrixXd::Zero(n, 1);
    Eigen::VectorXd q1_old = Eigen::MatrixXd::Zero(n, 1);
    Rcpp::Rcout << "q1 finished." << std::endl;
    Eigen::VectorXd q2_vec = Eigen::MatrixXd::Zero(n * (n - 1) / 2, 1);
    Eigen::VectorXd q2_old = Eigen::MatrixXd::Zero(n * (n - 1) / 2, 1); 
    Eigen::VectorXd q2_anchor = Eigen::MatrixXd::Zero(n * (n - 1) / 2, 1);
    Rcpp::Rcout << "q2 finished." << std::endl;
    Eigen::VectorXd q3_vec = Eigen::MatrixXd::Zero(p, 1);
    Eigen::VectorXd q3_old = Eigen::MatrixXd::Zero(p, 1);
    Rcpp::Rcout << "q3 finished." << std::endl;

    // variable and function for summarizing the result
    Rcpp::Function r_summary_iteration("RSAVS_Summary_Iteration");
    Rcpp::List current_iteration_info;
    Eigen::MatrixXd mu_improve_mat = Eigen::MatrixXd::Zero(lam1_len * lam2_len, n);    // the improved mu vector
    Eigen::MatrixXd group_num_mat = Eigen::MatrixXd::Constant(lam1_len, lam2_len, 1);    // the group number
    Eigen::MatrixXd active_num_mat = Eigen::MatrixXd::Constant(lam1_len, lam2_len, 0);    // the number of active covariates 

// prepare some values that will be needed in the algorithm
    // Generate the D matrix in sparse form
    Rcpp::Rcout << "D started." << std::endl;
    Eigen::SparseMatrix<double> d_mat = Generate_D_Matrix(n);
    Rcpp::Rcout << "D finished." << std::endl;
    // Generate the lefet matirx for updating beta. In this version, we pre-assume p < n
    Eigen::VectorXd beta_right = Eigen::MatrixXd::Zero(p, 1);
    Eigen::MatrixXd beta_left = r1 * x_mat.transpose() * x_mat + r3 * Eigen::MatrixXd::Identity(p, p);  
    Eigen::MatrixXd beta_tmp;    
    LDLT<Eigen::MatrixXd> beta_left_solver(p);
    LDLT<Eigen::MatrixXd> beta_left_solver_large_p(n);
    if(p > n){    // overwrite beta_left, and set the solver beta_left_solver_large_p if p > n
        beta_left = r1 * x_mat * x_mat.transpose() + r3 * Eigen::MatrixXd::Identity(n, n);
        beta_left_solver_large_p.compute(beta_left);
        beta_tmp = beta_left_solver_large_p.solve(x_mat);
        // beta_left_solver.compute(beta_left);    //    test when debugging large p
    }else{    // use the default solver beta_left_solver
        beta_left_solver.compute(beta_left);
    } 
    Rcpp::Rcout << "beta left finished." << std::endl;
    // Generate the left matrix for updating mu
    Eigen::VectorXd mu_right = Eigen::MatrixXd::Zero(n, 1);
    // Eigen::MatrixXd mu_left = r2 * d_mat.transpose() * d_mat + r1 * Eigen::MatrixXd::Identity(n, n);
    // use a fashion recommended by Eigen
    Eigen::MatrixXd mu_left = r1 * Eigen::MatrixXd::Identity(n, n);
    mu_left += r2 * d_mat.transpose() * d_mat;    
    LDLT<Eigen::MatrixXd> mu_left_solver(n);
    mu_left_solver.compute(mu_left);
    Rcpp::Rcout << "mu left finished." << std::endl;
// assign the correct updating functions    
    UpdateZ = UpdateZ_L1;    // default loss type is L1
    if(loss_type == "Huber"){
        UpdateZ = UpdateZ_Huber;
    }
    if(loss_type == "L2"){
        UpdateZ = UpdateZ_L2;
    }
    
    UpdateS = UpdateS_SCAD;    // default penalty type for s is SCAD
    if(p1_type == 'L'){
        UpdateS = UpdateS_Lasso;
    } else{
        if(p1_type == 'M'){
            UpdateS = UpdateS_MCP;
        }
    }
    
    UpdateW = UpdateW_SCAD;    // default penalty type for w is SCAD
    if(p2_type == 'L'){
        UpdateW = UpdateW_Lasso;
    } else{
        if(p2_type == 'M'){
            UpdateW = UpdateW_MCP;
        }
    }
    
// main algorithm
/*
    The algorithm iterates over the lambda matrix.
    The outer iteration is over lambda 1, the inner iteration is over lambda 2.
*/
    best_i = 0;
    best_j = 0;
    best_ind = 0;
    current_iteration_info = Rcpp::as<Rcpp::List>(r_summary_iteration(Rcpp::Named("y_vec", y_vec), 
                                                                      Rcpp::Named("x_mat", x_mat), 
                                                                      Rcpp::Named("mu_vec", mu_vec), 
                                                                      Rcpp::Named("beta_vec", beta_vec), 
                                                                      Rcpp::Named("s_vec", s_vec), 
                                                                      Rcpp::Named("w_vec", w_vec), 
                                                                      Rcpp::Named("loss_type", loss_type), 
                                                                      Rcpp::Named("loss_param", loss_param), 
                                                                      Rcpp::Named("phi", phi)));
    bic_mat(best_i, best_j) = Rcpp::as<double>(current_iteration_info["bic"]);
    mu_improve_mat.row(best_ind) = Rcpp::as<Eigen::VectorXd>(current_iteration_info["mu_vec"]);
    
    for(i = 0; i < lam1_len; i++){
        // find the appropriate starting point of j
        if(i == 0){
            j_start = 1;
        }else{
            j_start = 0;
        }
        
        for(j = j_start; j < lam2_len; j++){
            // assign lambda 1 and lambda 2
            lam1 = lam1_vec[i];
            lam2 = lam2_vec[j];
            p1_param[0] = lam1;
            p2_param[0] = lam2;
            
            // find (i_old, j_old), coordinate of initial values
            i_old = i - (j == 0);
            j_old = std::max(0, j - 1);
            
            // compute the index
            ind = Compute_Index(i, j, lam1_len, lam2_len);
            ind_old = Compute_Index(i_old, j_old, lam1_len, lam2_len);
            
            // prepare initial values
            beta_old = beta_mat.row(ind_old);
            mu_old = mu_mat.row(ind_old);
            z_old = z_mat.row(ind_old);
            w_old = w_mat.row(ind_old);
            q1_old = q1_mat.row(ind_old);
            q3_old = q3_mat.row(ind_old);
            if(j == 0){
                q2_old = q2_anchor;
                s_old = s_anchor;
            }
            // s_old = s_mat.row(ind_old);    // danger for large n
            // q2_old = q2_mat.row(ind_old);    // danger for large n
            
            // Rcpp::Rcout << "i = " << i << ", j = " << j << ", ind = " << ind << ". lam1 = " << lam1 << ", lam2 = " << lam2 << std::endl;
            // Rcpp::Rcout << "i_old = " << i_old << ", j_old = " << j_old << ", ind_old = " << ind_old << std::endl;
            // Rcpp::Rcout << "lam1 = " << lam1 << ", lam2 = " << lam2 << std::endl;
            // update steps
            k = 0;
            diff = 1.0;
            while((k < max_iter) && (diff > tol)){
                // update beta
                beta_right = r1 * x_mat.transpose() * (y_vec - mu_old - z_old) + r3 * w_old + x_mat.transpose() * q1_old - q3_old;
                if(p <= n){    // use the default solver
                    beta_vec = beta_left_solver.solve(beta_right);
                }else{    // use the large p solver
                    // beta_tmp = beta_left_solver.solve(x_mat);    // its compuation is done before
                    // beta_vec = beta_left_solver.solve(beta_right);    // test when debuggind large p
                    // beta_vec = 1.0 / r3 * (Eigen::MatrixXd::Identity(p, p) - r1 * x_mat.transpose() * beta_tmp) * beta_right;    // this one maybe slow
                    beta_vec = 1.0 / r3 * (beta_right - r1 * x_mat.transpose() * beta_tmp * beta_right);    // this one is a little faster.
                }
                
                // Rcpp::Rcout << "beta updated." << std::endl;
                // update mu
                mu_right = r1 * (y_vec - x_mat * beta_vec - z_old) + q1_old;
                mu_right += d_mat.transpose() * (r2 * s_old - q2_old);
                mu_vec = mu_left_solver.solve(mu_right);
                // Rcpp::Rcout << "mu updated." << std::endl;
                // update z
                z_vec = y_vec - mu_vec - x_mat * beta_vec + q1_old / r1;
                UpdateZ(z_vec, loss_param, r1);
                // Rcpp::Rcout << "z updated." << std::endl;
                // update s
                s_vec = d_mat * mu_vec + q2_old / r2;
                UpdateS(s_vec, p1_param, r2);
                // Rcpp::Rcout << "s updated." << std::endl;
                // update w
                w_vec = beta_vec + q3_old / r3;
                UpdateW(w_vec, p2_param, r3);
                // Rcpp::Rcout << "w updated." << std::endl;
                // update q1, q2 and q3
                q1_vec = q1_old + r1 * (y_vec - mu_vec - x_mat * beta_vec - z_vec);
                // Rcpp::Rcout << "q1 updated." << std::endl;
                q2_vec = q2_old + r2 * (d_mat * mu_vec - s_vec);
                // Rcpp::Rcout << "q2 updated." << std::endl;
                q3_vec = q3_old + r3 * (beta_vec - w_vec);
                // Rcpp::Rcout << "q3 updated." << std::endl;
                // compute residual and check for convergence
                k = k + 1;    // update the inner loop counter
                // primal residual
                res_part1 = (y_vec - mu_vec - x_mat * beta_vec - z_vec).squaredNorm();
                res_part2 = (d_mat * mu_vec - s_vec).squaredNorm();
                res_part3 = (beta_vec - w_vec).squaredNorm();
                primal_residual = res_part1 + res_part2 + res_part3;
                // Rcpp::Rcout << "primal updated." << std::endl;
                // dual residual
                res_part1 = (r1 * (z_vec - z_old) - r2 * d_mat.transpose() * (s_vec - s_old)).squaredNorm();
                res_part2 = (r1 * x_mat.transpose() * (z_vec - z_old) - r3 * (w_vec - w_old)).squaredNorm();
                dual_residual = res_part1 + res_part2;
                // Rcpp::Rcout << "dual updated." << std::endl;
                diff = std::max(primal_residual, dual_residual);
                
                // prepare for next inner iteration
                beta_old = beta_vec;
                mu_old = mu_vec;
                z_old = z_vec;
                s_old = s_vec;
                w_old = w_vec;
                q1_old = q1_vec;
                q2_old = q2_vec;
                q3_old = q3_vec;
            }
            
            // save the result
            beta_mat.row(ind) = beta_vec;
            mu_mat.row(ind) = mu_vec;
            z_mat.row(ind) = z_vec;
            w_mat.row(ind) = w_vec;
            q1_mat.row(ind) = q1_vec;
            q3_mat.row(ind) = q3_vec;
            s_old = s_vec;
            q2_old = q2_vec;
            if(j == 0){
                s_anchor = s_vec;
                q2_anchor = q2_vec;
            }
            // s_mat.row(ind) = s_vec;    // danger for large n
            // q2_mat.row(ind) = q2_vec;    // danger for large n
            
            k_mat(i, j) = k;
            // summary this result
            current_iteration_info = Rcpp::as<Rcpp::List>(r_summary_iteration(Rcpp::Named("y_vec", y_vec), 
                                                                              Rcpp::Named("x_mat", x_mat), 
                                                                              Rcpp::Named("mu_vec", mu_vec), 
                                                                              Rcpp::Named("beta_vec", beta_vec), 
                                                                              Rcpp::Named("s_vec", s_vec), 
                                                                              Rcpp::Named("w_vec", w_vec), 
                                                                              Rcpp::Named("loss_type", loss_type), 
                                                                              Rcpp::Named("loss_param", loss_param), 
                                                                              Rcpp::Named("phi", phi)));
            bic_mat(i, j) = Rcpp::as<double>(current_iteration_info["bic"]);
            mu_improve_mat.row(ind) = Rcpp::as<Eigen::VectorXd>(current_iteration_info["mu_vec"]);
            group_num_mat(i, j) = Rcpp::as<double>(current_iteration_info["group_num"]);
            active_num_mat(i, j) = Rcpp::as<double>(current_iteration_info["active_num"]);
            
            // check whether we should update s_best
            if(bic_mat(i, j) < bic_mat(best_i, best_j)){
                best_i = i; 
                best_j = j;
                best_ind = ind;
                s_best = s_vec;
            }
        }
    }

// create and return the results
    Rcpp::List res = Rcpp::List::create(
        Rcpp::Named("beta_mat", beta_mat), 
        Rcpp::Named("mu_mat", mu_mat), 
        Rcpp::Named("mu_improve_mat", mu_improve_mat),
        Rcpp::Named("z_mat", z_mat), 
        // Rcpp::Named("s_mat", s_mat),    // danger for large n
        Rcpp::Named("w_mat", w_mat),
        Rcpp::Named("q1_mat", q1_mat),
        Rcpp::Named("q3_mat", q3_mat),
        Rcpp::Named("bic_mat", bic_mat),
        Rcpp::Named("k_mat", k_mat),
        Rcpp::Named("group_num_mat", group_num_mat),
        Rcpp::Named("active_num_mat", active_num_mat), 
        Rcpp::Named("lam1_vec", lam1_vec),
        Rcpp::Named("lam2_vec", lam2_vec),
        Rcpp::Named("s_best", s_best),
        Rcpp::Named("best_ind", best_ind + 1),    // arrary starts from 1 in R 
        Rcpp::Named("best_i", best_i + 1),    // arrary starts from 1 in R
        Rcpp::Named("best_j", best_j + 1)    // arrary starts from 1 in R
    );
    return res;
}

// [[Rcpp::export]]
Rcpp::List RSAVS_LargeN_L2_Rcpp(const Eigen::MatrixXd x_mat, const Eigen::VectorXd y_vec, const int n, const int p, 
                             const char p1_type, Eigen::VectorXd p1_param, 
                             const char p2_type, Eigen::VectorXd p2_param, 
                             const Eigen::VectorXd lam1_vec, const Eigen::VectorXd lam2_vec,
                             const double r2, const double r3, 
                             const double phi, const double tol, const double max_iter){
/*
    This is the main ADMM algorithm for solving RSAVA when n is large.
    In this algorithm, we have lam1(for penalty 1, hence \mu) and lam2(for penalty 2, hence \beta) to iterate.
    Outside the function, it's designed as a matrix:
            +--- lam_2 ---+
            ||    -->     |
            ||    -->     |
            |+----------> |
          lam_1           |
            |             |
            |             |
            |             |
            +-------------+
    But in this function, it's converted to a vector in a column-major pattern.
    Refer to the function Compute_Index for the details.
    Args:
    Returns:
*/          

// prepare some variables
    int i, j, j_start, i_old, j_old, ind, ind_old, k, best_ind, best_i, best_j;
    /*
        i, j: iterator over lam1_vec and lam2_vec
        j_start: starting index of j. When i = 0, j starts from 1, otherwise j starts from 0
        i_old, j_old: the index of initial values for current (i, j)
        ind, ind_old: the index in the result matrix corresponding to coordinate (i, j), (i_old, j_old) of the lambda matrix.
        k: inner loop counter for the algorithm
        best_ind: the row index of the resulting matrix for the best solution
        best_i, best_j: the row and column index of the lambda matrix for the best solution
    */
    double lam1, lam2;
    double diff = 1.0;
    double primal_residual, dual_residual, res_part1, res_part2;
    const int lam1_len = lam1_vec.size();
    const int lam2_len = lam2_vec.size();
    Eigen::MatrixXd beta_mat = Eigen::MatrixXd::Zero(lam1_len * lam2_len, p);
    Rcpp::Function r_median("median");    // Access the R function median
    double y_median = Rcpp::as<double>(r_median(y_vec));    // For the initial value of mu vector
    Eigen::MatrixXd mu_mat = Eigen::MatrixXd::Constant(lam1_len * lam2_len, n, y_median);
    // Eigen::MatrixXd z_mat = Eigen::MatrixXd::Zero(lam1_len * lam2_len, n);
    // Eigen::MatrixXd s_mat = Eigen::MatrixXd::Zero(lam1_len * lam2_len, n * (n - 1) / 2);    // danger for large n
    Eigen::MatrixXd w_mat = Eigen::MatrixXd::Zero(lam1_len * lam2_len, p);
    // Eigen::MatrixXd q1_mat = Eigen::MatrixXd::Zero(lam1_len * lam2_len, n);
    // Eigen::MatrixXd q2_mat = Eigen::MatrixXd::Zero(lam1_len * lam2_len, n * (n - 1) / 2);    // danger for large n
    Eigen::MatrixXd q3_mat = Eigen::MatrixXd::Zero(lam1_len * lam2_len, p);
    Eigen::MatrixXd bic_mat = Eigen::MatrixXd::Zero(lam1_len, lam2_len);
    Eigen::MatrixXd k_mat = Eigen::MatrixXd::Zero(lam1_len, lam2_len);    // this matrix can tell us whether or not the algorithm converged at that lam1 and lam2
    // void (*UpdateZ)(VectorXd &, const Eigen::VectorXd &, const double &);
    void (*UpdateS)(VectorXd &, const Eigen::VectorXd &, const double &);
    void (*UpdateW)(VectorXd &, const Eigen::VectorXd &, const double &);
    
    // variables needed during algorithm
    Eigen::VectorXd beta_vec = Eigen::MatrixXd::Zero(p, 1); 
    Eigen::VectorXd beta_old = Eigen::MatrixXd::Zero(p, 1);
    
    Eigen::VectorXd mu_vec = Eigen::MatrixXd::Constant(n, 1, y_median); 
    Eigen::VectorXd mu_old = Eigen::MatrixXd::Constant(n, 1, y_median);
    // Eigen::VectorXd z_vec = Eigen::MatrixXd::Zero(n, 1);
    // Eigen::VectorXd z_old = Eigen::MatrixXd::Zero(n, 1);
    Eigen::VectorXd s_vec = Eigen::MatrixXd::Zero(n * (n - 1) / 2, 1);
    Eigen::VectorXd s_old = Eigen::MatrixXd::Zero(n * (n - 1) / 2, 1);
    Eigen::VectorXd s_anchor = Eigen::MatrixXd::Zero(n * (n - 1) / 2, 1);
    Eigen::VectorXd s_best = Eigen::MatrixXd::Zero(n * (n - 1) / 2, 1);
    Eigen::VectorXd w_vec = Eigen::MatrixXd::Zero(p, 1);
    Eigen::VectorXd w_old = Eigen::MatrixXd::Zero(p, 1);
    // Eigen::VectorXd q1_vec = Eigen::MatrixXd::Zero(n, 1);
    // Eigen::VectorXd q1_old = Eigen::MatrixXd::Zero(n, 1);
    Eigen::VectorXd q2_vec = Eigen::MatrixXd::Zero(n * (n - 1) / 2, 1);
    Eigen::VectorXd q2_old = Eigen::MatrixXd::Zero(n * (n - 1) / 2, 1); 
    Eigen::VectorXd q2_anchor = Eigen::MatrixXd::Zero(n * (n - 1) / 2, 1);
    Eigen::VectorXd q3_vec = Eigen::MatrixXd::Zero(p, 1);
    Eigen::VectorXd q3_old = Eigen::MatrixXd::Zero(p, 1);
    
    // variable and function for summarizing the result
    Rcpp::Function r_summary_iteration("RSAVS_Summary_Iteration");
    Rcpp::List current_iteration_info;
    Eigen::MatrixXd mu_improve_mat = Eigen::MatrixXd::Zero(lam1_len * lam2_len, n);    // the improved mu vector
    Eigen::MatrixXd group_num_mat = Eigen::MatrixXd::Constant(lam1_len, lam2_len, 1);    // the group number
    Eigen::MatrixXd active_num_mat = Eigen::MatrixXd::Constant(lam1_len, lam2_len, 0);    // the number of active covariates 
    
    // prepare some values that will be needed in the algorithm
    // Generate the D matrix in sparse form
    Eigen::SparseMatrix<double> d_mat = Generate_D_Matrix(n);
    
    // Generate the lefet matirx for updating beta. In this version, we pre-assume p < n
    Eigen::VectorXd beta_right = Eigen::MatrixXd::Zero(p, 1);
    Eigen::MatrixXd beta_left = 1.0 / n * x_mat.transpose() * x_mat + r3 / 2.0 * Eigen::MatrixXd::Identity(p, p);
    // beta_left = 1.0 / 2.0 * x_mat.transpose() * x_mat + r3 / 2.0 * Eigen::MatrixXd::Identity(p, p);    // use 1 / 2 \sum\loss, not 1 / n \sum\loss    
    LDLT<Eigen::MatrixXd> beta_left_solver(p);
    beta_left_solver.compute(beta_left);
    
    // Generate the left matrix for updating mu
    Eigen::VectorXd mu_right = Eigen::MatrixXd::Zero(n, 1);
    // Eigen::MatrixXd mu_left = r2 * d_mat.transpose() * d_mat + r1 * Eigen::MatrixXd::Identity(n, n);
    // use a fashion recommended by Eigen
    Eigen::MatrixXd mu_left = 1.0 / n * Eigen::MatrixXd::Identity(n, n);
    // mu_left = 1.0 / 2.0 * Eigen::MatrixXd::Identity(n, n);    // use 1 / 2 \sum\loss, not 1 / n \sum\loss
    mu_left += r2 / 2.0 * d_mat.transpose() * d_mat;    
    LDLT<Eigen::MatrixXd> mu_left_solver(n);
    mu_left_solver.compute(mu_left);
    
// assign the correct updating functions    
    // UpdateZ = UpdateZ_L1;    // default loss type is L1
    // if(loss_type == 'H'){
    //     UpdateZ = UpdateZ_Huber;
    // }
    
    UpdateS = UpdateS_SCAD;    // default penalty type for s is SCAD
    if(p1_type == 'L'){
        UpdateS = UpdateS_Lasso;
    } else{
        if(p1_type == 'M'){
            UpdateS = UpdateS_MCP;
        }
    }
    
    UpdateW = UpdateW_SCAD;    // default penalty type for w is SCAD
    if(p2_type == 'L'){
        UpdateW = UpdateW_Lasso;
    } else{
        if(p2_type == 'M'){
            UpdateW = UpdateW_MCP;
        }
    }
    
// main algorithm
/*
    The algorithm iterates over the lambda matrix.
    The outer iteration is over lambda 1, the inner iteration is over lambda 2.
*/

    best_i = 0;
    best_j = 0;
    best_ind = 0;
    current_iteration_info = Rcpp::as<Rcpp::List>(r_summary_iteration(Rcpp::Named("y_vec", y_vec), 
                                                                      Rcpp::Named("x_mat", x_mat), 
                                                                      Rcpp::Named("mu_vec", mu_vec), 
                                                                      Rcpp::Named("beta_vec", beta_vec), 
                                                                      Rcpp::Named("s_vec", s_vec), 
                                                                      Rcpp::Named("w_vec", w_vec), 
                                                                      Rcpp::Named("loss_type", "2"), 
                                                                      Rcpp::Named("loss_param", 0.0), 
                                                                      Rcpp::Named("phi", phi)));
    bic_mat(best_i, best_j) = Rcpp::as<double>(current_iteration_info["bic"]);
    mu_improve_mat.row(best_ind) = Rcpp::as<Eigen::VectorXd>(current_iteration_info["mu_vec"]);
    
    for(i = 0; i < lam1_len; i++){
        // find the appropriate starting point of j
        if(i == 0){
            j_start = 1;
        }else{
            j_start = 0;
        }
        
        for(j = j_start; j < lam2_len; j++){
            // assign lambda 1 and lambda 2
            lam1 = lam1_vec[i];
            lam2 = lam2_vec[j];
            p1_param[0] = lam1;
            p2_param[0] = lam2;
            
            // find (i_old, j_old), coordinate of initial values
            i_old = i - (j == 0);
            j_old = std::max(0, j - 1);
            
            // compute the index
            ind = Compute_Index(i, j, lam1_len, lam2_len);
            ind_old = Compute_Index(i_old, j_old, lam1_len, lam2_len);
            
            // prepare initial values
            beta_old = beta_mat.row(ind_old);
            mu_old = mu_mat.row(ind_old);
            // z_old = z_mat.row(ind_old);
            w_old = w_mat.row(ind_old);
            // q1_old = q1_mat.row(ind_old);
            q3_old = q3_mat.row(ind_old);
            if(j == 0){
                q2_old = q2_anchor;
                s_old = s_anchor;
            }
            // s_old = s_mat.row(ind_old);    // danger for large n
            // q2_old = q2_mat.row(ind_old);    // danger for large n
            
            // Rcpp::Rcout << "i = " << i << ", j = " << j << ", ind = " << ind << ". lam1 = " << lam1 << ", lam2 = " << lam2 << std::endl;
            // Rcpp::Rcout << "i_old = " << i_old << ", j_old = " << j_old << ", ind_old = " << ind_old << std::endl;
            // Rcpp::Rcout << "lam1 = " << lam1 << ", lam2 = " << lam2 << std::endl;
            // update steps
            k = 0;
            diff = 1.0;
            while((k < max_iter) && (diff > tol)){
                // update beta
                beta_right = x_mat.transpose() * (y_vec - mu_old) / n + r3 / 2.0 * w_old - q3_old / 2.0;    // 1 / n \sum\loss
                // beta_right = x_mat.transpose() * (y_vec - mu_old) / 2.0 + r3 / 2.0 * w_old - q3_old / 2.0;    // use 1 / 2 \sum\loss, not 1 / n \sum\loss
                beta_vec = beta_left_solver.solve(beta_right);
                
                // update mu
                mu_right = (y_vec - x_mat * beta_vec) / n;    // 1 / n \sum\loss
                // mu_right = (y_vec - x_mat * beta_vec) / 2.0;    // use 1 / 2 \sum\loss, not 1 / n \sum\loss
                mu_right += d_mat.transpose() * (r2 * s_old - q2_old) / 2.0;
                mu_vec = mu_left_solver.solve(mu_right);
                
                // update z
                // z_vec = y_vec - mu_vec - x_mat * beta_vec + q1_old / r1;
                // UpdateZ(z_vec, loss_param, r1);
                
                // update s
                s_vec = d_mat * mu_vec + q2_old / r2;
                UpdateS(s_vec, p1_param, r2);
                
                // update w
                w_vec = beta_vec + q3_old / r3;
                UpdateW(w_vec, p2_param, r3);
                
                // update q1, q2 and q3
                // q1_vec = q1_old + r1 * (y_vec - mu_vec - x_mat * beta_vec - z_vec);
                q2_vec = q2_old + r2 * (d_mat * mu_vec - s_vec);
                q3_vec = q3_old + r3 * (beta_vec - w_vec);
                
                // compute residual and check for convergence
                k = k + 1;    // update the inner loop counter
                // primal residual
                // res_part1 = (y_vec - mu_vec - x_mat * beta_vec - z_vec).squaredNorm();
                res_part1 = (d_mat * mu_vec - s_vec).squaredNorm();
                res_part2 = (beta_vec - w_vec).squaredNorm();
                primal_residual = res_part1 + res_part2;
                
                // dual residual
                res_part1 = (r2 * d_mat.transpose() * (s_vec - s_old)).squaredNorm();
                res_part2 = (r3 * (w_vec - w_old)).squaredNorm();
                dual_residual = res_part1 + res_part2;
                
                diff = std::max(primal_residual, dual_residual);
                
                // prepare for next inner iteration
                beta_old = beta_vec;
                mu_old = mu_vec;
                // z_old = z_vec;
                s_old = s_vec;
                w_old = w_vec;
                // q1_old = q1_vec;
                q2_old = q2_vec;
                q3_old = q3_vec;
            }
            
            // save the result
            beta_mat.row(ind) = beta_vec;
            mu_mat.row(ind) = mu_vec;
            // z_mat.row(ind) = z_vec;
            w_mat.row(ind) = w_vec;
            // q1_mat.row(ind) = q1_vec;
            q3_mat.row(ind) = q3_vec;
            s_old = s_vec;
            q2_old = q2_vec;
            if(j == 0){
                s_anchor = s_vec;
                q2_anchor = q2_vec;
            }
            // s_mat.row(ind) = s_vec;    // danger for large n
            // q2_mat.row(ind) = q2_vec;    // danger for large n
            
            k_mat(i, j) = k;
            // summary this result
            current_iteration_info = Rcpp::as<Rcpp::List>(r_summary_iteration(Rcpp::Named("y_vec", y_vec), 
                                                                              Rcpp::Named("x_mat", x_mat), 
                                                                              Rcpp::Named("mu_vec", mu_vec), 
                                                                              Rcpp::Named("beta_vec", beta_vec), 
                                                                              Rcpp::Named("s_vec", s_vec), 
                                                                              Rcpp::Named("w_vec", w_vec), 
                                                                              Rcpp::Named("loss_type", "2"), 
                                                                              Rcpp::Named("loss_param", 0.0), 
                                                                              Rcpp::Named("phi", phi)));
            bic_mat(i, j) = Rcpp::as<double>(current_iteration_info["bic"]);
            mu_improve_mat.row(ind) = Rcpp::as<Eigen::VectorXd>(current_iteration_info["mu_vec"]);
            group_num_mat(i, j) = Rcpp::as<double>(current_iteration_info["group_num"]);
            active_num_mat(i, j) = Rcpp::as<double>(current_iteration_info["active_num"]);
            
            // check whether we should update s_best
            if(bic_mat(i, j) < bic_mat(best_i, best_j)){
                best_i = i; 
                best_j = j;
                best_ind = ind;
                s_best = s_vec;
            }
        }
    }
 
// create and return the results
    Rcpp::List res = Rcpp::List::create(
        Rcpp::Named("beta_mat", beta_mat), 
        Rcpp::Named("mu_mat", mu_mat), 
        Rcpp::Named("mu_improve_mat", mu_improve_mat),
        // Rcpp::Named("z_mat", z_mat), 
        // Rcpp::Named("s_mat", s_mat),    // danger for large n
        Rcpp::Named("w_mat", w_mat),
        // Rcpp::Named("q1_mat", q1_mat),
        Rcpp::Named("q3_mat", q3_mat),
        Rcpp::Named("bic_mat", bic_mat),
        Rcpp::Named("k_mat", k_mat),
        Rcpp::Named("group_num_mat", group_num_mat),
        Rcpp::Named("active_num_mat", active_num_mat), 
        Rcpp::Named("lam1_vec", lam1_vec),
        Rcpp::Named("lam2_vec", lam2_vec), 
        Rcpp::Named("s_best", s_best),
        Rcpp::Named("best_ind", best_ind + 1),    // arrary starts from 1 in R 
        Rcpp::Named("best_i", best_i + 1),    // arrary starts from 1 in R
        Rcpp::Named("best_j", best_j + 1)    // arrary starts from 1 in R
    );
    return res;
}



Rcpp::List RSAVS_Solver_Cpp(const Eigen::VectorXd& y_vec, const Eigen::MatrixXd& x_mat, const int& n, const int& p, 
                            const std::string& l_type, const Eigen::VectorXd& l_param, 
                            const std::string& p1_type, const Eigen::VectorXd p1_param, 
                            const std::string& p2_type, const Eigen::VectorXd p2_param, 
                            const Eigen::VectorXd& const_r123, const Eigen::VectorXd& const_abc, 
                            const double& tol, const int& max_iter, 
                            const double& cd_tol, const int& cd_max_iter, 
                            const Rcpp::List& initial_values,
                            const double& phi){
    
    // Update functions for z, s and w
    void (*Update_Z)(Eigen::VectorXd &, const Eigen::VectorXd &, const double &);
    void (*Update_S)(Eigen::VectorXd &, const Eigen::VectorXd &, const double &);
    void (*Update_W)(Eigen::VectorXd &, const Eigen::VectorXd &, const double &);
    
    Update_Z = UpdateZ_L1;
    if(l_type == "L2"){
        Update_Z = UpdateZ_L2;
    }
    if(l_type == "Huber"){
        Update_Z = UpdateZ_Huber;
    }
    
    Update_S = UpdateW_Lasso;
    if(p1_type == "S"){
        Update_S = UpdateW_SCAD;
    }
    if(p1_type == "M"){
        Update_S = UpdateW_MCP;
    }
    
    Update_W = UpdateW_Lasso;
    if(p2_type == "S"){
        Update_W = UpdateW_SCAD;
    }
    if(p2_type == "M"){
        Update_W = UpdateW_MCP;
    }
    
    // prepare initial values for the algorithm
    Eigen::VectorXd beta_old = Rcpp::as<Eigen::VectorXd>(initial_values["beta_init"]);
    Eigen::VectorXd mu_old = Rcpp::as<Eigen::VectorXd>(initial_values["mu_init"]);
    Eigen::VectorXd z_old = Rcpp::as<Eigen::VectorXd>(initial_values["z_init"]);
    Eigen::VectorXd s_old = Rcpp::as<Eigen::VectorXd>(initial_values["s_init"]);
    Eigen::VectorXd w_old = Rcpp::as<Eigen::VectorXd>(initial_values["w_init"]);
    Eigen::VectorXd q1_old = Rcpp::as<Eigen::VectorXd>(initial_values["q1_init"]);
    Eigen::VectorXd q2_old = Rcpp::as<Eigen::VectorXd>(initial_values["q2_init"]);
    Eigen::VectorXd q3_old = Rcpp::as<Eigen::VectorXd>(initial_values["q3_init"]);
    
    // prepare variables for the algorithm
    Eigen::VectorXd beta_vec = Eigen::MatrixXd::Zero(p, 1);
    Eigen::VectorXd mu_vec = Eigen::MatrixXd::Zero(n, 1);
    Eigen::VectorXd z_vec = Eigen::MatrixXd::Zero(n, 1);
    Eigen::VectorXd s_vec = Eigen::MatrixXd::Zero(n * (n - 1) / 2, 1);
    Eigen::VectorXd w_vec = Eigen::MatrixXd::Zero(p, 1);
    Eigen::VectorXd q1_vec = Eigen::MatrixXd::Zero(n, 1);
    Eigen::VectorXd q2_vec = Eigen::MatrixXd::Zero(n * (n - 1) / 2, 1);
    Eigen::VectorXd q3_vec = Eigen::MatrixXd::Zero(p, 1);
    
    
