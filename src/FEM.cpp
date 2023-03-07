#include "FEM.hpp" 
#include <fstream>
#include <iostream>


Poisson::Poisson(string file_name, int deg, function<double(const vector<double>&)> func, int dim_){
    f = func;
    dim = dim_;
    degree = deg;
    load(file_name);
    
    construct();
}


void Poisson::load(string file_name){
    ifstream ifile(file_name, ios::in);
    if(!ifile.is_open()){
        cerr << "There was a problem opening the text file\n";
        exit(1);
    }

    ifile >> N;
    ifile >> M;

    double num;
    int v1, v2;
    for(int i = 0; i < N; i++){
        ifile >> num;
        coords.push_back(num);
    }

    for(int i = 0; i < M; i++){
        ifile >> v1;
        ifile >> v2; 
        edges.push_back(make_pair(v1, v2));
    }
}

void Poisson::construct(){
    MAT.assign(degree+1, vector<double>());
    for(int i = degree; i >= 0; i--){
        for(unsigned int j = 0; j < N - i ; j++)
            MAT[i].push_back( IntegrateLeft(degree*degree, j, j+i));
    }
    for(unsigned int i = 0; i < N; i++)
        F.push_back(IntegrateRight(5, i));
}


double Poisson::IntegrateLeft(int num, unsigned int i, unsigned int j){
    double res = 0;
    double eval, ell, sum;

    if(i == 0){
        ell = coords[1] - coords[0];
        sum = coords[1] + coords[0];
    }
    else if(j == N-1){
        ell = coords[N-1] - coords[N-2]; 
        sum = coords[N-1] + coords[N-2];
    }  
    else if(j == i){
        ell = coords[j+1] - coords[i-1];
        sum = coords[j+1] + coords[i-1];
    }
    else{
        ell = coords[j] - coords[i];
        sum = coords[j] + coords[i];
    }


    for(int k = 0; k < num; k++){
        eval = ell*LegRoots[num][k] + sum;
        eval = eval / 2;
        res += LegWeights[num][k] * ( dv(eval, i) * dv(eval, j) ); 
    }
        
    res = res * (ell / 2);
    return res;
}

double Poisson::IntegrateRight(int num, unsigned int i){
    double res = 0;
    double res_aux = 0;
    vector<double> eval(dim); 

    double ell, sum;

    if(i == 0){
        ell = coords[1] - coords[0];
        sum = coords[1] + coords[0];
    }
    else if(i == N-1){
        ell = coords[N-1] - coords[N-2]; 
        sum = coords[N-1] + coords[N-2];
    }  
    else{
        ell = coords[i] - coords[i-1];
        sum = coords[i] + coords[i-1];

        for(int k = 0; k < num; k++){
            eval[0] = ell*LegRoots[num][k] + sum;
            eval[0] = eval[0] / 2;
            res_aux += LegWeights[num][k] * ( f(eval) * v(eval[0], i) ); 
        }
        res_aux = res_aux * (ell / 2);

        ell = coords[i+1] - coords[i];
        sum = coords[i+1] + coords[i];
    }

    for(int k = 0; k < num; k++){
        eval[0] = ell*LegRoots[num][k] + sum;
        eval[0] = eval[0] / 2;
        res += LegWeights[num][k] * ( f(eval) * v(eval[0], i) ); 
    }
    res = res * (ell / 2);

    return (res + res_aux);
}

void Poisson::printSystem(){
    cout << "\nDiagonal of the stiffness matrix: \n";
    for(int i = degree; i >= 0; i--){
        cout << "\n";
        for(unsigned int j = 0; j < N - i ; j++)
            cout << MAT[i][j] / 0.05263157894736842 << "\t";
    }
    

    cout << "\nLoad vector: \n";
    cout << F[0] * (2 / 0.05263157894736842) << "\t";
    for(unsigned int i = 1; i < N-1; i++)
        cout << F[i] / 0.05263157894736842 << "\t";
    cout << F[N-1] * (2 / 0.05263157894736842) << "\t";
    cout << "\n";

}
