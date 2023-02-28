#include "FEM.hpp" 
#include <fstream>
#include <iostream>



Poisson::Poisson(string file_name, int degree){
    load(file_name);
    //construct(degree);
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
        cout << coords[i] << "\n";
    }

    for(int i = 0; i < M; i++){
        ifile >> v1;
        ifile >> v2; 
        edges.push_back(make_pair(v1, v2));
    }
}



void Poisson::construct(int deg){
    for(int i = deg; i >= 0; i--){
        for(int j = 0; j < N - i ; j++)
            MAT[i].push_back(GaussianQuadrature());
    }

    for(int i = 0; i < N; i++)
        F.push_back(GaussianQuadrature());
}



double Poisson::GaussianQuadrature(){
    
    return 0;
}