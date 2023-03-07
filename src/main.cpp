#include "FEM.hpp"
#include <iostream>
#include <vector>

using namespace std;

double func(const vector<double>& point){
    return 10;
}

int main(){
    Poisson Prueba("TEST/prueba.txt", 1, &func, 1);
    Prueba.printSystem();
    return 0;
}