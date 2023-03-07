#include <vector>
#include <string>
#include <math.h>
#include <functional>

using namespace std;

class Poisson{
    public:
    //Read mesh data, construct stiffness matrix and the forcing vector
        Poisson(string file_name, int deg, function<double(const vector<double>&)> func, int dim_);

    //Solve linear system
        vector<double> SolveDirichlet(string method);

    //Export stiffness matrix and load vector
        void printSystem();

    protected:
        vector<vector<double>> MAT; //Stiffness Matrix
        vector<double> F; //Forcing vector
        function<double( const vector<double>&)> f; // Initial condition function
        unsigned int N; //number of vertices
        unsigned int M; //number of edges
        int dim;
        int degree;
        double IntegrateLeft(int deg, unsigned int i, unsigned int j);
        double IntegrateRight(int deg, unsigned int i);


        // Table of weights and roots of Legendre's roots
        vector<vector<double>> LegRoots{
            {},
            {0.0}, 
            {-1.0/sqrt(3), 1.0/sqrt(3)},
            {-sqrt(3.0/5.0), 0, sqrt(3.0/5.0)},
            {-sqrt( (3.0 / 7.0) + (2.0/7.0)*sqrt(6.0 / 5.0)), -sqrt( (3.0 / 7.0) - (2.0/7.0)*sqrt(6.0/5.0)), sqrt( (3.0/7.0) - (2.0/7.0)*sqrt(6.0/5.0)), sqrt( (3.0/7.0) + (2.0/7.0)*sqrt(6.0/5.0)) },
            {-sqrt(5 + 2*sqrt(10.0/7.0)) / 3.0, -sqrt(5-2*sqrt(10.0/7.0)) / 3.0, 0, sqrt(5-2*sqrt(10.0/7.0)) / 3.0, sqrt(5+2*sqrt(10.0/7.0))/3.0}
        };
        vector<vector<double>> LegWeights{
            {},
            {2.0}, 
            {1.0, 1.0},
            {5.0/9.0, 8.0/9.0, 5.0/9.0},
            {(18.0 - sqrt(30.0)) / 36.0, (18.0 + sqrt(30.0)) / 36.0, (18.0 + sqrt(30.0)) / 36.0, (18.0 - sqrt(30.0)) / 36.0 },
            { (322.0 - 13*sqrt(70)) / 900.0, (322.0 + 13*sqrt(70)) / 900.0, 128.0 / 225.0 , (322.0 + 13*sqrt(70)) / 900.0 , (322 - 13*sqrt(70)) / 900.0 }
        };

        double v(double x, unsigned int i){
            if( (i > 0) && (coords[i-1] <= x && x <= coords[i]) )
                return (x - coords[i-1]) / (coords[i] - coords[i-1]);
            else if ( (i < N-1) && (coords[i] <= x && x <= coords[i+1]) )
                return (coords[i+1] - x) / (coords[i+1] - coords[i]);
            else 
                return 0;
        };

        double dv(double x, unsigned int i){
            if( (i > 0) && (coords[i-1] <= x && x <= coords[i]) )
                return 1.0;
            else if ( (i < N-1) && (coords[i] <= x && x <= coords[i+1]) )
                return -1.0;
            else 
                return 0;
        };

    //Mesh data
        vector<double> coords;
        vector<pair<unsigned int, unsigned int>> edges;

    
    private:
    //Load mesh
        void load(string file_name);

    //Construct system
        void construct();


};