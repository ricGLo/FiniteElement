#include <vector>
#include <string>

using namespace std;

class Poisson{
    public:
    //Read mesh data, construct stiffness matrix and the forcing vector
        Poisson(string file_name, int degree);

    //Solve linear system
        vector<double> Solve(string method);

    protected:
        vector<vector<double>> MAT; //Stiffness Matrix
        vector<double> F; //Forcing vector
        unsigned int N; //number of vertices
        unsigned int M; //number of edges
        double GaussianQuadrature();

    //Mesh data
        vector<double> coords;
        vector<pair<unsigned int, unsigned int>> edges;

    
    private:
    //Load mesh
        void load(string file_name);

    //Construct system
        void construct(int degree);


};