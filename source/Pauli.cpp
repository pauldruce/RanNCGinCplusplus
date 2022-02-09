#include "Pauli.h"

using namespace arma;
using namespace std;

Pauli::Pauli() :
Clifford(3,0)
{
   this->x = this->generators[1];
   this->y = -1 * this->generators[2];
   this->z = this->generators[0];
}

void Pauli::introduce()
{
   cout<< "My type is: ("<< p <<","<< q << ")\n";
   cout<<"My matrices will be"<< matrix_size<< "x" <<matrix_size<< ", and I have K0 dimension s ="<< s <<endl;
   
   // Print out the generators
   cout<< "The generators are:"<< endl;
   cout << "pauli.x = " << endl;
   this->x.print();
   cout << "pauli.y = "<< endl;
   this->y.print();
   cout << "pauli.z = "<< endl;
   this-> z.print();
   cout << endl;
   
   cout<< "The chirality operator is:"<< endl;
   chirality.print();
   
}
