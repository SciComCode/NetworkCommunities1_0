

#include "FuzzyCommunities.h"	      // Allowing access to the functionality in class FuzzyCommunities
using namespace NetworkCommunities;


int main(int argc, char** argv) {

  int c = 2;  // Number of communities


  //FuzzyCommunities* fc = new FuzzyCommunities("TestSoftX.edges", c);  // Creating object fc and reading network topology contained in file data.dat
 
  FuzzyCommunities* fc = new FuzzyCommunities("data.dat", c);  // Creating object fc and reading network topology contained in file data.dat
  fc->findCommunities();		// Calling method to determine memberships using the TRIBUNE algorithm.
  fc->writeU("results.out");		// Writing membership matrix to file results.out


}

