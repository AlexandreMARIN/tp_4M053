#include "include/Matrix.hpp"
#include <chrono>
#include <fstream>

using namespace std;
using namespace std::chrono;


Matrix create_A_N(int N){
  Matrix A_N(N);

  A_N.diag(2.0);
  A_N.diag(-1.0, 1);
  A_N.diag(-1.0, -1);

  return A_N;
}


int main(){

  Matrix A(0), B(0);
  high_resolution_clock::time_point begin, end;
  vector<int> sizes{5, 10, 20, 30, 60, 120};
  vector<int> LU_dur, Cholesky_dur;
  microseconds dur;

  LU_dur.reserve(sizes.size());
  Cholesky_dur.reserve(sizes.size());

  for(int s : sizes){
    A = create_A_N(s);
    B = A;
    begin = high_resolution_clock::now();
    A.decomp_LU();
    end = high_resolution_clock::now();
    dur = duration_cast<microseconds>(end - begin);
    LU_dur.push_back(dur.count());

    begin = high_resolution_clock::now();
    B.cholesky();
    end = high_resolution_clock::now();
    dur = duration_cast<microseconds>(end - begin);
    Cholesky_dur.push_back(dur.count());
  }

  //we write here a .json file
  ofstream file("LU.json");

  file << "{\n\t\"sizes\" : [" << sizes[0];
  for(vector<int>::size_type i = 1;i<sizes.size();i++){
    file << ", " << sizes[i];
  }
  file << "],\n\n\t\"LU_dur\" : [" << LU_dur[0];
  for(vector<int>::size_type i = 1;i<LU_dur.size();i++){
    file << ", " << LU_dur[i];
  }
  file << "],\n\n\t\"Cholesky_dur\" : [" << Cholesky_dur[0];
  for(vector<int>::size_type i = 1;i<Cholesky_dur.size();i++){
    file << ", " << Cholesky_dur[i];
  }
  file << "]\n\n}";

  return 0;
}
