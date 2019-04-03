#include "./include/Matrix.hpp"

using namespace std;

int main(){

  cout << "Give an integer : ";
  int n;
  cin >> n;
  if(cin.bad()){
    cout << "Bad argument\n";
    return 1;
  }


  COO L = COO::Laplacian(n);

  L.display();




  return 0;
}
