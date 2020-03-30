#include "aghMatrix.h"


// Parameter Constructor                                                                                                                                                      
template<typename T>
AGHMatrix<T>::AGHMatrix(const std::vector<std::vector<T>>& mat) 
{
  matrix.resize(mat.size());
  for (unsigned i = 0; i < mat.size(); i++) 
  {
    matrix[i].resize(mat[i].size());
    for(unsigned j = 0; j < mat[i].size(); j++)
    {
      matrix[i][j] = mat[i][j];
    }
  }
  rows = matrix.size();
  cols = matrix[1].size();
}

// Parameter Constructor                                                                                                                                                      
template<typename T>
AGHMatrix<T>::AGHMatrix(unsigned _rows, unsigned _cols, const T& _initial) 
{
  matrix.resize(_rows);
  for (unsigned i=0; i<matrix.size(); i++) 
  {
    matrix[i].resize(_cols, _initial);
  }
  rows = _rows;
  cols = _cols;
}

// Copy Constructor                                                                                                                                                           
template<typename T>
AGHMatrix<T>::AGHMatrix(const AGHMatrix<T>& rhs) 
{
  matrix = rhs.matrix;
  rows = rhs.get_rows();
  cols = rhs.get_cols();
}

// Get the number of rows of the matrix                                                                                                                                       
template<typename T>
unsigned AGHMatrix<T>::get_rows() const 
{
  return this->rows;
}

// Get the number of columns of the matrix                                                                                                                                    
template<typename T>
unsigned AGHMatrix<T>::get_cols() const 
{
  return this->cols;
}

// Assignment Operator                                                                                                                                                        
template<typename T>
AGHMatrix<T>& AGHMatrix<T>::operator=(const AGHMatrix<T>& rhs) 
{
  if (&rhs == this)
    return *this;

  unsigned new_rows = rhs.get_rows();
  unsigned new_cols = rhs.get_cols();

  matrix.resize(new_rows);
  for (unsigned i=0; i<matrix.size(); i++) 
  {
    matrix[i].resize(new_cols);
  }

  for (unsigned i=0; i<new_rows; i++) 
  {
    for (unsigned j=0; j<new_cols; j++) 
    {
      matrix[i][j] = rhs(i, j);
    }
  }
  rows = new_rows;
  cols = new_cols;

  return *this;
}

// Access the individual elements                                                                                                                                             
template<typename T>
T& AGHMatrix<T>::operator()(const unsigned& row, const unsigned& col) 
{
  return this->matrix[row][col];
}

// Access the individual elements (const)                                                                                                                                     
template<typename T>
const T& AGHMatrix<T>::operator()(const unsigned& row, const unsigned& col) const 
{
  return this->matrix[row][col];
}

// Addition of two matrices                                                                                                                                                   
template<typename T>
AGHMatrix<T> AGHMatrix<T>::operator+(const AGHMatrix<T>& rhs) 
{
  // Task 1 - implement addition of two matrices
  AGHMatrix<T> res(rows, cols, matrix[0][0]);
  if(cols != rhs.cols || rows!= rhs.rows){
      std::cout << "blad dodawania zle wymiary macierzy" << std::endl;
      exit(-1);
  }
  for(int i = 0; i < rows; i++){
      for(int j = 0; j < cols; j++){
          res.matrix[i][j] = (rhs.matrix[i][j] + matrix[i][j]);
      }
  }

  return res;
}

// Left multiplication of this matrix and another                                                                                                                              
template<typename T>
AGHMatrix<T> AGHMatrix<T>::operator*(const AGHMatrix<T>& rhs) 
{
  // Task 1 - implement multiplication of two matrices
  if(cols != rhs.rows){
      std::cout << "blad dodawania zle wymiary macierzy" << std::endl;
      exit(-1);
  }
  AGHMatrix<T> res(rows, rhs.get_cols());
  
  for(int i = 0; i < rows; i++){
      for(int j = 0; j < rhs.cols;j++){
          int sum = 0;
          for(int k = 0; k < cols; k++){
              sum += matrix[i][k] * rhs.matrix[k][j];
          }
          matrix[i][j] = sum;
      }
  }
  
}

// Printing matrix                                                                                                                        
template<typename T>
std::ostream& operator<<(std::ostream& stream, const AGHMatrix<T>& matrix) 
{
  for (int i=0; i<matrix.get_rows(); i++) 
  { 
    for (int j=0; j<matrix.get_cols(); j++) 
    {
        stream << matrix(i,j) << ", ";
    }
    stream << std::endl;
  }
  stream << std::endl;
  return stream;
}

//zad1
template<typename T>
bool AGHMatrix<T>::check_symmetric(void){
    if(cols != rows){
        return false;
    }
    for(int i = 0; i < rows; i++){
        for(int j = 0; j < cols; j++){
            if(matrix[i][j] != matrix[j][i]){
                return false;
            }
        }
    }
    return true;
}

//zad3
template<typename T>
void AGHMatrix<T>::swap(unsigned i, unsigned j, unsigned k, unsigned l){
    T temp;
    temp = matrix[i][j];
    matrix[i][j] = matrix[k][l];
    matrix[k][l] = temp;
}
template<typename T>
void AGHMatrix<T>::transpose(void){
    // if(rows != cols){
    //     std::cout << "macierz musi byc kwadratowa" << std::endl;
    //     exit(-1);
    // }
        for (int i = 0; i < rows; i++){
            for (int j = i+1; j < rows; j++){
                swap(i, j, j, i);               
            }
        }


        
}

void LU_double(AGHMatrix<double> &a, AGHMatrix<double> &l, AGHMatrix<double> &u){
  int n;
  n = a.get_cols();


  for(int k = 0; k < n; k++){
    l(k,k) = 1;
    for(int j = k; j < n; j++){
      double sum = 0;
      for(int s = 0; s < k; s++) sum += l(k,s)*u(s,j);

      u(k,j) = a(k,j) - sum;
      }for(int i = k; i < n; i++){

          double sum = 0;
          for(int s= 0; s < k; s++) sum += l(i,s)*u(s,k);

          l(i,k) = (a(i,k) - sum) / u(k,k);

      }
    
  }
}
double diagonal_det(AGHMatrix<double> &a){
  double det = 1;
  for(unsigned i = 0; i < a.get_cols(); i++){
      det *= a(i,i);
  }
  return det;
}


double LUdet(AGHMatrix<double> &a){
  double det;
  AGHMatrix<double> l(a.get_rows(), a.get_cols(), 0.0);
  AGHMatrix<double> u(a.get_rows(), a.get_cols(), 0.0);
  LU_double(a, l, u);

  det =  diagonal_det(u);
  return det;
}

void cholesky_decomposition(AGHMatrix<double> &a, AGHMatrix<double> &l){
  int n;
  n = a.get_cols();


  for(int k = 0; k < n; k++){
    double sum = 0;
    for(int s = 0; s < k; s++) sum += (l(k,s))*(l(k,s));
    l(k,k) = sqrt(a(k,k) - sum);

    for(int i = k + 1; i < n; i++){
      double sum = 0;
      for(int s = 0; s < k; s++) sum += l(i,s)*l(k,s);
      l(i,k) = (a(i,k) - sum)/l(k,k);
    }
  }


}

double Cholesky_det(AGHMatrix<double> &a){
  double det;
  AGHMatrix<double> l(a.get_rows(), a.get_cols(), 0.0);
  
  cholesky_decomposition(a, l);
  det =  diagonal_det(l);
  return det*det;
}

AGHMatrix<double> Gauss_alg(AGHMatrix<double> &a){
  double eps = 1e-12;
  int rows = a.get_rows();
  int cols = a.get_cols();
  for(int i = 0; i < rows; i++){ // dla każdego wiersza
    for(int j = i + 1; j < rows; j++){ // dla każdego większego wiersza
      if( fabs(a(i,i)) < eps){ // jeśli byłoby dzielenie przez zero
        std::cout << "błąd";
        exit(-1);
      }
      //obliczamy współczynnik przez który przemnożymy wszystkie kolumny wiersza, aby w "pierwszej" pojawiło się zero
      double m = -a(j,i)/ a(i,i);
      for(int k = i; k < cols; k++){
        a(j,k) += m * a(i,k); // dodajemy do całego wiersza i-ty wiersz przemnożony przez współczynnik
      }
    }

  }
  return a;
}

AGHMatrix<double> gauss_solve(AGHMatrix<double> &a){
  Gauss_alg(a);
  double eps = 1e-12;
  int n = a.get_rows();
  AGHMatrix<double> X(a.get_rows(), 1, 0.0);

  for(int i = n-1; i >= 0; i--){ // zaczynamy od konca
    double s = a(i,n);
    for(int j = n - 1; j >= i + 1; j--){ // liczenie rożnicy dla wszystkich wierszow w dol
      s -= a(i,j) * X(j,0);
    }
    if(fabs(a(i,i)) < eps) exit(-1); // na wypadek dzielenia przez zero
    X(i, 0) = s / a(i,i); // obliczanie wartosci niewiadomej
    
  }
  return X;

}


AGHMatrix<double> Jacoby(AGHMatrix<double> &a, int param){

  AGHMatrix<double> XX(a.get_rows(), 2, 0.0);
  AGHMatrix<double> N(a.get_rows(), 1, 0.0);
  AGHMatrix<double> M(a.get_rows(), a.get_rows(), 0.0);

  int n = a.get_rows();
  int m = a.get_cols();


  for(int i = 0; i < n; i++){
    N(i,0) = 1/a(i,i);
  }

for (int i=0; i< n; i++)
  for (int j=0; j< n; j++)
    if (i == j)
      M(i,j) = 0;
    else
      M(i,j) = - (a(i,j) * N(i,0));

  for(int k = 0; k < param; k++){
    for(int i = 0; i < n; i++){
      double sum = 0;
      for(int s = 0; s < n; s++){
        if(i != s){
          sum += a(i,s) * XX(s,0);
        }
      }
      XX(i,1) = (1.0/a(i,i)) * ( a(i,m-1) - sum);
    }
    for (int i = 0; i < n; i++)
      XX(i,0) = XX(i,1);
  }

  return XX;

}