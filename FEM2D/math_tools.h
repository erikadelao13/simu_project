#include <vector>
#include "math.h"
#include "stdlib.h"

using namespace std;

typedef vector<float> Vector;
typedef vector<Vector> Matrix;

void zeroes(Matrix &M,int n){
    for(int i=0;i<n;i++){
        vector<float> row(n,0.0);
        M.push_back(row);
    }
}

void zeroes(Matrix &M,int n,int m){
    for(int i=0;i<n;i++){
        vector<float> row(m,0.0);
        M.push_back(row);
    }
}

void zeroes(Vector &v,int n){
    for(int i=0;i<n;i++){
        v.push_back(0.0);
    }
}

void copyMatrix(Matrix A, Matrix &copy){
    zeroes(copy,A.size());
    for(int i=0;i<A.size();i++)
        for(int j=0;j<A.at(0).size();j++)
            copy.at(i).at(j) = A.at(i).at(j);
}

float calculateMember(int i,int j,int r,Matrix A,Matrix B){ 
    float member = 0;
    for(int k=0;k<r;k++) //la longitud que va a tener la fila i de A y la columna j de B
        member += A.at(i).at(k)*B.at(k).at(j); //se envian las dos matrices 
    return member;
}
//antes era matriz por numero real, hoy es matriz por matriz.
//recibo matriz a y matriz b, y recibo tres enteros. n es el numero de filas de a, m es el numero de columnas de b, de manera que r es el dato compartido, r corresponde tanto a las columnas de a, como a las columnas de b, que es para que se cumpla la condicion que menciono aca abajo
Matrix productMatrixMatrix(Matrix A,Matrix B,int n,int r,int m){ //recibo matriz A y matriz B que son las que voy a multiplicar, n es el numero de filas de a, m son las columnas de b, y r es el dato compartido que corresponde tanto a las columnas de a como a las filas de b, que es para que se cumpla la condicion que la columna de la primera tiene que ser igual al numero de filas de la segunda, de manera que se recorre utilizando de indices n y m  
    Matrix R;

    zeroes(R,n,m);
    //se hace un doble for para poder saber que posicion de la respuesta es la que se va a colocar.
    for(int i=0;i<n;i++) //fila i de la matriz a
        for(int j=0;j<m;j++) //columna j de la matriz b, todo esto se calcula en member
            R.at(i).at(j) = calculateMember(i,j,r,A,B); //calculatemember calcula el dato que va a ir en una casilla y devuelve un float y recibe los indices de la ecuacion que se va a calcular

    return R;
}

void productMatrixVector(Matrix A, Vector v, Vector &R){
    for(int f=0;f<A.size();f++){
        float cell = 0.0;
        for(int c=0;c<v.size();c++){
            cell += A.at(f).at(c)*v.at(c);
        }
        R.at(f) += cell;
    }
}

void productRealMatrix(float real,Matrix M,Matrix &R){
    zeroes(R,M.size());
    for(int i=0;i<M.size();i++)
        for(int j=0;j<M.at(0).size();j++)
            R.at(i).at(j) = real*M.at(i).at(j);
}

void getMinor(Matrix &M,int i, int j){
    //cout << "Calculando menor ("<<i+1<<","<<j+1<<")...\n";
    M.erase(M.begin()+i);
    for(int i=0;i<M.size();i++)
        M.at(i).erase(M.at(i).begin()+j);
}

float determinant(Matrix M){
    if(M.size() == 1) return M.at(0).at(0);
    else{
        float det=0.0;
        for(int i=0;i<M.at(0).size();i++){
            Matrix minor;
            copyMatrix(M,minor);
            getMinor(minor,0,i);
            det += pow(-1,i)*M.at(0).at(i)*determinant(minor);
        }
        return det;
    }
}

void cofactors(Matrix M, Matrix &Cof){
    zeroes(Cof,M.size());
    for(int i=0;i<M.size();i++){
        for(int j=0;j<M.at(0).size();j++){
            //cout << "Calculando cofactor ("<<i+1<<","<<j+1<<")...\n";
            Matrix minor;
            copyMatrix(M,minor);
            getMinor(minor,i,j);
            Cof.at(i).at(j) = pow(-1,i+j)*determinant(minor);
        }
    }
}

void transpose(Matrix M, Matrix &T){
    zeroes(T,M.at(0).size(),M.size());
    for(int i=0;i<M.size();i++)
        for(int j=0;j<M.at(0).size();j++)
            T.at(j).at(i) = M.at(i).at(j);
}

void inverseMatrix(Matrix M, Matrix &Minv){
    cout << "Iniciando calculo de inversa...\n";
    Matrix Cof, Adj;
    cout << "Calculo de determinante...\n";
    float det = determinant(M);
    if(det == 0) exit(EXIT_FAILURE);
    cout << "Iniciando calculo de cofactores...\n";
    cofactors(M,Cof);
    cout << "Calculo de adjunta...\n";
    transpose(Cof,Adj);
    cout << "Calculo de inversa...\n";
    productRealMatrix(1/det,Adj,Minv);
}
