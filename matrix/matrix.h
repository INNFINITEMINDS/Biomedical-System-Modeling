#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
using namespace std;

template<class T> class matrix{
    vector < vector <T> > data;
    int rows, cols;
    string majorType;
public:
    matrix (int rols, int getCols, string majorType);
    int getRows(){return rows;}
    int getCols(){return cols;}
    void setRows(int rows){rows = rows;}
    void setCols(int cols){cols = cols;}
    //set all the elements in matrix to the same number
    void setValue(double value);
    void addRow();
    void addCol();
    matrix<T> transpose();
    T dot(matrix<T>);
    T sum();
    // scalar operation
    void operator += (T);
    void operator -= (T);
    void operator /= (T);
    void operator *= (T);
    // indexing
    vector<T>& operator[](size_t i){return data[i];}
    T& operator()(size_t i, size_t j){return data[j][i];}
    // matrix operation
    matrix<T> operator + (matrix<T>&);
    matrix<T> operator - (matrix<T>&);
    matrix<T> operator * (matrix<T>&);
    // print the matrix dimensions, format, and contents
    friend ostream& operator << (ostream& out, matrix<T>& mat){
        ostringstream output;
        output << "dimension: " << mat.getRows() << " x " << mat.getCols() << endl;
        output << "format: column major" << endl;
        output << setprecision(4);
        for (int i = 0; i < mat.getRows(); i++){
            for (int j = 0; j < mat.getCols(); j++)
            output << mat.data[j][i] << "\t";
	    output << "\r" << endl;
        }
        out << output.str();
        return out;
    }
    matrix<T> forwardSub(matrix<T>&);
    matrix<T> backwardSub(matrix<T>&);
    matrix<T> pivot(int);
    void LU(matrix<T>&,matrix<T>&,matrix<T>&);
    T det();
};

//constructor
template<class T> matrix<T>::matrix (int rows, int cols, string majorType):rows(rows),cols(cols){
    if (majorType == "rows"){
        int temp = cols;
        cols = rows;
        rows = temp;
    }
    data.resize( cols );
    for( auto &q : data ) q.resize( rows );
}

//set all the elements in matrix to the same number
template<class T> void matrix<T>::setValue (double value){
    for (int i = 0; i < (*this).getRows(); i++)
    {
        for ( int j = 0; j < (*this).getCols(); j++ )
            (*this)[j][i] = value;
    }
}

// add a row or column to the end of matrix
template<class T> void matrix<T>::addRow(){
    for (size_t i = 0; i < this->data.size(); i++){
        this->data[i].resize(this->getRows() + 1);
    }
    this->setRows( this->getRows()+1 );
}
template<class T> void matrix<T>::addCol(){
    this->data.resize(this->getCols() + 1);
    this->data[this->getCols()].resize(this->getRows());
    this->setCols( this->getCols()+1 );
}

// scalar operations (+=, -=, /=, *=)
template<class T> void matrix<T>::operator += (T addend){
    for (int i = 0; i < (*this).getRows(); i++)
    {
        for ( int j = 0; j < (*this).getCols(); j++ )
            (*this)[j][i] += addend;
    }
}
template<class T> void matrix<T>::operator -= (T addend){
    for (int i = 0; i < (*this).getRows(); i++)
    {
        for ( int j = 0; j < (*this).getCols(); j++ )
            (*this)[j][i] -= addend;
    }
}
template<class T> void matrix<T>::operator /= (T addend){
    for (int i = 0; i < (*this).getRows(); i++)
    {
        for ( int j = 0; j < (*this).getCols(); j++ )
            (*this)[j][i] /= addend;
    }
}
template<class T> void matrix<T>::operator *= (T addend){
    for (int i = 0; i < (*this).getRows(); i++)
    {
        for ( int j = 0; j < (*this).getCols(); j++ )
            (*this)[j][i] *= addend;
    }
}

// matrix operations (+, -, *)
template<class T> matrix<T> matrix<T>::operator + (matrix<T> &M){
    matrix<T> new_matrix(this.getRows(),this.getCols(),"cols");
    for (int i = 0; i < new_matrix.getRows(); i++)
    {
        for ( int j = 0; j < new_matrix.getCols(); j++ )
            new_matrix[j][i] = (*this)[j][i] + M[j][i];
    }
    return new_matrix;
}
template<class T> matrix<T> matrix<T>::operator - (matrix<T> &M){
    matrix<T> new_matrix(this.getRows(),this.getCols(),"cols");
    for (int i = 0; i < new_matrix.getRows(); i++)
    {
        for ( int j = 0; j < new_matrix.getCols(); j++ )
            new_matrix[j][i] = (*this)[j][i] - M[j][i];
    }
    return new_matrix;
}
template<class T> matrix<T> matrix<T>::operator * (matrix<T> &M){
    matrix<T> new_matrix(this->getRows(),M.getCols(),"cols");
    for (int i = 0; i < new_matrix.getRows(); i++)
    {
        for ( int j = 0; j < new_matrix.getCols(); j++ ){
            T temp = 0;
            for( int k = 0; k < this->getCols(); k++ )
                temp += (*this)(i,k)*M(k,j);
            new_matrix(i,j) = temp;
        }
    }
    return new_matrix;
}

//matrix transpose
template<class T> matrix<T> matrix<T>::transpose (){
    matrix<T> matrix_trp((*this).getRows(),(*this).getCols(),"cols");
    for (int i = 0; i < matrix_trp.getRows(); i++)
    {
        for ( int j = 0; j < matrix_trp.getCols(); j++ )
            matrix_trp[j][i] = (*this)[i][j];
    }
    return matrix_trp;
}

//dot products
template<class T> T matrix<T>::dot (matrix<T> M){
    vector<T> M_line, data_line;
    double dproduct;
    for (int i = 0; i < M.getRows(); i++){
        for (int j = 0; j < M.getCols(); j++)
            M_line.push_back(M[j][i]);
    }
    for (int i = 0; i < (*this).getRows(); i++){
        for (int j = 0; j < (*this).getCols(); j++)
            data_line.push_back((*this)[j][i]);
    }
    for (int i = 0; i < data_line.size(); i++)
        dproduct += data_line[i] * M_line[i];
    return dproduct;
}

//additional functionality: summation of a matrix's elements
template<class T> T matrix<T>::sum(){
    T s;
    for (int i = 0; i < (*this).getRows(); i++){
        for (int j = 0; j < (*this).getCols(); j++)
            s += (*this)[j][i];
    }
    return s;
}

// forward substitution
template<class T> matrix<T> matrix<T>::forwardSub(matrix<T>& B){
    if( this->getRows() != B.getRows() ) throw "Dimension doesn't match.";
    matrix<T> solution(this->getRows(), 1, "");
    for( int i = 0; i < this->getRows(); ++i ){
        T A_const = 0;
        for( int j = 0; j < i; ++j )
            A_const += (*this)(i,j) * solution(j,0);
        solution(i,0) = (B(i,0) - A_const) / (*this)(i,i);
    }
    return solution;
}

// backward substitution
template<class T> matrix<T> matrix<T>::backwardSub(matrix<T>& B){
    if( this->getRows() != B.getRows() ) throw "Dimension doesn't match.";
    matrix<T> solution(this->getRows(), 1, "");
    for( int i = this->getRows()-1; i > -1; --i ){
        T A_const = 0;
        for( int j = this->getRows()-1; j > i; --j )
            A_const += (*this)(i,j) * solution(j,0);
        solution(i,0) = (B(i,0) - A_const) / (*this)(i,i);
    }
    return solution;
}

// pivoting
template<class T> matrix<T> matrix<T>::pivot( int pos ){
    // find the row with maximum pivot
    int max_row;
    T temp = 0;
    for( int i = pos; i < this->getRows(); ++i ){
        if( abs((*this)(i,pos)) > temp ){
            temp = abs((*this)(i,pos));
            max_row = i;
        }
    }

    // construct P(permutation matrix)
    matrix<T> P(this->getRows(),this->getRows(),"");
    for( int i = 0; i < this->getRows(); ++i )
        P(i,i) = 1;
    P(max_row,max_row) = 0; P(max_row,pos) = 1;
    P(pos,pos) = 0; P(pos,max_row) = 1;

    return P;
}

// LU decomposition
template<class T> void matrix<T>::LU(matrix<T>& L, matrix<T>& U, matrix<T>& P){
    // initialize
    for( int i = 0; i < L.getCols(); ++i )
        L(i,i) = 1;
    U = (*this) * L;
    P = L * L;

    for( int j = 0; j < this->getCols()-1; ++j ){
        // pivoting(change rows)
        matrix<T> P_temp(this->getRows(),this->getCols(),"");
        P_temp = U.pivot(j);
        U = P_temp * U;
        P = P_temp * P;

        // swap rows of L on previous columns if necessary
        if( j > 0 ){
            vector<T> swap_row;
            for( int k = j; k < P_temp.getCols(); ++k ){
                if( P_temp(k,k) != 1 )
                    swap_row.push_back(k);
            }
            if( !swap_row.empty() ){
                for( int k = 0; k < j; ++k ){
                    T temp = L(swap_row[0],k);
                    L(swap_row[0],k) = L(swap_row[1],k);
                    L(swap_row[1],k) = temp;
                }
            }
        }

        // elimination
        matrix<T> L_temp(this->getRows(),this->getCols(),"");
        for( int k = 0; k < L_temp.getCols(); ++k )
            L_temp(k,k) = 1;
        for( int i = j+1; i < L_temp.getRows(); ++i ){
            L_temp(i,j) = -U(i,j)/U(j,j);
            L(i,j) = -L_temp(i,j);
        }
        U = L_temp * U;
        // clear the error (small number from division)
        for( int i = j+1; i < U.getRows(); ++i )
            U(i,j) = 0;
    }
}

// determinant of the matrix
template<class T> T matrix<T>::det(){
    matrix<T> L(this->getRows(),this->getCols(),"");
    matrix<T> U(this->getRows(),this->getCols(),"");
    matrix<T> P(this->getRows(),this->getCols(),"");
    this->LU( L, U, P );

    int num_exchanges = 0;
    for( int i = 0; i < P.getCols(); ++i ){
        if( P(i,i) == 0 )
            num_exchanges++;
    }

    T det_U = 1;
    for( int i = 0; i < U.getCols(); ++i )
        det_U *= U(i,i);

    return pow(-1, num_exchanges)*det_U;
}
#endif
