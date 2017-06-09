#include<iostream>
#include"matrix.h"
using namespace std;

int main(){
    cout << "---- Test on given 5x5 matrix ----" << endl;
    matrix<double> mat(5,5,"");
    mat(0,0) = 1; mat(0,1) = 1; mat(0,2) = 2; mat(0,3) = 4; mat(0,4) = 8;
    mat(1,0) = 2; mat(1,1) = 2; mat(1,2) = 3; mat(1,3) = 5; mat(1,4) = 7;
    mat(2,0) = 10; mat(2,1) = 10; mat(2,2) = 7; mat(2,3) = 11; mat(2,4) = 13;
    mat(3,0) = 1; mat(3,1) = 5; mat(3,2) = 16; mat(3,3) = 1; mat(3,4) = 5;
    mat(4,0) = 100; mat(4,1) = 1; mat(4,2) = 2; mat(4,3) = 3; mat(4,4) = 4;

    matrix<double> mat_B(5,1,"");
    mat_B(0,0) = 2; mat_B(1,0) = 4; mat_B(2,0) = 6; mat_B(3,0) = 8; mat_B(4,0) = 10;

    matrix<double> L(5,5,"");
    matrix<double> U(5,5,"");
    matrix<double> P(5,5,"");
    mat.LU( L, U, P );
    cout << "Result:" << endl;
    cout << "L:" << endl << L << endl;
    cout << "U:" << endl << U << endl;
    cout << "P:" << endl << P << endl;

    matrix<double> Y(5,5,"");
    matrix<double> product(5,1,"");
    product = P * mat_B;
    Y = L.forwardSub(product);
    matrix<double> X(5,1,"");
    X = U.backwardSub(Y);

    matrix<double> test(5,1,"");
    cout << "multiply back to validate:" << endl;
    test = mat * X;
    cout << test << endl;

    // extra one
    cout << "---- Test on extra 3x3 matrix ----" << endl;
    matrix<double> extra(3,3,"");
    extra(0,0) = 1; extra(0,1) = 2; extra(0,2) = 3;
    extra(1,0) = 4; extra(1,1) = 5; extra(1,2) = 6;
    extra(2,0) = 7; extra(2,1) = 8; extra(2,2) = 9;

    matrix<double> extra_B(3,1,"");
    extra_B(0,0) = 10; extra_B(1,0) = 11; extra_B(2,0) = 12;

    cout << "coefficient matrix A:" << endl << extra << endl;
    cout << "constant matrix B:" << endl << extra_B << endl;

    matrix<double> Le(3,3,"");
    matrix<double> Ue(3,3,"");
    matrix<double> Pe(3,3,"");
    extra.LU( Le, Ue, Pe );
    cout << "Result:" << endl;
    cout << "L:" << endl << Le << endl;
    cout << "U:" << endl << Ue << endl;
    cout << "P:" << endl << Pe << endl;

    matrix<double> Ye(3,3,"");
    matrix<double> product_e(3,1,"");
    product_e = Pe * extra_B;
    Ye = Le.forwardSub(product_e);
    matrix<double> Xe(3,1,"");
    Xe = Ue.backwardSub(Ye);

    matrix<double> test_e(3,1,"");
    cout << "multiply back to validate:" << endl;
    test_e = extra * Xe;
    cout << test_e << endl;

    // determinant
    cout << "---- determinant test on 5x5 matrix ----" << endl;
    cout << "determinant: " << mat.det() << endl;

    return 0;
}
