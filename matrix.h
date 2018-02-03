#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <iostream>
#include <cstring>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <iomanip>
#include "stdio.h"

#define CellingDiv(a,b) (round((double)a/(double)b))
#define random_gen(min,max) (min+rand()%(max-min+1))
#define E(n) GenUnitaryMat(n)


using namespace std;

template<typename T> class Mat;
template<typename T1, typename T2> bool operator!=(Mat<T1> matA, Mat<T2> matB);
template<typename T1, typename T2> bool operator==(Mat<T1> matA, Mat<T2> matB);
template<typename T1, typename T2> Mat<double> operator*(T1 coff, Mat<T2> mat);
template<typename T1, typename T2> Mat<double> operator*(Mat<T1> mat, T2 coff);
Mat<double> GenUnitaryMat(int order);




template<typename T>
class Mat{
private:
    int rowNbr;
    int colNbr;
    vector<vector<T> > matrix;

public:
    Mat(){
        ;
    }

    Mat(int rowNbr, int colNbr){
        Build(rowNbr, colNbr);
    }

    int GetRowNbr(){
        return this->rowNbr;
    }

    int GetColNbr(){
        return this->colNbr;
    }

    void Resize(int rowNbr, int colNbr){
        Build(rowNbr, colNbr);
    }

    void Build(int rowNbr, int colNbr){
        this->rowNbr = rowNbr;
        this->colNbr = colNbr;
        matrix.resize(this->rowNbr);
        for(int i=0;i<rowNbr;i++){ //row scan
            matrix[i].resize(this->colNbr);
        }
    }

    void FillByRow(int rowIndex, vector<T> rowVec){
        matrix[rowIndex] = rowVec;
    }



    /* 求某行某列那个元素的代数余子式值
    */
    T GetCofactor(int rowIndex, int colIndex){
        T res = pow(-1,rowIndex+colIndex)*((this->GetSubDetMat(rowIndex,colIndex)).Determinant());
        return res;
    }


    /* calculate Determinant of a N*N matrix using function iteration method
    */
    double Determinant(){
        if(rowNbr != colNbr){
            return 0;
        }

        if(rowNbr == 2){
            return (matrix[0][0]*matrix[1][1] - matrix[0][1]*matrix[1][0]);
        }else{
            double sum = 0;
            for(int i=0;i<colNbr;i++){
                sum += pow(-1,i)*matrix[0][i]*(this->GetSubDetMat(0,i)).Determinant();
            }
            return sum;
        }
    }


    /*Get submatrix eliminating row&col according to rowIndex&colIndex
    */
    Mat<T> GetSubDetMat(int rowIndex, int colIndex){
        if(rowNbr != colNbr || rowNbr<=2 || colNbr <=2){
            return (*this);
        }

        Mat<T> subMat(rowNbr-1, colNbr-1);
        int subRowIndex, subColIndex;

        if(rowIndex == 0){
            subRowIndex = 1; // start from Row No.1
        }
        else{
            subRowIndex = 0; // start from Row No.0
        }

        if(colNbr == 0){
            subColIndex = 1;
        }else{
            subColIndex = 0;
        }

        for(int i=0;i<rowNbr-1;i++){
            for(int j=0;j<colNbr-1;j++){
                if(subRowIndex == rowIndex){
                    subRowIndex++;
                }
                if(subColIndex == colIndex){
                    subColIndex++;
                }
                subMat[i][j] = matrix[subRowIndex][subColIndex];
                subColIndex++;
            }
            subRowIndex++;
            if(colNbr == 0){
                subColIndex = 1;
            }else{
                subColIndex = 0;
            }
        }
        return subMat;
    }

    /* 求伴随矩阵
    */
    Mat<T> GetAdjointMat(){
        Mat<T> tempMat(rowNbr,colNbr);
        if(rowNbr > 2){
            for(int i=0;i<rowNbr;i++){
                for(int j=0;j<colNbr;j++){
                    tempMat[i][j] = this->GetCofactor(i,j);
                }
            }
        }else{
            tempMat[0][0] = matrix[1][1];
            tempMat[0][1] = -matrix[0][1];
            tempMat[1][0] = -matrix[1][0];
            tempMat[1][1] = matrix[0][0];
        }
        tempMat.Tranpose();
        return tempMat;
    }

    /* 求逆矩阵
    */
    Mat<double> GetReciprocalMat(){
        Mat<double> tempMat(rowNbr, colNbr);
        tempMat = this->GetDoubleFormat();
        double coff = 1/this->Determinant();
        tempMat = coff * (tempMat.GetAdjointMat());

        return tempMat;
    }

    void ShowContent(int precision=4){
        cout<<endl;
        for(int i=0;i<rowNbr;i++){
            for(int j=0;j<colNbr;j++){
                //matrix[i][j] = i+j;
                cout <<  setprecision(precision) << matrix[i][j] << "\t";
            }
            cout<<endl;
        }
    }

    /* 等号重载, 麻丹的, operator[]必须是成员函数
    */
    vector<T>& operator[](int const index){
        return matrix[index];
    }

    Mat<T> operator+(Mat<T> matB){
        Mat tempMat(rowNbr, colNbr);
        for(int i=0;i<rowNbr;i++){
            for(int j=0;j<colNbr;j++){
                tempMat[i][j] = matrix[i][j] + matB[i][j];
            }
        }
        return tempMat;
    }

    Mat<T> operator-(Mat<T> matB){
        Mat tempMat(rowNbr, colNbr);
        for(int i=0;i<rowNbr;i++){
            for(int j=0;j<colNbr;j++){
                tempMat[i][j] = matrix[i][j] - matB[i][j];
            }
        }
        return tempMat;
    }

    /* 等号重载, 麻丹的, operator=必须是成员函数
    */
    Mat<T>& operator=(Mat<T> matB){
        int rowNbr = matB.GetRowNbr();
        int colNbr = matB.GetColNbr();
        this->Resize(rowNbr, colNbr);
        for(int i=0;i<rowNbr;i++){ //matC row
            for(int j=0;j<colNbr;j++){ // matC col
                (*this)[i][j] = matB[i][j];
            }
        }

        return *this;
    }


    Mat<double> GetDoubleFormat(){
        Mat<double> tempMat(rowNbr, colNbr);
        for(int i=0;i<rowNbr;i++){
            for(int j=0; j<colNbr; j++){
                tempMat[j][i] = (double)matrix[j][i];
            }
        }
        return tempMat;
    }


    Mat<T>& Tranpose(){
        Mat tempMat(colNbr, rowNbr);

        for(int i=0;i<rowNbr;i++){
            for(int j=0; j<colNbr; j++){
                tempMat[j][i] = matrix[i][j];
            }
        }

        this->Resize(colNbr, rowNbr);
        *this = tempMat;

        return *this;
    }



};









/* 乘号重载，数字相乘，数字后置
*/
template<typename T1, typename T2>
Mat<double> operator*(Mat<T1> mat, T2 coff){
    int rowNbr = mat.GetRowNbr();
    int colNbr = mat.GetColNbr();
    Mat<double> tempMat(rowNbr, colNbr);
    for(int i=0;i<rowNbr;i++){ //matC row
        for(int j=0;j<colNbr;j++){ // matC col
            tempMat[i][j] = (double)(mat[i][j]*coff);
        }
    }
    return tempMat;
}



/* 乘号重载，数字相乘，数字前置
*/
template<typename T1, typename T2>
Mat<double> operator*(T1 coff, Mat<T2> mat){
    return mat*coff;
}



/* 乘号重载，矩阵相乘
*/
template<typename T1, typename T2>
Mat<double> operator*(Mat<T1> matA, Mat<T2> matB){
    int rowNbr = matA.GetRowNbr();
    int colNbr = matA.GetColNbr();
    Mat<double> tempMat(rowNbr, colNbr);
    for(int i=0;i<rowNbr;i++){ //matC row
        for(int j=0;j<colNbr;j++){ // matC col
            for(int u=0;u<colNbr;u++){ //matA col index
                tempMat[i][j] += matA[i][u] * matB[u][j];
            }
        }
    }
    return tempMat;
}
/****************************************************/
/** 另一种mat*mat算法
    v = zeros(10, 1);
    for i = 1:10
      for j = 1:10
        v(i) = v(i) + A(i, j) * x(j);
      end
    end
**/
/***************************************************/





Mat<double> GenUnitaryMat(int order){
    Mat<double> tempMat(order, order);
    for(int i=0;i<order;i++){
        tempMat[i][i] = 1;
    }
    return tempMat;
}


template<typename T1, typename T2>
bool operator==(Mat<T1> matA, Mat<T2> matB){
    if( matA.GetRowNbr() != matB.GetRowNbr()
     || matA.GetColNbr() != matB.GetColNbr()){
        return false;
    }

    for(int i=0;i<matA.GetRowNbr();i++){
        for(int j=0;j<matA.GetColNbr();j++){
            if( abs((double)matA[i][j]-(double)matB[i][j]) >= 0.000001 ){
                return false;
            }
        }
    }
    return true;
}


template<typename T1, typename T2>
bool operator!=(Mat<T1> matA, Mat<T2> matB){
    if( matA.GetRowNbr() != matB.GetRowNbr()
     || matA.GetColNbr() != matB.GetColNbr()){
        return true;
    }

    for(int i=0;i<matA.GetRowNbr();i++){
        for(int j=0;j<matA.GetColNbr();j++){
            if(  abs((double)matA[i][j]-(double)matB[i][j]) >= 0.000001  ){
                return true;
            }
        }
    }
    return false;
}

#endif // _MATRIX_H_
