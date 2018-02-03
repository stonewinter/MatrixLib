#include <iostream>
#include <vector>
#include "matrix.h"
#include <typeinfo>

using namespace std;

int main()
{
    //int a=500, b=500;

    vector<double> myVec;
    Mat<double> mat(5,5);

//    for(int i=0;i<a;i++){
//        for(int j=0;j<b;j++){
//            myVec.push_back(random_gen(0,999));
//        }
//        mat.FillByRow(i, myVec);
//    }
//    mat.ShowContent();
//
//    Mat<int> mat2(a,b);
//    for(int i=0;i<a;i++){
//        for(int j=0;j<b;j++){
//            myVec.push_back(random_gen(0,999));
//        }
//        mat2.FillByRow(i, myVec);
//    }
//    mat2.ShowContent();


//    myVec = {5,7,6,4,5};
//    mat.FillByRow(0, myVec);
//    myVec = {9,7,-10,-11,4};
//    mat.FillByRow(1, myVec);
//    myVec = {-3,2.5,43.7,43,8};
//    mat.FillByRow(2, myVec);
//    myVec = {-2,9.1,8,2,1};
//    mat.FillByRow(3, myVec);
//    myVec = {8,6,4,3,6};
//    mat.FillByRow(4, myVec);
//    myVec = {9,8,1,7,27,56};
//    mat.FillByRow(5, myVec);

//    mat.ShowContent();

//    int res = mat.Determinant();
//    cout<<"res = "<< res <<endl;


    vector<double> myVec2;
    Mat<double> mat2(3,3);
    myVec2 = {-11,10,3};
    mat2.FillByRow(0, myVec2);
    myVec2 = {4,-3,1};
    mat2.FillByRow(1, myVec2);
    myVec2 = {3,-1,1};
    mat2.FillByRow(2, myVec2);
    mat2.ShowContent();

    Mat<double> mat3;
    mat3 = mat2.GetReciprocalMat();
    mat3.ShowContent();

    //(mat2*mat3).ShowContent();

    if(mat3*mat2 == E(3)){
        cout<< "yyyyyyyy" << endl;
    }else{
        cout << "nnnnnnn" << endl;
    }

    return 0;
}
