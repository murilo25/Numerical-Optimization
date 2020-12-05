// Murilo Augusto Pinheiro
// 04-Dec-2020
// Solution for Diamond Kinetics - Dynamics Algorithm Implementation.pdf
// SteepestDescent.cpp : This file contains the 'main' function. Program execution begins and ends there.
// Header.h : This files contains all functions and includes used in this file.



#include "Header.h"
//#include "math_functions.h"

int main()
{
    std::vector<double> result;
    std::vector<double> ic;
    ic.push_back(1.2);
    ic.push_back(1.2);

    //result = optmization(9000, 0.00001, ic, STEEPEST_DESCENT);
    result = optmization(9000, 0.00001, ic, NEWTON);

    std::cout << "x1* = " << result[0] << "\nx2* = " << result[1] << "\ncost = " << result[2] << std::endl;

    return 0;
}


/*
{
    std::vector<double> line(DIMENSION + 1, 0);
    std::vector< std::vector<double> > A(DIMENSION, line);

    A[0][0] = 2;
    A[0][1] = 1;
    A[1][0] = 1;
    A[1][1] = 1;

    A[0][2] = 1;
    A[1][2] = 3;

    // Print input
    print(A);

    // Calculate solution
    std::vector<double> x(DIMENSION);
    x = gauss(A);

    // Print result
    std::cout << "Result:\t";
    for (int i = 0; i < DIMENSION; i++) {
        std::cout << x[i] << " ";
    }
    std::cout << std::endl;
}
*/