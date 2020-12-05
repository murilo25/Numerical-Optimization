// Murilo Augusto Pinheiro
// 04-Dec-2020
// Solution for Diamond Kinetics - Dynamics Algorithm Implementation.pdf
// SteepestDescent.cpp : This file contains the 'main' function. Program execution begins and ends there.
// Header.h : This files contains all functions and includes used in this file.



#include "Header.h"

int main()
{
    std::vector<double> result;
    std::vector<double> ic;
    ic.push_back(10);
    ic.push_back(-10);

    result = optmization(100, 0.00001, ic, STEEPEST_DESCENT);

    std::cout << "x1* = " << result[0] << "\nx2* = " << result[1] << "\ncost = " << result[2] << std::endl;

    return 0;
}