#pragma once

struct DecomposeMatrix
{
    double** l;
    double** u;
};

struct Matrix
{
    long n;
    long b;
    double* C;
    double** A;
    double* X;
};