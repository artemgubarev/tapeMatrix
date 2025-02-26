#pragma once

struct DecomposeMatrix
{
    double** l;
    double** u;
};

struct Matrix
{
    int32_t n;
    int32_t b;
    double* C;
    double** A;
    double* X;
};