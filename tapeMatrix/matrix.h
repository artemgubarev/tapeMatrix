#pragma once

struct DecomposeMatrix
{
    double** l;
    double** u;
};

struct Matrix
{
    std::uint32_t n;
    std::uint32_t b;
    double* C;
    double** A;
    double* X;
};