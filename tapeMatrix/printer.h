#pragma once
#include <stdio.h>

void print_1d(double* array, int n)
{
    for (int i = 0; i < n; i++)
    {
        printf("%lf ", array[i]);
    }
}

void print_2d(double** array, int n)
{
    for (int i = 0; i < n; i++) 
    {
        for (int j = 0; j < n; j++) 
        {
            printf("%lf ", array[i][j]);
        }
        printf("\n");
    }
}