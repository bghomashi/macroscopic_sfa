#pragma once

#include "maths.h"

dvector cos2_window(int N) {
    dvector w(N, 1.0);
    
    double A = Pi_2 / 100.0;
    int b = N - 101;

    for (int i = N - 101; i < N; i++) {
        double x = A * (i - b);
        w[i] = cos(x) * cos(x);
    }

    return w;
}