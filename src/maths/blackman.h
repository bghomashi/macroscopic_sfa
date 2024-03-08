#pragma once

#include "maths.h"

dvector blackman(int N) {
    dvector w(N, 1.0);
    int M = (N % 2 == 0 ? N / 2 : (N + 1) / 2);
    // for (int n = 0; n <= M - 1; n++)
    //     w[n] = 0.42 - 0.5 * cos(2 * Pi * n / (N - 1)) + 0.08 * cos(4 * Pi * n / (N - 1));
    for (int n = M; n <= N - 1; n++)
        w[n] = 0.42 - 0.5 * cos(2 * Pi * (N - 1 - n) / (N - 1)) + 0.08 * cos(4 * Pi * (N - 1 - n) / (N - 1));
    return w;
}