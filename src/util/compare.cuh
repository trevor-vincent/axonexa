#pragma once

__host__ __device__ bool doub_equal(double a, double b) {
    if (fabs(a - b) < EPSILON) {
        return true;
    } else {
        return false;
    }
}

__host__ __device__ bool real_equal(real a, real b, real eps) {
    if (fabs(a - b) < eps) {
        return true;
    } else {
        return false;
    }
}
