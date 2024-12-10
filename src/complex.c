#include <math.h>
#include "../include/complex.h"

const Cart CART_UNITY = {.real = 1, .imag = 0};
const Cart CART_ZERO = {.real = 0, .imag = 0};

const Polar POLAR_UNITY = {.mag = 1, .theta = 0};
const Polar POLAR_ZERO = {.mag = 0, .theta = 0};

Cart pol2cart(Polar in) {
    Cart out;
    out.real = in.mag * cos(in.theta);
    out.imag = in.mag * sin(in.theta);
    return out;
}

Polar cart2pol(Cart in) {
    Polar out;
    out.mag = sqrt(pow(in.real, 2) + pow(in.imag, 2));
    out.theta = atan2(in.imag, in.real);
    return out;
}

Cart add(Cart a, Cart b) {
    Cart out;
    out.real = a.real + b.real;
    out.imag = a.imag + b.imag;
    return out;
}

Cart subtract(Cart a, Cart b) {
    Cart out;
    out.real = a.real - b.real;
    out.imag = a.imag - b.imag;
    return out;
}

Cart multiply(Cart a, Cart b) {
    Polar out, a_polar, b_polar;
    a_polar = cart2pol(a); b_polar = cart2pol(b);
    out.mag = a_polar.mag * b_polar.mag;
    out.theta = a_polar.theta + b_polar.theta;
    return pol2cart(out);
}

Cart divide(Cart a, Cart b) {
    Polar out, a_polar, b_polar;
    a_polar = cart2pol(a); b_polar = cart2pol(b);
    out.mag = a_polar.mag / b_polar.mag;
    out.theta = a_polar.theta - b_polar.theta;
    return pol2cart(out);
}

Cart conjugate(Cart in) {
    Cart out;
    out.real = in.real;
    out.imag = in.imag * -1;
    return out;
}
