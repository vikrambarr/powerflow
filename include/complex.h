#ifndef COMPLEX_H
#define COMPLEX_H

typedef struct {
    double real;
    double imag;
} Cart;

typedef struct {
    double mag;
    double theta;
} Polar;

const Cart CART_UNITY;
const Cart CART_ZERO;

const Polar POLAR_UNITY;
const Polar POLAR_ZERO;

Cart pol2cart(Polar in);
Polar cart2pol(Cart in);

Cart add(Cart a, Cart b);
Cart subtract(Cart a, Cart b);
Cart multiply(Cart a, Cart b);
Cart divide(Cart a, Cart b);

Cart conjugate(Cart in);

#endif
