#ifndef POWERFLOW_H
#define POWERFLOW_H

#include "complex.h"

typedef struct {
    int num;
    int type;
    Polar phasor;
    Polar loadflow;
    Cart injected;
    double minq;
    double maxq;
} Node;

typedef struct {
    int from;
    int to;
    Cart z;
    double b_half;
} Line;

#endif
