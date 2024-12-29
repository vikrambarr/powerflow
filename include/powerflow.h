#ifndef POWERFLOW_H
#define POWERFLOW_H

#include "complex.h"

typedef struct Node {
    int num;
    int type;
    Polar phasor;
    Polar loadflow;
    Cart injected;
    double minq;
    double maxq;
} Node;

typedef struct Line {
    int from;
    int to;
    Cart impedance;
    double b_half;
} Line;

#endif
