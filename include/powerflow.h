#ifndef POWERFLOW_H
#define POWERFLOW_H

#include "complex.h"

typedef struct Node {
    int type; // 1: SLACK, 2: GENERATOR, 3: LOAD
    Polar voltage;
    Cart injected;
    double q_min;
    double q_max;
} Node;

typedef struct Line {
    int from;
    int to;
    Cart impedance;
    double b_half;
} Line;

#endif
