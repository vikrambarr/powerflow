#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "../include/powerflow.h"

Polar *create_loadflow_matrix(Node *nodes, size_t node_count) {
    Polar *loadflow = malloc(node_count * sizeof(Polar));
    for (int node = 0; node < node_count; node++) loadflow[node] = nodes[node].phasor;
    return loadflow;
}

Polar *create_admittance_matrix(Node *nodes, size_t node_count, Line *lines, size_t line_count) {
    Cart cart_admittance_matrix[node_count * node_count];
    memset(cart_admittance_matrix, 0, node_count * node_count * sizeof(Cart));

    for (int line = 0; line < line_count; line++) {
        int from = lines[line].from; int to = lines[line].to;
        Cart admittance = divide(CART_UNITY, lines[line].impedance);
        Cart self_admittance = {.real = admittance.real, .imag = admittance.imag + lines[line].b_half};

        // OFF-DIAGONAL
        cart_admittance_matrix[from * node_count + to] = subtract(cart_admittance_matrix[from * node_count + to], admittance);
        cart_admittance_matrix[to * node_count + from] = subtract(cart_admittance_matrix[to * node_count + from], admittance);

        // DIAGONAL
        cart_admittance_matrix[from * node_count + from] = add(cart_admittance_matrix[from * node_count + from], self_admittance);
        cart_admittance_matrix[to   * node_count + to  ] = add(cart_admittance_matrix[to * node_count + to], self_admittance);
    }

    Polar *admittance_matrix = malloc(node_count * node_count * sizeof(Polar));
    for (int y = 0; y < node_count * node_count; y++) admittance_matrix[y] = cart2pol(cart_admittance_matrix[y]);

    return admittance_matrix;
}

void populate_mismatch_matrix(double *mismatch_matrix, Polar *admittance_matrix, Polar *loadflow, Node *nodes, size_t node_count) {
    Cart node_power_balance;

    int real_eq_count = 0; int imag_eq_count = 0;
    double real_power_mismatch[node_count];
    double imag_power_mismatch[node_count];

    for (int r = 0; r < node_count; r++) {
        node_power_balance = CART_ZERO;
        for (int c = 0; c < node_count; c++) {
            node_power_balance.real += loadflow[r].mag * loadflow[c].mag * admittance_matrix[r * node_count + c].mag * cos(loadflow[r].theta - loadflow[c].theta - admittance_matrix[r * node_count + c].theta);
            node_power_balance.imag += loadflow[r].mag * loadflow[c].mag * admittance_matrix[r * node_count + c].mag * sin(loadflow[r].theta - loadflow[c].theta - admittance_matrix[r * node_count + c].theta);
        }

        if (nodes[r].type == 1) continue;
        real_power_mismatch[real_eq_count] = nodes[r].injected.real / 100 - node_power_balance.real;
        real_eq_count++;

        if (nodes[r].type == 2) continue;
        imag_power_mismatch[imag_eq_count] = nodes[r].injected.imag / 100 - node_power_balance.imag;
        imag_eq_count++;
    }

    for (int r_eq; r_eq < real_eq_count; r_eq++) {
        mismatch_matrix[r_eq] = real_power_mismatch[r_eq];
    }
    for (int i_eq; i_eq < imag_eq_count; i_eq++) {
        mismatch_matrix[real_eq_count + i_eq] = imag_power_mismatch[i_eq];
    }
}

void populate_jacobian_matrix(double *jacobian_matrix, Polar *admittance_matrix, Polar *loadflow, Node *nodes, size_t node_count) {
    int to_node_index[node_count];
    int real_eq_count = 0; int imag_eq_count = 0;

    for (int node = 0; node < node_count * 2; node++) {
        if (node < node_count) {
            if (nodes[node].type == 1) continue;
            to_node_index[real_eq_count] = node;
            real_eq_count++;
        } else {
            if (nodes[node - node_count].type == 1 || nodes[node - node_count].type == 2) continue;
            to_node_index[real_eq_count + imag_eq_count] = node - node_count;
            imag_eq_count++;
        }
    }

    for (int r = 0; r < real_eq_count; r++) {
        for (int c = 0; c < real_eq_count; c++) { // 1ST QUADRANT
            int k = r; int n = c;
            jacobian_matrix[(node_count + 1) * k + n] = loadflow[to_node_index[k]].mag * loadflow[to_node_index[n]].mag * admittance_matrix[(node_count + 1) * to_node_index[k] + to_node_index[n]].mag * sin(loadflow[to_node_index[k]].theta - loadflow[to_node_index[n]].theta - admittance_matrix[(node_count + 1) * to_node_index[k] + to_node_index[n]].theta);
            if (r != c) jacobian_matrix[(node_count + 1) * k + k] -= jacobian_matrix[(node_count + 1) * k + n];
        }
        for (int c = 0; c < imag_eq_count; c++) { // 2ND QUADRANT
            int k = r; int n = c + real_eq_count;
            jacobian_matrix[(node_count + 1) * k + n] = loadflow[to_node_index[k]].mag * admittance_matrix[(node_count + 1) * to_node_index[k] + to_node_index[n]].mag * cos(loadflow[to_node_index[k]].theta - loadflow[to_node_index[n]].theta - admittance_matrix[(node_count + 1) * to_node_index[k] + to_node_index[n]].theta);
            if (r != c) jacobian_matrix[(node_count + 1) * k + k] += jacobian_matrix[(node_count + 1) * k + n];
        }
    }

    for (int r = 0; r < imag_eq_count; r++) {
        for (int c = 0; c < real_eq_count; c++) { // 3RD QUADRANT
            int k = r + real_eq_count; int n = c;
            jacobian_matrix[(node_count + 1) * k + n] = -loadflow[to_node_index[k]].mag * loadflow[to_node_index[n]].mag * admittance_matrix[(node_count + 1) * to_node_index[k] + to_node_index[n]].mag * cos(loadflow[to_node_index[k]].theta - loadflow[to_node_index[n]].theta - admittance_matrix[(node_count + 1) * to_node_index[k] + to_node_index[n]].theta);
            if (r != c) jacobian_matrix[(node_count + 1) * k + k] -= jacobian_matrix[(node_count + 1) * k + n];
        }
        for (int c = 0; c < imag_eq_count; c++) { // 4TH QUADRANT
            int k = r + real_eq_count; int n = c + real_eq_count;
            jacobian_matrix[(node_count + 1) * k + n] = loadflow[to_node_index[k]].mag * admittance_matrix[(node_count + 1) * to_node_index[k] + to_node_index[n]].mag * sin(loadflow[to_node_index[k]].theta - loadflow[to_node_index[n]].theta - admittance_matrix[(node_count + 1) * to_node_index[k] + to_node_index[n]].theta);
            if (r != c) jacobian_matrix[(node_count + 1) * k + k] += jacobian_matrix[(node_count + 1) * k + n];
        }
    }
}
