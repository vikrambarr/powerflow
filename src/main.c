#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "../include/complex.h"
#include "../include/powerflow.h"

// DEFINE GRID
const Node nodes[3] = {
    {
        .type = 1,
        .voltage = {.mag = 1.0, .theta = 0.0},
        .injected = {.real = 0.0, .imag = 0.0},
    },
    {
        .type = 2,
        .voltage = {.mag = 1.0, .theta = 0.0},
        .injected = {.real = 0.5, .imag = -0.2},
    },
    {
        .type = 3,
        .voltage = {.mag = 1.0, .theta = 0.0},
        .injected = {.real = -1.7, .imag = -1.7},
    },
};

const Line lines[2] = {
    {
        .from = 0,
        .to = 1,
        .impedance = {.real = 0.1, .imag = 0.2},
        .b_half = 0.02,
    },
    {
        .from = 1,
        .to = 2,
        .impedance = {.real = 0.05, .imag = 0.2},
        .b_half = 0.02,
    },
};

const size_t node_count = 3;
const size_t line_count = 2;
const double tolerance = 1e-10;

int main() {
    // INITIALIZE
    Cart cart_admittance_matrix[node_count * node_count];
    Polar polar_admittance_matrix[node_count * node_count];

    Cart cart_power_flows[node_count * node_count];
    Polar polar_power_flows[node_count * node_count];
    Cart total_power_flows[node_count];

    double mismatch_matrix[node_count];
    double jacobian_matrix[(2 * node_count) * ((2 * node_count) - 1)];

    Polar loadflow[node_count];
    for (int node = 0; node < node_count; node++) loadflow[node] = nodes[node].voltage;

    // ADMITTANCE MATRIX
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
        cart_admittance_matrix[to   * node_count + to  ] = add(cart_admittance_matrix[to   * node_count + to  ], self_admittance);
    }

    for (int y = 0; y < node_count * node_count; y++) polar_admittance_matrix[y] = cart2pol(cart_admittance_matrix[y]);

    // ITERATION
    int max_iters = 100;
    double mismatch = tolerance;
    while (max_iters-- && mismatch >= tolerance) {

        // PRECALCULATE POWER FLOWS
        for (int k = 0; k < node_count; k++) {
            total_power_flows[k] = CART_ZERO;
            for (int n = 0; n < node_count; n++) {
                int node_index = k * node_count + n;

                polar_power_flows[node_index].mag = polar_admittance_matrix[node_index].mag * loadflow[k].mag * loadflow[n].mag;
                polar_power_flows[node_index].theta = loadflow[k].theta - loadflow[n].theta - polar_admittance_matrix[node_index].theta;

                cart_power_flows[node_index] = pol2cart(polar_power_flows[node_index]);
                total_power_flows[k] = add(total_power_flows[k], cart_power_flows[node_index]);
            }
        }

        // MISMATCH MATRIX
        int to_node_index[node_count];
        int real_eq_count = 0; int imag_eq_count = 0;

        for (int node = 0; node < node_count * 2; node++) {
            if (node < node_count) {
                if (nodes[node].type == 1) continue;
                to_node_index[real_eq_count] = node;
                mismatch_matrix[real_eq_count] = nodes[node].injected.real / 100 - total_power_flows[node].real;
                real_eq_count++;
            } else {
                if (nodes[node - node_count].type == 1 || nodes[node - node_count].type == 2) continue;
                to_node_index[real_eq_count + imag_eq_count] = node - node_count;
                mismatch_matrix[real_eq_count + imag_eq_count] = nodes[node - node_count].injected.imag / 100 - total_power_flows[node - node_count].imag;
                imag_eq_count++;
            }
        }

        int total_eq_count = real_eq_count + imag_eq_count;

        // JACOBIAN MATRIX
        for (int r = 0; r < real_eq_count; r++) {
            for (int c = 0; c < real_eq_count; c++) { // 1ST QUADRANT
                int k = r; int n = c;
                jacobian_matrix[(total_eq_count + 1) * k + n] = cart_power_flows[node_count * to_node_index[k] + to_node_index[n]].imag;
                if (r == c) jacobian_matrix[(total_eq_count + 1) * k + n] -= total_power_flows[to_node_index[k]].imag;
            }
            for (int c = 0; c < imag_eq_count; c++) { // 2ND QUADRANT
                int k = r; int n = c + real_eq_count;
                jacobian_matrix[(total_eq_count + 1) * k + n] = cart_power_flows[node_count * to_node_index[k] + to_node_index[n]].real / loadflow[to_node_index[n]].mag;
                if (r == c) jacobian_matrix[(total_eq_count + 1) * k + n] += total_power_flows[to_node_index[k]].real / loadflow[to_node_index[n]].mag;
            }
        }

        for (int r = 0; r < imag_eq_count; r++) {
            for (int c = 0; c < real_eq_count; c++) { // 3RD QUADRANT
                int k = r + real_eq_count; int n = c;
                jacobian_matrix[(total_eq_count + 1) * k + n] = -cart_power_flows[node_count * to_node_index[k] + to_node_index[n]].real;
                if (r == c) jacobian_matrix[(total_eq_count + 1) * k + n] += total_power_flows[to_node_index[k]].real;
            }
            for (int c = 0; c < imag_eq_count; c++) { // 4TH QUADRANT
                int k = r + real_eq_count; int n = c + real_eq_count;
                jacobian_matrix[(total_eq_count + 1) * k + n] = cart_power_flows[node_count * to_node_index[k] + to_node_index[n]].imag / loadflow[to_node_index[n]].mag;
                if (r == c) jacobian_matrix[(total_eq_count + 1) * k + n] += total_power_flows[to_node_index[k]].imag / loadflow[to_node_index[n]].mag;
            }
        }

        for (int solution = 0; solution < total_eq_count; solution++) {
            jacobian_matrix[(total_eq_count + 1) * solution + total_eq_count] = mismatch_matrix[solution];
        }

        // SOLVE LINALG
        for (int pivot = 0; pivot < total_eq_count - 1; pivot++) {
            for (int r = pivot + 1; r < total_eq_count; r++) {
                double ratio = jacobian_matrix[(total_eq_count + 1) * r + pivot] / jacobian_matrix[(total_eq_count + 1) * pivot + pivot];
                for (int c = pivot; c < total_eq_count + 1; c++) {
                    jacobian_matrix[(total_eq_count + 1) * r + c] -= jacobian_matrix[(total_eq_count + 1) * pivot + c] * ratio;
                }
            }
        }

        for (int pivot = total_eq_count - 1; pivot > 0; pivot--) {
            for (int r = pivot - 1; r >= 0; r--) {
                double ratio = jacobian_matrix[(total_eq_count + 1) * r + pivot] / jacobian_matrix[(total_eq_count + 1) * pivot + pivot];
                for (int c = pivot; c < total_eq_count + 1; c++) {
                    jacobian_matrix[(total_eq_count + 1) * r + c] -= jacobian_matrix[(total_eq_count + 1) * pivot + c] * ratio;
                }
            }
        }

        for (int pivot = 0; pivot < total_eq_count; pivot++) {
            jacobian_matrix[(total_eq_count + 1) * pivot + (total_eq_count)] /= jacobian_matrix[(total_eq_count + 1) * pivot + pivot];
        }

        // STEP VALUES & CALCULATE MISMATCH
        for (int eq = 0; eq < total_eq_count; eq++) {
            if (eq < real_eq_count) loadflow[to_node_index[eq]].theta += jacobian_matrix[(total_eq_count + 1) * eq + total_eq_count];
            else loadflow[to_node_index[eq]].mag += jacobian_matrix[(total_eq_count + 1) * eq + total_eq_count];

            mismatch = 0;
            if (fabs(jacobian_matrix[(total_eq_count + 1) * eq + total_eq_count]) > fabs(mismatch)) {
                mismatch = fabs(jacobian_matrix[(total_eq_count + 1) * eq + total_eq_count]);
            }
        }
    }

    for (int node = 0; node < node_count; node++) {
        printf("NODE %d:\nvoltage: %8.4lf ∠%8.4lf°\tpower: %8.4lf %8.4lfj\n\n", node + 1, loadflow[node].mag, loadflow[node].theta * 57.2958, total_power_flows[node].real * 100, total_power_flows[node].imag * 100);
    }

    for (int line = 0; line < line_count; line++) {
        Cart lineflow = subtract(cart_power_flows[node_count * lines[line].from + lines[line].to], cart_power_flows[node_count * lines[line].to + lines[line].from]);
        printf("LINE %d:\npower: %8.4lf %8.4lfj\n\n", line + 1, lineflow.real * 100, lineflow.imag * 100);
    }
}