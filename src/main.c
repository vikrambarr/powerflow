#include <stdio.h>
#include <math.h>
#include "../include/powerflow.h"

const int num_buses = 3;
const int num_lines = 2;

int main() {

    // setup inputs
    Node n0 = {0, 1, {1.0, 0.0}, {0.0, 0.0}, { 0.0,  0.0}};
    Node n1 = {1, 2, {1.0, 0.0}, {0.0, 0.0}, { 0.5, -0.2}};
    Node n2 = {2, 3, {1.0, 0.0}, {0.0, 0.0}, {-1.7, -1.7}};
    Node buses[num_buses] = {n0, n1, n2};

    Line l0 = {0, 1, {0.1, 0.2}, 0.02};
    Line l1 = {1, 2, {0.05, 0.2}, 0.02};
    Line lines[num_lines] = {l0, l1};

    // create admittance matrix
    Cart Y[num_buses][num_buses] = {{{.real = 0, .imag = 0}}};

    for (int b = 0; b < num_buses; b++) {

        // setup buses
        buses[b].loadflow = buses[b].phasor;

        for (int l = 0; l < num_lines; l++) {

            Line line = lines[l];
            int from = line.from; int to = line.to;
            Cart admittance = divide(CART_UNITY, line.z);
            Cart admittance_prime = {.real = admittance.real, .imag = admittance.imag + line.b_half};

            // off-diagonal
            if (b == 0) {
                Y[from][to] = subtract(Y[from][to], admittance);
                Y[to][from] = subtract(Y[to][from], admittance);
            }

            // diagonal
            if (line.from == b || line.to == b) {
                Y[b][b] = add(Y[b][b], admittance_prime);
            }
        }
    }

    printf("admittance matrix:\n");
    for (int y = 0; y < num_buses; y++) {
        for (int x = 0; x < num_buses; x++) {
            printf("%f, %fj\t", Y[x][y].real, Y[x][y].imag);
        }
        printf("\n");
    } printf("\n");

    // setup iteration
    double tol = 1;
    int iters = 10;
    double angles[num_buses][1] = {{0}};
    Cart S[num_buses][1] = {{{.real = 0, .imag = 0}}};
    double dP[num_buses - 1][1] = {{0}};
    double dQ[num_buses - 1][1] = {{0}};

    // iteration
    for (int iter = 0; iter < iters; iter++) {

        for (int i = 0; i < num_buses; i++) {
            S[i][0].real = 0; S[i][0].imag = 0;
            for (int k = 0; k < num_buses; k++) {
                Polar Y_ik = cart2pol(Y[i][k]);
                S[i][0].real += buses[i].loadflow.mag * buses[k].loadflow.mag * Y_ik.mag * cos(angles[i][0] - angles[k][0] - Y_ik.theta);
                S[i][0].imag += buses[i].loadflow.mag * buses[k].loadflow.mag * Y_ik.mag * sin(angles[i][0] - angles[k][0] - Y_ik.theta);
            }
        }

        printf("s matrix\n");
        for (int x = 0; x < num_buses; x++) {
            printf("%f %fj\n", S[x][0].real, S[x][0].imag);
        } printf("\n");

        // get dP and dQ
        // TODO: define BMVa (ratio between 1pu V and 1pu MVA) a bit better
        int num_dp = 0; int num_dq = 0;
        for (int i = 0; i < num_buses; i++) {
            if (buses[i].type != 1) {
                dP[num_dp][0] = buses[i].injected.real / 100 - S[i][0].real;
                num_dp++;
                if (buses[i].type != 2) {
                    dQ[num_dq][0] = buses[i].injected.imag / 100 - S[i][0].imag;
                    num_dq++;
                }
            }
        }

        // create dP dQ matrix
        double M[num_dp + num_dq][1];
        for (int dp_eq = 0; dp_eq < num_dp; dp_eq++) {
            M[dp_eq][0] = dP[dp_eq][0];
        }
        for (int dq_eq = 0; dq_eq < num_dq; dq_eq++) {
            M[dq_eq + num_dp][0] = dQ[dq_eq][0];
        }

        printf("pq matrix\n");
        for (int x = 0; x < num_dp + num_dq; x++) {
            printf("%f\n", M[x][0]);
        } printf("\n");

        // jacobian
        double J[num_dp + num_dq][num_dp + num_dq + 1];
        for (int x = 0; x < num_dp + num_dq; x++) {
            J[x][num_dp + num_dq] = M[x][0];
        }
        int y; int x;

        // J11
        y = 0; x = 0;
        for (int i = 0; i < num_buses; i++) {
            if (buses[i].type == 1) continue;
            for (int k = 0; k < num_buses; k++) {
                if (buses[k].type == 1) continue;
                J[y][x] = 0;
                if (y == x) {
                    for (int b = 0; b < num_buses; b++) {
                        if (b == i) continue;
                        Polar Y_ik = cart2pol(Y[i][b]);
                        J[y][x] -= buses[i].loadflow.mag * buses[b].loadflow.mag * Y_ik.mag * sin(angles[i][0] - angles[b][0] - Y_ik.theta);
                    }
                }
                else {
                    Polar Y_ik = cart2pol(Y[i][k]);
                    J[y][x] += buses[i].loadflow.mag * buses[k].loadflow.mag * Y_ik.mag * sin(angles[i][0] - angles[k][0] - Y_ik.theta);
                }
                x++;
            }
            y++; x = 0;
        }

        // J22
        y = num_dp; x = num_dp;
        for (int i = 0; i < num_buses; i++) {
            if (buses[i].type == 1 || buses[i].type == 2) continue;
            for (int k = 0; k < num_buses; k++) {
                if (buses[k].type == 1 || buses[k].type == 2) continue;
                J[y][x] = 0;
                if (y == x) {
                    Polar Y_ii = cart2pol(Y[i][k]);
                    J[y][x] -= 2 * buses[i].loadflow.mag * Y_ii.mag * sin(Y_ii.theta);
                    for (int b = 0; b < num_buses; b++) {
                        if (b == i) continue;
                        Polar Y_ik = cart2pol(Y[i][b]);
                        J[y][x] += buses[b].loadflow.mag * Y_ik.mag * sin(angles[i][0] - angles[b][0] - Y_ik.theta);
                    }
                }
                else {
                    Polar Y_ik = cart2pol(Y[i][k]);
                    J[y][x] += buses[i].loadflow.mag * Y_ik.mag * sin(angles[i][0] - angles[k][0] - Y_ik.theta);
                }
                x++;
            }
            y++; x = num_dp;
        }

        // J21
        y = num_dp; x = 0;
        for (int i = 0; i < num_buses; i++) {
            if (buses[i].type == 1 || buses[i].type == 2) continue;
            for (int k = 0; k < num_buses; k++) {
                if (buses[k].type == 1) continue;
                J[y][x] = 0;
                if (y == x) {
                    for (int b = 0; b < num_buses; b++) {
                        if (b == i) continue;
                        Polar Y_ik = cart2pol(Y[i][b]);
                        J[y][x] += buses[i].loadflow.mag * buses[b].loadflow.mag * Y_ik.mag * cos(angles[i][0] - angles[b][0] - Y_ik.theta);
                    }
                }
                else {
                    Polar Y_ik = cart2pol(Y[i][k]);
                    J[y][x] -= buses[i].loadflow.mag * buses[k].loadflow.mag * Y_ik.mag * cos(angles[i][0] - angles[k][0] - Y_ik.theta);
                }
                x++;
            }
            y++; x = 0;
        }

        // J12
        y = 0; x = num_dp;
        for (int i = 0; i < num_buses; i++) {
            if (buses[i].type == 1) continue;
            for (int k = 0; k < num_buses; k++) {
                if (buses[k].type == 1 || buses[k].type == 2) continue;
                J[y][x] = 0;
                if (y == x) {
                    Polar Y_ii = cart2pol(Y[i][k]);
                    J[y][x] += 2 * buses[i].loadflow.mag * Y_ii.mag * cos(Y_ii.theta);
                    for (int b = 0; b < num_buses; b++) {
                        if (b == i) continue;
                        Polar Y_ik = cart2pol(Y[i][b]);
                        J[y][x] += buses[b].loadflow.mag * Y_ik.mag * cos(angles[i][0] - angles[b][0] - Y_ik.theta);
                    }
                }
                else {
                    Polar Y_ik = cart2pol(Y[i][k]);
                    J[y][x] += buses[i].loadflow.mag * Y_ik.mag * cos(angles[i][0] - angles[k][0] - Y_ik.theta);
                }
                x++;
            }
            y++; x = num_dp;
        }

        printf("jacobian:\n");
        for (int y = 0; y < num_dp + num_dq; y++) {
            for (int x = 0; x < num_dp + num_dq + 1; x++) {
                printf("%f\t", J[y][x]);
            }
            printf("\n");
        } printf("\n");

        // linalg solve
        double solutions[num_dp + num_dq];
        int n = num_dp + num_dq;

        for (int col = 0; col < n; col++) {
            // pivoting for numerical stability
            int max_row = col;
            for (int i = col + 1; i < n; i++) {
                if (fabs(J[i][col]) > fabs(J[max_row][col])) {
                    max_row = i;
                }
            }
            // swap rows
            if (max_row != col) {
                for (int j = 0; j <= n; j++) {
                    double temp = J[col][j];
                    J[col][j] = J[max_row][j];
                    J[max_row][j] = temp;
                }
            }

            // eliminate below the pivot
            for (int row = col + 1; row < n; row++) {
                double factor = J[row][col] / J[col][col];
                for (int j = col; j <= n; j++) {
                    J[row][j] -= factor * J[col][j];
                }
            }
        }

        // back substitution
        for (int i = n - 1; i >= 0; i--) {
            solutions[i] = J[i][n];
            for (int j = i + 1; j < n; j++) {
                solutions[i] -= J[i][j] * solutions[j];
            }
            solutions[i] /= J[i][i];
        }

        printf("solutions\n");
        for (int x = 0; x < num_dp + num_dq; x++) {
            printf("%f\n", solutions[x]);
        } printf("\n");

        x = 0;
        for (int i = 0; i < num_buses; i++) {
            if (buses[i].type == 1) continue;
            angles[i][0] += solutions[x];
            buses[i].loadflow.theta = angles[i][0];
            x++;
        }

        for (int i = 0; i < num_buses; i++) {
            if (buses[i].type == 1 || buses[i].type == 2) continue;
            buses[i].loadflow.mag += solutions[x];
            x++;
        }

        int max_mismatch = 0;
        for (int i = 0; i < num_dp + num_dq; i++) {
            if (fabs(M[i][0]) > fabs(M[max_mismatch][0])) {
                max_mismatch = i;
            }
        }
        tol = fabs(M[max_mismatch][0]);
        printf("tol: %16.15lf\n", tol);
        printf("vlfs\n");
        for (int x = 0; x < num_buses; x++) {
            printf("%f\n", buses[x].loadflow.mag);
        } printf("\n");
        printf("angles\n");
        for (int x = 0; x < num_buses; x++) {
            printf("%f\n", angles[x][0]);
        } printf("\n");
    }

    // calculate powerflow


    return 0;
}