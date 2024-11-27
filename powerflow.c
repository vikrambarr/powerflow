#include <stdio.h>
#include <math.h>

struct Cart {
    double real;
    double imag;
};

struct Polar {
    double mag;
    double theta;
};

struct Node {
    int num;
    int type;
    double voltage;
    double theta;
    struct Cart injected;
    double minq;
    double maxq;
    double vlf;
    double thetalf;
};

struct Line {
    int from;
    int to;
    struct Cart z;
    double b_half;
};

struct Cart pol2cart(struct Polar in) {
    struct Cart out;
    out.real = in.mag * cos(in.theta);
    out.imag = in.mag * sin(in.theta);
    return out;
}

struct Polar cart2pol(struct Cart in) {
    struct Polar out;
    out.mag = sqrt(pow(in.real, 2) + pow(in.imag, 2));

    if (in.real != 0) out.theta = atan(in.imag / in.real);
    else out.theta = 0;

    if (in.real < 0) out.theta += M_PI;
    return out;
}

struct Cart conjugate(struct Cart in) {
    struct Cart out;
    out.real = in.real;
    out.imag = in.imag * -1;
    return out;
}

struct Cart multiply(struct Cart a, struct Cart b) {
    struct Polar out, a_polar, b_polar;
    a_polar = cart2pol(a); b_polar = cart2pol(b);
    out.mag = a_polar.mag * b_polar.mag;
    out.theta = a_polar.theta + b_polar.theta;
    return pol2cart(out);
}

struct Cart divide(struct Cart a, struct Cart b) {
    struct Polar out, a_polar, b_polar;
    a_polar = cart2pol(a); b_polar = cart2pol(b);
    out.mag = a_polar.mag / b_polar.mag;
    out.theta = a_polar.theta - b_polar.theta;
    return pol2cart(out);
}

struct Cart add(struct Cart a, struct Cart b) {
    struct Cart out;
    out.real = a.real + b.real;
    out.imag = a.imag + b.imag;
    return out;
}

struct Cart subtract(struct Cart a, struct Cart b) {
    struct Cart out;
    out.real = a.real - b.real;
    out.imag = a.imag - b.imag;
    return out;
}

const int num_buses = 3;
const int num_lines = 2;
struct Cart unity = {.real = 1, .imag = 0};
struct Cart nothing = {.real = 0, .imag = 0};

int main() {

    // setup inputs
    struct Node n0 = {0, 1, 1.0, 0.0, {0.0, 0.0}, 0.0, 0.0};
    struct Node n1 = {1, 2, 1.0, 0.0, {0.5, -0.2}, 0.0, 0.0};
    struct Node n2 = {2, 3, 1.0, 0.0, {-1.7, -1.7}, 0.0, 0.0};
    struct Node buses[num_buses] = {n0, n1, n2};

    struct Line l0 = {0, 1, {0.1, 0.2}, 0.02};
    struct Line l1 = {1, 2, {0.05, 0.2}, 0.02};
    struct Line lines[num_lines] = {l0, l1};

    // create admittance matrix
    struct Cart Y[num_buses][num_buses] = {{{.real = 0, .imag = 0}}};

    for (int b = 0; b < num_buses; b++) {

        // setup buses
        struct Node bus = buses[b];
        buses[b].vlf = buses[b].voltage; buses[b].thetalf = buses[b].theta;

        for (int l = 0; l < num_lines; l++) {

            struct Line line = lines[l];
            int from = line.from; int to = line.to;
            struct Cart admittance = divide(unity, line.z);
            struct Cart admittance_prime = {.real = admittance.real, .imag = admittance.imag + line.b_half};

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
    double tol = 0.00001;
    int iters = 1;
    double angles[num_buses][1] = {{0}};
    struct Cart S[num_buses][1] = {{{.real = 0, .imag = 0}}};
    struct Cart dS[num_buses][1] = {{{.real = 0, .imag = 0}}};;
    double dP[num_buses - 1][1] = {{0}};
    double dQ[num_buses - 1][1] = {{0}};

    // iteration
    for (int iter = 0; iter < iters; iter++) {
        for (int i = 0; i < num_buses; i++) {
            for (int k = 0; k < num_buses; k++) {
                struct Polar Y_ik = cart2pol(Y[i][k]);
                S[i][0].real += buses[i].vlf * buses[k].vlf * Y_ik.mag * cos(angles[i][0] - angles[k][0] - Y_ik.theta);
                S[i][0].imag += buses[i].vlf * buses[k].vlf * Y_ik.mag * sin(angles[i][0] - angles[k][0] - Y_ik.theta);
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
                        struct Polar Y_ik = cart2pol(Y[i][b]);
                        J[y][x] -= buses[i].vlf * buses[b].vlf * Y_ik.mag * sin(angles[i][0] - angles[b][0] - Y_ik.theta);
                    }
                }
                else {
                    struct Polar Y_ik = cart2pol(Y[i][k]);
                    J[y][x] += buses[i].vlf * buses[k].vlf * Y_ik.mag * sin(angles[i][0] - angles[k][0] - Y_ik.theta);
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
                    struct Polar Y_ii = cart2pol(Y[i][k]);
                    J[y][x] -= 2 * buses[i].vlf * Y_ii.mag * sin(Y_ii.theta);
                    for (int b = 0; b < num_buses; b++) {
                        if (b == i) continue;
                        struct Polar Y_ik = cart2pol(Y[i][b]);
                        J[y][x] += buses[b].vlf * Y_ik.mag * sin(angles[i][0] - angles[b][0] - Y_ik.theta);
                    }
                }
                else {
                    struct Polar Y_ik = cart2pol(Y[i][k]);
                    J[y][x] += buses[i].vlf * Y_ik.mag * sin(angles[i][0] - angles[k][0] - Y_ik.theta);
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
                        struct Polar Y_ik = cart2pol(Y[i][b]);
                        J[y][x] += buses[i].vlf * buses[b].vlf * Y_ik.mag * cos(angles[i][0] - angles[b][0] - Y_ik.theta);
                    }
                }
                else {
                    struct Polar Y_ik = cart2pol(Y[i][k]);
                    J[y][x] -= buses[i].vlf * buses[k].vlf * Y_ik.mag * cos(angles[i][0] - angles[k][0] - Y_ik.theta);
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
                    struct Polar Y_ii = cart2pol(Y[i][k]);
                    J[y][x] += 2 * buses[i].vlf * Y_ii.mag * cos(Y_ii.theta);
                    for (int b = 0; b < num_buses; b++) {
                        if (b == i) continue;
                        struct Polar Y_ik = cart2pol(Y[i][b]);
                        J[y][x] += buses[b].vlf * Y_ik.mag * cos(angles[i][0] - angles[b][0] - Y_ik.theta);
                    }
                }
                else {
                    struct Polar Y_ik = cart2pol(Y[i][k]);
                    J[y][x] += buses[i].vlf * Y_ik.mag * cos(angles[i][0] - angles[k][0] - Y_ik.theta);
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
        
    }

    // calculate powerflow


    return 0;
}