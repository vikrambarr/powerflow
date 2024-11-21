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
    struct Node n1 = {1, 3, 1.0, 0.0, {-0.2, 0.4}, 0.0, 0.0};
    struct Node n2 = {2, 3, 1.0, 0.0, {-0.7, -0.7}, 0.0, 0.0};
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
    double angles[num_buses][1];
    struct Cart S[num_buses][1];
    struct Cart dS[num_buses][1];
    double dP[num_buses - 1][1];
    double dQ[num_buses - 1][1];

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
            struct Node bus = buses[i];
            if (bus.type != 1) {
                dP[num_dp][0] = bus.injected.real / 100 - S[i][0].real;
                num_dp++;
                if (bus.type != 2) {
                    dQ[num_dq][0] = bus.injected.imag / 100 - S[i][0].imag;
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
        // TODO: code currently assumes that 1st node is slack

        double J[2 * (num_buses - 1)][2 * (num_buses - 1)];
        
        for (int i = 0; i < num_buses - 1; i++) {
            for (int k = 0; k < num_buses - 1; k++) {

                // set to zero
                J[i][k]                                 = 0;
                J[i + num_buses - 1][k]                 = 0;
                J[i][k + num_buses - 1]                 = 0;
                J[i + num_buses - 1][k + num_buses - 1] = 0;

                // diagonal
                if (i == k) {
                    struct Polar Y_ii = cart2pol(Y[i + 1][k + 1]);
                    J[i + num_buses - 1][k]                 += 2 * buses[i + 1].vlf * Y_ii.mag * cos(Y_ii.theta);
                    J[i + num_buses - 1][k + num_buses - 1] -= 2 * buses[i + 1].vlf * Y_ii.mag * sin(Y_ii.theta);
                    for (int b = 0; b < num_buses; b++) {
                        if (b == i + 1) continue;
                        struct Polar Y_ik = cart2pol(Y[i + 1][b]);
                        J[i][k]                                 -= buses[i + 1].vlf * buses[b].vlf * Y_ik.mag * sin(angles[i + 1][0] - angles[b][0] - Y_ik.theta);
                        J[i][k + num_buses - 1]                 += buses[i + 1].vlf * buses[b].vlf * Y_ik.mag * cos(angles[i + 1][0] - angles[b][0] - Y_ik.theta);
                        J[i + num_buses - 1][k]                 += buses[b].vlf * Y_ik.mag * cos(angles[i + 1][0] - angles[b][0] - Y_ik.theta);
                        J[i + num_buses - 1][k + num_buses - 1] += buses[b].vlf * Y_ik.mag * sin(angles[i + 1][0] - angles[b][0] - Y_ik.theta);
                    }
                }

                // off-diagonal
                else {
                    struct Polar Y_ik = cart2pol(Y[i + 1][k + 1]);
                    J[i][k]                                 += buses[i + 1].vlf * buses[k + 1].vlf * Y_ik.mag * sin(angles[i + 1][0] - angles[k + 1][0] - Y_ik.theta);
                    J[i][k + num_buses - 1]                 -= buses[i + 1].vlf * buses[k + 1].vlf * Y_ik.mag * cos(angles[i + 1][0] - angles[k + 1][0] - Y_ik.theta);
                    J[i + num_buses - 1][k]                 += buses[i + 1].vlf * Y_ik.mag * cos(angles[i + 1][0] - angles[k + 1][0] - Y_ik.theta);
                    J[i + num_buses - 1][k + num_buses - 1] += buses[i + 1].vlf * Y_ik.mag * sin(angles[i + 1][0] - angles[k + 1][0] - Y_ik.theta);
                }
            }
        }

        printf("jacobian:\n");
        for (int y = 0; y < 2 * (num_buses - 1); y++) {
            for (int x = 0; x < 2 * (num_buses - 1); x++) {
                printf("%f\t", J[x][y]);
            }
            printf("\n");
        } printf("\n");
    }

    // calculate powerflow


    return 0;
}