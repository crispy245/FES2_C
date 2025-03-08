#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// --------------------- Data Structures ---------------------

// Representation of a quadratic polynomial in GF(2) in n variables:
//   f(x) = constant + sum_{i} linear[i]*x_i + sum_{0 <= j < i < n} quad[i][j]*x_i*x_j
typedef struct {
    int constant;    // 0 or 1
    int *linear;     // array of length n (each entry 0 or 1)
    int **quad;      // 2D array: for i=0..n-1, quad[i][j] valid for j < i (others unused)
} Poly;

// The enumeration state used by fes_eval.
typedef struct {
    unsigned int i;        // current Gray code counter (from 0 up to (1<<n)-1)
    unsigned int y;        // aggregated constant (m-bit mask: bit i from poly i)
    unsigned int *d1;      // array of length n; each d1[k] is an m–bit mask
    unsigned int **d2;     // n x n array; d2[k][j] is an m–bit mask (only for k>j used)
} State;

// --------------------- Utility Functions ---------------------

// Return the index (0-indexed) of the lowest–set bit in x; if x == 0, return -1.
int bit1(unsigned int x) {
    if (x == 0) return -1;
    int pos = 0;
    while ((x & 1) == 0) {
        x >>= 1;
        pos++;
    }
    return pos;
}

// Return the index of the “second” (i.e. next) set bit, defined as bit1( x XOR (x & -x) ).
int bit2(unsigned int x) {
    unsigned int first = x & -x;
    unsigned int remainder = x ^ first;
    return bit1(remainder);
}

// Simple “getter” functions for the polynomial coefficients.
int get_constant(Poly* f) {
    return f->constant;
}
int get_linear(Poly* f, int k) {
    return f->linear[k];
}
int get_quad(Poly* f, int k, int j) {
    return f->quad[k][j];
}

// --------------------- Polynomial Functions ---------------------

// Evaluate polynomial f on vector x (of length n). All operations are mod 2.
int eval_poly(Poly* f, int *x, int n) {
    int result = f->constant;
    for (int i = 0; i < n; i++) {
        if (x[i])
            result ^= f->linear[i];
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < i; j++) {
            if (x[i] && x[j])
                result ^= f->quad[i][j];
        }
    }
    return result & 1;
}

// Create a new random quadratic polynomial in n variables.
// The constant, linear, and (lower‐triangle of) quadratic coefficients are chosen at random (0 or 1).
Poly* random_poly(int n) {
    Poly* poly = malloc(sizeof(Poly));
    poly->constant = rand() % 2;
    poly->linear = malloc(n * sizeof(int));
    for (int i = 0; i < n; i++) {
        poly->linear[i] = rand() % 2;
    }
    poly->quad = malloc(n * sizeof(int*));
    for (int i = 0; i < n; i++) {
        poly->quad[i] = malloc(n * sizeof(int));
        for (int j = 0; j < n; j++) {
            poly->quad[i][j] = 0;
        }
    }
    // Only fill the entries with j < i.
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < i; j++) {
            poly->quad[i][j] = rand() % 2;
        }
    }
    return poly;
}

// --------------------- FES Enumeration Functions ---------------------

// The "init" function builds the initial State from an array of m polynomials (polys)
// in n variables. (It “assembles” the bit–masks for the constant parts and the differentials.)
State init_state(Poly **polys, int m, int n) {
    State s;
    s.i = 0;
    s.y = 0;
    s.d1 = malloc(n * sizeof(unsigned int));
    s.d2 = malloc(n * sizeof(unsigned int*));
    for (int k = 0; k < n; k++) {
        s.d1[k] = 0;
        s.d2[k] = malloc(n * sizeof(unsigned int));
        for (int j = 0; j < n; j++) {
            s.d2[k][j] = 0;
        }
    }
    // In the original Python code the "X" are the generators (x0,x1,...,x_{n-1}).
    // For each polynomial (indexed by i), add its bits into y, d2, and d1.
    for (int i = 0; i < m; i++) {
        Poly* f = polys[i];
        s.y |= ((unsigned int)get_constant(f)) << i;
        for (int k = 0; k < n; k++) {
            for (int j = 0; j < k; j++) {
                s.d2[k][j] |= ((unsigned int)get_quad(f, k, j)) << i;
            }
        }
        s.d1[0] |= ((unsigned int)get_linear(f, 0)) << i;
        for (int k = 1; k < n; k++) {
            int bit = (((s.d2[k][k-1] >> i) & 1) ^ get_linear(f, k));
            s.d1[k] |= ((unsigned int)bit) << i;
        }
    }
    return s;
}

// The "step" function updates the State to move to the next point in the Gray code enumeration.
void step(State* s) {
    s->i = s->i + 1;
    int k1 = bit1(s->i);
    int k2 = bit2(s->i);
    if (k2 != -1) { // (in the original Python code: "if k2:" – note that 0 is falsy, so only nonzero k2 trigger an update)
        s->d1[k1] ^= s->d2[k2][k1];
    }
    s->y ^= s->d1[k1];
}

// fes_eval enumerates all 2^n points in GF(2)^n (using Gray code order)
// and returns those indices (converted via i^(i>>1)) where the “aggregate” value y is 0.
// (In our version we simply store each solution as an unsigned int whose binary
// expansion gives the solution vector.)
unsigned int* fes_eval(Poly **polys, int m, int n, int *nsol) {
    // Maximum number of points is 2^n.
    unsigned int *sols = malloc((1u << n) * sizeof(unsigned int));
    int sol_count = 0;
    State s = init_state(polys, m, n);
    // Enumerate over all points (note: the Python code uses "while s.i < 2^n - 1",
    // where 2^n is written as (1 << n)).
    while (s.i < ((1u << n) - 1)) {
        if (s.y == 0) {
            // Append the candidate solution: convert Gray code index to the corresponding binary vector
            sols[sol_count++] = s.i ^ (s.i >> 1);
        }
        step(&s);
    }
    // Free the memory allocated for the state arrays.
    free(s.d1);
    for (int k = 0; k < n; k++) {
        free(s.d2[k]);
    }
    free(s.d2);
    *nsol = sol_count;
    return sols;
}

// --------------------- Main ---------------------

int main(void) {
    int n = N;   // number of variables
    int m = M;  // number of polynomials

    srand((unsigned)time(NULL));

    // Generate a random "known solution" (an n–bit vector)
    int *sol = malloc(n * sizeof(int));
    printf("Known solution:\n");
    for (int i = 0; i < n; i++) {
        sol[i] = rand() % 2;
        printf("%d ", sol[i]);
    }
    printf("\n\n");

    // Create an array of m random quadratic polynomials in n variables.
    // Then “adjust” each so that the known solution is a root (i.e. f(sol) = 0).
    Poly **polys = malloc(m * sizeof(Poly*));
    for (int i = 0; i < m; i++) {
        polys[i] = random_poly(n);
        if (eval_poly(polys[i], sol, n) == 1) {
            // Flip the constant term so that f(sol) becomes 0.
            polys[i]->constant ^= 1;
        }
    }

    // Print the input polynomials.
    printf("Input polynomials:\n");
    for (int i = 0; i < m; i++) {
        printf("f%d = %d", i, polys[i]->constant);
        for (int j = 0; j < n; j++) {
            if (polys[i]->linear[j])
                printf(" + x%d", j);
        }
        for (int k = 0; k < n; k++) {
            for (int j = 0; j < k; j++) {
                if (polys[i]->quad[k][j])
                    printf(" + x%d*x%d", k, j);
            }
        }
        printf("\n");
    }
    printf("\nPerforming full evaluation...\n\n");

    // Evaluate (enumerate) all candidate solutions.
    int sol_count;
    unsigned int *sols = fes_eval(polys, m, n, &sol_count);

    // Print the found solutions (each printed as an n–bit binary string).
    printf("Found solutions:\n");
    for (int i = 0; i < sol_count; i++) {
        unsigned int sol_val = sols[i];
        for (int j = 0; j < n; j++) { // Iterate from 0 to n-1 (LSB to MSB)
            printf("%d", (sol_val >> j) & 1);
            printf(" ");
        }
        printf("\n");
    }

    // Free allocated memory.
    free(sols);
    free(sol);
    for (int i = 0; i < m; i++) {
        free(polys[i]->linear);
        for (int j = 0; j < n; j++) {
            free(polys[i]->quad[j]);
        }
        free(polys[i]->quad);
        free(polys[i]);
    }
    free(polys);

    return 0;
}
