def generate_unrolled_code(n, m, print_flag=True):
    """
    Generate fully unrolled and branchless C code for the polynomial evaluation
    with precomputed random values.
    
    Args:
        n: Number of variables
        m: Number of polynomials
        print_flag: Whether to include print statements in the generated code
    """
    code = []
    
    # Add headers
    code.append("""#include "../include/measure.h"
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>  // Include stdio.h for printf

// Define maximum sizes for arrays
#define MAX_N {}
#define MAX_M {}
#define MAX_SOLS (1 << MAX_N) // Maximum number of solutions (2^MAX_N)
""".format(n, m))

    # Add structs and utility functions, but omit the random function
    code.append("""
typedef struct {
  int constant;           // 0 or 1
  int linear[MAX_N];      // array of length n (each entry 0 or 1)
  int quad[MAX_N][MAX_N]; // 2D array: for i=0..n-1, quad[i][j]
} Poly;

// Getter functions for the polynomial coefficients.
int get_constant(Poly *f) { return f->constant; }
int get_linear(Poly *f, int k) { return f->linear[k]; }
int get_quad(Poly *f, int k, int j) { return f->quad[k][j]; }

// For bit manipulation
static const int MultiplyDeBruijnBitPosition[32] = {
    0,  1,  28, 2,  29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4,  8,
    31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6,  11, 5,  10, 9};

int bit1(unsigned int x) {
  unsigned int v = x;
  // If x is nonzero, compute the index; if zero, r is meaningless.
  unsigned int r = MultiplyDeBruijnBitPosition[((v & -v) * 0x077CB531U) >> 27];
  // Create a mask that is all 1's when x==0.
  unsigned int m = -(!x);
  //  x==0, return -1; otherwise return r. Who needs branches...
  return (int)((r & ~m) | ((unsigned int)-1 & m));
}

int bit2(unsigned int x) { return bit1(x ^ (x & -x)); }
""")

    # Generate precomputed random values for the polynomial
    import random
    random.seed() 

    def rand_xorshift_sim():
        """Simulate the rand_xorshift function to get predictable values"""
        return random.randint(0, 1)  # Simplified to return 0 or 1 directly
    
    # Generate a known solution
    known_solution = [rand_xorshift_sim() for _ in range(n)]
    
    # Generate polynomial constants, linear terms, and quadratic terms
    poly_constants = []
    poly_linear = []
    poly_quad = []
    
    for i in range(m):
        # Generate random polynomial
        constant = rand_xorshift_sim()
        linear_terms = [rand_xorshift_sim() for _ in range(n)]
        quad_terms = [[0 for _ in range(n)] for _ in range(n)]
        
        for k in range(n):
            for j in range(k):
                quad_terms[k][j] = rand_xorshift_sim()
        
        # Evaluate polynomial at the known solution to see if adjustment is needed
        result = constant
        for i_var, val in enumerate(known_solution):
            if val == 1:
                result ^= linear_terms[i_var]
        
        for k in range(n):
            for j in range(k):
                if known_solution[k] == 1 and known_solution[j] == 1:
                    result ^= quad_terms[k][j]
        
        # If evaluation is 1, flip constant term to make it 0
        if result == 1:
            constant ^= 1
        
        poly_constants.append(constant)
        poly_linear.append(linear_terms)
        poly_quad.append(quad_terms)
    
    # Generate preloaded polynomial code
    code.append("""
// Preloaded polynomials with precomputed random values
void init_preloaded_polys(Poly polys[MAX_M]) {""")
    
    for i in range(m):
        code.append(f"""
  // Initialize polynomial {i}
  polys[{i}].constant = {poly_constants[i]};""")
        
        for j in range(n):
            code.append(f"""
  polys[{i}].linear[{j}] = {poly_linear[i][j]};""")
        
        for k in range(n):
            for j in range(n):
                if k > j:
                    code.append(f"""
  polys[{i}].quad[{k}][{j}] = {poly_quad[i][k][j]};""")
                else:
                    code.append(f"""
  polys[{i}].quad[{k}][{j}] = 0;""")
    
    code.append("}\n")
    
    # Add eval_poly function
    code.append("""
// Evaluate polynomial f on vector x (of length n). All operations are mod 2.
int eval_poly(Poly *f, int *x, int n) {
  int result = f->constant;""")

    # Unroll eval_poly function
    eval_code = []
    for i in range(n):
        eval_code.append(f"  result ^= f->linear[{i}] & x[{i}];")
    
    for i in range(n):
        for j in range(i):
            eval_code.append(f"  result ^= f->quad[{i}][{j}] & x[{i}] & x[{j}];")
    
    code.append("\n".join(eval_code))
    code.append("  return result & 1;\n}\n")

    # Generate struct and init_state
    code.append("""
typedef struct {
  unsigned int i;
  unsigned int y;
  unsigned int d1[MAX_N];
  unsigned int d2[MAX_N][MAX_N];
} State;

void init_state(State *s, Poly polys[MAX_M], int m, int n) {
  s->i = 0;
  s->y = 0;
""")

    # Unroll initialization of d1 and d2
    for k in range(n):
        code.append(f"  s->d1[{k}] = 0;")
        for j in range(n):
            code.append(f"  s->d2[{k}][{j}] = 0;")
    
    code.append("\n  // Initialize arrays with polynomial values")
    for i in range(m):
        code.append(f"  s->y |= ((unsigned int)get_constant(&polys[{i}])) << {i};")
    
    for i in range(m):
        for k in range(n):
            for j in range(k):
                code.append(f"  s->d2[{k}][{j}] |= ((unsigned int)get_quad(&polys[{i}], {k}, {j})) << {i};")

    for i in range(m):
        code.append(f"  s->d1[0] |= ((unsigned int)get_linear(&polys[{i}], 0)) << {i};")
    
    for i in range(m):
        for k in range(1, n):
            code.append(f"""  {{
    int bit = (((s->d2[{k}][{k-1}] >> {i}) & 1) ^ get_linear(&polys[{i}], {k}));
    s->d1[{k}] |= ((unsigned int)bit) << {i};
  }}""")
    
    code.append("}\n")

    # Create step function
    code.append("""
void step(State *s) {
  s->i++;
  int k1 = bit1(s->i);
  int k2 = bit2(s->i);
  int is_k2_none = k2 == -1;

  // Branchless update to d1[k1]
  s->d1[k1] ^= s->d2[(is_k2_none + k2) & ~(!(k2 + 1))][k1];
  s->y ^= s->d1[k1];
}
""")

    # Create fully unrolled fes_eval
    code.append("""
void fes_eval(Poly polys[MAX_M], int m, int n, unsigned int sols[MAX_SOLS], int *nsol) {
  int sol_count = 0;
  State s;

  init_state(&s, polys, m, n);
""")
    
    # Fully unroll 
    for i in range(1, 2**n):
        code.append(f"""  // Step {i}
  {{
      volatile unsigned int mask = -(s.y == 0);
      volatile unsigned int gray_code = (s.i ^ (s.i >> 1));
      sols[sol_count] = (gray_code & mask) | (sols[sol_count] & ~mask);
      sol_count += (s.y == 0);
      step(&s);
  }}""")
    
    code.append("  *nsol = sol_count;\n}\n")

    # Main function with precomputed known solution
    sol_str = ", ".join(str(x) for x in known_solution[:n])  # Only use first n elements
    
    main_code = [f"""
int main(void) {{
  int n = {n};  // number of variables
  int m = {m}; // number of polynomials

  // Precomputed known solution
  int sol[MAX_N] = {{{sol_str}}};
"""]
    
    # Add print statements if requested
    if print_flag:
        main_code.append("""
  printf("Known solution:\\n");
  for (int i = 0; i < n; i++) {{
    printf("%d ", sol[i]);
  }}
  printf("\\n");
""")
    
    main_code.append("""
  // Create polynomial system with preloaded values
  Poly polys[MAX_M];
  init_preloaded_polys(polys);
""")

    # Add polynomial printing if requested
    if print_flag:
        main_code.append("""
  // Print the input polynomials.
  printf("\\nInput polynomials:\\n");
  for (int i = 0; i < m; i++) {{
    printf("f%d = %d", i, polys[i].constant);

    for (int j = 0; j < n; j++) {{
      if (polys[i].linear[j]) {{
        printf(" + x%d", j);
      }}
    }}

    for (int k = 0; k < n; k++) {{
      for (int j = 0; j < k; j++) {{
        if (polys[i].quad[k][j]) {{
          printf(" + x%d*x%d", k, j);
        }}
      }}
    }}

    printf("\\n");
  }}
""")

    main_code.append("""
  // Evaluate (enumerate) all candidate solutions.
  int sol_count;
  unsigned int sols[MAX_SOLS];

  ACTIVATE_COUNTER(0);
  fes_eval(polys, m, n, sols, &sol_count);
  DEACTIVATE_COUNTER(0);
""")

    # solution printing if requested
    if print_flag:
        main_code.append("""
  printf("\\nFound solution: ...\\n");
  // Print the found solutions (each printed as an nbit binary string).
  for (int i = 0; i < sol_count; i++) {{
    unsigned int sol_val = sols[i];
    for (int j = 0; j < n; j++) {{ // Iterate from 0 to n-1 (LSB to MSB)
      printf("%d ", (sol_val >> j) & 1);
    }}
    printf("\\n");
  }}
""")

    main_code.append("  return 0;\n}")
    code.append("".join(main_code))
    
    return "\n".join(code)

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Generate unrolled C code with precomputed random values')
    parser.add_argument('n', type=int, help='Number of variables')
    parser.add_argument('m', type=int, help='Number of polynomials')
    parser.add_argument('--print', dest='print_flag', action='store_true', 
                        help='Include print statements in generated code')
    parser.set_defaults(print_flag=False)
    
    args = parser.parse_args()
    
    unrolled_code = generate_unrolled_code(args.n, args.m, args.print_flag)
    print(unrolled_code)

if __name__ == '__main__':
    main()