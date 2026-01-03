/* Fabian Vargas Jimenez | Ingeniero Estadistico | 2026
 * Tesis: https://github.com/FabianVJ97/Tesis_Morfologia_Series_Sismicas?tab=readme-ov-file#tesis_morfologia_series_sismicas
 * LinkedIn: https://www.linkedin.com/in/fabi%C3%A1n-oscar-andr%C3%A9s-vargas-jim%C3%A9nez-4a85131a4/
 */
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
using namespace Rcpp;

inline double dist_manhattan_1d(double a, double b) {
  return std::abs(a - b);
}

inline double step_symmetric2(double g_left,
                              double g_up,
                              double g_diag,
                              double d_ij) {
  double c1 = g_left + d_ij;
  double c2 = g_up   + d_ij;
  double c3 = g_diag + 2.0*d_ij;
  return std::min(c1, std::min(c2, c3));
}

inline void sakoe_limits(int i, int n, int m, int w,
                         int &j_start, int &j_end) {
  j_start = std::max(1, i - w);
  j_end   = std::min(m, i + w);
}

// [[Rcpp::export]]
double dtw_sakoe_manhattan_symmetric2_cpp(const NumericVector& A,
                                          const NumericVector& B,
                                          int window_size) {
  int n = A.size();
  int m = B.size();
  
  if (window_size <= 0 || window_size > std::max(n, m)) {
    window_size = std::max(n, m);
  }
  
  const double INF = 1e20;
  NumericMatrix C(n + 1, m + 1);
  
//Inicializacion a infinito
  for (int i = 0; i <= n; ++i) {
    for (int j = 0; j <= m; ++j) {
      C(i, j) = INF;
    }
  }
  C(0, 0) = 0.0;
  
   // DP con banda Sakoe Chiba con patron symmetric2
  for (int i = 1; i <= n; ++i) {
    int j_start, j_end;
    sakoe_limits(i, n, m, window_size, j_start, j_end);
    
    for (int j = j_start; j <= j_end; ++j) {
      double d_ij = dist_manhattan_1d(A[i - 1], B[j - 1]);
      C(i, j) = step_symmetric2(
        C(i,   j - 1),  // izquierda
        C(i-1, j    ),  // arriba
        C(i-1, j - 1),  // diagonal
        d_ij
      );
    }
  }
  
  //Distancia DTW 
  return C(n, m);
}
