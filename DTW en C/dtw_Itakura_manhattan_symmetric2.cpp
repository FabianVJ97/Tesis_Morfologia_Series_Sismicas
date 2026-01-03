/* Fabian Vargas Jimenez | Ingeniero Estadistico | 2026
 * Tesis: https://github.com/FabianVJ97/Tesis_Morfologia_Series_Sismicas?tab=readme-ov-file#tesis_morfologia_series_sismicas
 * LinkedIn: https://www.linkedin.com/in/fabi%C3%A1n-oscar-andr%C3%A9s-vargas-jim%C3%A9nez-4a85131a4/
 */
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
using namespace Rcpp;

 //Distancia Manhattan
inline double dist_manhattan_1d(double a, double b) {
  return std::abs(a - b);
}

  //Patron symmetric2
inline double step_symmetric2(double g_left,
                              double g_up,
                              double g_diag,
                              double d_ij) {
  double c1 = g_left + d_ij;       //izquierda
  double c2 = g_up   + d_ij;       // arriba
  double c3 = g_diag + 2.0 * d_ij; // diagonal penalizada
  return std::min(c1, std::min(c2, c3));
}

// Limites tipo Itakura 
inline void itakura_limits(int i, int n, int m,
                           int &j_start, int &j_end) {
  // Evitar division por cero en series de longitud 1
  double t;
  if (n > 1) {
    t = static_cast<double>(i - 1) / static_cast<double>(n - 1);
  } else {
    t = 0.0;
  }

  //Limites en coordenadas normalizadas [0,1]
  double lower_norm = std::max(t / 2.0, 2.0 * t - 1.0);
  double upper_norm = std::min(2.0 * t, (t + 1.0) / 2.0);

  //Convertir a indices 1..m
  //m-1 porque j_norm = (j-1)/(m-1) aproximadamente
  int js = static_cast<int>(std::floor(lower_norm * (m - 1))) + 1;
  int je = static_cast<int>(std::ceil (upper_norm * (m - 1))) + 1;

  //Recortar a [1, m]
  js = std::max(1, std::min(m, js));
  je = std::max(1, std::min(m, je));

  // Si por alguna razon js > je, forzar un unico punto cercano a la diagonal
  if (js > je) {
    int mid = (n > 1)
      ? (static_cast<int>(std::round(t * (m - 1))) + 1)
      : 1;
    mid = std::max(1, std::min(m, mid));
    js = je = mid;
  }

  j_start = js;
  j_end   = je;
}

// [[Rcpp::export]]
double dtw_itakura_manhattan_symmetric2_cpp(const NumericVector& A,
                                            const NumericVector& B,
                                            int window_size) {
  int n = A.size();
  int m = B.size();
// 'window_size' se ignora en esta version Itakura,
// se mantiene solo por compatibilidad de firma.
  const double INF = 1e20;
  NumericMatrix C(n + 1, m + 1);

  //Inicializacion a infinito
  for (int i = 0; i <= n; ++i) {
    for (int j = 0; j <= m; ++j) {
      C(i, j) = INF;
    }
  }
  C(0, 0) = 0.0;

//DP con ventana Itakura y patron symmetric2
  for (int i = 1; i <= n; ++i) {
    int j_start, j_end;
    itakura_limits(i, n, m, j_start, j_end);

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
  // Distancia DTW (sin normalizar)
  return C(n, m);
}
