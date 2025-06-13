#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <time.h>


#define h 0.01 // Espaciado 
#define PI 3.141592653589793
#define N 1000 // Número de puntos espaciales
#define I _Complex_I // Definición de la unidad imaginaria
int T = 2000; // Número de pasos temporales

// Función para inicializar la función de onda Phi
void inicializar_phi(double complex *Phi, double n_ciclos) {
    double k_tilde = (2.0 * PI * n_ciclos) / N; // \tilde{k}_0 = 2π n_ciclos / N
    double x0 = (N * h) / 4.0;                 // x_0 = Nh / 4
    double sigma = (N * h) / 20.0;             // σ = Nh / 16

    for (int j = 1; j < N; j++) {
        double x = j * h; // Posición espacial
        Phi[j] = cexp(I * k_tilde * j) * cexp(-pow((x - x0), 2) / (2 * sigma * sigma)); // Φ_{j,0} = e^{i\tilde{k}_0 x} e^{-(x - x_0)^2 / (2σ^2)}
    }

    Phi[0] = 0.0;
    Phi[N] = 0.0;

    // Normalización de Phi
    double norma = 0.0;
    for (int j = 0; j <= N; j++) {
        norma += pow(cabs(Phi[j]), 2) * h; // Norma = Σ |Φ[j]|^2 * h
    }
    norma = sqrt(norma);

    for (int j = 1; j <N; j++) {
        Phi[j] /= norma; // Normalizamos Φ[j]
    }
}

// Función para inicializar el potencial V
void inicializar_potencial(double *V, double lambda, double k_tilde) {
    double V_j_h2 = lambda * k_tilde * k_tilde; // Altura del potencial
    int inicio = (2 * N) / 5;
    int fin = (3 * N) / 5;

    for (int j = 0; j <= N; j++) {
        if (j >= inicio && j <= fin) {
            V[j] = V_j_h2; // Dentro del rango del potencial
        } else {
            V[j] = 0.0; // Fuera del rango del potencial
        }
    }
}

// Función para encontrar el máximo global de un array PD de tamaño T
int maximo_global(double *PD, int T) {
    int idx_max = 0;
    double max = PD[0];
    for (int n = 1; n < T; n++) {
        if (PD[n] > max) {
            max = PD[n];
            idx_max = n;
        }
    }
    return idx_max;
}

// Valor esperado de j y su varianza
void valor_esperado_j_var(double complex *Phi, double *media, double *varianza) {
    double suma = 0.0, suma2 = 0.0;
    for (int j = 0; j <= N; j++) {
        double prob = pow(cabs(Phi[j]), 2) * h;
        suma += j * prob; // Valor esperado <j> = Σ j |Φ[j]|^2 * h
        suma2 += j * j * prob;
    }
    *media = suma;
    *varianza = suma2 - suma * suma;
}

// Valor esperado de la energía cinética, su varianza y <T^2> usando la cuarta derivada
void valor_esperado_energia_cinetica(double complex *Phi, double *media, double *varianza) {
    double suma = 0.0, suma_T2 = 0.0;
    // Segunda derivada para <T> y <T^2>
    for (int j = 2; j < N-1; j++) {
        // Laplaciano (segunda derivada central)
        double complex laplaciano = (Phi[j+1] - 2.0*Phi[j] + Phi[j-1]) / (h*h);
        double energia = -0.5 * creal(conj(Phi[j]) * laplaciano) * h;
        suma += energia;
        // Cuarta derivada central discreta para <T^2>
        double complex cuarta_derivada = (Phi[j-2] - 4.0*Phi[j-1] + 6.0*Phi[j] - 4.0*Phi[j+1] + Phi[j+2]) / pow(h, 4);
        // Operador cinético aplicado dos veces: (-1/2 d^2/dx^2)^2 = (1/4) d^4/dx^4
        double energia2 = 0.25 * creal(conj(Phi[j]) * cuarta_derivada) * h;
        suma_T2 += energia2;
    }
    *media = suma;
    *varianza = suma_T2 - suma * suma;
}

int main() {

    srand((unsigned)time(NULL));

    int m = 0; // Número de experimentos
    int m_T = 0;  // Número de transmisiones
    double PD[T]; // Guardar PD(n) en cada paso de tiempo
    int n_D = 0;  
    double PD_nD = 0.0;

    double complex Phi[N + 1], Phi_next[N + 1], chi[N + 1];
    double V[N + 1];
    complex double a[N + 1], b[N + 1], c[N + 1], d[N + 1];
    complex double gamma[N -1];

       // Parámetros iniciales
    int n_ciclos = N/16; // Número de ciclos
    double lambda = 0.5; // Altura del potencial
    double k_tilde = (2.0 * PI * n_ciclos) / N; 

  
    double tilde_s = 1.0 / (4.0 * k_tilde * k_tilde);

    FILE *f_pd = fopen("probabilidades_pd.txt", "w");
    if (f_pd == NULL) {
    fprintf(stderr, "Error al abrir el fichero para guardar PD.\n");
    return 1;
    } 
     FILE *f_j = fopen("j_esperado.dat", "w");
    if (f_j == NULL) {
        fprintf(stderr, "Error al abrir el fichero para guardar los datos.\n");
        return 1;
    }

     FILE *f_E = fopen("energia_esperada.dat", "w");
    if (f_E == NULL) {
        fprintf(stderr, "Error al abrir el fichero para guardar los datos.\n");
        return 1;
    }


    for (int experimento = 0; experimento <=m; experimento++) {

    a[0] = 0.0;
    b[0] = 0.0;
    c[0] = 0.0; 
    d[0] = 0.0; 

    // Inicializamos la función de onda y el potencial
    inicializar_phi(Phi, n_ciclos);
    inicializar_potencial(V, lambda, k_tilde);

     // Abrimos el fichero para guardar las normas
     FILE *f_normas = fopen("normas.txt", "w");
     if (f_normas == NULL) {
         fprintf(stderr, "Error al abrir el fichero para guardar las normas.\n");
         return 1;
     }

     // Abrimos el fichero para guardar los datos de la función de onda
    FILE *f_datos = fopen("schrodinger.dat", "w");
    if (f_datos == NULL) {
        fprintf(stderr, "Error al abrir el fichero para guardar los datos.\n");
        return 1;
    }

   
      for (int j = 1; j < N; j++) {
        a[j] = 1.0;
        b[j] = -2.0 + 2.0 * I / tilde_s - V[j];
        c[j] = 1.0;
      }

      // Calculamos alpha fuera del bucle temporal
      complex double alpha[N];
      alpha[N-1] = 0.0;

        for (int j = N-1; j > 0; j--) {
            gamma[j] = b[j] + c[j] * alpha[j];
            alpha[j-1] = -a[j] / gamma[j];
        }

        int j;

    // Iteración en el tiempo
    for (int n = 0; n < T; n++) {

         // Guardar datos para animacion.py: x, |Phi|^2, V
        for (j = 0; j <= N; j++) {
            double x = j * h;
            double densidad = pow(cabs(Phi[j]), 2);
            fprintf(f_datos, "%.8f,%.8f,%.8f\n", x, densidad, V[j]);
        }
        fprintf(f_datos, "\n"); // Línea en blanco entre bloques

        // Construimos los coeficientes del sistema tridiagonal
        for (int j = 1; j < N; j++) {
            d[j] = 4.0 * I / tilde_s * Phi[j];
        }

        // Resolvemos el sistema tridiagonal para chi
        complex double beta[N]; // Beta aún depende del tiempo
        beta[N - 1] = 0.0 + 0.0*I;      // β_{N-1,n} = 0
        // Orden descendente
        for (int j = N - 1; j > 0; j--) {
            beta[j - 1] = (d[j] - c[j] * beta[j]) / gamma[j];
        }

        // Orden ascendente para calcular χ
        chi[0] = 0.0;
        chi[N] = 0.0;
        for (int j = 0; j < N - 1; j++) {
            chi[j + 1] = alpha[j]*chi[j] + beta[j];
        }

        // Calculamos Phi_next
        for (int j = 1; j < N; j++) {
            Phi_next[j] = chi[j] - Phi[j];
        }

        // Actualizamos Phi
        for (int j = 1; j < N; j++) {
            Phi[j] = Phi_next[j];
        }

    
        double norma = 0.0;
        for (int j = 1; j <N; j++) {
            norma += pow(cabs(Phi[j]), 2) * h; // Norma = Σ |Φ[j]|^2 * h
        }

       
        fprintf(f_normas, "Paso %d: Norma = %.10f\n", n, norma);
        
        if(experimento==0) {
            double pd = 0.0;
            for (int j = (4*N)/5; j <= N; j++) {
                pd += pow(cabs(Phi[j]), 2) * h;
            }
            PD[n] = pd;

             fprintf(f_pd, "%d %.20f\n", n, PD[n]);

                double media_j, var_j;
   
                 double media_E, var_E;
            valor_esperado_j_var(Phi, &media_j, &var_j);
            valor_esperado_energia_cinetica(Phi, &media_E, &var_E);

        fprintf(f_j, "%d %.10f %.10f\n", n, media_j, sqrt(var_j));
        fprintf(f_E, "%d %.10f %.10f\n", n, media_E, sqrt(var_E));

        }
   
    }
        fclose(f_normas);
        fclose(f_datos);

        // Buscar el primer máximo global de PD(n)
        if(experimento==0) {
            n_D = maximo_global(PD, T);
          printf("Máximo global de PD(n) en n = %d\n", n_D);

           PD_nD = PD[n_D];
             printf("PD(n_D) = %.20f\n", PD_nD);
        }
       
       T=n_D;
       int detectada;
       if(experimento>0) {
        double p = (double)rand() / ((double)RAND_MAX + 1.0);
       // printf("Número aleatorio p = %.20f\n", p);
             if (p < PD_nD) {
                 m_T++;
                 detectada=1;
            }else 
            {
                 detectada=0;
             }
     }

}
    fclose(f_pd);
    fclose(f_j);
    fclose(f_E);
    double K = (double)m_T / m;
    printf("Número de transmisiones m_T = %d\n", m_T);
    printf("Coeficiente de transmisión K = %.6f\n", K);

    return 0;
}




