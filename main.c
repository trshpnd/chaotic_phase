#include "pendulo_caotico.h"
#include <time.h>

const double x_0 = 1.00;
const double tv_0 = 0.00;
const double transiente = 800.00;

	//// 1 arquivo com t, x, v
	//// outro para poincare (x, v)
	//// fazer 2 vetores de arquivos de tamanho 7 (1 vetor p/ os 2 primeiros graficos, outro p/ poincaré)
	//// criar arquivos com loop

typedef struct Pendulo{
	double x[7];			// vetor para os valores de x de 0.4 a 1.0
	double v[7];			// vetor para os valores de v de 0.4 a 1.0
	FILE *phase_arr[7];		// vetor para os arquivos dos dados dos espaços de fase de 0.4 a 1.0
	FILE *poincare_arr[7];	// vetor para os arquivos dos dados das seções de poincaré de 0.4 a 1.0
} Pendulo;

Pendulo Pendulo_default = {
	.x = x_0,
	.v = tv_0
};

int main(void) {
	clock_t tStart = clock();
	printf("Criando arquivos... ");

	double t_final = transiente + 20000.00;
	double t_phase = 5000.00;
	double periodo = (2*PI)/w;
	double f = 0.4;
	int tam_vetor = 7;

	Pendulo p = Pendulo_default;

	gerador_arquivos(p.phase_arr, p.poincare_arr, tam_vetor); 	// Cria um arquivo de saída com o respectivo filename
	printf("Pronto.\n");										// em cada posição do vetor do struct
	 
	double t;
	printf("Plotando pontos... ");

	for(t = tv_0; t <= t_final + h; t = t + h){

		for(int i = 0; i < tam_vetor; i++){ //loop prints

			if(t > transiente){
				if(t <= t_phase){
					fprintf(p.phase_arr[i], "%f %f %f\n", t-transiente, p.x[i], p.v[i]);//fprint phase_0.x.dat (t, x, v)
				}
				if(fabs(fmod(t, periodo)) <= 0.01){		// apenas se for múltiplo do periodo
					fprintf(p.poincare_arr[i], "%f %f\n", p.x[i], p.v[i]);				//fprint poincare_0.x.dat (x, v)
				}
			}
		}

		f = 0.4;	//reinicia f
		
		for(int i = 0; i < 7; i++){ 					//loop rk: atualizando x e v a cada f (0.4 ~ 1.0)
			runge_kutta4(t, &p.x[i], &p.v[i], f);
			f = f + 0.1;
		}
	}

	printf("Pronto.\n");
	printf("Processo finalizado em %f segundos.\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
	
	return 0;
}