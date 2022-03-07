#include <stdio.h>
#include <math.h>

#define PI M_PI

const double h = 0.01;
const double c = 0.05;
const double w = 0.7;

double dx(double v){	//dx/dt
	return v;
}

double dv(double t, double x, double v, double f){	//dv/dt
	return -c * v -sin(x) + f * (cos(w * t));
}

void runge_kutta4(double t, double *x, double *v, double f){	//rk4	
	double k1_x = dx(*v);
	double k1_v = dv(t, *x, *v, f);

	double k2_x = dx(*v + (k1_v * h)/2);
	double k2_v = dv(t + (h/2), *x + ((k1_x/2) * h), *v + ((k1_v/2) * h), f);

	double k3_x = dx(*v + (k2_v * h)/2);
	double k3_v = dv(t + (h/2), *x + ((k2_x/2) * h), *v + ((k2_v/2) * h), f);

	double k4_x = dx(*v + (k3_v * h));
	double k4_v = dv(t + h, *x + (k3_x * h), *v + (k3_v * h), f);

	*x = *x + (((k1_x + (2 * k2_x) + (2 * k3_x) + k4_x) * h) / 6);
	*v = *v + (((k1_v + (2 * k2_v) + (2 * k3_v) + k4_v) * h) / 6);

	if(fabs(*x) > PI){ // ajuste de -pi a pi

		if(*x > 0){
			*x = *x - 2*PI;
		}

		else{
			*x = *x + 2*PI;
		}
	}
}

void gerador_arquivos(FILE **arr_1, FILE **arr_2, int numfiles){
	float f = 0.4;
	for (int i = 0; i < numfiles; i++)
	{
		char filename[20];

    	sprintf(filename, "phase_%.1f.dat", f);
    	arr_1[i] = fopen(filename, "w");

		sprintf(filename, "poincare_%.1f.dat", f);
		arr_2[i] = fopen(filename, "w");

		f = f + 0.1;
	}
}