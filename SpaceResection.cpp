#define _CRT_SECURE_NO_WARNINGS


#include<stdio.h>
#include<iostream>
#include<malloc.h>
#include<stdlib.h>
#include<math.h>
#include<conio.h>
#include<ctime>

using namespace std;

#define N 4 
#define M 4 



void RotationMatrix(double omega, double phi, double kappa, double* a, double* b, double* c);

void ErrorEquation(double omega, double kappa, double *A, double *l, double *a, double *b, double *c, double f, double *x, double *y, double *x0, double *y0, double *z0);

double *inv(double *a, int n);

void Multiply(double *A, double *B, double *C, int rows1, int cols1, int rows2, int cols2);

void Transpose(double *A, double AT[], int m, int n);

/*-----------------------------------------------主函数--------------------------------------------*/
int main()    
{

	
	double  x[N], y[N], x0[N], y0[N], z0[N], X[N], Y[N], Z[N];
	x[0] = -86.15; y[0] = -68.99; X[0] = 36589.41; Y[0] = 25273.32; Z[0] = 2195.17;
	x[1] = -53.40; y[1] = 82.21;  X[1] = 37631.08; Y[1] = 31324.51; Z[1] = 728.69;
	x[2] = -14.78; y[2] = -76.63;  X[2] = 39100.97; Y[2] = 24934.98; Z[2] = 2386.50;
	x[3] = 10.46; y[3] = 64.43;  X[3] = 40426.54; Y[3] = 30319.81; Z[3] = 757.31;
	int i, j = 0;
	printf("控制点坐标：\n");
	for (i = 0; i < N; i++)
		printf("%lf    %lf    %lf   %lf   %lf\n", x[i], y[i], X[i], Y[i], Z[i]);

	
	int m = 50000;
	double f = 153.24;
	
	printf("\n摄影机主距和摄影比例尺分母：f = %f, m = %d\n", f, m);
	f = f / 1000.0;


	double phi, omega, kappa, Xs0, Ys0, Zs0;


	double Sum1 = 0, Sum2 = 0, Sum3 = 0;
	for (i = 0; i < N; i++)
	{
		x[i] /= 1000.0;
		y[i] /= 1000.0;
		Sum1 += X[i];
		Sum2 += Y[i];
		Sum3 += Z[i];
	}
	Xs0 = Sum1 / N;
	Ys0 = Sum2 / N;                          
	Zs0 = m * f + Sum3 / N;
	phi = omega = kappa = 0;
	printf("\n外方位元素初始值：Xs0 = %f,Ys0 = %f, Zs0 = %f\nphi = %f, omega = %f, kappa = %f", Xs0, Ys0, Zs0, phi, omega, kappa);

	double m0;      
	while (j < M)	
	{
		
		double a[3], b[3], c[3];
		RotationMatrix(omega, phi, kappa, a, b, c);

		for (i = 0; i < N; i++)
		{
			x0[i] = -f * (a[0] * (X[i] - Xs0) + b[0] * (Y[i] - Ys0) + c[0] * (Z[i] - Zs0)) / (a[2] * (X[i] - Xs0) + b[2] * (Y[i] - Ys0) + c[2] * (Z[i] - Zs0));  
			y0[i] = -f * (a[1] * (X[i] - Xs0) + b[1] * (Y[i] - Ys0) + c[1] * (Z[i] - Zs0)) / (a[2] * (X[i] - Xs0) + b[2] * (Y[i] - Ys0) + c[2] * (Z[i] - Zs0));
			z0[i] = a[2] * (X[i] - Xs0) + b[2] * (Y[i] - Ys0) + c[2] * (Z[i] - Zs0);
		}

		double A[2 * N * 6], l[2 * N];
		
		ErrorEquation(omega, kappa, A, l, a, b, c, f, x, y, x0, y0, z0);
		
		
		double AT[6 * 2 * N], ATA[6 * 6], *InvofATA = NULL, ATl[6], V[6];
		Transpose(A, AT, 2 * N, 6);
		Multiply(AT, A, ATA, 6, 2 * N, 2 * N, 6); 

		
		InvofATA = inv(ATA, 6);

		
		Multiply(AT, l, ATl, 6, 2 * N, 2 * N, 1);
		Multiply(InvofATA, ATl, V, 6, 6, 6, 1);

		
		Xs0 += V[0];
		Ys0 += V[1];
		Zs0 += V[2];
		phi += V[3];
		omega += V[4];
		kappa += V[5];
		printf("\n第%d轮迭代：外方位元素为：\n", j);
		printf("Xs=%.5lf,Ys=%.5lf,Zs=%.5lf,t=%.5lf,w=%.5lf,k=%.5lf\n", Xs0, Ys0, Zs0, phi, omega, kappa);

		
		double s = 0, AX[8], v[8];
		Multiply(A, V, AX, 8, 6, 6, 1); 
		for (i = 0; i < 8; i++)
		{
			v[i] = AX[i] - l[i];
			s += v[i] * v[i];
		}
		m0 = sqrt(s / 2);
		j++;
	}

	printf("\n\n\nXs=%.5lf,Ys=%.5lf,Zs=%.5lf,t=%.5lf,w=%.5lf,k=%.5lf\n", Xs0, Ys0, Zs0, phi, omega, kappa);
	printf("\n中误差为：%f\n", m0);
	system("pause");
	return 0;
}



void RotationMatrix(double omega, double phi, double kappa, double* a, double* b, double* c)
{
	a[0] = cos(phi)*cos(kappa) - sin(phi)*sin(omega)*sin(kappa);
	a[1] = -cos(phi)*sin(kappa) - sin(phi)*sin(omega)*cos(kappa);
	a[2] = -sin(phi)*cos(omega);
	b[0] = cos(omega)*sin(kappa);
	b[1] = cos(omega)*cos(kappa);           
	b[2] = -sin(omega);
	c[0] = sin(phi)*cos(kappa) + cos(phi)*sin(omega)*sin(kappa);
	c[1] = -sin(phi)*sin(kappa) + cos(phi)*sin(omega)*cos(kappa);
	c[2] = cos(phi)*cos(omega);
}


void ErrorEquation(double omega, double kappa, double *A, double *l, double *a, double *b, double *c, double f, double *x, double *y, double *x0, double *y0, double *z0)
{
	for (int i = 0; i < N; i++)
	{
		A[i * 12 + 0] = (a[0] * f + a[2] * x[i]) / z0[i];
		A[i * 12 + 1] = (b[0] * f + b[2] * x[i]) / z0[i];
		A[i * 12 + 2] = (c[0] * f + c[2] * x[i]) / z0[i];
		A[i * 12 + 3] = y[i] * sin(omega) - (x[i] * (x[i] * cos(kappa) - y[i] * sin(kappa)) / f + f * cos(kappa))*cos(omega);
		A[i * 12 + 4] = -f * sin(kappa) - x[i] * (x[i] * sin(kappa) + y[i] * cos(kappa)) / f;
		A[i * 12 + 5] = y[i];
		A[i * 12 + 6] = (a[1] * f + a[2] * y[i]) / z0[i];
		A[i * 12 + 7] = (b[1] * f + b[2] * y[i]) / z0[i];
		A[i * 12 + 8] = (c[1] * f + c[2] * y[i]) / z0[i];
		A[i * 12 + 9] = -x[i] * sin(omega) - (y[i] * (x[i] * cos(kappa) - y[i] * sin(kappa)) / f - f * sin(kappa))*cos(omega);
		A[i * 12 + 10] = -f * cos(kappa) - y[i] * (x[i] * sin(kappa) + y[i] * cos(kappa)) / f;
		A[i * 12 + 11] = -x[i];
		l[i * 2 + 0] = x[i] - x0[i];
		l[i * 2 + 1] = y[i] - y0[i];
	}
}


void Transpose(double *A, double AT[], int m, int n)    //计算矩阵的转置
{
	int i, j;
	for (i = 0; i < m; i++)
		for (j = 0; j < n; j++)
			AT[j*m + i] = A[i*n + j];
}


void Multiply(double *A, double *B, double *C, int rows1, int cols1, int rows2, int cols2)  //计算两矩阵相乘
{
	int i, j, l, u;
	if (cols1 != rows2)
	{
		printf("error!cannot do the multiplication.\n");
		return;
	}
	for (i = 0; i < rows1; i++)
		for (j = 0; j < cols2; j++)
		{
			u = i * cols2 + j;
			C[u] = 0.0;
			for (l = 0; l < cols1; l++)
				C[u] += A[i*cols1 + l] * B[l*cols2 + j];
		}
	return;
}


double *inv(double *a, int n)       
{                                 
	int *is, *js, i, j, k, l, u, v;
	double d, p;
	is = (int*)malloc(n * sizeof(int));
	js = (int*)malloc(n * sizeof(int));
	for (k = 0; k <= n - 1; k++)
	{
		d = 0.0;
		for (i = k; i < n; i++)
			for (j = k; j < n; j++)
			{
				l = i * n + j;
				p = fabs(a[l]);
				if (p > d)
				{
					d = p;
					is[k] = i;
					js[k] = j;
				}
			}
		if (d + 1.0 == 1.0)
		{
			free(is);
			free(js);
			printf("error not inv\n");
			return NULL;
		}
		if (is[k] != k)
			for (j = 0; j < n; j++)
			{
				u = k * n + j;
				v = is[k] * n + j;
				p = a[u];
				a[u] = a[v];
				a[v] = p;
			}
		if (js[k] != k)
			for (i = 0; i < n; i++)
			{
				u = i * n + k;
				v = i * n + js[k];
				p = a[u];
				a[u] = a[v];
				a[v] = p;
			}
		l = k * n + k;
		a[l] = 1.0 / a[l];
		for (j = 0; j < n; j++)
			if (j != k)
			{
				u = k * n + j;
				a[u] = a[u] * a[l];
			}
		for (i = 0; i < n; i++)
			if (i != k)
				for (j = 0; j < n; j++)
					if (j != k)
					{
						u = i * n + j;
						a[u] = a[u] - a[i*n + k] * a[k*n + j];
					}
		for (i = 0; i < n; i++)
			if (i != k)
			{
				u = i * n + k;
				a[u] = -a[u] * a[l];
			}
	}
	for (k = n - 1; k >= 0; k--)
	{
		if (js[k] != k)
			for (j = 0; j <= n - 1; j++)
			{
				u = k * n + j;
				v = js[k] * n + j;
				p = a[u];
				a[u] = a[v];
				a[v] = p;
			}
		if (is[k] != k)
			for (i = 0; i < n; i++)
			{
				u = i * n + k;
				v = i * n + is[k];
				p = a[u];
				a[u] = a[v];
				a[v] = p;
			}
	}
	free(is);
	free(js);
	return a;
}