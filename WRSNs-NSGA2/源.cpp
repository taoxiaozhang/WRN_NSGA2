#define _CRT_SECURE_NO_DEPRECATE
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include <iostream>  
#include <string.h> 
#include <algorithm>  


#define maxgeneration 500
#define popsize 50
#define Nmax 25                                  //reader充电Locaton数量
#define Node_number 100                         //需要充电Node的数量
#define chromlength Nmax*4
#define Pi 3.14
#define PC 0.7
#define PM 0.1
#define area 100                                //检测区域范围
#define MaxChargeTime 20
using namespace std;
const double eps = 1e-8;


struct Point
{
	double x, y;
};

Point p[5000];

struct individual
{
	double chrom[chromlength];
	double X[Nmax][2];
	double T[Nmax];
	double S[Nmax];
	double f1;
	double f2;
	double f3;
	double f4;
	double fitness;
	int irank;
	double idistance;

};


struct individual population[popsize];
struct individual Cpopulation[popsize];
struct individual Rt[2 * popsize];
struct individual Bt[10 * popsize];
struct individual Best_population[maxgeneration];
Point c;
double r;

double dist(Point A, Point B)                                             //求AB两点的距离
{
	return sqrt((A.x - B.x)*(A.x - B.x) + (A.y - B.y)*(A.y - B.y));
}

/***返回三角形的外心 */
Point circumcenter(Point A, Point B, Point C)
{
	Point ret;
	double a1 = B.x - A.x, b1 = B.y - A.y, c1 = (a1*a1 + b1*b1) / 2;
	double a2 = C.x - A.x, b2 = C.y - A.y, c2 = (a2*a2 + b2*b2) / 2;
	double d = a1*b2 - a2*b1;
	ret.x = A.x + (c1*b2 - c2*b1) / d;
	ret.y = A.y + (a1*c2 - a2*c1) / d;
	return ret;
}

/***c为圆心，r为半径 */      //结果存在c和r中
void min_cover_circle(Point *p, int n, Point &c, double &r)
{
	random_shuffle(p, p + n);
	c = p[0]; r = 0;
	for (int i = 1; i<n; i++)
	{
		if (dist(p[i], c)>r + eps)   //第一个点  
		{
			c = p[i]; r = 0;
			for (int j = 0; j<i; j++)
				if (dist(p[j], c)>r + eps)  //第二个点  
				{
					c.x = (p[i].x + p[j].x) / 2;
					c.y = (p[i].y + p[j].y) / 2;
					r = dist(p[j], c);
					for (int k = 0; k<j; k++)
						if (dist(p[k], c)>r + eps)  //第三个点  
						{   //求外接圆圆心，三点必不共线  
							c = circumcenter(p[i], p[j], p[k]);
							r = dist(p[i], c);
						}
				}
		}
	}
}

void Initialize(int NodeNumber)                       //初始化程序
{
	int i, j;
	int n = NodeNumber;

	for (i = 0; i < n; i++)                       //随机生成n个Node(需要充电的点)
	{
		p[i].x = rand() % (area * 10 + 1) / 10.0;
		p[i].y = rand() % (area * 10 + 1) / 10.0;

	}
	min_cover_circle(p, n, c, r);                     //求n个点的最小覆盖圆
	for (i = 0; i<popsize; i++)
	{
		for (j = 0; j<chromlength; j = j + 4)               //初始化位置范围为最小覆盖圆c内，半径为0到r
		{
			int index;
			index = j * 4;
			double radius = r*(rand() % 1001 / 1000.0);
			double angle = 2 * Pi*(rand() % 1001 / 1000.0);
			population[i].chrom[j] = c.x + radius*cos(angle);
			population[i].chrom[j + 1] = c.y + radius*sin(angle);
			population[i].chrom[j + 2] = MaxChargeTime*(rand() % 1001 / 1000.0);
			population[i].chrom[j + 3] = (rand() % 1001 / 1000.0);

		}


		int t = 0;

		while (t < Nmax)
		{
			population[i].X[t][0] = population[i].chrom[t * 4];
			population[i].X[t][1] = population[i].chrom[t * 4 + 1];
			population[i].T[t] = population[i].chrom[t * 4 + 2];
			population[i].S[t] = population[i].chrom[t * 4 + 3];
			t++;
		}
	}


}


void Evaluate()
{
	int i, j, m;
	int NodeNumber = 100;         //区域中需要充电节点的数量,max1000

	for (i = 0; i < popsize; i++)
	{
		double E[Nmax][1000] = { 0 };
		double Total_E[1000] = { 0 };
		double d;
		double Reader_power;
		Point Rp;
		for (j = 0; j < Nmax; j++)                          //求节点充电能量E与Total_E
		{
			Rp.x = population[i].X[j][0];
			Rp.y = population[i].X[j][1];

			for (m = 0; m < NodeNumber; m++)
			{

				d = dist(Rp, p[m]);
				Reader_power = 36 / ((d + 30)*(d + 30));
				E[j][m] = Reader_power*population[i].T[j];
				Total_E[m] = Total_E[m] + E[j][m];

			}

		}

		for (j = 0; j < NodeNumber; j++)
		{
			if (Total_E[j] < 2)
			{
				double P[1000] = { 0 };
				double D[1000] = { 0 };
				int flag[Nmax] = { 0 };
				for (m = 0; m < Nmax; m++)
					flag[m] = m;


				for (m = 0; m < Nmax; m++)
				{
					Rp.x = population[i].X[m][0];
					Rp.y = population[i].X[m][1];
					D[m] = dist(p[j], Rp);
					d = D[m];
					P[m] = 36 / ((d + 30)*(d + 30));


				}

				for (m = 1; m<Nmax; m++)
					for (int n = 0; n < Nmax - m; n++)
					{
						if (P[n] < P[n + 1])
						{
							double p = P[n];
							P[n] = P[n + 1];
							P[n + 1] = p;
							int index = flag[n];
							flag[n] = flag[n + 1];
							flag[n + 1] = index;
						}
					}

				double sum = 0;
				for (m = 0; m < Nmax; m++)
					sum = sum + P[m];

				double delta_t = (2 - Total_E[j]) / sum;
				for (m = 0; m < Nmax; m++)
					population[i].T[flag[m]] = population[i].T[flag[m]] + delta_t;

				for (m = 0; m <NodeNumber; m++)
				{
					Total_E[m] = 0;
					for (int n = 0; n < Nmax; n++)
						E[n][m] = 0;
				}

				for (m = 0; m < Nmax; m++)                          //求节点充电能量E与Total_E
				{
					Rp.x = population[i].X[m][0];
					Rp.y = population[i].X[m][1];

					for (int n = 0; n < NodeNumber; n++)
					{

						d = dist(Rp, p[n]);
						Reader_power = 36 / ((d + 30)*(d + 30));
						E[m][n] = Reader_power*population[i].T[m];
						Total_E[n] = Total_E[n] + E[m][n];

					}

				}

			}

		}

		double exchange_E[1000] = { 0 };

		for (j = 1; j < Nmax; j++)
		{
			for (m = 0; m < Nmax - j; m++)
				if (population[i].S[m] > population[i].S[m + 1])
				{

					for (int mm = 0; mm < NodeNumber; mm++)
					{
						exchange_E[mm] = E[m][mm];
						E[m][mm] = E[m + 1][mm];
						E[m + 1][mm] = exchange_E[mm];
					}

					double x1 = population[i].X[m][0];
					double y1 = population[i].X[m][1];
					double t1 = population[i].T[m];
					double s1 = population[i].S[m];


					population[i].X[m][0] = population[i].X[m + 1][0];
					population[i].X[m + 1][0] = x1;
					population[i].X[m][1] = population[i].X[m + 1][1];
					population[i].X[m + 1][1] = y1;
					population[i].T[m] = population[i].T[m + 1];
					population[i].T[m + 1] = t1;
					population[i].S[m] = population[i].S[m + 1];
					population[i].S[m + 1] = s1;


				}
		}



		//计算f1,f2,f3,,f4,fitness
		double Total_t = 0;
		double Total_Psi = 0;
		double Total_wastePower = 0;
		double Total_distance = 0;
		for (j = 0; j < Nmax; j++)
		{
			Total_t = Total_t + population[i].T[j];
		}
		for (j = 0; j < NodeNumber; j++)
		{
			Total_wastePower = Total_wastePower + (Total_E[j] - 2);
		}
		double distance;
		Point P1;
		Point P2;
		for (j = 0; j < Nmax - 1; j++)
		{
			P1.x = population[i].X[j][0];
			P1.y = population[i].X[j][1];
			P2.x = population[i].X[j + 1][0];
			P2.y = population[i].X[j + 1][1];
			distance = dist(P1, P2);
			Total_distance = Total_distance + distance;
		}
		population[i].f1 = Total_t;
		population[i].f3 = Total_wastePower;
		population[i].f4 = Total_distance;
		int Charge_flag1[1000] = { 0 };
		int Charge_flag2[1000] = { 0 };
		double Max_chargeNumber = 0;
		double T_E[1000] = { 0 };
		double fullnode[Nmax] = { 0 };
		for (j = 0; j < Nmax; j++)
		{
			for (m = 0; m < NodeNumber; m++)
			{
				T_E[m] = T_E[m] + E[j][m];
				if (T_E[m]>2)
					Charge_flag2[m] = 1;
				else
					Charge_flag2[m] = 0;
			}
			int count = 0;
			for (m = 0; m < NodeNumber; m++)
			{
				count = count + (Charge_flag2[m] - Charge_flag1[m]);
			}
			fullnode[j] = count;
			if (count > Max_chargeNumber)
				Max_chargeNumber = count;

			for (m = 0; m < NodeNumber; m++)
				Charge_flag1[m] = Charge_flag2[m];
		}

		population[i].f2 = Max_chargeNumber;

	}

	double f1max = population[0].f1;
	double f2max = population[0].f2;
	double f3max = population[0].f3;
	double f4max = population[0].f4;
	for (i = 0; i < popsize; i++)
	{
		if (population[i].f1 > f1max)
			f1max = population[i].f1;
		if (population[i].f2 > f2max)
			f2max = population[i].f2;
		if (population[i].f3 > f3max)
			f3max = population[i].f3;
		if (population[i].f4 > f4max)
			f4max = population[i].f4;
	}

	for (i = 0; i < popsize; i++)
		population[i].fitness = 1 - (0.6*(population[i].f1 / f1max) + 0.3*(population[i].f2 / f2max) + 0.1*(population[i].f3 / f3max) + 0 * (population[i].f4 / f4max));


}

void Crossover()
{
	int i, j;
	double p;
	for (i = 0; i < popsize; i++)
	{
		Cpopulation[i] = population[i];
	}

	for (i = 0; i < popsize; i++)
	{
		p = rand() % 1001 / 1000.0;
		if (p < PC)
		{
			int Temp = (int)rand() % popsize;

			Point P1, P2;
			for (j = 0; j < Nmax; j++)
			{
				double D[Nmax] = { 0 };
				P1.x = population[i].X[j][0];
				P1.y = population[i].X[j][1];
				for (int m = 0; m < Nmax; m++)
				{
					P2.x = population[Temp].X[m][0];
					P2.y = population[Temp].X[m][1];
					D[m] = dist(P1, P2);
				}
				int flag = 0;
				double minD = D[0];
				for (int m = 0; m<Nmax; m++)
					if (minD > D[m])
					{
						flag = m;
						minD = D[m];
					}
				Cpopulation[i].X[j][0] = (population[Temp].X[flag][0] + population[i].X[j][0]) / 2;
				Cpopulation[i].X[j][1] = (population[Temp].X[flag][1] + population[i].X[j][1]) / 2;
				Cpopulation[i].T[j] = (population[Temp].T[flag] + population[i].T[j]) / 2;
				Cpopulation[i].S[j] = (population[Temp].S[flag] + population[i].S[j]) / 2;
			}

		}
	}

}

void Mutation()
{
	int i, j;
	double p;
	for (i = 0; i < popsize; i++)
	{
		for (j = 0; j < Nmax; j++)
		{
			p = rand() % 1001 / 1000.0;
			if (p <= PM)
			{
				double p1 = rand() % 1001 / 1000.0;
				double flag = 1;
				if (p1 > 0.5)
					flag = -1;

				Cpopulation[i].X[j][0] = (1 + flag*0.2)*Cpopulation[i].X[j][0];
				Cpopulation[i].X[j][1] = (1 + flag*0.2)*Cpopulation[i].X[j][1];
				Cpopulation[i].T[j] = (1 + flag*0.2)*Cpopulation[i].T[j];
				Cpopulation[i].S[j] = (1 + flag*0.2)*Cpopulation[i].S[j];

			}
		}
	}
}

void Fitness_Sort()
{
	int i, j;
	for (i = 1; i<popsize - 1; i++)
		for (j = 0; j < popsize - i; j++)
		{
			if (population[j].fitness < population[j + 1].fitness)
			{
				individual p = population[j];
				population[j] = population[j + 1];
				population[j + 1] = p;

			}
		}
}
void Selection()
{
	int i, j;
	double p;
	for (i = 0; i < popsize; i++)
		Rt[i] = population[i];

	for (i = 0; i < popsize; i++)
		population[i] = Cpopulation[i];


	Evaluate();
	Fitness_Sort();

	for (i = 0; i < popsize; i++)
		Rt[i + popsize] = population[i];


	for (i = 0; i<2 * popsize - 1; i++)
		for (j = 0; j < 2 * popsize - 1 - i; j++)
		{
			if (Rt[j].fitness <Rt[j + 1].fitness)
			{
				individual p = Rt[j];
				Rt[j] = Rt[j + 1];
				Rt[j + 1] = p;

			}
		}

	for (i = 0; i < popsize; i++)
		population[i] = Rt[i];

}

void OutputData()
{
	FILE *fp;
	fp = fopen("Point.txt", "w");
	for (int i = 0; i < 100; i++)
		fprintf(fp, "%f,%f\n", p[i].x, p[i].y);
	fclose(fp);

	FILE *fc;
	fc = fopen("Cicle.txt", "w");
	fprintf(fc, "%f,%f,%f", c.x, c.y, r);
	fclose(fc);

	fp = fopen("Fitness.txt", "w");
	for (int i = 0; i < maxgeneration; i++)
		fprintf(fp, "%f\n", Best_population[i].fitness);
	fclose(fp);
	fp = fopen("F.txt", "w");
	for (int i = 0; i < maxgeneration; i++)
		fprintf(fp, "%f,%f,%f,%f\n", Best_population[i].f1, Best_population[i].f2, Best_population[i].f3, Best_population[i].f4);
	fclose(fp);

	fp = fopen("Location.txt", "w");
	for (int i = 0; i<Nmax; i++)
		fprintf(fp, "%f,%f\n", population[0].X[i][0], population[0].X[i][1]);
	fclose(fp);
}


void Fastnondominatedsort()
{
	int i, j;
	int point = 0;
	int np[2 * popsize] = { 0 };
	int Fi[2 * popsize] = { 0 };
	int Sp[2 * popsize][2 * popsize] = { 0 };
	for (i = 0; i<2 * popsize; i++)
	{
		for (j = 0; j<2 * popsize; j++)
		{
			if (Rt[i].f1<Rt[j].f1 || Rt[i].f2<Rt[j].f2)
				Sp[i][j] = 1;
			else if (Rt[i].f1>=Rt[j].f1 || Rt[i].f2>=Rt[j].f2)
				np[i] = np[i] + 1;
		}
		if (np[i] == 0)
		{
			Rt[i].irank = 1;
			Fi[point] = i;
			point++;
		}
	}
	int t = 1;
	int prank = t;
	while (point != 0)
	{
		int Q[2 * popsize] = { 0 };
		int Qpoint = 0;
		for (i = 0; i<point; i++)
			for (j = 0; j<2 * popsize; j++)
			{
				if (Sp[Fi[i]][j] == 1)
				{
					np[j] = np[j] - 1;
					if (np[j] == 0)
					{
						Rt[j].irank = t + 1;
						Q[Qpoint] = j;
						Qpoint++;
					}
				}
			}

		t++;
		prank = t;
		for (i = 0; i<Qpoint; i++)
			Fi[i] = Q[i];
		point = Qpoint;

	}
}
void Crowdingdistancesort()
{
	int i, j;
	double fmin, fmax;
	struct individual temp;
	for (i = 0; i<2 * popsize; i++)
		Rt[i].idistance = 0;

	for (i = 0; i < 2 * popsize; i++)
		for (j = 0; j < 2 * popsize - i - 1; j++)
		{
			if (Rt[j].f1 > Rt[j + 1].f1)
			{
				temp = Rt[j];
				Rt[j] = Rt[j + 1];
				Rt[j + 1] = temp;
			}
		}
	Rt[0].idistance = 1000000000;
	Rt[2 * popsize - 1].idistance = 1000000000;
	fmin = Rt[i].f1;
	fmax = Rt[2 * popsize - 1].f2;
	if (fmin != fmax)
	{
		for (i = 1; i < 2 * popsize - 1; i++)
			Rt[i].idistance = Rt[i].idistance + (Rt[i + 1].f1 - Rt[i - 1].f1) / (fmax - fmin);
	}

	for (i = 0; i < 2 * popsize; i++)
		for (j = 0; j < 2 * popsize - i - 1; j++)
		{
			if (Rt[j].f2 > Rt[j + 1].f2)
			{
				temp = Rt[j];
				Rt[j] = Rt[j + 1];
				Rt[j + 1] = temp;
			}
		}
	Rt[0].idistance = 1000000000;
	Rt[2 * popsize - 1].idistance = 1000000000;
	fmin = Rt[i].f2;
	fmax = Rt[2 * popsize - 1].f2;
	if (fmin != fmax)
	{
		for (i = 1; i < 2 * popsize - 1; i++)
			Rt[i].idistance = Rt[i].idistance + (Rt[i + 1].f2 - Rt[i - 1].f2) / (fmax - fmin);
	}
}
void NSGA2()
{
	int i;
	for (i = 0; i < popsize; i++)
		Rt[i] = population[i];

	for (i = 0; i < popsize; i++)
		Rt[i + popsize] = Cpopulation[i];

	Fastnondominatedsort();
	Crowdingdistancesort();
	int a = 1;

}
void main()
{

	srand(time(NULL));
	Initialize(100);
	Evaluate();
	Fitness_Sort();
	//记录最好解的初始布局位置
	FILE *fp;
	fp = fopen("Location0.txt", "w");
	for (int i = 0; i<Nmax; i++)
	fprintf(fp, "%f,%f\n", population[0].X[i][0], population[0].X[i][1]);
	fclose(fp);

	int t = 0;
	while (t < maxgeneration)
	{

		Crossover();
		Mutation();

		NSGA2();
		//Selection();

		t = t + 1;
		printf("Generation: %d\n", t);

		getchar();
	}
}