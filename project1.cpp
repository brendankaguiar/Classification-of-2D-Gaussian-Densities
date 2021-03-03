/*
Code written and edited by Julia Adamczyk and Brendan Aguiar
Version: 0.2

Version History edits
0.1 Built ranf and box_muller function -Brendan Aguiar
0.2 Edited ranf and box muller function -Julia Adamczyk/Brendan Aguiar
0.2 Assigned mean and covariance matrices and started setA -Julia Adamczyk/Brendan Aguiar
*/

#include <iostream>
#include <math.h>

using namespace std;
float ranf(float m);
float box_muller(float m, float s);


int main() {

	float mean1[2] = { 1, 1 };
	float cov1[2] = { 1 , 1 };
	float mean2[2] = { 4, 4 };
	float cov2[2] = { 1 , 1 };
	float setA[2][2000];
	
	for (int i = 0; i < 2000; i++)
	{
		setA[0][i] = box_muller(mean1[0], cov1[0]);
		setA[1][i] = box_muller(mean1[1], cov1[1]);
	}
	
	for (int i = 0; i < 2000; i++)
	{
		cout << "x = " << setA[0][i] << " y = " << setA[1][i] << endl;
	}

    return 0;
}

float ranf(float m) {

	return (m * rand() / (float)RAND_MAX);
}

/* boxmuller.c           Implements the Polar form of the Box-Muller
						 Transformation
					  (c) Copyright 1994, Everett F. Carter Jr.
						  Permission is granted by the author to use
			  this software for any application provided this
			  copyright notice is preserved.
*/

float box_muller(float m, float s)	/* normal random variate generator */
{				        /* mean m, standard deviation s */
	float x1, x2, w, y1;
	static float y2;
	static int use_last = 0;

	if (use_last)		        /* use value from previous call */
	{
		y1 = y2;
		use_last = 0;
	}
	else
	{
		do {
			x1 = 2.0 * ranf(1) - 1.0;
			x2 = 2.0 * ranf(1) - 1.0;
			w = x1 * x1 + x2 * x2;
		} while (w >= 1.0);

		w = sqrt((-2.0 * log(w)) / w);
		y1 = x1 * w;
		y2 = x2 * w;
		use_last = 1;
	}

	return(m + y1 * s);
}
