/*
Code written and edited by Julia Adamczyk and Brendan Aguiar
Version: 0.4

Version History edits
0.1 Built ranf and box_muller function -Brendan Aguiar
0.2 Edited ranf and box muller function -Julia Adamczyk/Brendan Aguiar
0.2 Assigned mean and covariance matrices and started setA -Julia Adamczyk/Brendan Aguiar
0.3 added the method to calculate probability - Julia Adamczyk
0.4 Started Classify and Case1 Function - Brendan Aguiar
*/
#include <iostream>
#include <math.h>
#include <list> 

#define PLACEHOLDER 1
using namespace std;
float ranf(float m);
float box_muller(float m, float s);
void generateDistr(list <float> set[], float m1[], float s1[][2], float m2[], float s2[][2]);
void printDistr(list <float> setA[]);
void classify(list <float> set[], float m1[], float s1[][2], float m2[], float s2[][2], float r[]);
float case1(list<float>::iterator i1, list<float>::iterator i2, float m[], float s, float prior);

int main() {
	//Assignments
	float cov1[2][2] = { {1, 0},
						{0, 1} };
	float cov2[2][2] = { {1, 0},
						{0, 1} };
	float mean1[2] = { 1, 1 };
	float mean2[2] = { 4, 4 };
	list <float> setA[2];
	float errorRates[4]; //r[0] = class 1 missclassification, r[1] = class 2 missclassification,
					// r[2] = Total missclassification, r[3] = Bhattacharyya bound
	

	//Set A generation
	generateDistr(setA, mean1, cov1, mean2, cov2);
	
	//Classifcation
	classify(setA, mean1, cov1, mean2, cov2, errorRates);
	
	//Report and log generation
	//printDistr(setA);

    return 0;
}

/*ranf() written by George Bebis*/
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

/*
Description: Method to generate gaussian distribution.
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

/*
Description: Generates gaussian sample distribution using box muller and provided means and covariances.
*/
void generateDistr(list <float> set[], float m1[], float s1[][2], float m2[], float s2[][2])
{
	for (int i = 1; i <= 40000; i++) //sample generation for class w1
	{
		set[0].push_back(box_muller(m1[0], s1[0][0])); //pushes x value 
		set[1].push_back(box_muller(m1[1], s1[1][1])); //pushes y value
	}
	for (int i = 40001; i <= 200000; i++)//sample generation for class w2
	{
		set[0].push_back(box_muller(m2[0], s2[0][0])); //pushes x value
		set[1].push_back(box_muller(m2[1], s2[1][1])); //pushes y value
	}
}

/*
Description: Prints sample set to console. Could be modified to print to file.
*/
void printDistr(list <float> set[]) 
{
	int it3 = 1;
	list<float>::iterator it2 = set[1].begin();
	for (list<float>::iterator it = set[0].begin(); it != set[0].end(); ++it)
	{
		cout << "n = " << it3 << ", x = " << *it << ", y = " << *it2 << endl;
		it3++;
		++it2;
	}
}

/*
Description: Checks covariances of classes to decide determinant case. Generates P(w1/x) and P(w2/x).
Uses logic from minimum error rate: Decides w1 if P(w1/x) > P(w2/x), else decides w2. 
Increments misclassifcation observations for class 1.
*/
void classify(list <float> set[], float m1[], float s1[][2], float m2[], float s2[][2], float r[])
{
	int miss1 = 0;//missclassification incrementor
	int miss2 = 0;//missclassification incrementor
	float prior1 = .3;//probability of class 1 40,000/200,000
	float prior2 = .7;//probability of class 1 40,000/200,000

	list<float>::iterator it1 = set[0].begin();
	list<float>::iterator it2 = set[1].begin();
	if (s1[0][0] == s2[0][0] && s1[1][1] == s2[1][1]) //Case I where covariances are equal
	{
		float s = s1[1][1];
		for (int i = 1; i <= 40000; i++)//w1 samples
		{
			float g1 = case1(it1, it2, m1, s, prior1);//g1(x) = P(w1/x)
			float g2 = case1(it1, it2, m2, s, prior2);//g2(x) = P(w2/x)
			if (g1 <= g2) //if w1 is missclassied
				miss1++; //increments missclassification rate
			++it1;
			++it2;
		}
		for (int i = 40000; i <= 200000; i++)// w1 samples
		{
			float g1 = PLACEHOLDER;//g1(x) = P(w1/x)
			float g2 = PLACEHOLDER;//g2(x) = P(w2/x)
			if (g1 <= g2) //if w1 is missclassied
				miss2++; //increments missclassification rate
		}
	}
	
	else //Case III where covariances are unequal
	{

	}

	//Error Calculations
	r[0] = miss1 / 40000;
	r[1] = miss2 / 160000;
	r[2] = (miss1 + miss2) / 200000;
	r[3] = PLACEHOLDER;// Bhattacharyya bound
}

/*
Description: Returns the discriminant of case I where covariances of class 1 and class 2 are equal.
*/
float case1(list<float>::iterator i1, list<float>::iterator i2, float m[], float s, float prior)
{
	float e = PLACEHOLDER;//euclideanDistance() or  ||x-m||^2
	return log2(prior) + (e / (2 * s));
}