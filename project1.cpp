/*
Code written and edited by Julia Adamczyk and Brendan Aguiar
Version: 0.8

Version History edits
0.1 Built ranf and box_muller function -Brendan Aguiar
0.2 Edited ranf and box muller function -Julia Adamczyk/Brendan Aguiar
0.2 Assigned mean and covariance matrices and started setA -Julia Adamczyk/Brendan Aguiar
0.3 added the method to calculate probability - Julia Adamczyk
0.4 Started Classify and Case1 Function - Brendan Aguiar
0.5 Completed case1, case3, and euclidean discriminant functions (need a review) - Julia Adamczyk
0.6 Reviewed discriminant functions. Cleaned up descriptions - Brendan Aguiar
0.7 Completed classify function for case3, provide calculations for error rates, added global floats for sizes - Julia Adamczyk
0.8 Implemented Bhattacharyya Bound, printErrorReport, classifyEuclidean, menu, and switch case - Brendan Aguiar
*/
//Library and namespace inclusions
#include <iostream>
#include <math.h>
#include <list>
#include <fstream>
using namespace std;

//Gloabl assignments
const float CLASS_1_SIZE = 60000.0;
const float CLASS_2_SIZE = 140000.0;
const float TOTAL_SIZE = 200000.0;

//Function declarations
float ranf(float m);
float box_muller(float m, float s);
void generateDistr(list <float> set[], float m1[], float s1[][2], float m2[], float s2[][2]);
void printDistr(list <float> setA[]);
void classify(list <float> set[], float m1[], float s1[][2], float m2[], float s2[][2], float r[]);
float case1(list<float>::iterator i1, list<float>::iterator i2, float m[], float s, float prior);
float case3(list<float>::iterator i1, list<float>::iterator i2, float m[], float s[][2], float prior);
float euclidean(list<float>::iterator i1, list<float>::iterator i2, float m[]);
float determinant_of_diagonal(float mat[][2]);
float BatBound(float m1[], float s1[][2], float m2[], float s2[][2]);
void printMenu();
void classifyEuclidean(list <float> set[], float m1[], float s1[][2], float m2[], float s2[][2], float r[]);
void printErrorReport(float r[], string dataset);

int main() {

	//Assignments for setA
	list <float> setA[2];
	float cov1[2][2] = { {1, 0},
						{0, 1} };
	float cov2[2][2] = { {1, 0},
						{0, 1} };
	float mean1[2] = { 1, 1 };
	float mean2[2] = { 4, 4 };
	float errorRates1[4]; //r[0] = class 1 missclassification, r[1] = class 2 missclassification,
					// r[2] = Total missclassification, r[3] = Bhattacharyya bound

	//Assignments for setB
	list <float> setB[2];
	float cov3[2][2] = { {1, 0},
						{0, 1} };
	float cov4[2][2] = { {4, 0},
						{0, 8} };
	float mean3[2] = { 1, 1 };
	float mean4[2] = { 4, 4 };
	float errorRates2[4];
	
	//Program Control
	bool again = true;
	int switch_on;
	while (again)
	{
		printMenu();
		cin >> switch_on;
		switch (switch_on)
		{
		case 1:
			//Set A generation
			generateDistr(setA, mean1, cov1, mean2, cov2);
			break;
		case 2:
			//Set B generation
			generateDistr(setB, mean3, cov3, mean4, cov4);
			break;
		case 3:
			//Print Set A to text file
			printDistr(setA);
			break;
		case 4:
			//Print Set B to text file
			printDistr(setB);
			break;
		case 5:
			//Set A classification using discriminant cases
			classify(setA, mean1, cov1, mean2, cov2, errorRates1);
			break;
		case 6:
			//Set A and Set B classification using discriminant cases
			classify(setB, mean3, cov3, mean4, cov4, errorRates2);
			break;
		case 7:
			//Set A classification using euclidean distances
			classifyEuclidean(setA, mean1, cov1, mean2, cov2, errorRates1);
			break;
		case 8:
			//Set B classification using euclidean distances
			classifyEuclidean(setB, mean3, cov3, mean4, cov4, errorRates2);
			break;
		case 9:
			//Print Set A distribution 
			printDistr(setA);
			break;
		case 10:
			//Print Set B distribution 
			printDistr(setB);
			break;
		case 11:
			//Print Error Report Set A
			printErrorReport(errorRates1, "Set A");
			break;
		case 12:
			//Print Error Report Set B
			printErrorReport(errorRates2, "Set B");
			break;
		case 13:
			again = false;
			break;
		default:
			break;
		}
	}
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
			x1 = 2.0f * ranf(1) - 1.0f;
			x2 = 2.0f * ranf(1) - 1.0f;
			w = x1 * x1 + x2 * x2;
		} while (w >= 1.0f);

		w = sqrt((-2.0f * log(w)) / w);
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
	for (int i = 1; i <= 60000; i++) //sample generation for class w1
	{
		set[0].push_back(box_muller(m1[0], s1[0][0])); //pushes x value
		set[1].push_back(box_muller(m1[1], s1[1][1])); //pushes y value
	}
	for (int i = 60001; i <= 200000; i++)//sample generation for class w2
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
	cout << "Enter file name: " << endl;
	string filename;
	cin >> filename;
	ofstream fout;
	fout.open(filename);
	cout << "Writing to file now..." << endl;
	list<float>::iterator it2 = set[1].begin(); //iterator for y values
	for (list<float>::iterator it = set[0].begin(); it != set[0].end(); ++it)// initialize loop with iterator for x values
	{
		fout << *it << " " << *it2 << endl;
		++it2;
	}
	cout << "File successfully written." << endl;
	fout.close();
}

/*
Description: Checks covariances of classes to decide discriminant case. Generates P(w1/x) and P(w2/x).
Uses logic from minimum error rate: Decides w1 if P(w1/x) > P(w2/x), else decides w2.
Increments misclassifcation observations for class 1.
*/
void classify(list <float> set[], float m1[], float s1[][2], float m2[], float s2[][2], float r[])
{
	int miss1 = 0;//missclassification incrementor
	int miss2 = 0;//missclassification incrementor
	float prior1 = .3f;//probability of class 1 40,000/200,000
	float prior2 = .7f;//probability of class 1 40,000/200,000

	//decision boundary points when g1(x) == g2(x)
	//for the future we need to store them somewhere to plot
	//i don't know if it works cause it comes up as 0 so we can't really do the decision boundary like that
	float count_decision_boundary = 0;

	list<float>::iterator it1 = set[0].begin();
	list<float>::iterator it2 = set[1].begin();
	if (s1[0][0] == s2[0][0] && s1[1][1] == s2[1][1]) //Case I where covariances are equal
	{
		float s = s1[1][1];
		for (int i = 1; i <= 60000; i++)//w1 samples should hold more weight in g1
		{
			float g1 = case1(it1, it2, m1, s, prior1);//g1(x) = P(w1/x)
			float g2 = case1(it1, it2, m2, s, prior2);//g2(x) = P(w2/x)
			if (g1 >= g2) //if w1 is missclassied
				miss1++; //increments missclassification rate
			if (g1 == g2)
				count_decision_boundary++;
			++it1;
			++it2;
		}
		for (int i = 60001; i <= 200000; i++)// w2 samples should hold more weight in g2
		{
			float g1 = case1(it1, it2, m1, s, prior1);//g1(x) = P(w1/x)
			float g2 = case1(it1, it2, m2, s, prior2);
			if (g1 < g2) //if w2 is missclassified
				miss2++; //increments missclassification rate
			if (g1 == g2)
				count_decision_boundary++;
			++it1;
			++it2;
		}
		cout << "Number of samples of class 1 that are missclassified: " << miss1 << endl;
		cout << "Number of samples of class 2 that are missclassified: " << miss2 << endl;
		cout << "Decision boundary points:" << count_decision_boundary << endl;
	}
	else //Case III where covariances are unequal
	{
		for (int i = 1; i <= 40000; i++)//w1 samples should hold more weight in g1
		{
			float g1 = case3(it1, it2, m1, s1, prior1);//g1(x) = P(w1/x)
			float g2 = case3(it1, it2, m2, s2, prior2);//g2(x) = P(w2/x)
			if (g1 >= g2) //if w1 is missclassied
				miss1++; //increments missclassification rate
			if (g1 == g2)
				count_decision_boundary++;
			++it1;
			++it2;
		}
		for (int i = 40001; i <= 200000; i++)// w2 samples should hold more weight in g2
		{
			float g1 = case3(it1, it2, m1, s1, prior1);//g1(x) = P(w1/x)
			float g2 = case3(it1, it2, m2, s2, prior2);
			if (g1 < g2) //if w2 is missclassied
				miss2++; //increments missclassification rate
			if (g1 == g2)
				count_decision_boundary++;
			++it1;
			++it2;
		}
		cout << "Number of samples of class 1 that are missclassified: " << miss1 << endl;
		cout << "Number of samples of class 2 that are missclassified: " << miss2 << endl;
		cout << "Decision boundary points:" << count_decision_boundary << endl;
		//Error Calculations
		cout << "Error rates for problem 2:" << endl;
	}
	//Error Calculations
	cout << "Error rates for problem 1:" << endl;
	r[0] = miss1 / CLASS_1_SIZE;
	r[1] = miss2 / CLASS_2_SIZE;
	r[2] = (miss1 + miss2) / TOTAL_SIZE;
	r[3] = BatBound(m1, s1, m1, s2);// Bhattacharyya bound
	cout << r[0] << " " << r[1] << " " << r[2] << " " << r[3] << endl;
}

/*
Description: Returns the discriminant of case I where covariances of class 1 and class 2 are equal.
*/
float case1(list<float>::iterator i1, list<float>::iterator i2, float m[], float s, float prior)
{
	float e = euclidean(i1, i2, m); // ||x-m||^2
	return log(prior) + (e / (2 * s * s));
}

float case3(list<float>::iterator i1, list<float>::iterator i2, float m[], float s[][2], float prior)
{
	float inverse_sx = 1 / s[0][0];//should equal s[1][1]
	float inverse_sy = 1 / s[1][1];//should equal s[0][0]
	float determinant = determinant_of_diagonal(s);
	float addend1 = (*i1) * (-0.5f * inverse_sx) * (*i1) + (*i2) * (-0.5f * inverse_sy) * (*i2); //(x^t)*Wi*x
	float addend2 = (inverse_sx * m[0] * (*i1)) + (inverse_sy * m[1] * (*i2)); //(wi)^t*x
	float addend3 = -0.5f * ((m[0] * m[0] * inverse_sx) + (m[1] * m[1] * inverse_sy) - 0.5f * log2(determinant) + log2(prior)); //wi0
	return (addend1 + addend2 + addend3);
}

/*
Description: Returns the euclidean distance.
*/
float euclidean(list<float>::iterator i1, list<float>::iterator i2, float m[])
{
	return (((*i1) - m[0]) * ((*i1) - m[0]) + ((*i2) - m[1]) * ((*i2) - m[1])); //euclideanDistance() or  ||x-m||^2
}

/*
Description: Returns the determinant of the passed in diagonal matrix.
*/
float determinant_of_diagonal(float mat[][2])
{
	float det = 1;
	for (int i = 0; i < 2; i++) {
		det *= mat[i][i];
	}
	return det;
}

/*
Description: Calculates and returns the Bhattacharyya bound.
*/
float BatBound(float m1[], float s1[][2], float m2[], float s2[][2])
{
	float addend1 = pow(m1[0] - m2[0], 2) * (1 / (.5f * (s1[0][0] + s2[0][0]))) + (pow(m1[1] - m2[1], 2) * (1 / (.5f * (s1[1][1] + s2[1][1]))));
	addend1 = addend1 * .125f;
	float mat1[2][2] = { {.5f * (s1[0][0] + s2[0][0]), 0},
		{0, .5f * (s1[1][1] + s2[1][1])} };
	float det1 = determinant_of_diagonal(mat1);
	float det2 = determinant_of_diagonal(s1);
	float det3 = determinant_of_diagonal(s2);
	float addend2 = .5f * log2(det1 / (pow(det2, .5f) * pow(det3, .5f)));
	return addend1 + addend2;
}

/*
Description: Prints menu for switch cases.
*/
void printMenu() {
	cout << "Select from the following choices...\n1. Generate SetA \n2. Generate SetB \n3. Print SetA \n4. Print SetB";
	cout << "\n5. Classify SetA \n6. Classify SetB \n7. Classify(Euclidean) SetA \n8. Classify(Euclidean) SetB \n9. Print SetA \n10. Print SetB";
	cout << "\n11. Print Error Report SetA \n12. Print Error Report SetB \n13. Quit Program";
	cout << "\n\nNote: Generate data before you classify it. Classify data before you generate an error report.\n\n";
	cout << "Choice: ";
}

/*
Description: Classifies a dataset using only euclidean distance.
*/
void classifyEuclidean(list <float> set[], float m1[], float s1[][2], float m2[], float s2[][2], float r[])
{
	int miss1 = 0;//missclassification incrementor
	int miss2 = 0;//missclassification incrementor
	float prior1 = .3f;//probability of class 1 40,000/200,000
	float prior2 = .7f;//probability of class 1 40,000/200,000

	list<float>::iterator it1 = set[0].begin();
	list<float>::iterator it2 = set[1].begin();
	if (s1[0][0] == s2[0][0] && s1[1][1] == s2[1][1]) //Case I where covariances are equal
	{
		float s = s1[1][1];
		for (int i = 1; i <= 60000; i++)//w1 samples should hold more weight in g1
		{
			float g1 = euclidean(it1, it2, m1);//g1(x) = P(w1/x)
			float g2 = euclidean(it1, it2, m2);//g2(x) = P(w2/x)
			if (g1 >= g2) //if w1 is missclassied
				miss1++; //increments missclassification rate
			if (g1 == g2)

			++it1;
			++it2;
		}
		for (int i = 60001; i <= 200000; i++)// w2 samples should hold more weight in g2
		{
			float g1 = euclidean(it1, it2, m1);//g1(x) = P(w1/x)
			float g2 = euclidean(it1, it2, m2);//g2(x) = P(w2/x)
			if (g1 < g2) //if w2 is missclassied
				miss2++; //increments missclassification rate
			++it1;
			++it2;
		}
		cout << "Number of samples of class 1 that are missclassified: " << miss1 << endl;
		cout << "Number of samples of class 2 that are missclassified: " << miss2 << endl;
	}
	//Error Calculations
	cout << "Error rates for problem 1:" << endl;
	r[0] = miss1 / CLASS_1_SIZE;
	r[1] = miss2 / CLASS_2_SIZE;
	r[2] = (miss1 + miss2) / TOTAL_SIZE;
	r[3] = BatBound(m1, s1, m1, s2);// Bhattacharyya bound
	cout << r[0] << " " << r[1] << " " << r[2] << " " << r[3] << endl;
}

/*
Description: Print out the error report for the passed in data set to a file.
*/
void printErrorReport(float r[], string dataset)
{
	//r[0] = class 1 missclassification, r[1] = class 2 missclassification,
					// r[2] = Total missclassification, r[3] = Bhattacharyya bound
	cout << "Enter file name: " << endl;
	string filename;
	cin >> filename;
	ofstream fout;
	fout.open(filename);
	cout << "Writing to file now..." << endl;
	fout << "Error Report for " << dataset << ":" << endl;
	fout << "Class 1 Missclassification Rate: " << r[0] << endl;
	fout << "Class 2 Missclassification Rate: " << r[1] << endl;
	fout << "Total Missclassification Rate: " << r[2] << endl;
	fout << "Bhattacharyya bound: " << r[3] << endl;
	cout << "File successfully written." << endl;
	fout.close();
}