#include <iostream>
#include <math.h>
using namespace std;
double ranf(double m);
void boxMuller(int m, double &y1, double &y2);
int main() {
    double x1, x2, y1, y2;
    bool again = 1;
    double m = 1;
    while (again)
    {
        boxMuller(m, y1, y2);
        cout << y1 << " <-y1 y2-> " << y2 << endl << "again? 0 for no.\n";
        cin >> again;
    }
    return 0;
}

double ranf(double m) {

	return (m * rand() / (double)RAND_MAX);
}

void boxMuller(int m, double &y1, double &y2) {
    double w, x1, x2;
    do {
        x1 = 2.0 * ranf(m) - 1.0;
        x2 = 2.0 * ranf(m) - 1.0;
        w = x1 * x1 + x2 * x2;
    } while (w >= 1.0);

    w = sqrt((-2.0 * log2(w)) / w);
    y1 = x1 * w;
    y2 = x2 * w;
}