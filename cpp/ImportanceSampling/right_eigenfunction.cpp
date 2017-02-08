#include <algorithm>
#include <armadillo>
#include <cmath>
#include <map>

using namespace std;
using namespace arma;

struct garch21_param
{
    double omega;
    double a1;
    double a2;
    double b1;
};

										  double right_eigenfunction(struct garch21_param &param, unsigned )
{
    
}

int main(int argc, char **argv)
{
    struct garch21_param param = {
	9.225747e-06, 8.835834e-02, 9.685783e-02, 6.543018e-01
    };
    unsigned n;
    double alpha = 4.47103;

    // cout << "alpha[1]: ";
    // cin >> param.a1;

    // cout << "alpha[2]: ";
    // cin >> param.a2;

    // cout << "beta[1]: ";
    // cin >> param.b1;

    // cout << "tail index: ";
    // cin >> alpha;
    
    cout << "Number of evaluation points: ";
    cin >> n;

    double ang_max = atan(1/param.a1);
    vec eigenfun(n);
    vec angles(n);
    vec tangents(n);
    vec cosines(n);
    mat P(n, n);
    cx_vec eigenval(n);
    cx_mat eigenmat(n, n);
    
    
    double d = ang_max/(n + 1);
    for (unsigned i = 0; i < n; i++) {
	angles[i] = d * (i + 1);
	tangents[i] = tan(angles[i]);
	cosines[i] = cos(angles[i]);
    }

    double c1 = pow(pow(param.a2, 2) + pow(param.b1, 2), 0.5);
    double ang1 = atan(param.b1 / param.a2);
    for (unsigned i = 0; i < n; i++) {
	for (unsigned j = 0; j < n; j++) {
	    double t1 = sin(angles[i] + ang1);
	    double t2 = sqrt(tangents[i] * param.a2 + param.b1);
	    double t3 = (param.a2 * tangents[i] + param.b1);
	    t3 *= tangents[j] / (1 - tangents[j] * param.a1) / 2;
	    t3 = exp(-t3);

	    double t5 = sqrt(2 * M_PI * tangents[j]);
	    double t6 = pow(cosines[j], alpha + 2);
	    double t7 = pow(1 - tangents[j] * param.a1, alpha + 1.5);
	    P(i, j) = pow(c1 * t1, alpha) * t2 * t3 / t5 / t6 / t7;
	}
    }
    cout << P << endl;
    eig_gen(eigenval, eigenmat, P);
    cout << "eigenvalues" << endl;
    cout << eigenval << endl;

    cout << "eigenvectors" << endl;
    cout << eigenmat << endl;


    return 0;
}
