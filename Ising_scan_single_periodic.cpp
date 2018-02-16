#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <string>
#include <random>

#define L 41 // L is odd, consider L*L data qubits
double q;
double nlogp; // p parameter given as -log p
double T; // temperature
long long M = 5000000000; // number of steps
int seed = 0;
std::string str;

int *lattice[L - 1][(L + 1) / 2][2];
int spins[L * L]; // spin number is less than L^2
int coup[L - 1][L + 1]; // ferro (1) or antiferro (-1)
int location[L * L][3]; // row and column indices for each spin
// and left (0) or right (1) or both (-1)
int num; // number of spins

std::default_random_engine eng;
std::uniform_real_distribution<double> unirnd(0.0, 1.0);

void lattice_structure()
{
	int i, j;
	num = 0;
	// spins
	for (i = 0; i < (L - 1) / 2; i++)
	{
		for (j = 0; j < (L + 1) / 2; j++)
			if (unirnd(eng) < q)
			{
				lattice[2 * i][j][0] = &spins[num];
				lattice[2 * i][j][1] = &spins[num];
				location[num][0] = 2 * i;
				location[num][1] = j;
				location[num][2] = -1;
				num++;
			}
			else
			{
				lattice[2 * i][j][0] = &spins[num];
				location[num][0] = 2 * i;
				location[num][1] = j;
				location[num][2] = 0;
				num++;
				lattice[2 * i][j][1] = &spins[num];
				location[num][0] = 2 * i;
				location[num][1] = j;
				location[num][2] = 1;
				num++;
			}
		for (j = 0; j < (L + 1) / 2; j++)
			if (unirnd(eng) < q)
			{
				lattice[2 * i + 1][j][0] = &spins[num];
				lattice[2 * i + 1][j][1] = &spins[num];
				location[num][0] = 2 * i + 1;
				location[num][1] = j;
				location[num][2] = -1;
				num++;
			}
			else
			{
				lattice[2 * i + 1][j][0] = &spins[num];
				location[num][0] = 2 * i + 1;
				location[num][1] = j;
				location[num][2] = 0;
				num++;
				lattice[2 * i + 1][j][1] = &spins[num];
				location[num][0] = 2 * i + 1;
				location[num][1] = j;
				location[num][2] = 1;
				num++;
			}
	}

	// bonds. 1 for ferro, -1 for antiferro
	double p = std::exp(-nlogp);
	for (i = 0; i < L - 1; i++)
		for (j = 0; j < L + 1; j++)
			if (unirnd(eng) < p)
				coup[i][j] = -1;
			else
				coup[i][j] = 1;
}

void initialize_spins()
{
	for (int i = 0; i < num; i++)
		spins[i] = 1;
}

void measure(double &E, double &mag)
{
	int i, j;
	long k;
	E = 0;
	E += (*lattice[0][0][1]) * (*lattice[1][0][0]) * coup[0][0]
		+ (*lattice[0][0][0]) * (*lattice[1][(L - 1) / 2][1]) * coup[0][L];
	for (j = 1; j <= (L - 1) / 2; j++)
	{
		E += (*lattice[1][j - 1][1]) * (*lattice[0][j][0])
			* coup[0][2 * j - 1];
		E += (*lattice[0][j][1]) * (*lattice[1][j][0]) * coup[0][2 * j];
		E += (*lattice[L - 2][j - 1][1]) * (*lattice[0][j][0])
			* coup[L - 2][2 * j - 1];
		E += (*lattice[0][j][1]) * (*lattice[L - 2][j][0])
			* coup[L - 2][2 * j];
	}
	for (i = 1; i < (L - 1) / 2; i++)
	{
		E += (*lattice[2 * i][0][1]) * (*lattice[2 * i - 1][0][0])
			* coup[2 * i - 1][0];
		E += (*lattice[2 * i][0][0])
			* (*lattice[2 * i - 1][(L - 1) / 2][1]) * coup[2 * i - 1][L];
		for (j = 1; j <= (L - 1) / 2; j++)
		{
			E += (*lattice[2 * i - 1][j - 1][1]) * (*lattice[2 * i][j][0])
				* coup[2 * i - 1][2 * j - 1];
			E += (*lattice[2 * i][j][1]) * (*lattice[2 * i - 1][j][0])
				* coup[2 * i - 1][2 * j];
		}
		E += (*lattice[2 * i][0][1]) * (*lattice[2 * i + 1][0][0])
			* coup[2 * i][0];
		E += (*lattice[2 * i][0][0])
			* (*lattice[2 * i + 1][(L - 1) / 2][1]) * coup[2 * i][L];
		for (j = 1; j <= (L - 1) / 2; j++)
		{
			E += (*lattice[2 * i + 1][j - 1][1]) * (*lattice[2 * i][j][0])
				* coup[2 * i][2 * j - 1];
			E += (*lattice[2 * i][j][1]) * (*lattice[2 * i + 1][j][0])
				* coup[2 * i][2 * j];
		}
	}
	E = -E;

	mag = 0;
	for (k = 0; k < num; k++)
		mag += spins[k];
}

void diff(long id, double &dE, double &dmag)
{
	int i = location[id][0], j = location[id][1], k = location[id][2],
		lasti = (i - 1 + L - 1) % (L - 1), nexti = (i + 1) % (L - 1),
		lastj = (j - 1 + (L + 1) / 2) % ((L + 1) / 2),
		nextj = (j + 1) % ((L + 1) / 2);
	bool left = (k == 0) || (k == -1), right = (k == 1) || (k == -1);
	dE = 0;
	if (i % 2 == 0)
	{
		if (left)
		{
			dE += coup[lasti][2 * lastj + 1] * (*lattice[i][j][0])
				* (*lattice[lasti][lastj][1]);
			dE += coup[i][2 * lastj + 1] * (*lattice[i][j][0])
				* (*lattice[nexti][lastj][1]);
		}
		if (right)
		{
			dE += coup[lasti][2 * j] * (*lattice[i][j][1])
				* (*lattice[lasti][j][0]);
			dE += coup[i][2 * j] * (*lattice[i][j][1])
				* (*lattice[nexti][j][0]);
		}
	}
	else
	{
		if (left)
		{
			dE += coup[lasti][2 * j] * (*lattice[i][j][0])
				* (*lattice[lasti][j][1]);
			dE += coup[i][2 * j] * (*lattice[i][j][0])
				* (*lattice[nexti][j][1]);
		}
		if (right)
		{
			dE += coup[lasti][2 * j + 1] * (*lattice[i][j][1])
				* (*lattice[lasti][nextj][0]);
			dE += coup[i][2 * j + 1] * (*lattice[i][j][1])
				* (*lattice[nexti][nextj][0]);
		}
	}
	dE *= 2;

	dmag = -2 * spins[id];
}

void show_average(long long step, double E_sum, double E2_sum,
	double mag_sum, double mag2_sum, double mag4_sum)
{
	using namespace std;
	E_sum /= step;
	E2_sum /= step;
	mag_sum /= step;
	mag2_sum /= step;
	mag4_sum /= step;
	ofstream file(str, ios::app);
	file << "q = " << q << ", p = " << exp(-nlogp) << ", L = " << L << endl;
	file << "Number of spins is: " << num << endl;
	file << "Current number of steps = " << step << endl;
	file << "At temperature T = " << T << endl;
	file << "Average energy E = " << E_sum << endl;
	file << "Specific heat C = " << (E2_sum - E_sum * E_sum) / T / T * num
		<< endl;
	file << "Average magnetization M = " << mag_sum << endl;
	file << "Binder cumulant U = " << 1 - mag4_sum / 3 / mag2_sum / mag2_sum
		<< endl;
	file.close();
}

void monte_carlo()
{
	long id;
    long long k = 1;
	double dE, dmag, E, mag, E_sum = 0, E2_sum = 0,
		mag_sum = 0, mag2_sum = 0, mag4_sum = 0;
	lattice_structure();
	initialize_spins();
	measure(E, mag);
	std::uniform_int_distribution<long> randint(0, num - 1);
	for (long long i = 1; i < M; i++)
	{
		id = randint(eng);
		diff(id, dE, dmag);
		if (unirnd(eng) < std::exp(-dE / T))
		{
			spins[id] = -spins[id];
			E += dE;
			mag += dmag;
		}
		if (i >= M / 10)
		{
			double temp;
			temp = E / num;
			E_sum += temp;
			E2_sum += temp * temp;
			temp = std::abs(mag) / num;
			mag_sum += temp;
			temp *= temp;
			mag2_sum += temp;
			mag4_sum += temp * temp;
			if (i + 1 - M / 10 == k)
			{
				show_average(k, E_sum, E2_sum, mag_sum, mag2_sum, mag4_sum);
				k *= 2;
			}
		}
	}
}

int main(int argc, char *argv[])
{
	using namespace std;
	if (argc < 3)
	{
		cout << "Not enough parameters!" << endl;
		return(-1);
	}
	q = atof(argv[1]);
	nlogp = atof(argv[2]);
	if (nlogp <= log(2.0))
	{
		cout << "Error! -log p must be greater than log 2." << endl;
		return(-1);
	}
	T = 2.0 / (log(1 - exp(-nlogp)) + nlogp);
	if (argc > 3)
		M = atoll(argv[3]);
	if (argc > 4)
		seed = atoi(argv[4]);
	str = "L" + to_string(L) + "_q" + argv[1] + "_nlogp" + argv[2]
		+ "_M" + to_string(M) + "_seed" + to_string(seed) + ".txt";
	eng.seed(seed);
	monte_carlo();
	return(0);
}