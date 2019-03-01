#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <string>
#include <random>
#include <iomanip>

struct Bond2Type
{
	int x, y;
	double J; // ferro (+) or antiferro (-) for each 2-body bond
	bool conn; // connectivity for 2-body bond in a cluster
};
struct Bond4Type
{
	int x, y, z, w;
	double J; // ferro (+) or antiferro (-) for each 4-body bond
	bool conn; // connectivity for 4-body bond in a cluster
};

int V, E, G; // number of spins, 2-body interaction, 4-body interaction
int *spins; // spin up (1) or down (-1) for each spin
Bond2Type *bond2; // list of 2-body bonds
Bond4Type *bond4; // list of 4-body bonds
int *deg2; // degree of each node for 2-body interaction
int *deg4; // degree of each node for 4-body interaction
Bond2Type ***A2; // adjacency list for 2-body, pointing to entries of bond2
Bond4Type ***A4; // adjacency list for 4-body, pointing to entries of bond4
int n; // number of clusters
int *cluster; // cluster label for each spin, ranging from 1 to n
int *cluster_size; // size of each cluster, index from 1 to n

double T = 1; // temperature set to 1
long M = 100000000; // number of steps
long long seed = 0; // random seed

std::string str = "output_for_", str_graph = "input.txt";

std::mt19937_64 eng;
std::uniform_real_distribution<double> unirnd(0.0, 1.0);

void lattice_structure()
{
	using namespace std;

	// load graph structure
	/* File in the following format:
	V E G
	i j J // E rows of 2-body bonds
	i j k l J // G rows of 4-body bonds
	*/
	ifstream file(str_graph);
	file >> V >> E >> G;
	spins = new int[V];
	cluster = new int[V];
	cluster_size = new int[V + 1]; // at most V clusters
	bond2 = new Bond2Type[E];
	bond4 = new Bond4Type[G];
	deg2 = new int[V];
	deg4 = new int[V];
	for (int i = 0; i < V; i++)
	{
		deg2[i] = 0;
		deg4[i] = 0;
	}

	// load bonds and compute deg2 and deg4
	for (int i = 0; i < E; i++)
	{
		file >> bond2[i].x >> bond2[i].y >> bond2[i].J;
		deg2[bond2[i].x]++;
		deg2[bond2[i].y]++;
	}
	for (int i = 0; i < G; i++)
	{
		file >> bond4[i].x >> bond4[i].y >> bond4[i].z >> bond4[i].w >> bond4[i].J;
		deg4[bond4[i].x]++;
		deg4[bond4[i].y]++;
		deg4[bond4[i].z]++;
		deg4[bond4[i].w]++;
	}

	// allocate memory for A2 and A4
	A2 = new Bond2Type**[V];
	for (int i = 0; i < V; i++)
		A2[i] = new Bond2Type*[deg2[i]];
	A4 = new Bond4Type**[V];
	for (int i = 0; i < V; i++)
		A4[i] = new Bond4Type*[deg4[i]];
	// convert bond list into adjacency list form
	for (int i = 0; i < V; i++)
	{
		deg2[i] = 0;
		deg4[i] = 0;
	}
	for (int i = 0; i < E; i++)
	{
		A2[bond2[i].x][deg2[bond2[i].x]] = bond2 + i;
		A2[bond2[i].y][deg2[bond2[i].y]] = bond2 + i;
		deg2[bond2[i].x]++;
		deg2[bond2[i].y]++;
	}
	for (int i = 0; i < G; i++)
	{
		A4[bond4[i].x][deg4[bond4[i].x]] = bond4 + i;
		A4[bond4[i].y][deg4[bond4[i].y]] = bond4 + i;
		A4[bond4[i].z][deg4[bond4[i].z]] = bond4 + i;
		A4[bond4[i].w][deg4[bond4[i].w]] = bond4 + i;
		deg4[bond4[i].x]++;
		deg4[bond4[i].y]++;
		deg4[bond4[i].z]++;
		deg4[bond4[i].w]++;
	}
}

void dealloc_data()
{
	delete[] spins;
	delete[] cluster;
	delete[] cluster_size;
	delete[] bond2;
	delete[] bond4;
	delete[] deg2;
	delete[] deg4;
	for (int i = 0; i < V; i++)
	{
		delete[] A2[i];
		delete[] A4[i];
	}
	delete[] A2;
	delete[] A4;
}

void initialize_spins()
{
    for (int i = 0; i < V; i++)
        spins[i] = 1;
}

void measure(double &En, double &mag)
{
	int i, j;
	long k;
	En = 0;
	for (int i = 0; i < E; i++)
		En -= bond2[i].J * spins[bond2[i].x] * spins[bond2[i].y];
	for (int i = 0; i < G; i++)
		En -= bond4[i].J * spins[bond4[i].x] * spins[bond4[i].y]
			* spins[bond4[i].z] * spins[bond4[i].w];

	mag = 0;
	for (k = 0; k < V; k++)
		mag += spins[k];
}

void diff(long id, double &dE, double &dmag)
{
	dE = 0;
	for (int i = 0; i < deg2[id]; i++)
		dE += 2 * A2[id][i]->J * spins[A2[id][i]->x] * spins[A2[id][i]->y];
	for (int i = 0; i < deg4[id]; i++)
		dE += 2 * A4[id][i]->J * spins[A4[id][i]->x] * spins[A4[id][i]->y]
		* spins[A4[id][i]->z] * spins[A4[id][i]->w];

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
	ofstream file(str + ".txt", ios::app);
	file << "Number of spins is: " << V << endl;
	file << "Current number of steps = " << step << endl;
	file << "At temperature T = " << T << endl;
	file << "Average energy E = " << E_sum << endl;
	file << "Specific heat C = " << (E2_sum - E_sum * E_sum) / T / T * V
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
	double dE, dmag, En, mag, E_sum = 0, E2_sum = 0,
		mag_sum = 0, mag2_sum = 0, mag4_sum = 0;
	lattice_structure();
	initialize_spins();
	measure(En, mag);
	std::uniform_int_distribution<long> randint(0, V - 1);
	for (long long i = 1; i < M; i++)
	{
		id = randint(eng);
		diff(id, dE, dmag);
		if (unirnd(eng) < std::exp(-dE / T))
		{
			spins[id] = -spins[id];
			En += dE;
			mag += dmag;
		}
		if (i >= M / 10)
		{
			double temp;
			temp = En / V;
			E_sum += temp;
			E2_sum += temp * temp;
			temp = std::abs(mag) / V;
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
	dealloc_data();
}

int main(int argc, char *argv[])
{
	using namespace std;
	if (argc > 1)
		M = atol(argv[1]);
	if (argc > 2)
		seed = atoi(argv[2]);
	if (argc > 3)
		str_graph.assign(argv[3]);
	
	str = str + str_graph + "_M" + to_string(M) + "_seed" + to_string(seed);
	eng.seed(seed);
	monte_carlo();
    return(0);
}