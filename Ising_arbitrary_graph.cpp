#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <string>
#include <random>

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
long M = 10000000; // number of steps
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

void construct_graph()
{
    // if two spins are in low-E state, connect them with probability p
    // otherwise they are disconnected
    for (int i = 0; i < E; i++)
		if (unirnd(eng) < std::exp(-2 * bond2[i].J
			* spins[bond2[i].x] * spins[bond2[i].y] / T))
			bond2[i].conn = 0;
		else
			bond2[i].conn = 1;
	// if four spins are in low-E state, connect them with probability p
	// otherwise they are disconnected
	for (int i = 0; i < G; i++)
		if (unirnd(eng) < std::exp(-2 * bond4[i].J
			* spins[bond4[i].x] * spins[bond4[i].y]
			* spins[bond4[i].z] * spins[bond4[i].w] / T))
			bond4[i].conn = 0;
		else
			bond4[i].conn = 1;
}

void BFS(int s0)
{
	// Breadth first search to find a cluster starting from s0
	int *queue = new int[V];
	queue[0] = s0;
	n++;
	cluster[s0] = n;

	int *head = queue, *tail = queue;
	while (tail - head >= 0)
	{
		int s = *head;
		head++;
		for (int i = 0; i < deg2[s]; i++)
			if (A2[s][i]->conn)
			{
				if (cluster[A2[s][i]->x] == 0)
				{
					tail++;
					*tail = A2[s][i]->x;
					cluster[A2[s][i]->x] = n;
				}
				if (cluster[A2[s][i]->y] == 0)
				{
					tail++;
					*tail = A2[s][i]->y;
					cluster[A2[s][i]->y] = n;
				}
			}
		for (int i = 0; i < deg4[s]; i++)
			if (A4[s][i]->conn)
			{
				if (cluster[A4[s][i]->x] == 0)
				{
					tail++;
					*tail = A4[s][i]->x;
					cluster[A4[s][i]->x] = n;
				}
				if (cluster[A4[s][i]->y] == 0)
				{
					tail++;
					*tail = A4[s][i]->y;
					cluster[A4[s][i]->y] = n;
				}
				if (cluster[A4[s][i]->z] == 0)
				{
					tail++;
					*tail = A4[s][i]->z;
					cluster[A4[s][i]->z] = n;
				}
				if (cluster[A4[s][i]->w] == 0)
				{
					tail++;
					*tail = A4[s][i]->w;
					cluster[A4[s][i]->w] = n;
				}
			}
	}

	delete[] queue;
}

void find_clusters()
{
	// initialize
	n = 0;
	for (int i = 0; i < V; i++)
		cluster[i] = 0; // 0 means unused

	for (int i = 0; i < V; i++)
		if (cluster[i] == 0)
			BFS(i);

	for (int i = 1; i <= V; i++)
		cluster_size[i] = 0;
	for (int i = 0; i < V; i++)
		cluster_size[cluster[i]] += spins[i];
}

void flip_clusters()
{
    bool *flip = new bool[n + 1];
    for (int i = 1; i <= n; i++)
        if (unirnd(eng) < 0.5)
            flip[i] = 1;
        else
            flip[i] = 0;
    for (int i = 0; i < V; i++)
        if (flip[cluster[i]])
            spins[i] = -spins[i];
	delete[] flip;
}

void measure(double &mag2, double &mag4)
{
	double temp;
    
	mag2 = 0;
	mag4 = 0;
	for (int k = 1; k <= n; k++)
	{
		temp = (double)cluster_size[k] * cluster_size[k];
		mag2 += temp;
		mag4 -= 2 * temp * temp;
	}
	mag4 += 3 * mag2 * mag2;
}

void export_data(double *mag2, double *mag4)
{
    using namespace std;
	ofstream file;
    file.open("mag2.bin", ios::out | ios::binary);
    file.write((char *)mag2, sizeof(double) * M);
	file.close();
    file.open("mag4.bin", ios::out | ios::binary);
    file.write((char *)mag4, sizeof(double) * M);
	file.close();
}

double mean(double *arr, long size)
{
	double v = 0;
	for (long i = 0; i < size; i++)
		v += arr[i];
	v /= size;
	return(v);
}

double std_dev(double *arr, long size)
{
	double v = 0, u = mean(arr, size);
	for (long i = 0; i < size; i++)
		v += (arr[i] - u) * (arr[i] - u);
	v /= size;
	return(std::sqrt(v));
}

void auto_correlation(double *mag2, double *mag4)
{
	// This function will destroy the values inside mag2 and mag4
	// thus must be placed at the end of the program.
	using namespace std;
	int i, l = 0;
	long step = 2;
	while (step <= M)
	{
		l++;
		step *= 2;
	}
	step /= 2;
	double *cor1 = new double[l], *cor2 = new double[l];
	mag2 = mag2 + M - step;
	mag4 = mag4 + M - step;
	for (i = 0; i < l; i++)
	{
		cor1[i] = std_dev(mag2, step) / sqrt(step - 1);
		cor2[i] = std_dev(mag4, step) / sqrt(step - 1);
		step /= 2;
		for (long j = 0; j < step; j++)
		{
			mag2[j] = (mag2[2 * j] + mag2[2 * j + 1]) / 2;
			mag4[j] = (mag4[2 * j] + mag4[2 * j + 1]) / 2;
		}
	}
	ofstream file("auto_cor_" + str + ".txt");
	for (i = 0; i < l; i++)
		file << cor1[i] << ' ';
	file << endl;
	for (i = 0; i < l; i++)
		file << cor2[i] << ' ';
	file << endl;
	file.close();
	delete[] cor1;
	delete[] cor2;
}

void data_analysis(double *mag2, double *mag4)
{
    using namespace std;
    double mag2_av = 0, mag4_av = 0;
    for (long i = M / 10; i < M; i++) // skip initial 1/10 steps to thermalize
    {
        mag2_av += mag2[i];
        mag4_av += mag4[i];
    }
    mag2_av /= (M - M / 10);
    mag4_av /= (M - M / 10);
	ofstream file(str + ".txt");
    file << "Number of spins is: " << V << endl;
    file << "At temperature T = " << T << endl;
    file << "Binder cumulant U = " << 1 - mag4_av / 3 / mag2_av / mag2_av
        << endl;
	file.close();
}

void monte_carlo()
{
    double *mag2 = new double[M], *mag4 = new double[M];
    lattice_structure();
    initialize_spins();
    for (long i = 0; i < M; i++)
    {
        construct_graph();
        find_clusters();
        measure(mag2[i], mag4[i]);
        flip_clusters();
    }
    export_data(mag2, mag4);
    data_analysis(mag2, mag4);
	auto_correlation(mag2, mag4);
	delete[] mag2;
    delete[] mag4;
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