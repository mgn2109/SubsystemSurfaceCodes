#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <string>
#include <random>

#define L 11 // L is odd, consider L*L data qubits
double q;
double nlogp; // p parameter given as -log p
double T; // temperature
long M = 1000000; // number of steps
int seed = 0;
std::string str;

int *lattice[L - 1][(L + 1) / 2][2];
int spins[L * L]; // spin number is less than L^2
int coup[L - 1][L + 1]; // ferro (1) or antiferro (-1)
bool conn[L - 1][L + 1]; // cluster connectivity
long num; // number of spins

std::mt19937_64 eng;
std::uniform_real_distribution<double> unirnd(0.0, 1.0);

void lattice_structure()
{
	int i, j;
	num = 0;
	// spins
	for (i = 0; i < L - 1; i++)
	{
		for (j = 0; j < (L + 1) / 2; j++)
			if (unirnd(eng) < q)
			{
				lattice[i][j][0] = &spins[num];
				lattice[i][j][1] = &spins[num];
				num++;
			}
			else
			{
				lattice[i][j][0] = &spins[num];
				num++;
				lattice[i][j][1] = &spins[num];
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

void construct_graph(double p0)
{
    int i, j;
    // if two adjacent spins are parallel, connect them with probability p
    // otherwise they are disconnected
    for (i = 0; i < (L - 1) / 2; i++)
    {
        for (j = 0; j <= (L - 1) / 2; j++)
        {
            if (((*lattice[2 * i][j][1]) * (*lattice[2 * i + 1][j][0])
                == coup[2 * i][2 * j]) && unirnd(eng) < p0)
                conn[2 * i][2 * j] = 1;
            else
                conn[2 * i][2 * j] = 0;
            if (((*lattice[2 * i + 1][j][1])
                * (*lattice[2 * i][(j + 1) % ((L + 1) / 2)][0])
                == coup[2 * i][2 * j + 1]) && unirnd(eng) < p0)
                conn[2 * i][2 * j + 1] = 1;
            else
                conn[2 * i][2 * j + 1] = 0;
            
            if (((*lattice[(2 * i + 2) % (L - 1)][j][1])
                * (*lattice[2 * i + 1][j][0])
                == coup[2 * i + 1][2 * j]) && unirnd(eng) < p0)
                conn[2 * i + 1][2 * j] = 1;
            else
                conn[2 * i + 1][2 * j] = 0;
            if (((*lattice[2 * i + 1][j][1]) * 
            (*lattice[(2 * i + 2) % (L - 1)][(j + 1) % ((L + 1) / 2)][0])
                == coup[2 * i + 1][2 * j + 1]) && unirnd(eng) < p0)
                conn[2 * i + 1][2 * j + 1] = 1;
            else
                conn[2 * i + 1][2 * j + 1] = 0;
        }
    }
}

int find_root(long x, long *root)
{
    if (root[x] == x)
        return(x);
    else
    {
        root[x] = find_root(root[x], root);
        return(root[x]);
    }
}

void find_clusters(long &nl, long *label, long *size)
{
    int i, j, col = (L + 1) / 2;
	long k;
    nl = 0; // current label;
    long *root = new long[L * L];
	for (k = 0; k < L * L; k++)
		root[k] = k;
    // first 2 rows of sites
    for (j = 0; j < col; j++)
        if (lattice[0][j][0] == lattice[0][j][1])
            label[lattice[0][j][0] - spins] = nl++;
        else
        {
            label[lattice[0][j][0] - spins] = nl++;
            label[lattice[0][j][1] - spins] = nl++;
        }
    for (j = 0; j < col; j++)
        if (lattice[1][j][0] == lattice[1][j][1])
        {
            if (conn[0][2 * j] && conn[0][2 * j + 1])
            {
                label[lattice[1][j][0] - spins] =
                    find_root(label[lattice[0][j][1] - spins], root);
                if (label[lattice[1][j][0] - spins] !=
                    find_root(label[lattice[0][(j + 1) % col][0] - spins],
                        root))
                    root[root[label[lattice[0][(j + 1) % col][0] - spins]]] =
                        label[lattice[1][j][0] - spins];
                continue;
            }
            if (conn[0][2 * j])
            {
                label[lattice[1][j][0] - spins] =
                    find_root(label[lattice[0][j][1] - spins], root);
                continue;
            }
            if (conn[0][2 * j + 1])
            {
                label[lattice[1][j][0] - spins] =
                    find_root(label[lattice[0][(j + 1) % col][0] - spins],
                        root);
                continue;
            }
            label[lattice[1][j][0] - spins] = nl++;
        }
        else
        {
            if (conn[0][2 * j])
                label[lattice[1][j][0] - spins] =
                    find_root(label[lattice[0][j][1] - spins], root);
            else
                label[lattice[1][j][0] - spins] = nl++;
            if (conn[0][2 * j + 1])
                label[lattice[1][j][1] - spins] =
                    find_root(label[lattice[0][(j + 1) % col][0] - spins],
                        root);
            else
                label[lattice[1][j][1] - spins] = nl++;
        }
    // remaining rows
    for (i = 1; i < (L - 1) / 2; i++)
    {
        // even rows
        for (j = 0; j < col; j++)
            if (lattice[2 * i][j][0] == lattice[2 * i][j][1])
            {
                if (conn[2 * i - 1][(2 * j + L) % (L + 1)]
                    && conn[2 * i - 1][2 * j])
                {
                    label[lattice[2 * i][j][0] - spins] =
                        find_root(label[lattice[2 * i - 1]
                            [(j - 1 + col) % col][1] - spins], root);
                    if (label[lattice[2 * i][j][0] - spins] !=
                        find_root(label[lattice[2 * i - 1][j][0] - spins],
                                    root))
                        root[root[label[lattice[2 * i - 1][j][0] - spins]]] =
                            label[lattice[2 * i][j][0] - spins];
                    continue;
                }
                if (conn[2 * i - 1][(2 * j + L) % (L + 1)])
                {
                    label[lattice[2 * i][j][0] - spins] =
                        find_root(label[lattice[2 * i - 1]
                            [(j - 1 + col) % col][1] - spins], root);
                    continue;
                }
                if (conn[2 * i - 1][2 * j])
                {
                    label[lattice[2 * i][j][0] - spins] =
                        find_root(label[lattice[2 * i - 1][j][0] - spins],
                                    root);
                    continue;
                }
                label[lattice[2 * i][j][0] - spins] = nl++;
            }
            else
            {
                if (conn[2 * i - 1][(2 * j + L) % (L + 1)])
                    label[lattice[2 * i][j][0] - spins] =
                        find_root(label[lattice[2 * i - 1]
                            [(j - 1 + col) % col][1] - spins], root);
                else
                    label[lattice[2 * i][j][0] - spins] = nl++;
                if (conn[2 * i - 1][2 * j])
                    label[lattice[2 * i][j][1] - spins] =
                        find_root(label[lattice[2 * i - 1][j][0] - spins],
                                    root);
                else
                    label[lattice[2 * i][j][1] - spins] = nl++;
            }
        // odd rows
        for (j = 0; j < col; j++)
            if (lattice[2 * i + 1][j][0] == lattice[2 * i + 1][j][1])
            {
                if (conn[2 * i][2 * j] && conn[2 * i][2 * j + 1])
                {
                    label[lattice[2 * i + 1][j][0] - spins] =
                        find_root(label[lattice[2 * i][j][1] - spins], root);
                    if (label[lattice[2 * i + 1][j][0] - spins] !=
                        find_root(label[lattice[2 * i][(j + 1) % col][0]
                            - spins], root))
                        root[root[label[lattice[2 * i][(j + 1) % col][0]
                            - spins]]]
                                = label[lattice[2 * i + 1][j][0] - spins];
                    continue;
                }
                if (conn[2 * i][2 * j])
                {
                    label[lattice[2 * i + 1][j][0] - spins] =
                        find_root(label[lattice[2 * i][j][1] - spins], root);
                    continue;
                }
                if (conn[2 * i][2 * j + 1])
                {
                    label[lattice[2 * i + 1][j][0] - spins] =
                        find_root(label[lattice[2 * i][(j + 1) % col][0]
                            - spins], root);
                    continue;
                }
                label[lattice[2 * i + 1][j][0] - spins] = nl++;
            }
            else
            {
                if (conn[2 * i][2 * j])
                    label[lattice[2 * i + 1][j][0] - spins] =
                        find_root(label[lattice[2 * i][j][1] - spins], root);
                else
                    label[lattice[2 * i + 1][j][0] - spins] = nl++;
                if (conn[2 * i][2 * j + 1])
                    label[lattice[2 * i + 1][j][1] - spins] =
                        find_root(label[lattice[2 * i][(j + 1) % col][0]
                            - spins], root);
                else
                    label[lattice[2 * i + 1][j][1] - spins] = nl++;
            }
    }
    // connect first and last rows
	for (j = 0; j < col; j++)
	{
		if (conn[L - 2][(2 * j + L) % (L + 1)])
			if (find_root(label[lattice[0][j][0] - spins], root) !=
					find_root(label[lattice[L - 2][(j - 1 + col) % col][1]
					- spins], root))
				root[root[label[lattice[L - 2][(j - 1 + col) % col][1]
					- spins]]] = root[label[lattice[0][j][0] - spins]];
		if (conn[L - 2][2 * j])
			if (find_root(label[lattice[0][j][1] - spins], root) !=
					find_root(label[lattice[L - 2][j][0] - spins], root))
				root[root[label[lattice[L - 2][j][0] - spins]]] =
					root[label[lattice[0][j][1] - spins]];
	}
    
    for (k = 0; k < num; k++)
        label[k] = find_root(label[k], root);
    for (k = 0; k < nl; k++)
        root[k] = -1;
    nl = 0;
    for (k = 0; k < num; k++)
        if (root[label[k]] == -1)
        {
            root[label[k]] = nl;
            label[k] = nl++;
            size[label[k]] = spins[k];
        }
        else
        {
            label[k] = root[label[k]];
            size[label[k]]+= spins[k];
        }
    delete[] root;
}

void flip_clusters(long nl, long *label)
{
    long i;
    bool flip[L * L];
    for (i = 0; i < nl; i++)
        if (unirnd(eng) < 0.5)
            flip[i] = 1;
        else
            flip[i] = 0;
    for (i = 0; i < num; i++)
        if (flip[label[i]])
            spins[i] = -spins[i];
}

long long tonumber(long nl, long lab1, long lab2)
{
    using namespace std;
    return((long long)nl * min(lab1, lab2) + max(lab1, lab2));
}

void measure(long nl, long *size, double &mag2, double &mag4)
{
    long k;
	double temp;
    
	mag2 = 0;
	mag4 = 0;
	for (k = 0; k < nl; k++)
	{
		temp = (double)size[k] * size[k];
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

void data_analysis(double *mag2, double *mag4)
{
    using namespace std;
    long i;
    double mag2_av = 0, mag4_av = 0;
    for (i = M / 10; i < M; i++) // skip initial 1/10 steps to thermalize
    {
        mag2_av += mag2[i];
        mag4_av += mag4[i];
    }
    mag2_av /= (M - M / 10);
    mag4_av /= (M - M / 10);
	ofstream file(str);
    file << "Number of spins is: " << num << endl;
    file << "At temperature T = " << T << endl;
    file << "Binder cumulant U = " << 1 - mag4_av / 3 / mag2_av / mag2_av
        << endl;
	file.close();
}

void monte_carlo()
{
    long i;
    long nl; // number of labels
	long *label = new long[L * L], *size = new long[L * L];
    double *mag2 = new double[M], *mag4 = new double[M];
    lattice_structure();
    initialize_spins();
    for (i = 0; i < M; i++)
    {
        construct_graph(1 - std::exp(-2.0 / T));
        find_clusters(nl, label, size);
        measure(nl, size, mag2[i], mag4[i]);
        flip_clusters(nl, label);
    }
    export_data(mag2, mag4);
    data_analysis(mag2, mag4);
    delete[] label;
    delete[] size;
	delete[] mag2;
    delete[] mag4;
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
		M = atol(argv[3]);
	if (argc > 4)
		seed = atoi(argv[4]);
	str = "L" + to_string(L) + "_q" + argv[1] + "_nlogp" + argv[2]
		+ "_M" + to_string(M) + "_seed" + to_string(seed) + ".txt";
	eng.seed(seed);
    monte_carlo();
    return(0);
}