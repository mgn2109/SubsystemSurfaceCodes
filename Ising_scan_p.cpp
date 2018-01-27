#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <string>
#include <random>

#define L 101 // L is odd, consider L*L data qubits
double q;
double nlogp; // p parameter given as -log p
double T; // temperature
long M = 100000; // number of steps
int seed = 0;
std::string str;

int *lattice[L - 1][(L + 1) / 2][2];
int spins[L * L]; // spin number is less than L^2
int coup[L - 2][L]; // ferro (1) or antiferro (-1)
bool conn[L - 2][L]; // cluster connectivity
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
		lattice[2 * i][0][0] = &spins[num];
		lattice[2 * i][0][1] = &spins[num];
		num++;
		for (j = 1; j < (L + 1) / 2; j++)
			if (unirnd(eng) < q)
			{
				lattice[2 * i][j][0] = &spins[num];
				lattice[2 * i][j][1] = &spins[num];
				num++;
			}
			else
			{
				lattice[2 * i][j][0] = &spins[num];
				num++;
				lattice[2 * i][j][1] = &spins[num];
				num++;
			}
		for (j = 0; j < (L - 1) / 2; j++)
			if (unirnd(eng) < q)
			{
				lattice[2 * i + 1][j][0] = &spins[num];
				lattice[2 * i + 1][j][1] = &spins[num];
				num++;
			}
			else
			{
				lattice[2 * i + 1][j][0] = &spins[num];
				num++;
				lattice[2 * i + 1][j][1] = &spins[num];
				num++;
			}
		lattice[2 * i + 1][(L - 1) / 2][0] = &spins[num];
		lattice[2 * i + 1][(L - 1) / 2][1] = &spins[num];
		num++;
	}

	// bonds. 1 for ferro, -1 for antiferro
	double p = std::exp(-nlogp);
	for (i = 0; i < L - 2; i++)
		for (j = 0; j < L; j++)
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
    if (((*lattice[0][0][1]) * (*lattice[1][0][0]) == coup[0][0])
		&& unirnd(eng) < p0)
        conn[0][0] = 1;
    else
        conn[0][0] = 0;
    for (j = 1; j <= (L - 1) / 2; j++)
    {
        if (((*lattice[1][j - 1][1]) * (*lattice[0][j][0]) ==
			coup[0][2 * j - 1]) && unirnd(eng) < p0)
            conn[0][2 * j - 1] = 1;
        else
            conn[0][2 * j - 1] = 0;
        if (((*lattice[0][j][1]) * (*lattice[1][j][0]) == coup[0][2 * j])
			&& unirnd(eng) < p0)
            conn[0][2 * j] = 1;
        else
            conn[0][2 * j] = 0;
    }
    for (i = 1; i < (L - 1) / 2; i++)
    {
        if (((*lattice[2 * i][0][1]) * (*lattice[2 * i - 1][0][0])
			== coup[2 * i - 1][0]) && unirnd(eng) < p0)
            conn[2 * i - 1][0] = 1;
        else
            conn[2 * i - 1][0] = 0;
        for (j = 1; j <= (L - 1) / 2; j++)
        {
            if (((*lattice[2 * i - 1][j - 1][1]) * (*lattice[2 * i][j][0])
                == coup[2 * i - 1][2 * j - 1]) && unirnd(eng) < p0)
                conn[2 * i - 1][2 * j - 1] = 1;
            else
                conn[2 * i - 1][2 * j - 1] = 0;
            if (((*lattice[2 * i][j][1]) * (*lattice[2 * i - 1][j][0])
                == coup[2 * i - 1][2 * j]) && unirnd(eng) < p0)
                conn[2 * i - 1][2 * j] = 1;
            else
                conn[2 * i - 1][2 * j] = 0;
        }
        if (((*lattice[2 * i][0][1]) * (*lattice[2 * i + 1][0][0])
                == coup[2 * i][0]) && unirnd(eng) < p0)
            conn[2 * i][0] = 1;
        else
            conn[2 * i][0] = 0;
        for (j = 1; j <= (L - 1) / 2; j++)
        {
            if (((*lattice[2 * i + 1][j - 1][1]) * (*lattice[2 * i][j][0])
                == coup[2 * i][2 * j - 1]) && unirnd(eng) < p0)
                conn[2 * i][2 * j - 1] = 1;
            else
                conn[2 * i][2 * j - 1] = 0;
            if (((*lattice[2 * i][j][1]) * (*lattice[2 * i + 1][j][0])
                == coup[2 * i][2 * j]) && unirnd(eng) < p0)
                conn[2 * i][2 * j] = 1;
            else
                conn[2 * i][2 * j] = 0;
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

void find_clusters(long &nl, long *label, long *size, long *H)
{
    int i, j;
	long k;
    nl = 0; // current label;
    long *root = new long[L * L];
	for (k = 0; k < L * L; k++)
		root[k] = k;
    // first 2 rows of sites
    label[lattice[0][0][0] - spins] = nl++; // first column
    for (j = 1; j < (L + 1) / 2; j++)
        if (lattice[0][j][0] == lattice[0][j][1])
            label[lattice[0][j][0] - spins] = nl++;
        else
        {
            label[lattice[0][j][0] - spins] = nl++;
            label[lattice[0][j][1] - spins] = nl++;
        }
    for (j = 0; j < (L - 1) / 2; j++)
        if (lattice[1][j][0] == lattice[1][j][1])
        {
            if (conn[0][2 * j] && conn[0][2 * j + 1])
            {
                label[lattice[1][j][0] - spins] =
                    find_root(label[lattice[0][j][1] - spins], root);
                if (label[lattice[1][j][0] - spins] !=
                        find_root(label[lattice[0][j + 1][0] - spins], root))
                    root[root[label[lattice[0][j + 1][0] - spins]]] =
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
                    find_root(label[lattice[0][j + 1][0] - spins], root);
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
                    find_root(label[lattice[0][j + 1][0] - spins], root);
            else
                label[lattice[1][j][1] - spins] = nl++;
        }
    if (conn[0][L - 1]) // last column
        label[lattice[1][(L - 1) / 2][0] - spins] =
            find_root(label[lattice[0][(L - 1) / 2][1] - spins], root);
    else
        label[lattice[1][(L - 1) / 2][0] - spins] = nl++;
    // remaining rows
    for (i = 1; i < (L - 1) / 2; i++)
    {
        // even rows
        if (conn[2 * i - 1][0]) // first column
            label[lattice[2 * i][0][1] - spins] =
                find_root(label[lattice[2 * i - 1][0][0] - spins], root);
        else
            label[lattice[2 * i][0][1] - spins] = nl++;
        for (j = 1; j < (L + 1) / 2; j++)
            if (lattice[2 * i][j][0] == lattice[2 * i][j][1])
            {
                if (conn[2 * i - 1][2 * j - 1] && conn[2 * i - 1][2 * j])
                {
                    label[lattice[2 * i][j][0] - spins] =
                        find_root(label[lattice[2 * i - 1][j - 1][1]
                                    - spins], root);
                    if (label[lattice[2 * i][j][0] - spins] !=
                        find_root(label[lattice[2 * i - 1][j][0] - spins],
                                    root))
                        root[root[label[lattice[2 * i - 1][j][0] - spins]]] =
                            label[lattice[2 * i][j][0] - spins];
                    continue;
                }
                if (conn[2 * i - 1][2 * j - 1])
                {
                    label[lattice[2 * i][j][0] - spins] =
                        find_root(label[lattice[2 * i - 1][j - 1][1]
                                    - spins], root);
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
                if (conn[2 * i - 1][2 * j - 1])
                    label[lattice[2 * i][j][0] - spins] =
                        find_root(label[lattice[2 * i - 1][j - 1][1]
                                    - spins], root);
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
        for (j = 0; j < (L - 1) / 2; j++)
            if (lattice[2 * i + 1][j][0] == lattice[2 * i + 1][j][1])
            {
                if (conn[2 * i][2 * j] && conn[2 * i][2 * j + 1])
                {
                    label[lattice[2 * i + 1][j][0] - spins] =
                        find_root(label[lattice[2 * i][j][1] - spins], root);
                    if (label[lattice[2 * i + 1][j][0] - spins] !=
                        find_root(label[lattice[2 * i][j + 1][0] - spins],
                                    root))
                        root[root[label[lattice[2 * i][j + 1][0] - spins]]] =
                            label[lattice[2 * i + 1][j][0] - spins];
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
                        find_root(label[lattice[2 * i][j + 1][0] - spins],
                                    root);
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
                        find_root(label[lattice[2 * i][j + 1][0] - spins],
                                    root);
                else
                    label[lattice[2 * i + 1][j][1] - spins] = nl++;
            }
        if (conn[2 * i][L - 1]) // last column
            label[lattice[2 * i + 1][(L - 1) / 2][0] - spins] =
                find_root(label[lattice[2 * i][(L - 1) / 2][1] - spins],
                            root);
        else
            label[lattice[2 * i + 1][(L - 1) / 2][0] - spins] = nl++;
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
            size[label[k]] = 1;
        }
        else
        {
            label[k] = root[label[k]];
            size[label[k]]++;
        }
    delete[] root;
    
    // count H in each cluster
    for (k = 0; k < nl; k++)
        H[k] = 0;
    H[label[lattice[0][0][1] - spins]] += *lattice[0][0][1];
    for (j = 1; j < (L + 1) / 2; j++)
    {
        H[label[lattice[0][j][0] - spins]] += *lattice[0][j][0];
        H[label[lattice[0][j][1] - spins]] += *lattice[0][j][1];
    }
    for (j = 0; j < (L - 1) / 2; j++)
    {
        H[label[lattice[L - 2][j][0] - spins]] += *lattice[L - 2][j][0];
        H[label[lattice[L - 2][j][1] - spins]] += *lattice[L - 2][j][1];
    }
    H[label[lattice[L - 2][(L - 1) / 2][0] - spins]] +=
                                        *lattice[L - 2][(L - 1) / 2][0];
}

void flip_clusters(double p, long nl, long *label, long *H)
{
    long i;
    bool flip[L * L];
    double temp;
    for (i = 0; i < nl; i++)
    {
        if (H[i] > 0)
            temp = std::pow(p, H[i]);
        else
            temp = 1;
        if (unirnd(eng) < 0.5 * temp)
            flip[i] = 1;
        else
            flip[i] = 0;
    }
    for (i = 0; i < num; i++)
        if (flip[label[i]])
            spins[i] = -spins[i];
}

long long tonumber(long nl, long lab1, long lab2)
{
    using namespace std;
    return((long long)nl * min(lab1, lab2) + max(lab1, lab2));
}

void measure(double &E, double &mag)
{
    int i, j;
    long k;
    E = 0;
    E += (*lattice[0][0][1]) * (*lattice[1][0][0]) * coup[0][0];
    E += *lattice[0][0][1]; // H field
    for (j = 1; j <= (L - 1) / 2; j++)
    {
        E += (*lattice[1][j - 1][1]) * (*lattice[0][j][0])
			* coup[0][2 * j - 1];
        E += (*lattice[0][j][1]) * (*lattice[1][j][0]) * coup[0][2 * j];
        E += *lattice[0][j][0] + *lattice[0][j][1]; // H field
    }
    for (i = 1; i < (L - 1) / 2; i++)
    {
        E += (*lattice[2 * i][0][1]) * (*lattice[2 * i - 1][0][0])
			* coup[2 * i - 1][0];
        for (j = 1; j <= (L - 1) / 2; j++)
        {
            E += (*lattice[2 * i - 1][j - 1][1]) * (*lattice[2 * i][j][0])
				* coup[2 * i - 1][2 * j - 1];
            E += (*lattice[2 * i][j][1]) * (*lattice[2 * i - 1][j][0])
				* coup[2 * i - 1][2 * j];
        }
        E += (*lattice[2 * i][0][1]) * (*lattice[2 * i + 1][0][0])
			* coup[2 * i][0];
        for (j = 1; j <= (L - 1) / 2; j++)
        {
            E += (*lattice[2 * i + 1][j - 1][1]) * (*lattice[2 * i][j][0])
				* coup[2 * i][2 * j - 1];
            E += (*lattice[2 * i][j][1]) * (*lattice[2 * i + 1][j][0])
				* coup[2 * i][2 * j];
        }
    }
    // H field
    for (j = 0; j < (L - 1) / 2; j++)
        E += *lattice[L - 2][j][0] + *lattice[L - 2][j][1];
    E += *lattice[L - 2][(L - 1) / 2][0];
    E = -E / num;
    
    mag = 0;
    for (k = 0; k < num; k++)
        mag += spins[k];
    mag /= num;
}

void export_data(double *E, double *mag)
{
    using namespace std;
	ofstream file;
    file.open("E.bin", ios::out | ios::binary);
    file.write((char *)E, sizeof(double) * M);
	file.close();
    file.open("mag.bin", ios::out | ios::binary);
    file.write((char *)mag, sizeof(double) * M);
	file.close();
}

void data_analysis(double *E, double *mag)
{
    using namespace std;
    long i;
    double E_av = 0, E2_av = 0, mag_av = 0, mag2_av = 0, mag4_av = 0;
    for (i = M / 10; i < M; i++) // skip initial 1/10 steps to thermalize
    {
        E_av += E[i];
        E2_av += E[i] * E[i];
        mag_av += abs(mag[i]);
        double temp = mag[i] * mag[i];
        mag2_av += temp;
        mag4_av += temp * temp;
    }
    E_av /= (M - M / 10);
    E2_av /= (M - M / 10);
    mag_av /= (M - M / 10);
    mag2_av /= (M - M / 10);
    mag4_av /= (M - M / 10);
	ofstream file(str);
    file << "Number of spins is: " << num << endl;
    file << "At temperature T = " << T << endl;
    file << "Average energy E = " << E_av << endl;
    file << "Specific heat C = " << (E2_av - E_av * E_av) / T / T * num
        << endl;
    file << "Average magnetization M = " << mag_av << endl;
    file << "Binder cumulant U = " << 1 - mag4_av / 3 / mag2_av / mag2_av
        << endl;
	file.close();
}

void monte_carlo()
{
    long i;
    long nl; // number of labels
    long *label = new long[L * L], *size = new long[L * L],
        *H = new long[L * L]; // spins subjected to magnetic field in a
                              // cluster; a field of 2 counts for 2 spins
    double *E = new double[M], *mag = new double[M];
    lattice_structure();
    initialize_spins();
    for (i = 0; i < M; i++)
    {
        construct_graph(1 - std::exp(-2.0 / T));
        find_clusters(nl, label, size, H);
        measure(E[i], mag[i]);
        flip_clusters(std::exp(-2.0 / T), nl, label, H);
    }
    //export_data(E, mag);
    data_analysis(E, mag);
    delete[] label;
    delete[] size;
    delete[] H;
	delete[] E;
    delete[] mag;
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