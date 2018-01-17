#include <fstream>
#include <cstdlib>
#include <cmath>
#include <random>

#define L 11 // L is odd, consider L*L data qubits
#define q 0.5
#define M 1000 // number of steps
#define T 1 // temperature

int *lattice[L - 1][(L + 1) / 2][2];
int spins[L * L]; // spin number is less than L^2
int num; // number of spins

std::default_random_engine eng(0);
std::uniform_real_distribution<double> unirnd(0.0, 1.0);

void lattice_structure()
{
    int i, j;
    num = 0;
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
}

void initialize_spins()
{
    for (int i = 0; i < num; i++)
        spins[i] = 1;
}

void construct_graph(double p, bool conn[L - 2][L])
{
    int i, j;
    // if two adjacent spins are parallel, connect them with probability p
    // otherwise they are disconnected
    if (*lattice[0][0][1] == *lattice[1][0][0] && unirnd(eng) < p)
        conn[0][0] = 1;
    else
        conn[0][0] = 0;
    for (j = 1; j <= (L - 1) / 2; j++)
    {
        if (*lattice[1][j - 1][1] == *lattice[0][j][0] && unirnd(eng) < p)
            conn[0][2 * j - 1] = 1;
        else
            conn[0][2 * j - 1] = 0;
        if (*lattice[0][j][1] == *lattice[1][j][0] && unirnd(eng) < p)
            conn[0][2 * j] = 1;
        else
            conn[0][2 * j] = 0;
    }
    for (i = 1; i < (L - 1) / 2; i++)
    {
        if (*lattice[2 * i][0][1] == *lattice[2 * i - 1][0][0]
            && unirnd(eng) < p)
            conn[2 * i - 1][0] = 1;
        else
            conn[2 * i - 1][0] = 0;
        for (j = 1; j <= (L - 1) / 2; j++)
        {
            if (*lattice[2 * i - 1][j - 1][1] == *lattice[2 * i][j][0]
                && unirnd(eng) < p)
                conn[2 * i - 1][2 * j - 1] = 1;
            else
                conn[2 * i - 1][2 * j - 1] = 0;
            if (*lattice[2 * i][j][1] == *lattice[2 * i - 1][j][0]
                && unirnd(eng) < p)
                conn[2 * i - 1][2 * j] = 1;
            else
                conn[2 * i - 1][2 * j] = 0;
        }
        if (*lattice[2 * i][0][1] == *lattice[2 * i + 1][0][0]
                && unirnd(eng) < p)
            conn[2 * i][0] = 1;
        else
            conn[2 * i][0] = 0;
        for (j = 1; j <= (L - 1) / 2; j++)
        {
            if (*lattice[2 * i + 1][j - 1][1] == *lattice[2 * i][j][0]
                && unirnd(eng) < p)
                conn[2 * i][2 * j - 1] = 1;
            else
                conn[2 * i][2 * j - 1] = 0;
            if (*lattice[2 * i][j][1] == *lattice[2 * i + 1][j][0]
                && unirnd(eng) < p)
                conn[2 * i][2 * j] = 1;
            else
                conn[2 * i][2 * j] = 0;
        }
    }
}

int find_root(long x, long root[L * L])
{
    if (root[x] == x)
        return(x);
    else
    {
        root[x] = find_root(root[x], root);
        return(root[x]);
    }
}

void find_clusters(bool conn[L - 2][L], long &nl, long label[L * L],
                    long size[L * L])
{
    int i, j;
	long k;
    nl = 0; // current label;
    long root[L * L];
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
}

void flip_clusters(long nl, long label[L * L])
{
    long i;
    bool flip[L * L];
    for (i = 0; i < nl; i++)
        if (unirnd(eng) < 0.5)
            flip[i] = 0;
        else
            flip[i] = 1;
    for (i = 0; i < num; i++)
        if (flip[label[i]])
            spins[i] = -spins[i];
}

void measure(long nl, long label[L * L], long size[L * L],
        double &mag, double &mag2, double &mag4)
{
    long i, temp;
    mag = 0;
    for (i = 0; i < num; i++)
        mag += spins[i];
    mag /= num;
    mag2 = 0;
    mag4 = 0;
    for (i = 0; i < nl; i++)
    {
        temp = size[i] * size[i];
        mag2 += temp;
        mag4 += temp * temp;
    }
    mag4 = 3 * mag2 * mag2 - 2 * mag4;
    temp = num * num;
    mag2 /= temp;
    mag4 /= temp * temp;
}

void export_data(double mag[M], double mag2[M], double mag4[M])
{
    using namespace std;
	ofstream file;
    file.open("mag.bin", ios::out | ios::binary);
    file.write((char *)mag, sizeof(double) * M);
	file.close();
    file.open("mag2.bin", ios::out | ios::binary);
    file.write((char *)mag2, sizeof(double) * M);
	file.close();
    file.open("mag4.bin", ios::out | ios::binary);
    file.write((char *)mag4, sizeof(double) * M);
	file.close();
}

void monte_carlo()
{
    long i;
    long nl; // number of labels
    bool conn[L - 2][L];
    long label[L * L], size[L * L];
    double mag[M], mag2[M], mag4[M];
    lattice_structure();
    initialize_spins();
    for (i = 0; i < M; i++)
    {
        construct_graph(1 - std::exp(-2 / T), conn);
        find_clusters(conn, nl, label, size);
        measure(nl, label, size, mag[i], mag2[i], mag4[i]);
        flip_clusters(nl, label);
    }
    export_data(mag, mag2, mag4);
}

int main()
{
    monte_carlo();
    return(0);
}