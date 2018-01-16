#include <fstream>
#include <cstdlib>
#include <cmath>
#include <random>

#define L 11 // L is odd, consider L*L data qubits
#deffine q 0.5
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
        for (j = 0; j < (L - 1) / 2; j++)
        {
            if (unirnd(eng) < q)
            {
                lattice[2 * i][j + 1][0] = &spins[num];
                lattice[2 * i][j + 1][1] = &spins[num];
                num++;
            }
            else
            {
                lattice[2 * i][j + 1][0] = &spins[num];
                num++;
                lattice[2 * i][j + 1][1] = &spins[num];
                num++;
            }
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
        root[x] = find_root(x, root);
        return(root[x]);
    }
}

void find_clusters(bool conn[L - 2][L], long &nl, long label[L * L])
{
    int i, j;
    nl = 0; // current label;
    long root[L * L];
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
    long k;
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
        }
        else
            label[k] = root[label[k]];
}

void flip_cluster()
{
    
}

void monte_carlo(double E[M], double mag[M])
{
    long i;
    long nl; // number of labels
    bool conn[L - 2][L];
    long label[L * L];
    double E[M], E2[M], mag[M], mag2[M], mag4[M];
    lattice_structure();
    initialize_spins();
    for (i = 0; i < M; i++)
    {
        construct_graph(1 - std::exp(-2 / T), conn);
        find_cluster(conn, nl, label);
    }
}

void data_analysis(double E[M], double mag[M])
{
    
}

int main()
{
    monte_carlo(E, mag);
    
}