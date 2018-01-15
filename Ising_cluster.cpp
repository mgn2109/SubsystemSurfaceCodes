#include <fstream>
#include <cstdlib>
#include <random>

#define L 11 // L is odd, consider L*L data qubits
#deffine q 0.5
#define M 1000 // number of steps
#define T 1 // temperature

int *lattice[L - 1][(L + 1) / 2][2];
int spins[L * L]; // spin number is less than L^2

std::default_random_engine eng(0);
std::uniform_real_distribution<double> unirnd(0.0, 1.0);

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

int root_label(long x, long *label)
{
    
}

void find_clusters(bool conn[L - 2][L])
{
    
}

void flip_cluster()
{
    
}

void monte_carlo(double E[M], double mag[M])
{
    
}

void data_analysis(double E[M], double mag[M])
{
    
}

int main()
{
    double E[M], mag[M];
    monte_carlo(E, mag);
    
}