#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <random>

#define L 11 // L is odd, consider L*L data qubits
#define q 0.5
#define M 100000 // number of steps
#define T 1 // temperature

int *lattice[L - 1][(L + 1) / 2][2]; //lattice of qubits
int spins[L * L]; // spins, spin number is less than L^2
bool conn[L - 2][L]; // cluster connectivity
int num; // number of spins

std::default_random_engine eng(0);
std::uniform_real_distribution<double> unirnd(0.0, 1.0);


//picks out which qubits are the same, and which are different, 
//according to the lattice connectivity
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
			if (unirnd(eng) < q) //with probability q, spins are same
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


//initializes all spins in the up state
void initialize_spins()
{
    for (int i = 0; i < num; i++)
        spins[i] = 1;
}

/*
void construct_graph(double p)
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
*/


/*
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
*/

/*
void find_clusters(long &nl, long *label, long *size)
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
}
*/


/*
void flip_clusters(long nl, long *label)
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
*/

//single spin-flip function
void flip_spin()
{
	int i = (int) unirnd(eng) * (L-1);
	int j = (int) unirnd(eng) * (L+1)/2;
	int k = (int) unirnd(eng) * 2;

	int delta_E;

	int candidate_spin = *lattice[i][j][k];

	if(i != 0 && i != (L-2) && j != 0 && j != ((L+1)/2 - 1)) //Checks if chosen spin lies in bulk.
	{	
		if (i%2 == 0)
		{
			if(lattice[i][j][k] == lattice[i][j][!k]) //Checks if interaction is 4-local.
			{
				delta_E = (-1) * ((-1) * candidate_spin - candidate_spin) * 
				(*lattice[i-1][j][1] + *lattice[i-1][j+1][0] + 
				 *lattice[i+1][j][1] + *lattice[i+1][j+1][0]);
		
			else
			{
				delta_E = (-1) * ((-1) * candidate_spin - candidate_spin) * 
				(*lattice[i-1][j + k][!k] +
				*lattice[i+1][j + k][!k]);
			}
		}
		else
		{
			if(lattice[i][j][k] == lattice[i][j][!k]) //Checks if interaction is 4-local.
			{
				delta_E = (-1) * ((-1) * candidate_spin - candidate_spin) * 
				(*lattice[i-1][j-1][1] + *lattice[i-1][j][0] + 
				 *lattice[i+1][j-1][1] + *lattice[i+1][j][0]);
		
			else
			{
				delta_E = (-1) * ((-1) * candidate_spin - candidate_spin) * 
				(*lattice[i-1][j - (!k)][!k] +
				*lattice[i+1][j - (!k)][!k]);
			}
		}

	}
		
	
	else if(i == 0 && j == 0) //top left corner case
	{
		if(lattice[i][j][k] == lattice[i][j][!k]) //Checks if interaction is 4-local.
		{
			delta_E = (-1) * ((-1) * candidate_spin - candidate_spin) * 
			(2 + 
			 *lattice[i+1][j][1] + *lattice[i+1][j+1][0]);
		}
		else
		{
			delta_E = (-1) * ((-1) * candidate_spin - candidate_spin) * 
			(1 +
			*lattice[i+1][j+k][!k]);
		}
	}

	else if(i == 0 && j == ((L+1)/2)-1) //top right corner case
	{
		delta_E = (-1) * ((-1) * candidate_spin - candidate_spin) * 
			(1 +
			*lattice[i+1][j][1]);
	}

	else if(i == L-2 && j == 0) //bottom left corner case
	{
		delta_E = (-1) * ((-1) * candidate_spin - candidate_spin) * 
			(1 +
			*lattice[i-1][j][0]);
	}

	else if(i == L-2 && j == ((L+1)/2) - 1) //bottom right corner case
	{
		if(lattice[i][j][k] == lattice[i][j][!k]) //Checks if interaction is 4-local.
		{
			delta_E = (-1) * ((-1) * candidate_spin - candidate_spin) * 
			(2 + 
			 *lattice[i-1][j-1][1] + *lattice[i-1][j][0]);
		}
		else
		{
			delta_E = (-1) * ((-1) * candidate_spin - candidate_spin) * 
			(1 +
			*lattice[i-1][j][0]);
		}
	}

	else if(i == 0 && j != 0 && j != (L+1)/2 - 1) //top border case
	{
		if(lattice[i][j][k] == lattice[i][j][!k]) //Checks if interaction is 4-local.
		{
			delta_E = (-1) * ((-1) * candidate_spin - candidate_spin) * 
			(2 + 
			 *lattice[i+1][j][1] + *lattice[i+1][j+1][0]);

		else
		{
			delta_E = (-1) * ((-1) * candidate_spin - candidate_spin) * 
			(1 +
			*lattice[i+1][j + k][!k]);
		}
	}

	else if(i == L-2 && j != 0 && j != (L+1)/2 - 1) //bottom border case
	{
		if(lattice[i][j][k] == lattice[i][j][!k]) //Checks if interaction is 4-local.
		{
			delta_E = (-1) * ((-1) * candidate_spin - candidate_spin) * 
			(*lattice[i-1][j][1] + *lattice[i-1][j+1][0] + 
			 2);
		
		else
		{
			delta_E = (-1) * ((-1) * candidate_spin - candidate_spin) * 
			(*lattice[i-1][j + k][!k] +
			1);
		}
	}

	else if(i != 0 && i != L-2 && j == 0) //left border case
	{
		delta_E = (-1) * ((-1) * candidate_spin - candidate_spin) * 
			(*lattice[i-1][0][0] +
			*lattice[i+1][0][0]);
	}

	else if(i != 0 && i != L-2 && j == ((L+1)/2) - 1) //right border case
	{
		delta_E = (-1) * ((-1) * candidate_spin - candidate_spin) * 
			(*lattice[i-1][j][1] +
			*lattice[i+1][j][1]);
	}

	
	if(unirnd(eng) < std::exp(-delta_E / T)) //flip spin according to local Boltzmann
	{
		*lattice[i][j][k] = (-1) * (*lattice[i][j][k]);
	}


}


//data function
long long tonumber(long nl, long lab1, long lab2)
{
    using namespace std;
    return((long long)nl * min(lab1, lab2) + max(lab1, lab2));
}


//
void measure(long nl, long *label, double &E, double &E2)
{
    int i, j;
    long k, nb = 0;
    long long *bond = new long long[(L - 2) * L];
    E = 0;
    if (label[lattice[0][0][1] - spins] == label[lattice[1][0][0] - spins])
        E++;
    else
        bond[nb++] = tonumber(nl, label[lattice[0][0][1] - spins],
                                label[lattice[1][0][0] - spins]);
    for (j = 1; j <= (L - 1) / 2; j++)
    {
        if (label[lattice[1][j - 1][1] - spins] ==
                label[lattice[0][j][0] - spins])
            E++;
        else
            bond[nb++] = tonumber(nl, label[lattice[1][j - 1][1] - spins],
                                    label[lattice[0][j][0] - spins]);
        if (label[lattice[0][j][1] - spins] ==
                label[lattice[1][j][0] - spins])
            E++;
        else
            bond[nb++] = tonumber(nl, label[lattice[0][j][1] - spins],
                                    label[lattice[1][j][0] - spins]);
    }
    for (i = 1; i < (L - 1) / 2; i++)
    {
        if (label[lattice[2 * i][0][1] - spins] ==
                label[lattice[2 * i - 1][0][0] - spins])
            E++;
        else
            bond[nb++] = tonumber(nl, label[lattice[2 * i][0][1] - spins],
                                    label[lattice[2 * i - 1][0][0] - spins]);
        for (j = 1; j <= (L - 1) / 2; j++)
        {
            if (label[lattice[2 * i - 1][j - 1][1] - spins] ==
                    label[lattice[2 * i][j][0] - spins])
                E++;
            else
                bond[nb++] = tonumber(nl,
                        label[lattice[2 * i - 1][j - 1][1] - spins],
                        label[lattice[2 * i][j][0] - spins]);
            if (label[lattice[2 * i][j][1] - spins] ==
                    label[lattice[2 * i - 1][j][0] - spins])
                E++;
            else
                bond[nb++] = tonumber(nl,
                        label[lattice[2 * i][j][1] - spins],
                        label[lattice[2 * i - 1][j][0] - spins]);
        }
        if (label[lattice[2 * i][0][1] - spins] ==
                label[lattice[2 * i + 1][0][0] - spins])
            E++;
        else
            bond[nb++] = tonumber(nl, label[lattice[2 * i][0][1] - spins],
                                    label[lattice[2 * i + 1][0][0] - spins]);
        for (j = 1; j <= (L - 1) / 2; j++)
        {
            if (label[lattice[2 * i + 1][j - 1][1] - spins] ==
                    label[lattice[2 * i][j][0] - spins])
                E++;
            else
                bond[nb++] = tonumber(nl,
                        label[lattice[2 * i + 1][j - 1][1] - spins],
                        label[lattice[2 * i][j][0] - spins]);
            if (label[lattice[2 * i][j][1] - spins] ==
                    label[lattice[2 * i + 1][j][0] - spins])
                E++;
            else
                bond[nb++] = tonumber(nl,
                        label[lattice[2 * i][j][1] - spins],
                        label[lattice[2 * i + 1][j][0] - spins]);
        }
    }
    E = -E / num;
    std::sort(bond, bond + nb);
    long long last = -1;
    long count = 0;
    E2 = 0;
    for (k = 0; k < nb; k++)
        if (bond[k] == last)
            count++;
        else
        {
            E2 += count * count;
            last = bond[k];
            count = 1;
        }
    E2 += count * count;
    E2 = E2 / num / num  + E * E;
    delete[] bond;
}


//writes out data
void export_data(double *E, double *E2)
{
    using namespace std;
	ofstream file;
    file.open("E.bin", ios::out | ios::binary);
    file.write((char *)E, sizeof(double) * M);
	file.close();
    file.open("E2.bin", ios::out | ios::binary);
    file.write((char *)E2, sizeof(double) * M);
	file.close();
}


//chews data
void data_analysis(double *E, double *E2)
{
    using namespace std;
    long i;
    double E_av = 0, E2_av = 0;
    for (i = M / 10; i < M; i++) // skip initial 1/10 steps to thermalize
    {
        E_av += E[i];
        E2_av += E2[i];
    }
    E_av /= (M - M / 10);
    E2_av /= (M - M / 10);
    cout << "Number of spins is: " << num << endl;
    cout << "At temperature T = " << T << endl;
    cout << "Average energy E = " << E_av << endl;
    cout << "Specific heat C = " << (E2_av - E_av * E_av) / T / T * num
        << endl;
}



void monte_carlo()
{
    long i;
    long nl; // number of labels
    long *label = new long[L * L], *size = new long[L * L];
    double *E = new double[M], *E2 = new double[M];
    lattice_structure();
    initialize_spins();
    for (i = 0; i < M; i++)
    {
        construct_graph(1 - std::exp(-2.0 / T));
        find_clusters(nl, label, size);
        measure(nl, label, E[i], E2[i]);
        flip_clusters(nl, label);
    }
    //export_data(E, E2);
    data_analysis(E, E2);
    delete[] label;
    delete[] size;
	delete[] E;
	delete[] E2;
}

int main()
{
    monte_carlo();
    return(0);
}