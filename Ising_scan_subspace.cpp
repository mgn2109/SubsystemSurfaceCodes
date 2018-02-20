#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <string>
#include <random>
#include <iomanip>



#define L 11 //L is odd, consider L*L data qubits.
double q = 1 ; //probability of minimal X-type plaquette.
double nlogp = 3; //p parameter given as -ln(p).
double T; //temperature
long long M = 1000; //number of steps
int seed = 1; //seed number
std::string str; //output name


int lattice_template[L - 1][L - 1]; //Template that gives connectivity and spin structure of lattice.
int *gauge_generators[L - 1][L]; //Indexes the vertical gauge generators as pointers to spins.
int mySpins[L * L]; //(Over-)indexes the spins.
int myCoup[L][L]; //Determines the interactions, ferro- or anti-, one per qubit.
long counter;  //number of mySpins
int myLocation[L * L][2][2]; //holds locations of spins in terms of their "upper left corner qubit"
                             //and "upper right corner qubit."


std::mt19937_64 eng;
std::uniform_real_distribution<double> unirnd(0.0, 1.0);





void generate_lattice_template()  //Template that describes spin structure and connectivity.
{
	int i, j;

	for (i = 0; i < L - 1; i++)
	{
		for (j = 0; j < L - 1; j++)
		{
			if (unirnd(eng) < q)
				lattice_template[i][j] = 1;
			else
				lattice_template[i][j] = 0;
		}
	}



	/*// Test.

	std::cout << std::endl << "generate_lattice_template" << std::endl << std::endl;

	for (int i = 0; i < L - 1; i++)
	{
		for (int j = 0; j < L - 1; j++)
		{
			std::cout << lattice_template[i][j] << " ";
		}
		std::cout << std::endl;
	}

	std::cout << std::endl << std::endl;

	// End Test.*/

}





void assign_generators_to_spins() //Points the gauge generators to their corresponding spins.
{                                 //Also saves the locations of the spins in the lattice.
	counter = 0;

	int i, j;

	for (i = 0; i < L - 1; i++)
	{

		myLocation[counter][0][0] = i; //The first spin in any row always begins
		myLocation[counter][0][1] = 0; //with the first qubit in that row

		for (j = 0; j < L; j++)
		{
			gauge_generators[i][j] = &mySpins[counter];

			if (j < L - 1)
			{
				if (lattice_template[i][j] == 0)
				{
					myLocation[counter][1][0] = i; //Saves the end location of this spin.
					myLocation[counter][1][1] = j;

					counter++;

					myLocation[counter][0][0] = i;  //Saves the beginning location of this spin.
					myLocation[counter][0][1] = j + 1;

				}
			}
			else
			{
				myLocation[counter][1][0] = i; //The last spin in any row will end at the last qubit.
				myLocation[counter][1][1] = j;

				counter++;
			}
		}
	}



	/*//Test.

	std::cout << std::endl << "assign_generators_to_spins" << std::endl << std::endl;

	for (int i = 0; i < counter; i++)
	{
		std::cout << "(" << myLocation[i][0][0] << "," << myLocation[i][0][1] 
		<< ") ---> (" << myLocation[i][1][0] << "," << myLocation[i][1][1] << ")" << std::endl;
		std::cout << &mySpins[i] << std::endl << std::endl;
	}

	std::cout << "-------------------------------------" << std::endl << std::endl;

	for (int i = 0; i < L-1; i++)
	{
		for (int j = 0; j < L; j++)
		{
			std::cout << "(" << i << "," << j << "): " << gauge_generators[i][j] << std::endl;
		}
	}

	std::cout << std::endl;

	//End Test.*/
}





void initialize_couplings() //Initialize frustrated couplings: one per qubit.
{
	double p = std::exp(-nlogp);
	int i,j;

	for (i = 0; i < L; i++)
	{
		for (j = 0; j < L; j++)
		{

			if (unirnd(eng) < p)
				myCoup[i][j] = -1;  //1 for ferromagnetic, -1 for anti-ferro.  So p=0 gives full ferro.
			else
				myCoup[i][j] = 1;
		}
	}



	/*//Test.

	std::cout << "initialize_couplings" << std::endl;

	for (int i = 0; i < L; i++)
	{
		
		std::cout << std::endl;
		
		for (int j = 0; j < L; j++)
		{
			std::cout << myCoup[i][j] << " "; 
		}
	}

	std::cout << std::endl;

	//End Test.*/
}





void initialize_mySpins() //Initializes mySpins.
{
	for (int i = 0; i < counter; i++)
		mySpins[i] = 1;
}





void myMeasure(double &E, double &mag)
{
	int i, j;
	long k;
	E = 0;

	for(i = 1; i < L - 1; i++)  //These are bulk contributions
	{
		for(j = 0; j < L; j++)  
		{
			E += (*gauge_generators[i - 1][j]) * (*gauge_generators[i][j]) * myCoup[i][j];
			//Configuration energy in the bulk.  Fix the (i,j) bulk qubit.  There's a coupling
			//associated by myCoup[i][j].  Multiply this by the connected spins, whose values 
			//are given by pointers from their constituent gauge-generator pointers.
		}
	}

	for(j = 0; j < L; j++)
	{
		E += (*gauge_generators[0][j]) * myCoup[0][j] //Top boundary contribution.
		   + (*gauge_generators[L - 2][j]) * myCoup[L-1][j]; //Bottom boundary contribution.
	}

	E = -E; //Correct energy sign.

	mag = 0;  //Add up magnetization.  We have now measured the magnetization and 
			  //energy of the configuration.

	for (k = 0; k < counter; k++)
		mag += mySpins[k];



	/*//Test.

	std::cout << "myMeasure" << std::endl;
	std::cout << "mag is " << mag << " and E is " << E << std::endl << std::endl;

	*///End Test.
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
	file << "q = " << q << ", p = " << setprecision(16) << exp(-nlogp)
        << ", L = " << L << endl;
	file << "Number of spins is: " << counter << endl;
	file << "Current number of steps = " << step << endl;
	file << "At temperature T = " << T << endl;
	file << "Average energy E = " << E_sum << endl;
	file << "Specific heat C = " << (E2_sum - E_sum * E_sum) / T / T * counter
		<< endl;
	file << "Average magnetization M = " << mag_sum << endl;
	file << "Binder cumulant U = " << 1 - mag4_sum / 3 / mag2_sum / mag2_sum
		<< endl << endl;
	file.close();
}





void myDiff(long id, double &dE, double &dmag)
{
	
	dmag = -2 * mySpins[id];  //Change in magnetization.

	dE = 0;  //We compute the local energy.

	int j;
	int startQubit[2];
	int endQubit[2];

	startQubit[0] = myLocation[id][0][0];
	startQubit[1] = myLocation[id][0][1];

	endQubit[0] = myLocation[id][1][0];
	endQubit[1] = myLocation[id][1][1];

	if (startQubit[0] == 0) //Upper boundary case.
	{
		for (j = startQubit[1]; j <= endQubit[1]; j++)
		{
			dE += mySpins[id] * myCoup[startQubit[0]][j]; //Adds up the upper-boundary interactions

			dE += mySpins[id] * myCoup[startQubit[0] + 1][j] * (*gauge_generators[startQubit[0] + 1][1]);
			//Adds up the bottom interactions as the spin (mySpins[id]) times
			//the coupling interactions (myCoup[1 qubit down][j]) times
			//the adjacent spin, whose value is pointed to by the gauge_generator on column
			//j and row startQubit + 1.
		}


	}
	else if (startQubit[0] == L-2) //Bottom boundary case.
	{
		for (j = startQubit[1]; j <= endQubit[1]; j++)
		{
			dE += mySpins[id] * myCoup[startQubit[0]][j] * (*gauge_generators[startQubit[0] - 1][j]);
			// Adds up the top interactions

			dE += mySpins[id] * myCoup[startQubit[0] + 1][j];
			// Adds up the bottom-boundary interactions
		}
	}

	else //Bulk case
	{
		for(j = startQubit[1]; j <= endQubit[1]; j++)
		{
			dE += mySpins[id] * myCoup[startQubit[0]][j] * (*gauge_generators[startQubit[0] - 1][j]);
			//Adds up the top interactions

			dE += mySpins[id] * myCoup[startQubit[0] + 1][j] * (*gauge_generators[startQubit[0] + 1][j]);
			//Adds up the bottom interactions
		}
	}

	dE *= 2; //difference is double the local energy.




	/*//Test

	std::cout << std::endl << "myDiff" << std::endl << std::endl;
	std::cout << "id is: " << id << " and dE is: " << dE;
	std::cout << std::endl << std::endl;

	*///End Test.
}





void my_monte_carlo()
{
	long id;
    long long k = 1;
	double dE, dmag, E, mag, E_sum = 0, E2_sum = 0,
		mag_sum = 0, mag2_sum = 0, mag4_sum = 0;

	generate_lattice_template();     //Create lattice template.
	assign_generators_to_spins();  //According to lattice template, point generators to spins

	initialize_mySpins(); //Initialize all mySpins to up.
	initialize_couplings(); //Initialize the interactions.
	

	//We now have the spins and couplings initialized, where the interaction is determined
	//for each spin by the vertically adjacent gauge generators, which point to their
	//associated spins.

	myMeasure(E, mag);  //Update E and mag according to lattice.

	std::uniform_int_distribution<long> randint(0, counter - 1);  //Creates an RNG over the spins.

	for (long long i = 1; i < M; i++)
	{

		id = randint(eng); //id indexes the candidate spin to be flipped.

		myDiff(id, dE, dmag); //this should compute the change in energy dE and change in
							  //magnetization

		if (unirnd(eng) < std::exp(-dE / T))
		{
			mySpins[id] = -mySpins[id];
			E += dE;
			mag += dmag;
		}
		
		if (i >= M / 10)
		{
			double temp;
			temp = E / counter;
			E_sum += temp;
			E2_sum += temp * temp;
			temp = std::abs(mag) / counter;
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
	q = atof(argv[1]);     //Pass in q first.
	nlogp = atof(argv[2]);     //Pass in -ln(p)=ln(1/p) second.
	if (nlogp <= log(2.0))
	{
		cout << "Error! -log p must be greater than log 2." << endl;
		return(-1);
	}
	T = 2.0 / (log(1 - exp(-nlogp)) + nlogp);     //Fixes T according to the Nishimori line. 
	if (argc > 3)
		M = atoll(argv[3]);     //Fixes the number of time steps.
	if (argc > 4)
		seed = atoi(argv[4]);     //Fixes the initial random seed.
	str = "L_" + to_string(L) + "___q_" + argv[1] + "___nlogp_" + argv[2]     //String (for naming) as system size, q, ln(1/p), number of time steps, and seed.
		+ "___M_" + to_string(M) + "___seed_" + to_string(seed) + "_.txt";
	eng.seed(seed);     //Passes seed to RNG.
	my_monte_carlo();     //Runs Monte Carlo experiment.
	return(0);
}