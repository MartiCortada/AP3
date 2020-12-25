/* ------------------- METAHEURISTIC: SIMULATED ANNEALING ------------------- */
// AUTHORS: Jofre Poch Soler, Martí Cortada Garcia

/* ------------------------- ALL REQUIRED LIBRARIES ------------------------- */
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <utility>
#include <vector>

using namespace std;

/* ------------------ GLOBAL VARIABLES AND DATA STRUCTURES ------------------ */
typedef std::chrono::high_resolution_clock Time;
typedef std::chrono::duration<float> fsec;

int W, N; // W = roll width, N = total number of orders
auto start = Time::now(); // variable used for computing the time til we reach the solution

// matrix useful to determine whether a new piece can be placed in some positions or not
using UsedCanvas = vector<vector<bool>>;

// Data structure that stores information for the output
struct PiecePlace {
    pair<int, int> l; // Coordinates of the top-left point
    pair<int, int> r; // Coordinates of the bottom-right point
};

// Vector of positions (coordinates) of pieces placed in the roll, needed for the output
using Coordinates = vector<PiecePlace>;

// Data structure that represents a roll (piece of dimensions "p x q")
struct Roll {
    int p;
    int q;
};

// Data structure that stores information of the optimal result
struct OptimalResult {
    int L; // optimal length
    Coordinates coord; // optimal vector of coordinates
};

OptimalResult opt_res; // optimal result

vector<Roll> opt_s; // it will store our initial solution

double T = 100; // Temperature parameter (useful for simulated annealing)

/* ------------------------------- FUNCTIONS ------------------------------- */
/* Given a vector of coordinates coord and the length L of our output configuration,
it prints the solution in the way we are asked (whose correctnes will be checked
with an alternative program) and sends it in a file called output */
void write_output(string output, const Coordinates& coord, int L)
{
    // duration time manipulation
    auto end = Time::now(); // variable used for computing the time til we reach the solution
    fsec duration = end - start;

    // writing on the output file
    ofstream file;
    file.setf(ios::fixed);
    file.precision(1);
    file.open(output);
    file << duration.count() << endl;
    file << L << endl;
    for (int c = 0; c < int(coord.size()); ++c) {
        file << coord[c].l.second << " " << coord[c].l.first << "   ";
        file << coord[c].r.second << " " << coord[c].r.first << endl;
    }
    file << endl;
    file.close();
}

/* It defines the probability of accepting a worsening move using Boltzmann distribution.
Given the length of a configuration l_i, the optimal length l and a defined temperature T,
we calculate the probability as follows: exp(-(l_i - l)/(T)) */
double get_probability(int l_i, int l)
{
    if (l_i == l)
        ++l_i;
    // Let's clarify what the following equation defines: the farther l_i  is from l,
    // the lower the probability of accepting a move (given a fixed T).
    double probability = exp(-(l_i - l) / T);
    return probability;
}

/* It return the new updated temperature after each iteration using a geometric law
of parameter alpha (defined for us), computing the following: T_{k+1} = alpha * T_{k},
where T_{k} is a function of the current temperature.*/
void update_temperature()
{
    double alpha = 0.999;
    T *= alpha;
}

/* Given a solution (vector of pieces), returns another solution belonging to its
neighborhood. A vector of pieces s' is considered a neighbour of s if either swapping
two pieces of the vector or inverting the position of one, s' = s.
*/
vector<Roll> random_neighbour(const vector<Roll>& initial_solution)
{
    vector<Roll> neighbour_solution = initial_solution;

    int pos1 = rand() % (initial_solution.size());
    int pos2 = rand() % (initial_solution.size());

    // invert the coordinates of the roll situated in that position
    if (rand() % 2) {
        if (initial_solution[pos1].q <= W) {
            neighbour_solution[pos1].p = initial_solution[pos1].q;
            neighbour_solution[pos1].q = initial_solution[pos1].p;
        }
        if (initial_solution[pos2].q <= W) {
            neighbour_solution[pos2].p = initial_solution[pos2].q;
            neighbour_solution[pos2].q = initial_solution[pos2].p;
        }
    }
    // we swap the rolls situated in the positions of the random numbers
    else {
        neighbour_solution[pos1] = initial_solution[pos2];
        neighbour_solution[pos2] = initial_solution[pos1];
    }

    return neighbour_solution;
}


/* Given a boolean matrix B and 4 integers representing the dimensions
of a piece "p x q" and the positions of the matrix (i,j), returs whether it is
possible or not to place that piece in the matrix. */
bool legal(const UsedCanvas& B, int p, int q, int i, int j)
{
    if (j + p > W or B[i][j])
        return false; // if it is already filled or exceeds the matrix dimensions
    // if B[i][j] == false, which means that place is not filled, then iterate
    // for all position that would occupy the peace and check whether it's empty or not
    for (int r = i; r < i + q; ++r) {
        for (int s = j; s < j + p; ++s) {
            if (B[r][s])
                return false;
        }
    }
    return true; // We can place the piece!
}

/* Given a vector of rolls and a length, returns a solution in the output's format: a length L
and a vector of coordinates representing the top-left and the low-right cell occupied for every
piece. */
OptimalResult get_solution(const vector<Roll>& rolls, int max_length)
{
    UsedCanvas B(max_length, vector<bool>(W, false)); // worst case -> dimension "max_length x W"
    Coordinates coord(N); // will be filled by the coordinates of our output
    int l = 0; // length at the beginning
    bool placed; // will became true when a piece is placed
    // coordinates i and j to know where has been placed the last piece
    int coord_i = 0;
    int coord_j = 0;

    for (int idx = 0; idx < int(rolls.size()); ++idx) {
        placed = false;
        while (not placed) {
            // the roll doesn't have enought space on this row
            if (W - coord_j < rolls[idx].p) {
                ++coord_i;
                coord_j = 0;
            }
            // the current cell is already occupied
            else if (B[coord_i][coord_j]) {
                // if it's the last cell of the row, we jump into the next one
                if (coord_j == W - 1) {
                    ++coord_i;
                    coord_j = 0;
                } else {
                    ++coord_j;
                }
            } else {
                // check if all the cells that will occupy the roll are available
                if (legal(B, rolls[idx].p, rolls[idx].q, coord_i, coord_j)) {
                    for (int i = 0; i < rolls[idx].q; ++i) {
                        for (int j = 0; j < rolls[idx].p; ++j) {
                            B[coord_i + i][coord_j + j] = true;
                        }
                    }
                    l = max(l, coord_i + rolls[idx].q);
                    placed = true;
                    coord[idx].l = { coord_i, coord_j };
                    coord[idx].r = { coord_i + rolls[idx].q - 1, coord_j + rolls[idx].p - 1 };

                    // the coordinates are set to the top-right corner of the placed roll
                    coord_j = coord_j + rolls[idx].p - 1;
                } else {
                    ++coord_j;
                }
            }
        }
    }
    return { l, coord };
}

/* It performs the typical structure of a simulated annealing algorithm, updating
a given temperatura after each iteration untill termination conditions are met.
In other words, we try to approximatie the global optimum of a given function. */
void simulated_annealing(string output, int max_length)
{
    int k = 0;
    while (T > 0.0001) { // almost 0
        vector<Roll> s1 = random_neighbour(opt_s);
        OptimalResult S1 = get_solution(s1, max_length);

        if (S1.L < opt_res.L) {
            opt_res = S1;
            write_output(output, opt_res.coord, opt_res.L);
        } else {
            double prob = get_probability(S1.L, opt_res.L);

            if (((double)rand() / (RAND_MAX)) <= prob) {
                opt_res = S1;
                write_output(output, opt_res.coord, opt_res.L);
            }
        }

        update_temperature();
        ++k;
    }
}

/* It generates an initial solution randomly as follows: first select randomly
which piece will be placed first and then (again randomly) if we store it with the
dimensions p x q or q x p, and so on and so forth. */
vector<Roll> generate_initial_solution(vector<Roll>& rolls)
{
    vector<Roll> initial_solution(int(rolls.size()));

    // select randomly which piece we choose at each step
    int idx;
    int idx1 = 0;
    for (int i = int(rolls.size()) - 1; i > 0; --i) {
        idx = (rand() % i) + 1; // random value in [0, ..., i]
        initial_solution[idx1] = rolls[idx];
        rolls[idx] = rolls[int(rolls.size()) - 1];
        rolls.pop_back();
        ++idx1;
    }
    initial_solution[idx1] = rolls[0];

    // select randomly, for each position, if we place the piece as p x q or q x p
    int p, q;
    for (int i = 0; i < int(initial_solution.size()); ++i) {
        // 0 o 1 with prob 1/2: if 0 we get p x q, if 1 we get q x p
        if (rand() % 2) { // we change p x q -> q x p
            if (initial_solution[i].q <= W) {
                q = initial_solution[i].q;
                p = initial_solution[i].p;
                initial_solution[i].p = q;
                initial_solution[i].q = p;
            }
        }
    }

    return initial_solution;
}

/* Main loop:
    it reads from a file all the data we need and store this data in a vector.
    Then, we send this data on a simulated annealing (metaheuristic) algorithm
    in order to find a solution (not necessarily the optimal)!
*/
int main(int argc, char* argv[])
{
    // read input from a file
    ifstream f(argv[1]);
    f >> W >> N;

    vector<Roll> rolls; // data structure that stores our initial configuration
    int N_copy = N; // careful! We do not want to change N since it is global variable!
    int n, p, q; // n = number of orders of a specific type of dimensions "p x q"
    int max_length; // the L of the worst possible combination, used to define the lenght of
        // the boolean matrix (canvas) and the first optimal lenght

    max_length = 0;
    while (N_copy > 0) {
        f >> n >> p >> q;
        for (int i = 0; i < n; ++i)
            rolls.push_back({ p, q });
        max_length += n * max(p, q);
        N_copy -= n;
    }
    f.close();
    ++max_length;
    opt_res.L = max_length; // max_length at the beginning (worst case)

    srand(time(NULL)); // each time the random selection will be different!
    vector<Roll> initial_solution = generate_initial_solution(rolls);

    opt_res = get_solution(initial_solution, max_length); // optimal result
    opt_s = initial_solution; // initial solution

    simulated_annealing(argv[2], max_length); // simulated annealing algorithm
}
