/* ----------------------------- METAHEURISTIC ----------------------------- */
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

OptimalResult opt_res; // optimal result (global variabla) -> it will change while iterating

vector<Roll> opt_s;

double T = 30;

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
    if (l_i == l) ++l_i;
    // Let's clarify what the following equation do: the farther l_i  is from l,
    // the lower the probability of accepting a move (given a fixed T).
    double probability = exp(-(l_i - l) / (T));
    return probability;
}

/* It return the new updated temperature after each iteration using a geometric law
of parameter alpha defined for us, computing the following: T_{k+1} = alpha * T_{k},
where T_{k} is a function of the current temperature with the iteration counter k.*/
void update_temperature()
{   
    double alpha = 0.99;
    T *= alpha;
}

/* Given a solution (vector of pieces), returns another solution belonging to its
neighborhood. A vector of pieces s' is considered a neighbour of s if either swapping
two pieces of the vector or inverting the position of one, s' = s.
*/
vector<Roll> random_neighbour(const vector<Roll>& initial_solution)
{
    vector<Roll> neighbour_solution = initial_solution;

    // invert the coordinates of the roll situated in that position
    if (rand() % 2) {
        int pos = rand() % (initial_solution.size());
        neighbour_solution[pos].p = initial_solution[pos].q;
        neighbour_solution[pos].q = initial_solution[pos].p;
    }
    // we swap the rolls situated in the positions of the random numbers
    else {
        // two random numbers between 0 and the number of pieces
        int pos1 = rand() % (initial_solution.size());
        int pos2 = rand() % (initial_solution.size());
        neighbour_solution[pos1] = initial_solution[pos2];
        neighbour_solution[pos2] = initial_solution[pos1];
    }

    return neighbour_solution;
}

void print_solution(const vector<Roll>& s)
{
    for (auto roll : s)
        cout << roll.p << " " << roll.q << endl;
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

    for (int idx = 0; idx < rolls.size(); ++idx) {
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

void simulated_annealing(string output, int max_length)
{
    int k = 0;
    while (T > 0.01) { // almost 0
        vector<Roll> s1 = random_neighbour(opt_s);
        OptimalResult S1 = get_solution(s1, max_length);
         cout << "Kklk" << endl;

        if (S1.L < opt_res.L) {
            opt_res = S1;
            write_output(output, opt_res.coord, opt_res.L);
        }
        else {
            double prob = get_probability(S1.L, opt_res.L);
            cout << "prob: " << prob << endl;

            if (((double) rand() / (RAND_MAX)) <= prob) {
                opt_res = S1;
                write_output(output, opt_res.coord, opt_res.L);
            }
        }

        update_temperature();
        cout << "T: " << T << endl;
        cout << "k: "<< k << endl;
        ++k;
    }
}

int main(int argc, char* argv[])
{
    // read input from a file
    ifstream f(argv[1]);
    f >> W >> N;

    vector<Roll> initial_solution; // data structure that stores our initial configuration
    int N_copy = N; // careful! We do not want to change N since it is global variable!
    int n, p, q; // n = number of orders of a specific type of dimensions "p x q"
    int max_length; // the L of the worst possible combination, used to define the lenght of
        // the boolean matrix (canvas) and the first optimal lenght

    max_length = 0;
    while (N_copy > 0) {
        f >> n >> p >> q;
        for (int i = 0; i < n; ++i)
            initial_solution.push_back({ p, q });
        max_length += n * max(p, q);
        N_copy -= n;
    }
    f.close();
    ++max_length;
    opt_res.L = max_length; // max_length at the beginning (worst case)

    opt_res = get_solution(initial_solution, max_length);
    opt_s = initial_solution;

    srand(time(NULL));

    simulated_annealing(argv[2], max_length);
}
