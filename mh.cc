/* ----------------------------- METAHEURISTIC ----------------------------- */
// AUTHORS: Jofre Poch Soler, Martí Cortada Garcia

/* ------------------------- ALL REQUIRED LIBRARIES ------------------------- */
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <utility>
#include <vector>
#include <ctime>
#include <cstdlib>

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

/* ------------------------------- FUNCTIONS ------------------------------- */

/* Given a solution (vector of pieces), returns another solution belonging to its
neighborhood. A vector of pieces s' is considered a neighbour of s if either swapping
two pieces of the vector or inverting the position of one, s' = s.
*/
vector<Roll> random_neighbour(const vector<Roll>& initial_solution)
{   
    vector<Roll> neighbour_solution = initial_solution;
    srand((unsigned) time(0));
    
    // two random numbers between 0 and the number of pieces
    int pos1 = rand() % (initial_solution.size());
    int pos2 = rand() % (initial_solution.size());
    cout << pos1 << " " << pos2 << endl << endl;
    
    // if the two random numbers are equal, we just invert the coordinates
    // of the roll situated in that position
    if (pos1 == pos2) {
        neighbour_solution[pos1].p = initial_solution[pos1].q;
        neighbour_solution[pos1].q = initial_solution[pos1].p;
    }
    // we swap the rolls situated in the positions of the random numbers
    else {
        neighbour_solution[pos1] = initial_solution[pos2];
        neighbour_solution[pos2] = initial_solution[pos1];
    }

    return neighbour_solution;
}

void print_solution(const vector<Roll>& s) 
{
    for (auto roll : s) cout << roll.p << " " << roll.q << endl;
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

    vector<Roll> neighbour = random_neighbour(initial_solution);
    print_solution(initial_solution);
    cout << endl;
    print_solution(neighbour);

    // it will store the coordinates of every roll of the partial solutions
    // Coordinates coord(N);

    // boolean matrix that represents the canvas. The 1-cells mean that a roll has been placed there
    // and the 0-cells mean that that position is free.
    // It is initially defined with size W x max_length, to cover all the possible partial solutions
    // UsedCanvas B(max_length, vector<bool>(W, false));

    // it stores the coordinates of the left-upper corner of the last roll placed in the canvas
    // pair<int, int> current_coordinates = { 0, 0 };

    // boolean vector that shows which rolls have been used and which have not for every partial solution
    // vector<bool> rolls_b(N, false);
    // generate_canvas(argv[2], B, rolls, rolls_b, coord, current_coordinates, 0, 0);
}
