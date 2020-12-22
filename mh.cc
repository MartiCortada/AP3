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
