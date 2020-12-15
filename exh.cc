/* ---------------------- EXHAUSTIVE SEARCH ALGORITHM ---------------------- */
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
/* Given a vector of coordinates coord and the length L (both declared as global),
it prints the solution in the way we are asked (whose correctnes will be checked
with an alternative program) and sends it in a file called output.txt */
void write_output(string output)
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
    file << opt_res.L << endl;
    for (int c = 0; c < int(opt_res.coord.size()); ++c) {
        file << opt_res.coord[c].l.second << " " << opt_res.coord[c].l.first << "   ";
        file << opt_res.coord[c].r.second << " " << opt_res.coord[c].r.first << endl;
    }
    file << endl;
    file.close();
}

/* Given a boolean matrix B and 4 integers representing the dimensions
of a piece "p x q" and the positions of the matrix (i,j), returs whether it is
possible or not to place that piece in the matrix. */
bool legal(UsedCanvas& B, int p, int q, int i, int j)
{
    if (j + p > W or B[i][j])
        return false; // if it is already filled or exceeds the matrix dimensions
    // if B[i][j] == false, which means that place is not filled, then iterate
    // for all position that would ocupate the piece and check whether it's empty or not
    for (int r = i; r < i + q; ++r) {
        for (int s = j; s < j + p; ++s) {
            if (B[r][s])
                return false;
        }
    }
    return true; // We can place the piece!
}

/* Given a boolean matrix B, this function sets to True o False the cells occupied by the roll of size
p x q starting from the cell i, j. The input parameter 'option' specifies if the roll {p, q} has to be
marked ('1') or unmarked ('0') on B.
Important to take into account:
    We will always call mark_or_unmark_roll() after checking it is legal, so there won't be illegal
    positions of the matrix accessed. */
void mark_or_unmark_roll(UsedCanvas& B, int p, int q, pair<int, int>& current_coordinates, pair<int, int> initial_coordinates, bool option)
{
    int i = current_coordinates.first;
    int j = current_coordinates.second;
    if (option) {
        for (int r = i; r < i + q; ++r) {
            for (int s = j; s < j + p; ++s) {
                B[r][s] = true;
            }
        }
    } else { // option == False
        for (int r = i; r < i + q; ++r) {
            for (int s = j; s < j + p; ++s)
                B[r][s] = false;
        }
        // if the roll is unmarked on B, the variable current_coordinates
        // takes de value of the coordinates of the last placed roll
        current_coordinates = initial_coordinates;
    }
}

/* Function that places a roll into the first available position of the matrix B. A loop iterates from
left to right and from top to bottom until the roll is placed. Once a 'legal' position is found, all the
variables (current_coordinates, coord, current_L) are updated and the cells occupied by the
new piece are marked as true. */
void place_roll(const Roll& roll, UsedCanvas& B, Coordinates& coord, int& current_L,
    int k, pair<int, int>& current_coordinates, pair<int, int> initial_coordinates)
{
    bool placed = false; // is set to True when the roll is placed
    for (int i = 0; not placed; ++i) {
        for (int j = 0; j < W and not placed; ++j) {
            // check if the current position i, j of the matrix B and the rest of the cells that would
            // occupy the piece are available (False)
            if (legal(B, roll.p, roll.q, i, j)) {
                current_coordinates = { i, j };
                // the cells of B that occupie the placed roll are set to True
                mark_or_unmark_roll(B, roll.p, roll.q, current_coordinates, initial_coordinates, 1);
                // the coordinates of the partial solution are filled with the top-left and
                // the bottom right position of the roll on the matrix
                coord[k].l = current_coordinates;
                coord[k].r = { i + roll.q - 1, j + roll.p - 1 };
                // current_L is updated only if the lower position of the placed roll is larger than it
                current_L = max(current_L, i + roll.q);
                placed = true;
            }
        }
    }
}

/* It generates all possible combinations of pieces (dimensions "p x q" & "q x p") and checks whether we
can obtain or not an optimal solution recursively and using external funcions (check place_roll() &
mark_or_unmark_roll() for more information). At each iteration we use the current lenght current_L to
determine if we are generating an optimal solution or not, that is, it is our decision variable. */
void generate_canvas(string output, UsedCanvas& B, const vector<Roll>& rolls, vector<bool>& rolls_b, Coordinates& coord,
    pair<int, int> current_coordinates, int current_L, int k)
{
    if (current_L >= opt_res.L)
        return; // optimitzation
    if (current_L < opt_res.L) {
        if (k == N) {
            opt_res.L = current_L;
            opt_res.coord = coord;
            write_output(output); // we reach a possible optimal solution
        } else {
            int initial_L = current_L;
            pair<int, int> initial_coordinates = current_coordinates;
            for (int s = 0; s < N; ++s) {
                if (not rolls_b[s]) {
                    rolls_b[s] = true;

                    // generating normal orientation
                    place_roll(rolls[s], B, coord, current_L, k, current_coordinates, initial_coordinates);
                    generate_canvas(output, B, rolls, rolls_b, coord, current_coordinates, current_L, k + 1);
                    current_L = initial_L;
                    mark_or_unmark_roll(B, rolls[s].p, rolls[s].q, current_coordinates, initial_coordinates, 0);

                    // generating inverted orientation
                    if (rolls[s].q <= W and rolls[s].q != rolls[s].p) {
                        place_roll({ rolls[s].q, rolls[s].p }, B, coord, current_L, k, current_coordinates, initial_coordinates);
                        generate_canvas(output, B, rolls, rolls_b, coord, current_coordinates, current_L, k + 1);
                        current_L = initial_L;
                        mark_or_unmark_roll(B, rolls[s].q, rolls[s].p, current_coordinates, initial_coordinates, 0);
                    }
                    rolls_b[s] = false;
                }
            }
        }
    }
}

/*
Main loop:
   it reads from a file all the data we need and store this data sorted in a vector.
   Then, we send this data on an exhaustive search algorithm in order to find a solution
   that has to be the optimal!
*/
int main(int argc, char* argv[])
{
    // read input from a file
    ifstream f(argv[1]);
    f >> W >> N;

    vector<Roll> rolls; // data structure that stores all information about pieces
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

    // boolean matrix that represents the canvas. The 1-cells mean that a roll has been placed there
    // and the 0-cells mean that that position is free.
    // It is initially defined with size W x max_length, to cover all the possible partial solutions
    UsedCanvas B(max_length, vector<bool>(W, false));

    // it will store the coordinates of every roll of the partial solutions
    Coordinates coord(N);

    // it stores the coordinates of the left-upper corner of the last roll placed in the canvas
    pair<int, int> current_coordinates = { 0, 0 };

    // boolean vector that shows which rolls have been used and which have not for every partial solution
    vector<bool> rolls_b(N, false);
    generate_canvas(argv[2], B, rolls, rolls_b, coord, current_coordinates, 0, 0);
}
