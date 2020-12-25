/* ---------------------------- GREEDY ALGORITHM ---------------------------- */
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
int l; // length of the configuration

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

/* ------------------------------- FUNCTIONS ------------------------------- */
// Criteria used to sort a vector of <Roll> elements in descending order of area
bool criteria(const Roll& r1, const Roll& r2)
{
    return r1.p * r1.q > r2.p * r2.q;
}

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
        file << coord[c].l.first << " " << coord[c].l.second << "   ";
        file << coord[c].r.first << " " << coord[c].r.second << endl;
    }
    file << endl;
    file.close();
}

/* Given a boolean matrix B and 4 integers representing the dimensions
of a piece "p x q" and the positions of the matrix (i,j), returs whether it is
possible or not to place that piece in the matrix. */
bool legal(const UsedCanvas& B, int p, int q, int i, int j)
{
    if (j + p > W or B[i][j])
        return false; // if it is already filled or exceeds the matrix dimensions
    // if B[i][j] == false, which means that place is not filled, then iterate
    // for all position that would ocupate the peace and check whether it's empty or not
    for (int r = i; r < i + q; ++r) {
        for (int s = j; s < j + p; ++s) {
            if (B[r][s])
                return false;
        }
    }
    return true; // We can place the piece!
}

/* Given a boolean matrix B and 4 integers representing the dimensions
of a piece "p x q" and the positions of the matrix (i,j), returns the same
matrix B, but with the postitions {i,...,i+q-1} & {j,...,j+p-1} filled, that
is, all these positions become true so no more peaces can be placed there.
Important to take into account:
    We will always call place_roll() after checking it is legal! */
UsedCanvas place_roll(UsedCanvas& B, int p, int q, int i, int j)
{
    for (int r = i; r < i + q; ++r) {
        for (int s = j; s < j + p; ++s) {
            B[r][s] = true;
        }
    }
    return B;
}

/*
Greedy algorithm:
    Given the sorted vector of <Roll> elements, a boolean matrix B and a vector of coordinates
    we analyze iterating for each piece of rolls if fits in B starting from the top-left corner,
    that is i = j = 0. If rolls[k] fits there, we place it and keep iterating for other pieces.
    We find whether a peace can be placed or not with the auxiliar function legal() and, if so,
    we place the piece in B with the function place_roll(). We save the coordinates of each piece
    storing the results in the position coord[k] of coord.
Important to take into account:
    We have to analyze two different possibilities. The piece can be placed as "p x q" or
    "q x p". So two different implementations will be considered and take the better of them!
*/
void greedy(string output, const vector<Roll>& rolls, UsedCanvas& B, Coordinates& coord)
{
    bool place1, place2; // will became true when a piece can be placed
    int l1, l2; //  lengths of the 2 different possibilites we will have (we will compare them)
    int i_copy, j_copy, i1_copy, j1_copy; // to keep the position (i,j) where the piece will be placed
    for (int k = 0; k < N; ++k) {
        place1 = place2 = false;
        l1 = l2 = l;
        i_copy = j_copy = i1_copy = j1_copy = 0;
        for (int i = 0; not place1 or not place2; ++i) {
            for (int j = 0; j < W and (not place1 or not place2); ++j) {
                if (not place1 and legal(B, rolls[k].p, rolls[k].q, i, j)) {
                    if (i < l)
                        l1 = max(i + rolls[k].q - 1, l);
                    else
                        l1 = i + rolls[k].q - 1; // case i >= l
                    // we find where to put the piece and save a copy of the coordinates
                    place1 = true;
                    i_copy = i;
                    j_copy = j;
                }
                // case where the inverted roll cannot be placed because it's length is larger than W
                if (rolls[k].q > W)
                    place2 = true;
                if (not place2 and legal(B, rolls[k].q, rolls[k].p, i, j)) {
                    if (i < l)
                        l2 = max(i + rolls[k].p - 1, l);
                    else
                        l2 = i + rolls[k].p - 1; // case i >= l
                    // we find where to put the piece and save a copy of the coordinates
                    place2 = true;
                    i1_copy = i;
                    j1_copy = j;
                }
            }
        }
        if (rolls[k].q > W)
            l = l1; // if the inverted roll cannot be placed, we take l1
        else
            l = min(l1, l2); // update l taking the minimum one
        // {place1 = place2 = true} -> 2 possibilites -> we proceed taking into acccount the length
        if (l1 < l2 or rolls[k].q > W) {
            B = place_roll(B, rolls[k].p, rolls[k].q, i_copy, j_copy);
            coord[k].l = { j_copy, i_copy };
            coord[k].r = { j_copy + rolls[k].p - 1, i_copy + rolls[k].q - 1 };
        } else { // l2 >= l1
            B = place_roll(B, rolls[k].q, rolls[k].p, i1_copy, j1_copy);
            coord[k].l = { j1_copy, i1_copy };
            coord[k].r = { j1_copy + rolls[k].q - 1, i1_copy + rolls[k].p - 1 };
        }
    }
    write_output(output, coord, l + 1);
}

/* We define a boolean matrix and the output coordinates vector that will be send
   in our algorithm along with the sorted <Roll> elements vector */
void greedy(string output, const vector<Roll>& rolls, int max_length)
{
    UsedCanvas B(max_length, vector<bool>(W, false)); // worst case -> dimension "max_length x W"
    Coordinates coord(N); // will be filled by the coordinates of our output
    l = 0; // length at the beginning (we recall it was declared as a global variable)
    greedy(output, rolls, B, coord); // greedy algorithm
}

/*
Main loop:
   it reads from a file all the data we need and store this data sorted in a vector as
   criteria specifies. Then, we send this data on a greedy algorithm in order to
   find a solution that, maybe, will not be the optimal!
*/
int main(int argc, char* argv[])
{
    // read input from a file
    ifstream f(argv[1]);
    f >> W >> N;

    vector<Roll> rolls; // data structure that stores all information about pieces
    int N_copy = N; // careful! We do not want to change N since it is global variable!
    int n, p, q; // n = number of orders of a specific type of dimensions "p x q"
    int max_length = 0;
    while (N_copy > 0) {
        f >> n >> p >> q;
        for (int i = 0; i < n; ++i)
            rolls.push_back({ p, q });
        max_length += n * max(p, q);
        N_copy -= n;
    }
    f.close();

    // we sort the <roll> elements vector using criteria as a way of sorting
    sort(rolls.begin(), rolls.end(), criteria); // check criteria() for more information!
    greedy(argv[2], rolls, max_length); // greedy algorithm
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
        if (initial_solution[pos].q <= W) {
            neighbour_solution[pos].p = initial_solution[pos].q;
            neighbour_solution[pos].q = initial_solution[pos].p;
        }
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


/* Given a solution (vector of pieces), returns another solution belonging to its
neighborhood. A vector of pieces s' is considered a neighbour of s if either swapping
two pieces of the vector or inverting the position of one, s' = s.
*/
vector<Roll> random_neighbour(const vector<Roll>& initial_solution)
{
    vector<Roll> neighbour_solution = initial_solution;

    // two random numbers between 0 and the number of pieces
    int pos1 = rand() % (initial_solution.size());
    int pos2 = rand() % (initial_solution.size());
    
     if (initial_solution[pos1].q <= W) {
        neighbour_solution[pos1].p = initial_solution[pos1].q;
        neighbour_solution[pos1].q = initial_solution[pos1].p;
    }
    neighbour_solution[pos1] = initial_solution[pos2];
    neighbour_solution[pos2] = initial_solution[pos1];
    

    return neighbour_solution;
}
