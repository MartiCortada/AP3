#include <iostream>
#include <fstream>
#include <assert.h>
#include <map>
#include <vector>

using namespace std;

typedef pair<int, int> P;
typedef map<P, int>    M;
typedef pair<P, P>     PP;
typedef vector<PP>     VPP;
typedef vector<VPP>    VVPP;

const PP UNDEF({-1, -1}, {-1, -1});

ifstream in;
int  W, L;
M    n;
VPP  sol;
VVPP board;


ostream& operator << (ostream& out, const P& z) {
  out << "(" << z.first << ", " << z.second << ")";
  return out;
}


bool within_bounds(const P& point) {
  return
    0 <= point.first  and point.first  < W and
    0 <= point.second and point.second < L;
}


bool solution_ok() {

  board = VVPP(L, VPP(W, UNDEF));

  int found_L = 0;
  
  for (auto tl_br : sol) {

    P tl = tl_br.first;
    if (not within_bounds(tl)) {
      cout << "Error: position " << tl << " is out of bounds" << endl;
      return false;
    }

    P br = tl_br.second;
    if (not within_bounds(br)) {
      cout << "Error: position " << br << " is out of bounds" << endl;
      return false;
    }

    found_L = max(found_L,br.second + 1);
    int p = br.first  - tl.first  + 1;
    int q = br.second - tl.second + 1;

    if (p <= 0 or q <= 0) {
      cout << "Error: top-left corner "  << tl
           << " and bottom-right corner " << br
           << " do not define a valid rectangle" << endl;
      return false;
    }

    P d = (p <= q ? P(p, q) : P(q, p));
    if (n.count(d) == 0) {
      cout << "Error: rectangle of dimensions " << d.first << "x" << d.second
           << " defined by top-left corner " << tl
           << " and bottom-right corner "    << br
           << " does not match any in input data" << endl;
      return false;
    }

    --n[d];

    if (n[d] < 0) {
      cout << "Error: too many rectangles of dimensions "
           << d.first << "x" << d.second << endl;
      return false;
    }

    for (int i = tl.second; i <= br.second; ++i)
      for (int j = tl.first; j <= br.first; ++j) {
	if (board[i][j] != UNDEF) {
          P tl2 = board[i][j].first;
          P br2 = board[i][j].second;
	  cout << "Error:"
	       << " rectangle defined by top-left corner "          << tl
               << " and bottom-right corner "                       << br
	       << " overlaps rectangle defined by top-left corner " << tl2
               << " and bottom-right corner "                       << br2
	       << " at position " << P(j, i) << endl;
          return false;
        }
	board[i][j] = {tl, br};
      }
  }

  if (found_L != L) {
    cout << "Error:"
         << " Solution file indicates L = " << L
         << " but the rectangles determine that L should be " << found_L << endl;
    return false;
  }
  
  for (auto di_ni : n) {
    auto di = di_ni.first;
    auto ni = di_ni.second;
    if (ni > 0) {
      cout << "Error: " << ni << " rectangle"
           << (ni == 1 ? "" : "s")
           << " of dimensions "
           << di.first << "x" << di.second
           << " missing" << endl;
      return false;
    }
  }
  return true;
}


void read_instance(const char* file) {
  ifstream in(file);
  int N;
  in >> W >> N;
  int ni, pi, qi;
  while (N != 0) {
    in >> ni >> pi >> qi;
    N -= ni;
    n[{pi, qi}] = ni;
  }
}


void read_solution(const char* file) {
  ifstream in(file);
  double time;
  in >> time;
  in >> L;
  int xtl, ytl, xbr, ybr;
  while (in >> xtl >> ytl >> xbr >> ybr)
    sol.push_back({{xtl, ytl}, {xbr, ybr}});
}


int main(int argc, char** argv) {

  // Write help message.
  if (argc == 1) {
    cout << "Makes a sanity check of a solution" << endl;
    cout << "Usage: " << argv[0] << " INPUT_FILE OUTPUT_FILE" << endl;
    exit(0);
  }

  assert(argc == 3);

  read_instance(argv[1]);
  read_solution(argv[2]);

  if (solution_ok()) {cout << "OK"     << endl; exit(0); }
  else               {cout << "FAILED" << endl; exit(1); }
}
