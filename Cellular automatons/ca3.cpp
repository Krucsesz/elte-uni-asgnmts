#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <cstdlib>

using namespace std;

const int L = 4096; 

void box_count(const vector<vector<int>>& grid, const string& filename) {
    ofstream file(filename);
    file << "s N\n"; // s = dobozméret, N = dobozok száma
    
    // dobozméret + 2 hatvány
    for (int s = 1; s <= L / 2; s *= 2) {
        int count = 0;
        for (int y = 0; y < L; y += s) {
            for (int x = 0; x < L; x += s) {
                bool has_point = false;
                // élők keresése
                for (int dy = 0; dy < s && y + dy < L; ++dy) {
                    for (int dx = 0; dx < s && x + dx < L; ++dx) {
                        if (grid[y + dy][x + dx] > 0) {
                            has_point = true;
                            break;
                        }
                    }
                    if (has_point) break;
                }
                if (has_point) count++;
            }
        }
        file << s << " " << count << "\n";
    }
    file.close();
}

// 1. 1D sejtautomata
void run_1d_ca() {
    vector<vector<int>> grid(L, vector<int>(L, 0));
    grid[0][L/2] = 1; // Középről indul 1 sejt
    for (int t = 1; t < L; ++t) {
        for (int x = 1; x < L - 1; ++x) {
            grid[t][x] = grid[t-1][x-1] ^ grid[t-1][x+1];
        }
    }
    box_count(grid, "box_1d.data");
    cout << "1D automaton box counting done." << endl;
}

// 2. erdőtűz
void run_forest_fire() {
    vector<vector<int>> grid(L, vector<int>(L, 0));
    double p = 0.05, f = 0.001; // fa növése és villámcsapás

    for(int i=0; i<L; ++i) for(int j=0; j<L; ++j) grid[i][j] = (rand()%100 < 50) ? 1 : 0;

    for (int step = 0; step < 200; ++step) {
        vector<vector<int>> next_grid = grid;
        for (int y = 0; y < L; ++y) {
            for (int x = 0; x < L; ++x) {
                if (grid[y][x] == 2) {
                    next_grid[y][x] = 0;
                } else if (grid[y][x] == 1) {
                    bool fire_neighbor = false;
                    int dx[] = {0,0,1,-1}, dy[] = {1,-1,0,0};
                    for(int k=0; k<4; ++k) {
                        int nx = x+dx[k], ny = y+dy[k];
                        if(nx>=0 && nx<L && ny>=0 && ny<L && grid[ny][nx]==2) fire_neighbor = true;
                    }
                    if (fire_neighbor || (rand() / (double)RAND_MAX) < f) next_grid[y][x] = 2;
                } else {
                    if ((rand() / (double)RAND_MAX) < p) next_grid[y][x] = 1;
                }
            }
        }
        grid = next_grid;
    }

    // csak a tüzet nézzük
    vector<vector<int>> fire_grid(L, vector<int>(L, 0));
    for(int i=0; i<L; ++i) for(int j=0; j<L; ++j) if(grid[i][j]==2) fire_grid[i][j]=1;

    box_count(fire_grid, "box_forest.data");
    cout << "Forest fire box counting done." << endl;
}

// 3. életjáték
void run_gol() {
    vector<vector<int>> grid(L, vector<int>(L, 0));
    for(int i=0; i<L; ++i) for(int j=0; j<L; ++j) grid[i][j] = (rand()%100 < 20) ? 1 : 0;

    for (int step = 0; step < 100; ++step) {
        vector<vector<int>> next_grid = grid;
        for (int y = 0; y < L; ++y) {
            for (int x = 0; x < L; ++x) {
                int neighbors = 0;
                for (int dy = -1; dy <= 1; ++dy) {
                    for (int dx = -1; dx <= 1; ++dx) {
                        if (dx == 0 && dy == 0) continue;
                        int nx = (x + dx + L) % L;
                        int ny = (y + dy + L) % L;
                        neighbors += grid[ny][nx];
                    }
                }
                if (grid[y][x] == 1 && (neighbors < 2 || neighbors > 3)) next_grid[y][x] = 0;
                else if (grid[y][x] == 0 && neighbors == 3) next_grid[y][x] = 1;
            }
        }
        grid = next_grid;
    }
    box_count(grid, "box_gol.data");
    cout << "Game of Life box counting done." << endl;
}

int main() {
    srand(42);
    cout << "Simulations starting..." << endl;
    
    run_1d_ca();
    run_forest_fire();
    run_gol();
    
    cout << "All data saved to files." << endl;
    return 0;
}