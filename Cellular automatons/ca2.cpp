#include <iostream>
#include <vector>
#include <queue>
#include <cstdlib>
#include <fstream>
#include <map>

using namespace std;

const int L = 50;
int grid[L][L] = {0};

// homokszem leejtése és lavina kezelése
int drop_sand() {
    int x = rand() % L;
    int y = rand() % L;
    grid[y][x]++;
    
    if (grid[y][x] < 4) return 0;

    int avalanche_size = 0;
    queue<pair<int, int>> q;
    q.push({x, y});

    while (!q.empty()) {
        auto [cx, cy] = q.front();
        q.pop();

        if (grid[cy][cx] >= 4) {
            grid[cy][cx] -= 4;
            avalanche_size++;

            int dx[] = {0, 0, 1, -1};
            int dy[] = {1, -1, 0, 0};
            for (int i = 0; i < 4; i++) {
                int nx = cx + dx[i];
                int ny = cy + dy[i];
                // Nyílt peremfeltétel
                if (nx >= 0 && nx < L && ny >= 0 && ny < L) {
                    grid[ny][nx]++;
                    if (grid[ny][nx] == 4) { 
                        q.push({nx, ny});
                    }
                }
            }
        }
    }
    return avalanche_size;
}

int main() {
    srand(42); // fix seed a reprodukálhatóságért
    map<int, int> avalanche_counts;
    
    int burn_in = 20000;
    int total_drops = 200000;

    cout << "Building sandpile..." << endl;
    for (int i = 0; i < burn_in; i++) drop_sand();

    cout << "Measuring avalanches..." << endl;
    for (int i = 0; i < total_drops; i++) {
        int size = drop_sand();
        if (size > 0) {
            avalanche_counts[size]++;
        }
    }

    ofstream file("datas/sandpile.data");
    file << "Size Count\n";
    for (auto const& [size, count] : avalanche_counts) {
        file << size << " " << count << "\n";
    }
    file.close();

    cout << "Done, results saved to datas/sandpile.data" << endl;
    return 0;
}