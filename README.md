# Orthogonal-Hamiltonian-path

My interest in space-filling curves led me to take a closer look at higher-order iterations of the 3D Hilbert and Moore curves. While examining their structure, I noticed that they often contain collinear segments. This led me to explore the possibility of constructing a path through all the points of a 3D grid in which the direction changes strictly by 90 degrees at each step — that is, with no consecutive collinear edges.

In essence, the goal was to find a Hamiltonian path that satisfies an additional constraint: a right-angle turn at every step. It turned out that many such paths are possible, which naturally raised the next question: how can such a path be defined procedurally for any even-sized grid? (For grids with odd dimensions, such solutions do not exist.)

Initial attempts involved recursive algorithms with backtracking and greedy search methods. However, these approaches either failed to achieve full coverage or became computationally expensive on large grids.

Eventually, I developed a deterministic algorithm capable of generating such a path efficiently on grids of any even size, with guaranteed full coverage and strict 90-degree turns.

The implementation was done in Houdini by VEX code. The entire process is handled in a single Detail Wrangle: first generating the grid points, then computing the path, and finally building a polyline. The grid resolution is defined by a promoted "count" parameter. The code is provided below.
#
<p align="center">
<img src="https://github.com/user-attachments/assets/a8041424-aca0-4e11-83e3-62e339d76cb4">
<img src="https://github.com/user-attachments/assets/4cee8947-780c-4f81-b7ea-b94f12136501">
</p>

#
```c
int count = chi("count");
int grid_size = count - (count % 2); 

int getPositionIndex(vector pos; int grid_size) {
    int x = int(pos.x);
    int y = int(pos.y); 
    int z = int(pos.z);
    
    if (x < 0 || x >= grid_size || y < 0 || y >= grid_size || z < 0 || z >= grid_size) 
        return -1;
    
    return x * grid_size * grid_size + y * grid_size + z;
}

int point_indices[];
resize(point_indices, grid_size * grid_size * grid_size);

int idx = 0;
for(int i = 0; i < grid_size; i++) {
    for(int j = 0; j < grid_size; j++) {
        for(int k = 0; k < grid_size; k++) {
            int pt_idx = addpoint(0, set(i, j, k));
            point_indices[idx] = pt_idx;
            idx++;
        }
    }
}

int npts = len(point_indices);
int lvl = 0;
int path[] = array();
vector pos = {1, 0, 0};
vector prev_dir = {0, 1, 0};
vector dirs[] = {{0,1,0}, {0,-1,0}, {0,0,1}, {0,0,-1}, {1,0,0}, {-1,0,0}};

int current_idx = getPositionIndex(pos, grid_size);
if (current_idx == -1) {
    warning("Starting position is invalid!");
} else {
    append(path, point_indices[current_idx]);

    int visited[];
    resize(visited, npts);
    visited[current_idx] = 1;

    while (len(path) < npts) {
        
        int x = current_idx / (grid_size * grid_size);
        int y = (current_idx / grid_size) % grid_size;
        int z = current_idx % grid_size;
        pos = set(x, y, z);
        
        int check = 0;
         
        for (int i = 0; i < 2 && !check; i++) {
            foreach (vector d; dirs) {
                
                vector nbr_pos = pos + d;
                int nbrpt = getPositionIndex(nbr_pos, grid_size);
                
                if (nbrpt == -1 || visited[nbrpt] || lvl - nbr_pos.y <= -2) 
                    continue;
                
                int condition_met = 0;
                if (i == 0 && (nbr_pos.x == 0 || nbr_pos.z == grid_size - 1 || nbr_pos.x == grid_size - 1)) {              
                    condition_met = 1;
                } else if (i == 1) {   
                    condition_met = 1;
                }
                
                if (condition_met) {
                    current_idx = nbrpt;
                    prev_dir = d;
                    append(path, point_indices[current_idx]);
                    visited[current_idx] = 1;
                    check = 1;
                    break;
                }
            }
        }
        
        if (!check) {
            lvl += 2;
            if (lvl > cbrt(npts)) break;
        }
    }

    if (len(path) > 1) {
        addprim(0, "polyline", path);
    }
}
```
