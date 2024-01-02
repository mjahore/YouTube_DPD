// Allocate these arrays in the main code.
int *cell_head, *cell_list;

// These are the indices used to calculate the locations of all 
// 14 neighboring cells when we process interactions.
int nix[] = {0,-1,-1,-1,0,0,-1,1,-1, 0, 1,-1,0,1};
int niy[] = {0, 0,-1, 1,1,0, 0,0,-1,-1,-1, 1,1,1};
int niz[] = {0, 0, 0, 0,0,1, 1,1, 1, 1, 1, 1,1,1};
