#include <Grid/Grid.h>


int main (int argc, char **argv) {

    // Initializing Grid library environment.
    Grid::Grid_init(&argc,&argv);








    // finalising Grid environment.
    Grid::Grid_finalize();

}
