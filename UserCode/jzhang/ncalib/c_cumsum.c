/*
c_cumsum.c

*/

void c_cumsum (unsigned char* array, int m, int n, int l, int axis) {

    int i, j, k;
    int index = 0 ;

    if (axis == 2) {
        for (i = 0; i < m; i++) {
            for (j = 0; j < n; j++) {
                for (k = 1; k < l; k++ ) {
                    index = i * n * l + j * l + k;
                    array[index] += array[index - 1];
                }
            }
        }
    }
    else if (axis == 1) {
        for (i = 0; i < m; i++) {
            for (k = 0; k < l; k++ ) {
                for (j = 1; j < n; j++) {
                    index = i * n * l + j * l + k;
                    array[index] += array[index - l];
                }
            }
        }
    }
    else if (axis == 0) {
        for (j = 0; j < n; j++) {
            for (k = 0; k < l; k++ ) {
                for (i = 1; i < m; i++) {
                    index = i * n * l + j * l + k;
                    array[index] += array[index - n * l];
                }
            }
        }
    }

    return;
}
