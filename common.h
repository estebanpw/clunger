#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include "structs.h"

void terror(const char * err_msg);
char buffered_fgetc(char *buffer, uint64_t *pos, uint64_t *read, FILE *f);
char * get_mem_from_pool(uint64_t bytes);
void traverseClusters(Cluster * head, FILE * out);