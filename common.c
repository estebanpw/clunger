#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include "structs.h"

void terror(const char * err_msg){
	fprintf(stderr, "\nERR: %s\n", err_msg);
	exit(-1);
}

/* Buffered reader to read char by char
 *  @param buffer: An already allocated buffer of chars
 *  @param pos: Current position of last char read in buffer (should be initialized to READBUF+1
 *  @param read: Number of chars read at iteration
 *	@param f: File descriptor from which to read
 *	@return: The read char
 */

char buffered_fgetc(char *buffer, uint64_t *pos, uint64_t *read, FILE *f) {
    if (*pos >= READBUF) {
        *pos = 0;
        memset(buffer, 0, READBUF);
        *read = fread(buffer, 1, READBUF, f);
    }
    *pos = *pos + 1;
    return buffer[*pos-1];
}

/* Memory pool: allocates space once to later on give pointer positions
   	@param bytes: Requested bytes
   	@return: A void pointer that points to a position with the requested bytes
*/

char * get_mem_from_pool(int64_t bytes){
	static uint64_t bytes_used = 0;
	static char * mem;
	if(bytes_used == 0){
		mem = (char *) malloc(INIT_CLUSTERS * sizeof(Cluster));
		if(mem == NULL) terror("Could not allocate requested memory block");
	}
	if(bytes == -1){
		//Free
		free(mem-bytes_used);
		return NULL;
	}

	bytes_used += bytes;
	if(bytes_used > INIT_CLUSTERS*sizeof(Cluster)) terror("Increase memory for clusters");
	mem += bytes;

	return (mem-bytes);
}

/* Traverses the list of clusters showing the prototypes and the repetitions
	@head:	The head cluster

*/

void traverseClusters(Cluster * head, FILE * out){
	fprintf(out, "======TRAVERSING CLUSTERS=====");
	l_item * l;
	Cluster * c = head;
	while(c != NULL){
		fprintf(out, "\n>%s@", c->prototype);
		
		l = c->reps;
		while(l != NULL){
			fprintf(out, "(%"PRIu64", %"PRIu64")", l->pos, l->seq);
			l = l->next;
		}
		
		c = c->next;
	}
	fprintf(out, "\n======DONE=====\n");
}