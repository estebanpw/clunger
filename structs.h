#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>

#define MAX_K 32
#define INIT_CLUSTERS 100000
#define READBUF 100000
#define MIN_SCORE 28
#define IGAP -5
#define EGAP -2

//Not used at the moment
//#pragma pack(push, 1)

typedef struct node_list{
	uint64_t pos;	//The position in the file
	uint64_t seq;	//The sequence to which the kmer belongs to
	struct node_list * next;	//The next repetition
} l_item;

typedef struct c_cluster{
	char prototype[MAX_K];	//The initial prototype
	l_item * reps;	//A pointer to the linked list holding the repetitions
	struct c_cluster * next;
} Cluster;

struct cell{
	int64_t score;
	uint64_t xe;
	uint64_t ye;
	uint64_t xs;
	uint64_t ys;
	uint64_t igaps;
	uint64_t egaps;
	uint64_t ident;
};