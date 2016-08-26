#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "common.h"
#include "NWscore2rows.h"

int main(int argc, char ** av){
	
	if(argc != 4) terror("USE: buildcluster <input.fasta> <cluster.out> <seed_size>");
	
	
	//Database to read kmers from and output dictionary
	FILE * database, * dout;
	
	//Open database
	database = fopen64(av[1], "rt");
	if(database == NULL) terror("Could not open database");
	
	//Out database
	dout = fopen64(av[2], "wb");
	if(dout == NULL) terror("Could not open output database dictionary");
	
	//Get size of db
	fseeko64(database, 0L, SEEK_END);
	uint64_t dbByteSize = ftello64(database);
	fseeko64(database, 0L, SEEK_SET);

	//Get size of seeds
	uint16_t seedSize = (uint16_t) atoi(av[3]);	

	//Create initial cluster and a head pointer
	Cluster * clu = (Cluster *) get_mem_from_pool(sizeof(Cluster));
	if(clu == NULL) terror("Could not allocate initial clusters");
	Cluster * c_head = clu;
	Cluster * previous; //To traverse the list


	//Allocate rows for the alignments
	struct cell * r0 = (struct cell *) malloc(seedSize * sizeof(struct cell));
	struct cell * r1 = (struct cell *) malloc(seedSize * sizeof(struct cell));
	struct cell * mc = (struct cell *) malloc(seedSize * sizeof(struct cell));
	struct cell align_result;
	if(r0 == NULL || r1 == NULL || mc == NULL) printf("Could not allocate alignment rows");

	//Variables to read kmers
	char c = 'N'; //Char to read characters

	char wf[seedSize], wr[seedSize];
	wf[0]='\0';
	wr[0]='\0';
	
	//Variables to account for positions
	int strandF = 1, strandR = 0;
	uint64_t pos = 0, crrSeqL = 0, protoNum = 0; //Absolute coordinates; current sequence length; Number of used prototypes
	int64_t seqN = -1, alignScore; //Sequence Number; alignment score for a word against a prototype
	
	//Print info
	fprintf(stdout, "[INFO] Creating cluster prototypes using seed size=%"PRIu16"\n", seedSize);

	//Variables to read from buffer
	uint64_t posBuffer = READBUF + 1, tReadBuffer = 0;
	char *readBuffer = (char *) malloc(READBUF * sizeof(char));
	if (readBuffer == NULL) terror("Could not allocate memory for reading buffer");

	c = buffered_fgetc(readBuffer, &posBuffer, &tReadBuffer, database);
    while (!feof(database) || (feof(database) && posBuffer < tReadBuffer)){
		// Check if it's a special line
        if (!isupper(toupper(c))) { // Comment, empty or quality (+) line
            if (c == '>') { // Comment line
				c = buffered_fgetc(readBuffer, &posBuffer, &tReadBuffer, database);
                while (c != '\n') c = buffered_fgetc(readBuffer, &posBuffer, &tReadBuffer, database); //Avoid comment line
                seqN++; // New sequence
                crrSeqL = 0; // Reset buffered sequence length
                pos++; // Absolute coordinates: add one for the "*"
            }
	    	c = buffered_fgetc(readBuffer, &posBuffer, &tReadBuffer, database); //First char
        	continue;
        }

        //See if sequence has to be broken
        if(c == 'A' || c == 'C' || c == 'G' || c == 'T') crrSeqL++; else crrSeqL = 0;
        
        //Shift left and right here
        memcpy(wf, wf+1, seedSize-1); //Copy all characters one position to the left
        wf[seedSize-1] = c;
        wf[seedSize] = '\0';
        
        //TODO right shift
        
        pos++;

        if (crrSeqL >= (uint64_t) seedSize) { // Full well formed sequence
            if (strandF) {

            	//We have a word to align

            	if(protoNum == 0){ //First prototype
            		clu->prototype[seedSize] = '\0';
            		strncpy(clu->prototype, wf, seedSize);
            		clu->reps = (l_item *) get_mem_from_pool(sizeof(l_item));
            		clu->reps->pos = pos;
            		clu->reps->seq = seqN;
            		clu->next = NULL;
            		protoNum++;
            	}else{
            		//Get node references
					clu = c_head;
	            	previous = c_head;
	            	alignScore = 0;
	            	//While there are still prototypes and we have not got enough score to join a prototype
					while(alignScore < MIN_SCORE && clu != NULL){
						if(clu->prototype[0] != '\0'){
            				//If we can align it (first prototype)
							//ALIGN HERE
							align_result = NWscore2rows(wf, 0, seedSize, clu->prototype, 0, seedSize, IGAP, EGAP, mc, r0, r1);
							//printf("Aligned:\n%s\n%s\nWith score :%"PRId64" ident: %"PRIu64" len: %"PRIu64"\n", wf, clu->prototype, align_result.score, align_result.ident, align_result.xe-align_result.xs);
							//getchar();
            				if(alignScore > MIN_SCORE){
            					//Add to current
            					l_item * aux_rep = clu->reps->next;
            					clu->reps->next = (l_item *) get_mem_from_pool(sizeof(l_item));
            					clu->reps->pos = pos;
            					clu->reps->seq = seqN;
            					clu->reps->next = aux_rep;
            					break;
            				}

            			}
            			//Go to next node
            			previous = clu;
            			clu = clu->next;
	            	}

	            	if(clu == NULL){
	            		//Add new prototype
	            		clu = (Cluster *) get_mem_from_pool(sizeof(Cluster));
	            		strncpy(clu->prototype, wf, seedSize);
	            		clu->prototype[seedSize]='\0';
	            		previous->next = clu;
	            		clu->reps = (l_item *) get_mem_from_pool(sizeof(l_item));
	            		clu->reps->pos = pos;
	            		clu->reps->seq = seqN;
	            		clu->next = NULL;
	            		protoNum++;

	            	}
            	}


            }
            
        }
        long double computed = (long double)pos*100/dbByteSize;
        if(pos % 1000 == 0){
        	fprintf(stdout, "Processed %"PRIu64" bytes from %"PRIu64", that is, a %Le\n", pos, dbByteSize, computed);
        }
		c = buffered_fgetc(readBuffer, &posBuffer, &tReadBuffer, database);


	}
	fprintf(stdout, "[INFO] Sequence of length %"PRIu64" has %"PRIu64" mers of size k=%d\n", pos, pos-seedSize, seedSize);
	traverseClusters(c_head, dout);
	get_mem_from_pool(-1); //Deallocate
	free(r0);
	free(r1);
	free(mc);
	fclose(database);
	fclose(dout);
}