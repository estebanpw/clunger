/*********

File		NWscore2rows.c
Author		EPW <estebanpw@uma.es>
Description	Calculates a Needleman-Wunsch scores matrix using two rows and stores the number of identities and gaps, but does not retrieve the alignment.
		It is mainly inteded for use with reads vs genomes.

INPUT		<char * X>	Sequence X translated to numbers stored as one byte each letter
		<uint64_t Xstart>	Position in X sequence to start the alignment
		<uint64_t Xend>		Position in X sequence to end the alignment
		<char * Y>	Sequence Y translated to numbers stored as one byte each letter
		<uint64_t Xstart>	Position in Y sequence to start the alignment
		<uint64_t Xend>		Position in Y sequence to end the alignment
		<int iGap>		Penalty to open gap
		<int eGap>		Penalty to extend gap
		<int **PAM>		Matrix storing the match and miss rewards and penalties for each combination of letters
		<struct cell * ...>	Cell structs allocated from outside this function
		
RETURNS
		<struct cell bc>	Bottom cell with the best score

**********/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <unistd.h>
#include <inttypes.h>
#include <math.h>
#include <string.h>
#include "structs.h"


#define max(a,b)    (((a)>=(b)) ? (a):(b))
#define min(x,y)    (((x) < (y)) ? (x) : (y))


void goodPrint(int64_t a){ if(a>=0) printf("  %"PRId64"",a); else printf(" %"PRId64"",a);}
static int PAM[5][5]={5,-3,-3,-3,0,-3,5,-3,-3,0,-3,-3,5-3,0,-3,-3,-3,5,0,0,0,0,0,0};
int valOfNucl(char c){
	if(c=='A') return 0;
	if(c=='C') return 1;
	if(c=='G') return 2;
	if(c=='T') return 3;
	return 4;
}
/*

Calculates NW table with two rows and stores a cellpath of scores, identities, gaps and starting and ending positions

*/

struct cell NWscore2rows (char * X, uint64_t Xstart, uint64_t Xend, char * Y, uint64_t Ystart, uint64_t Yend, int iGap, int eGap, struct cell * mc, struct cell * f0, struct cell * f1){
    
    uint64_t i,j,k,iCounter=0,jCounter,kCounter,currEgapR,currEgapG;
	int64_t scoreDiagonal,scoreLeft,scoreRight,score;
	uint64_t offset = Yend-Ystart;
	struct cell * faux;
	int percentage=0;
	
    struct cell mf;
    

    if(mc == NULL || f0 == NULL || f1 == NULL){
    	printf("Could not allocate memory.\n");
    	exit(-1);
    }
    
    
    mf.score = INT64_MIN;
    

    //First row. iCounter serves as counter from zero
    //printf("..0%%");
    for(i=Ystart;i<Yend;i++){
    	f0[iCounter].score = PAM[valOfNucl(X[Xstart])][valOfNucl(Y[i])];
    	f0[iCounter].igaps = 0;
    	f0[iCounter].egaps = 0;
    	f0[iCounter].ident = (((f0[iCounter].score) > (0)) ? (1) : (0));
    	f0[iCounter].xs = Xstart;
    	f0[iCounter].ys = i;
    	f0[iCounter].xe = Xstart;
    	f0[iCounter].ye = i;

    	//Set every column max
    	mc[iCounter] = f0[iCounter];
    	//goodPrint(f0[iCounter].score);
    	
    	iCounter++;
	}
	
	//Set row max
	mf = f0[0];

	
	iCounter=1;

	//Go through full matrix with 2 rows
	for(i=Xstart+1;i<Xend;i++){
		//Fill first rowcell
		//printf("ROW: %"PRId64" ||",i);
		//printf("%.2f \n", ((float)i/Xend));
		
		f1[0].score = PAM[valOfNucl(X[i])][valOfNucl(Y[Ystart])];
		f1[0].xs = i;
		f1[0].ys = Ystart;
		f1[0].ident = (((PAM[valOfNucl(X[i])][valOfNucl(Y[Ystart])]) > (0)) ? (1) : (0));
		f1[0].igaps = 0;
		f1[0].egaps = 0;		
		f1[0].xe = i;
		f1[0].ye = Ystart;

		mf = f0[0];

		//goodPrint(f1[0].score);
		jCounter=1;
		for(j=Ystart+1;j<Yend;j++){
			//Check if max in row has changed
			if(jCounter > 1 && mf.score <= f0[jCounter-2].score){
				mf = f0[jCounter-2];
				mf.xe = i-1;
				mf.ye = j-2;
			}
			
			score = PAM[valOfNucl(X[i])][valOfNucl(Y[j])];
			scoreDiagonal = f0[jCounter-1].score + score;
			if(jCounter>1){
				scoreLeft = mf.score + iGap + (j - (mf.ye+1))*eGap + score;
				currEgapR = (j - (mf.ye+1));
				}else{
					scoreLeft = INT_MIN;
				}
				
			if(iCounter>1){
				scoreRight = mc[jCounter-1].score + iGap + (i - (mc[j-1].xe+1))*eGap + score;
				currEgapG =  (i - (mc[j-1].xe+1));
				}else{
					scoreRight = INT_MIN;
				}
			
			//Choose maximum
			//f1[jCounter] = max(max(scoreDiagonal,scoreLeft),scoreRight);
			
			if(scoreDiagonal >= max(scoreLeft, scoreRight)){
				//Diagonal
				f1[jCounter] = f0[jCounter-1];
				f1[jCounter].score = scoreDiagonal;
				if(PAM[valOfNucl(X[i])][valOfNucl(Y[j])] > 0) f1[jCounter].ident += 1;
								
			}else if(scoreRight >= scoreLeft){
				//Gap in genome
				f1[jCounter] = mc[jCounter-1];
				f1[jCounter].score = scoreRight;
				f1[jCounter].igaps += 1;
				f1[jCounter].egaps += currEgapG;
				
			}else{
				//Gap in read
				f1[jCounter] = mf;
				f1[jCounter].score = scoreLeft;
				f1[jCounter].igaps += 1;
				f1[jCounter].egaps += currEgapR;
			}

			//Update movement
			f1[jCounter].xe = i;
			f1[jCounter].ye = j;
			//goodPrint(f1[jCounter].score);
			jCounter++;
		}
		


		kCounter=0;
		for(k=Ystart;k<Yend;k++){
			//Update column maximum at j
			if(mc[kCounter].score <= f0[kCounter].score){
				mc[kCounter] = f0[kCounter];
				mc[kCounter].xe = i-1;
				mc[kCounter].ye = kCounter;
			}

			kCounter++;
		}
		//Switch rows
		
		faux = f0;
		f0 = f1;
		f1 = faux;
		
		iCounter++;
	}



    int64_t bestScore=f0[0].score, bestId = 0;
    for(k=Ystart;k<Yend;k++){
    	//printf("start (%"PRIu64",%"PRIu64") end (%"PRIu64",%"PRIu64") score [%"PRId64"] gaps [%"PRIu64"] ident [%"PRIu64"]\n", f0[k].xs, f0[k].ys, f0[k].xe, f0[k].ye, f0[k].score, f0[k].gaps, f0[k].ident);
    	if(f0[k].score >= bestScore){
    		bestScore = f0[k].score;
    		//printf("start (%"PRIu64",%"PRIu64") end (%"PRIu64",%"PRIu64") score [%"PRId64"] gaps [%"PRIu64"] ident [%"PRIu64"]\n", f0[k].xs, f0[k].ys, f0[k].xe, f0[k].ye, f0[k].score, f0[k].gaps, f0[k].ident);
    		bestId = k;
    	}
    }
	


    //FORCING GLOBAL ALIGNMENT
    //return f0[Yend-1];
    
    return f0[bestId];
}


/*
gcc NWscore2rows.c  -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -O3 -o nwalign
./nwalign x y 8 0 pamDNAident.txt 0
*/

/*
int main(int ac, char **av){
	
	FILE *f, *g;
	
	if(ac != 7) terror("USE: ./nwalign seqx seqy igap egap pamDNAident.txt show");

	int iGap = - atoi(av[3]);
	int eGap = - atoi(av[4]);
	int show = atoi(av[6]);


	char Ax[MAXALF],Ay[MAXALF];
	int Tx[MAXALF],Ty[MAXALF];
	int lAx,lAy;
	int **PAM;

	PAM = LoadScores(av[5], 1, Ax,Ay,&lAx,&lAy);
	Traductor(Ax,Tx);


	uint64_t offset = 1000;
	struct cell * mc = (struct cell *) malloc(offset * sizeof(struct cell));
    struct cell * f0 = (struct cell *) malloc(offset * sizeof(struct cell));
    struct cell * f1 = (struct cell *) malloc(offset * sizeof(struct cell));


	char X[50] = "TTCTAGTACACATAGAG"; //GENOME
	char Y[50] = "TTCTAGTT"; //READ
	uint64_t xlen = strlen(X);
	uint64_t ylen = strlen(Y);
	SeqToNum(Y, Tx, Y, strlen(Y));
	SeqToNum(X, Tx, X, strlen(X));
	


   	struct cell bc = NWscore2rows(X, 0, xlen, Y, 0, ylen, iGap, eGap, PAM, mc, f0, f1);
   	
   	printf("start (%"PRIu64",%"PRIu64") end (%"PRIu64",%"PRIu64") score [%"PRId64"] gaps [%"PRIu64"] ident [%"PRIu64"]\n", bc.xs, bc.ys, bc.xe, bc.ye, bc.score, bc.gaps, bc.ident);

   	free(mc);
	free(f0);
	free(f1);
	
	return 0;	
}

*/
