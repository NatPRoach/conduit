
#include "default.h"
#include "seq_util.h"



/** randomizes seq[] by shuffling, and places the result in randseq[];
 if randseq[] and seq[] are distinct, seq[] is left unchanged */
void shuffle_seq(int len,
		char seq[],
		char randseq[])
{
  int i,j;
  char c;

  for (i=0;i<len;i++)  {
    j=rand()%len; /* CHOOSE RANDOM POSITION j */
    c=seq[i]; /* SWAP POSITIONS i AND j */
    randseq[i]=seq[j];
    randseq[j]=c;
  }

  return;
}






/**@memo TRANSLATE FROM ASCII LETTERS TO COMPACTED NUMBERICAL INDEX: 
    index_symbols(seq[i].length,seq[i].sequence,seq[i].sequence,
		  m->nsymbol,m->symbol);
*/
/** converts characters in seq[] to the INDEX of the matching character in
 symbols[], and returns the result in out[] */
void index_symbols(int nseq,char seq[],char out[],
		   int nsymbs,char symbols[])
{
  int i,j,k;
  LOOP (i,nseq) {
    k=nsymbs-1; /* DEFAULT: UNMATCHABLE SYMBOL */
    LOOP (j,nsymbs) {  /* FIND MATCHING SYMBOL */
      if (symbols[j]==seq[i]) {  /* FOUND IT! */
	k=j;
	break;
      }
    }
    out[i]=k; /* SAVE THE TRANSLATED CODE */
  }
  return;
}






int *Score_matrix_row=NULL;

int best_match_qsort_cmp(const void *void_a,const void *void_b)
{
  int *a=(int *)void_a,*b=(int *)void_b;

  if (Score_matrix_row[*a]>Score_matrix_row[*b])
    return -1;
  else if (Score_matrix_row[*a]<Score_matrix_row[*b])
    return 1;
  else /* EQUAL */
    return 0;
}




#ifdef SOURCE_EXCLUDED
char DNA_symbols[1024];
float DNA_rescale_score;
#endif

/** reads an alignment scoring matrix in the pam format */
int read_score_matrix(char filename[],ResidueScoreMatrix_T *m)
{
  int i,j,k,nsymb=0,found_symbol_line=0,isymb;
  char line[1024],dna_codes[256];
  FILE *ifile;

   /* GAP PENALTY DEFAULTS */
  m->gap_penalty_set[0][0]=m->gap_penalty_set[1][0]=12; /*SAVE PENALTIES*/
  m->gap_penalty_set[0][1]=m->gap_penalty_set[1][1]=2;
  m->gap_penalty_set[0][2]=m->gap_penalty_set[1][2]=0;
  m->trunc_gap_length = TRUNCATE_GAP_LENGTH;
  m->decay_gap_length = DECAY_GAP_LENGTH;
  
  ifile=fopen(filename,"r");
  if (!ifile) {
    WARN_MSG(USERR,(ERRTXT,"Can't open alignment matrix from %s\n",filename),"$Revision: 1.2.2.2 $");
    return -2; /* FAILED TO FIND FILE TO READ */
  }

  while (fgets(line,1023,ifile)) {
    if ('#'==line[0] || '\n'==line[0]) /* SKIP COMMENT OR BLANK LINES */
      continue;
    
    else if (1==sscanf(line,"GAP-TRUNCATION-LENGTH=%d",&i)) {
      m->trunc_gap_length = i;
    }
    
    else if (1==sscanf(line,"GAP-DECAY-LENGTH=%d",&i)) {
      m->decay_gap_length = i;
    }
    
    else if (3==sscanf(line,"GAP-PENALTIES=%d %d %d",&i,&j,&k)) {
      m->gap_penalty_set[0][0]=m->gap_penalty_set[1][0]=i; /*SAVE PENALTIES*/
      m->gap_penalty_set[0][1]=m->gap_penalty_set[1][1]=j;
      m->gap_penalty_set[0][2]=m->gap_penalty_set[1][2]=k;
    }

    else if (3==sscanf(line,"GAP-PENALTIES-X=%d %d %d",&i,&j,&k)) {
      m->gap_penalty_set[1][0]=i; /*SAVE PENALTIES ONLY FOR X DIRECTION*/
      m->gap_penalty_set[1][1]=j;
      m->gap_penalty_set[1][2]=k;
    }

#ifdef SOURCE_EXCLUDED
    else if (1==sscanf(line,"DNACODES=%99s",dna_codes)) { /* READ DNACODES*/
      strcpy(DNA_symbols,dna_codes);/*SYMBOLS COUNTED AS DNA FOR AUTORECOG*/
    }

    else if (1==sscanf(line,"DNASCALE=%f",&DNA_rescale_score))
      continue;
#endif
    
    else if (!found_symbol_line) { /* READ THIS LINE AS LIST OF SEQ SYMBOLS*/
      for (i=0;'\0'!=line[i];i++)
	if (!isspace(line[i])) /* IGNORE WHITESPACE */
	  m->symbol[nsymb++]=line[i]; /* SAVE TO LIST OF SYMBOLS */
      found_symbol_line=1; /* SET FLAG SO WE NOW READ MATRIX SCORE VALUES */
    }
    
    else { /* READ SCORING MATRIX LINES */
      found_symbol_line=0; /* DEFAULT: FAILED TO FIND MATCHING SYMBOL IN LIST*/
      LOOP (isymb,nsymb) /* FIND MATCH TO THIS SYMBOL */
	if (m->symbol[isymb]==line[0]) {
	  found_symbol_line=1; /* SIGNAL THAT WE SUCCESFULLY FOUND MATCH */
	  j=1; /* SKIP FIRST CHARACTER: OUR SEQUENCE SYMBOL */
	  LOOPF (i,nsymb) { /* READ ALL THE SCORE VALUES ON THIS LINE */
	    if (1==sscanf(line+j,"%d%n",&(m->score[isymb][i]),&k))
	      j+=k; /* ADVANCE THE READING POSITION */
	    else { /* MISSING SCORE DATA: ERROR! */
	      IF_GUARD(1,5.23,(ERRTXT,"Missing score value for pair %c:%c",
			      m->symbol[isymb],m->symbol[i]),TRAP)
		;
	      fclose(ifile); /* CLOSE OUR STREAM */
	      return -1;
	    }
	  }
	  break;
	}
      IF_GUARD(!found_symbol_line,1.5,(ERRTXT,"Missing or unknown sequence symbol: %c",line[0]),TRAP) { /* ERROR: AN INVALID SYMBOL, NOT IN LIST */
	fclose(ifile); /* CLOSE OUR STREAM */
	return -1;
      }
    }
  }
  fclose(ifile);

  /* CONSTRUCT GAP PENALTY ARRAYS FROM GAP PARAMETERS: */
  m->max_gap_length = m->trunc_gap_length + m->decay_gap_length;
  CALLOC (m->gap_penalty_x, m->max_gap_length+2, LPOScore_T);
  CALLOC (m->gap_penalty_y, m->max_gap_length+2, LPOScore_T);

  /*** GAP OPENING PENALTY @ L=0->1 */
  m->gap_penalty_x[0] = m->gap_penalty_set[0][0];
  m->gap_penalty_y[0] = m->gap_penalty_set[1][0];
  
  /*** 1st AFFINE EXTENSION PENALTY (A1) @ L=1->2,2->3,...T-1->T */
  for (i=1;i<m->trunc_gap_length;i++) {
    m->gap_penalty_x[i] = m->gap_penalty_set[0][1];
    m->gap_penalty_y[i] = m->gap_penalty_set[1][1];
  }
  
  /*** DECAYING EXTENSION PENALTY (A1-->A2; skipped if D=0) @ L=T->T+1,...T+D-1->T+D */
  for (i=0;i<m->decay_gap_length;i++) {
    double dec_x = (m->gap_penalty_set[0][1] - m->gap_penalty_set[0][2]) / ((double)(m->decay_gap_length + 1));
    double dec_y = (m->gap_penalty_set[1][1] - m->gap_penalty_set[1][2]) / ((double)(m->decay_gap_length + 1));
    m->gap_penalty_x[i+m->trunc_gap_length] = m->gap_penalty_set[0][1] - (i+1) * dec_x;
    m->gap_penalty_y[i+m->trunc_gap_length] = m->gap_penalty_set[1][1] - (i+1) * dec_y;
  }
  
  /*** 2nd AFFINE EXTENSION PENALTY (A2) @ L>=T+D */
  m->gap_penalty_x[m->max_gap_length] = m->gap_penalty_set[0][2];
  m->gap_penalty_y[m->max_gap_length] = m->gap_penalty_set[1][2];
  
  m->gap_penalty_x[m->max_gap_length+1] = 0;  /* DON'T REMOVE THIS!... SPECIAL STATE USED IN align_lpo. */
  m->gap_penalty_y[m->max_gap_length+1] = 0;  /* DON'T REMOVE THIS!... SPECIAL STATE USED IN align_lpo. */
  
  
  LOOPF (i,nsymb) {
    Score_matrix_row= m->score[i]; /* ROW TO USE FOR SORTING best_match */
    LOOP (j,nsymb)
      m->best_match[i][j] = j;
    qsort(m->best_match[i],nsymb,sizeof(m->best_match[0][0]),
	  best_match_qsort_cmp);
#ifdef SOURCE_EXCLUDED
    printf("%c SORT",m->symbol[i]); /* TEST: PRINT OUT SORTED TABLE */
    LOOPF (j,nsymb)
      printf("\t%c:%d",m->symbol[m->best_match[i][j]],
	     m->score[i][m->best_match[i][j]]);
    printf("\n");
#endif
  }

  m->symbol[nsymb]='\0'; /* TERMINATE THE SYMBOL STRING */
  m->nsymbol=nsymb;
  return nsymb;
}


/*
New function default_score_matrix below written by Nathan Roach 03 / 2020;
Labeled as new code in accordance with the GPLv2 license in which poaV2 is distributed
*/
void default_score_matrix(ResidueScoreMatrix_T *m){
  int i;
  m->symbol[ 0] = 'A';
  m->symbol[ 1] = 'T';
  m->symbol[ 2] = 'G';
  m->symbol[ 3] = 'C';
  m->symbol[ 4] = 'U';
  m->symbol[ 5] = 'S';
  m->symbol[ 6] = 'W';
  m->symbol[ 7] = 'R';
  m->symbol[ 8] = 'Y';
  m->symbol[ 9] = 'K';
  m->symbol[10] = 'M';
  m->symbol[11] = 'B';
  m->symbol[12] = 'V';
  m->symbol[13] = 'H';
  m->symbol[14] = 'D';
  m->symbol[15] = 'N';
  m->symbol[16] = '\0'; /* TERMINATE THE SYMBOL STRING */
  m->nsymbol = 16;
  m->score[ 0][ 0] =  0; m->score[ 0][ 1] = -8; m->score[ 0][ 2] = -8; m->score[ 0][ 3] = -8; m->score[ 0][ 4] = -8; m->score[ 0][ 5] = -4; m->score[ 0][ 6] =  1; m->score[ 0][ 7] =  1; m->score[ 0][ 8] = -4; m->score[ 0][ 9] = -4; m->score[ 0][10] =  1; m->score[ 0][11] = -4; m->score[ 0][12] = -1; m->score[ 0][13] = -1; m->score[ 0][14] = -1; m->score[ 0][15] = -2;
  m->score[ 1][ 0] = -8; m->score[ 1][ 1] =  3; m->score[ 1][ 2] = -8; m->score[ 1][ 3] = -8; m->score[ 1][ 4] =  3; m->score[ 1][ 5] = -4; m->score[ 1][ 6] =  1; m->score[ 1][ 7] = -4; m->score[ 1][ 8] =  1; m->score[ 1][ 9] =  1; m->score[ 1][10] = -4; m->score[ 1][11] = -1; m->score[ 1][12] = -4; m->score[ 1][13] = -1; m->score[ 1][14] = -1; m->score[ 1][15] = -2;
  m->score[ 2][ 0] = -8; m->score[ 2][ 1] = -8; m->score[ 2][ 2] =  3; m->score[ 2][ 3] = -8; m->score[ 2][ 4] = -8; m->score[ 2][ 5] =  1; m->score[ 2][ 6] = -4; m->score[ 2][ 7] =  1; m->score[ 2][ 8] = -4; m->score[ 2][ 9] =  1; m->score[ 2][10] = -4; m->score[ 2][11] = -1; m->score[ 2][12] = -1; m->score[ 2][13] = -4; m->score[ 2][14] = -1; m->score[ 2][15] = -2;
  m->score[ 3][ 0] = -8; m->score[ 3][ 1] = -8; m->score[ 3][ 2] = -8; m->score[ 3][ 3] =  3; m->score[ 3][ 4] = -8; m->score[ 3][ 5] =  1; m->score[ 3][ 6] = -4; m->score[ 3][ 7] = -4; m->score[ 3][ 8] =  1; m->score[ 3][ 9] = -4; m->score[ 3][10] =  1; m->score[ 3][11] = -1; m->score[ 3][12] = -1; m->score[ 3][13] = -1; m->score[ 3][14] = -4; m->score[ 3][15] = -2;
  m->score[ 4][ 0] = -8; m->score[ 4][ 1] =  3; m->score[ 4][ 2] = -8; m->score[ 4][ 3] = -8; m->score[ 4][ 4] =  3; m->score[ 4][ 5] = -4; m->score[ 4][ 6] = -4; m->score[ 4][ 7] = -4; m->score[ 4][ 8] = -4; m->score[ 4][ 9] = -4; m->score[ 4][10] = -4; m->score[ 4][11] = -4; m->score[ 4][12] = -4; m->score[ 4][13] = -4; m->score[ 4][14] = -4; m->score[ 4][15] = -4;
  m->score[ 5][ 0] = -4; m->score[ 5][ 1] = -4; m->score[ 5][ 2] =  1; m->score[ 5][ 3] =  1; m->score[ 5][ 4] = -4; m->score[ 5][ 5] = -1; m->score[ 5][ 6] = -4; m->score[ 5][ 7] = -2; m->score[ 5][ 8] = -2; m->score[ 5][ 9] = -2; m->score[ 5][10] = -2; m->score[ 5][11] = -1; m->score[ 5][12] = -1; m->score[ 5][13] = -3; m->score[ 5][14] = -3; m->score[ 5][15] = -1;
  m->score[ 6][ 0] =  1; m->score[ 6][ 1] =  1; m->score[ 6][ 2] = -4; m->score[ 6][ 3] = -4; m->score[ 6][ 4] = -4; m->score[ 6][ 5] = -4; m->score[ 6][ 6] = -1; m->score[ 6][ 7] = -2; m->score[ 6][ 8] = -2; m->score[ 6][ 9] = -2; m->score[ 6][10] = -2; m->score[ 6][11] = -3; m->score[ 6][12] = -3; m->score[ 6][13] = -1; m->score[ 6][14] = -1; m->score[ 6][15] = -1;
  m->score[ 7][ 0] =  1; m->score[ 7][ 1] = -4; m->score[ 7][ 2] =  1; m->score[ 7][ 3] = -4; m->score[ 7][ 4] = -4; m->score[ 7][ 5] = -2; m->score[ 7][ 6] = -2; m->score[ 7][ 7] = -1; m->score[ 7][ 8] = -4; m->score[ 7][ 9] = -2; m->score[ 7][10] = -2; m->score[ 7][11] = -3; m->score[ 7][12] = -1; m->score[ 7][13] = -3; m->score[ 7][14] = -1; m->score[ 7][15] = -1;
  m->score[ 8][ 0] = -4; m->score[ 8][ 1] =  1; m->score[ 8][ 2] = -4; m->score[ 8][ 3] =  1; m->score[ 8][ 4] = -4; m->score[ 8][ 5] = -2; m->score[ 8][ 6] = -2; m->score[ 8][ 7] = -4; m->score[ 8][ 8] = -1; m->score[ 8][ 9] = -2; m->score[ 8][10] = -2; m->score[ 8][11] = -1; m->score[ 8][12] = -3; m->score[ 8][13] = -1; m->score[ 8][14] = -3; m->score[ 8][15] = -1;
  m->score[ 9][ 0] = -4; m->score[ 9][ 1] =  1; m->score[ 9][ 2] =  1; m->score[ 9][ 3] = -4; m->score[ 9][ 4] = -4; m->score[ 9][ 5] = -2; m->score[ 9][ 6] = -2; m->score[ 9][ 7] = -2; m->score[ 9][ 8] = -2; m->score[ 9][ 9] = -1; m->score[ 9][10] = -4; m->score[ 9][11] = -1; m->score[ 9][12] = -3; m->score[ 9][13] = -3; m->score[ 9][14] = -1; m->score[ 9][15] = -1;
  m->score[10][ 0] =  1; m->score[10][ 1] = -4; m->score[10][ 2] = -4; m->score[10][ 3] =  1; m->score[10][ 4] = -4; m->score[10][ 5] = -2; m->score[10][ 6] = -2; m->score[10][ 7] = -2; m->score[10][ 8] = -2; m->score[10][ 9] = -4; m->score[10][10] = -1; m->score[10][11] = -3; m->score[10][12] = -1; m->score[10][13] = -1; m->score[10][14] = -3; m->score[10][15] = -1;
  m->score[11][ 0] = -4; m->score[11][ 1] = -1; m->score[11][ 2] = -1; m->score[11][ 3] = -1; m->score[11][ 4] = -4; m->score[11][ 5] = -1; m->score[11][ 6] = -3; m->score[11][ 7] = -3; m->score[11][ 8] = -1; m->score[11][ 9] = -1; m->score[11][10] = -3; m->score[11][11] = -1; m->score[11][12] = -2; m->score[11][13] = -2; m->score[11][14] = -2; m->score[11][15] = -1;
  m->score[12][ 0] = -1; m->score[12][ 1] = -4; m->score[12][ 2] = -1; m->score[12][ 3] = -1; m->score[12][ 4] = -4; m->score[12][ 5] = -1; m->score[12][ 6] = -3; m->score[12][ 7] = -1; m->score[12][ 8] = -3; m->score[12][ 9] = -3; m->score[12][10] = -1; m->score[12][11] = -2; m->score[12][12] = -1; m->score[12][13] = -2; m->score[12][14] = -2; m->score[12][15] = -1;
  m->score[13][ 0] = -1; m->score[13][ 1] = -1; m->score[13][ 2] = -4; m->score[13][ 3] = -1; m->score[13][ 4] = -4; m->score[13][ 5] = -3; m->score[13][ 6] = -1; m->score[13][ 7] = -3; m->score[13][ 8] = -1; m->score[13][ 9] = -3; m->score[13][10] = -1; m->score[13][11] = -2; m->score[13][12] = -2; m->score[13][13] = -1; m->score[13][14] = -2; m->score[13][15] = -1;
  m->score[14][ 0] = -1; m->score[14][ 1] = -1; m->score[14][ 2] = -1; m->score[14][ 3] = -4; m->score[14][ 4] = -4; m->score[14][ 5] = -3; m->score[14][ 6] = -1; m->score[14][ 7] = -1; m->score[14][ 8] = -3; m->score[14][ 9] = -1; m->score[14][10] = -3; m->score[14][11] = -2; m->score[14][12] = -2; m->score[14][13] = -2; m->score[14][14] = -1; m->score[14][15] = -1;
  m->score[15][ 0] = -2; m->score[15][ 1] = -2; m->score[15][ 2] = -2; m->score[15][ 3] = -2; m->score[15][ 4] = -4; m->score[15][ 5] = -1; m->score[15][ 6] = -1; m->score[15][ 7] = -1; m->score[15][ 8] = -1; m->score[15][ 9] = -1; m->score[15][10] = -1; m->score[15][11] = -1; m->score[15][12] = -1; m->score[15][13] = -1; m->score[15][14] = -1; m->score[15][15] = -1;

  m->best_match[ 0][ 0] =  0; m->best_match[ 0][ 1] = 10; m->best_match[ 0][ 2] =  6; m->best_match[ 0][ 3] =  7; m->best_match[ 0][ 4] = 13; m->best_match[ 0][ 5] = 14; m->best_match[ 0][ 6] = 12; m->best_match[ 0][ 7] = 15; m->best_match[ 0][ 8] =  9; m->best_match[ 0][ 9] =  5; m->best_match[ 0][10] = 11; m->best_match[ 0][11] =  8; m->best_match[ 0][12] =  1; m->best_match[ 0][13] =  4; m->best_match[ 0][14] = 3; m->best_match[ 0][15] = 2;
  m->best_match[ 1][ 0] =  4; m->best_match[ 1][ 1] =  1; m->best_match[ 1][ 2] =  9; m->best_match[ 1][ 3] =  6; m->best_match[ 1][ 4] =  8; m->best_match[ 1][ 5] = 11; m->best_match[ 1][ 6] = 13; m->best_match[ 1][ 7] = 14; m->best_match[ 1][ 8] = 15; m->best_match[ 1][ 9] = 12; m->best_match[ 1][10] = 10; m->best_match[ 1][11] =  5; m->best_match[ 1][12] =  7; m->best_match[ 1][13] =  3; m->best_match[ 1][14] = 2; m->best_match[ 1][15] = 0;
  m->best_match[ 2][ 0] =  2; m->best_match[ 2][ 1] =  5; m->best_match[ 2][ 2] =  7; m->best_match[ 2][ 3] =  9; m->best_match[ 2][ 4] = 12; m->best_match[ 2][ 5] = 14; m->best_match[ 2][ 6] = 11; m->best_match[ 2][ 7] = 15; m->best_match[ 2][ 8] =  8; m->best_match[ 2][ 9] =  6; m->best_match[ 2][10] = 10; m->best_match[ 2][11] = 13; m->best_match[ 2][12] =  4; m->best_match[ 2][13] =  1; m->best_match[ 2][14] = 3; m->best_match[ 2][15] = 0;
  m->best_match[ 3][ 0] =  3; m->best_match[ 3][ 1] =  8; m->best_match[ 3][ 2] =  5; m->best_match[ 3][ 3] = 10; m->best_match[ 3][ 4] = 13; m->best_match[ 3][ 5] = 12; m->best_match[ 3][ 6] = 11; m->best_match[ 3][ 7] = 15; m->best_match[ 3][ 8] =  6; m->best_match[ 3][ 9] =  9; m->best_match[ 3][10] = 14; m->best_match[ 3][11] =  7; m->best_match[ 3][12] =  0; m->best_match[ 3][13] =  4; m->best_match[ 3][14] = 2; m->best_match[ 3][15] = 1;
  m->best_match[ 4][ 0] =  4; m->best_match[ 4][ 1] =  1; m->best_match[ 4][ 2] = 15; m->best_match[ 4][ 3] = 12; m->best_match[ 4][ 4] = 13; m->best_match[ 4][ 5] = 14; m->best_match[ 4][ 6] =  5; m->best_match[ 4][ 7] =  6; m->best_match[ 4][ 8] =  7; m->best_match[ 4][ 9] =  8; m->best_match[ 4][10] =  9; m->best_match[ 4][11] = 10; m->best_match[ 4][12] = 11; m->best_match[ 4][13] =  3; m->best_match[ 4][14] = 2; m->best_match[ 4][15] = 0;
  m->best_match[ 5][ 0] =  2; m->best_match[ 5][ 1] =  3; m->best_match[ 5][ 2] = 11; m->best_match[ 5][ 3] = 15; m->best_match[ 5][ 4] = 12; m->best_match[ 5][ 5] =  5; m->best_match[ 5][ 6] =  8; m->best_match[ 5][ 7] =  7; m->best_match[ 5][ 8] =  9; m->best_match[ 5][ 9] = 10; m->best_match[ 5][10] = 13; m->best_match[ 5][11] = 14; m->best_match[ 5][12] =  1; m->best_match[ 5][13] =  6; m->best_match[ 5][14] = 4; m->best_match[ 5][15] = 0;
  m->best_match[ 6][ 0] =  0; m->best_match[ 6][ 1] =  1; m->best_match[ 6][ 2] = 15; m->best_match[ 6][ 3] =  6; m->best_match[ 6][ 4] = 13; m->best_match[ 6][ 5] = 14; m->best_match[ 6][ 6] =  7; m->best_match[ 6][ 7] =  8; m->best_match[ 6][ 8] =  9; m->best_match[ 6][ 9] = 10; m->best_match[ 6][10] = 12; m->best_match[ 6][11] = 11; m->best_match[ 6][12] =  5; m->best_match[ 6][13] =  3; m->best_match[ 6][14] = 4; m->best_match[ 6][15] = 2;
  m->best_match[ 7][ 0] =  2; m->best_match[ 7][ 1] =  0; m->best_match[ 7][ 2] = 15; m->best_match[ 7][ 3] =  7; m->best_match[ 7][ 4] = 12; m->best_match[ 7][ 5] = 14; m->best_match[ 7][ 6] =  6; m->best_match[ 7][ 7] =  9; m->best_match[ 7][ 8] = 10; m->best_match[ 7][ 9] =  5; m->best_match[ 7][10] = 11; m->best_match[ 7][11] = 13; m->best_match[ 7][12] =  4; m->best_match[ 7][13] =  1; m->best_match[ 7][14] = 8; m->best_match[ 7][15] = 3;
  m->best_match[ 8][ 0] =  3; m->best_match[ 8][ 1] =  1; m->best_match[ 8][ 2] = 15; m->best_match[ 8][ 3] =  8; m->best_match[ 8][ 4] = 11; m->best_match[ 8][ 5] = 13; m->best_match[ 8][ 6] =  6; m->best_match[ 8][ 7] =  9; m->best_match[ 8][ 8] = 10; m->best_match[ 8][ 9] =  5; m->best_match[ 8][10] = 14; m->best_match[ 8][11] = 12; m->best_match[ 8][12] =  4; m->best_match[ 8][13] =  7; m->best_match[ 8][14] = 0; m->best_match[ 8][15] = 2;
  m->best_match[ 9][ 0] =  2; m->best_match[ 9][ 1] =  1; m->best_match[ 9][ 2] = 15; m->best_match[ 9][ 3] = 11; m->best_match[ 9][ 4] =  9; m->best_match[ 9][ 5] = 14; m->best_match[ 9][ 6] =  8; m->best_match[ 9][ 7] =  5; m->best_match[ 9][ 8] =  6; m->best_match[ 9][ 9] =  7; m->best_match[ 9][10] = 12; m->best_match[ 9][11] = 13; m->best_match[ 9][12] = 10; m->best_match[ 9][13] =  0; m->best_match[ 9][14] = 4; m->best_match[ 9][15] = 3;
  m->best_match[10][ 0] =  3; m->best_match[10][ 1] =  0; m->best_match[10][ 2] = 15; m->best_match[10][ 3] = 10; m->best_match[10][ 4] = 12; m->best_match[10][ 5] = 13; m->best_match[10][ 6] =  6; m->best_match[10][ 7] =  7; m->best_match[10][ 8] =  8; m->best_match[10][ 9] =  5; m->best_match[10][10] = 11; m->best_match[10][11] = 14; m->best_match[10][12] =  1; m->best_match[10][13] =  2; m->best_match[10][14] = 4; m->best_match[10][15] = 9;
  m->best_match[11][ 0] = 15; m->best_match[11][ 1] =  1; m->best_match[11][ 2] =  2; m->best_match[11][ 3] =  3; m->best_match[11][ 4] =  5; m->best_match[11][ 5] =  8; m->best_match[11][ 6] =  9; m->best_match[11][ 7] = 11; m->best_match[11][ 8] = 13; m->best_match[11][ 9] = 12; m->best_match[11][10] = 14; m->best_match[11][11] =  7; m->best_match[11][12] = 10; m->best_match[11][13] =  6; m->best_match[11][14] = 0; m->best_match[11][15] = 4;
  m->best_match[12][ 0] =  0; m->best_match[12][ 1] =  2; m->best_match[12][ 2] =  3; m->best_match[12][ 3] =  5; m->best_match[12][ 4] =  7; m->best_match[12][ 5] = 10; m->best_match[12][ 6] = 12; m->best_match[12][ 7] = 15; m->best_match[12][ 8] = 13; m->best_match[12][ 9] = 11; m->best_match[12][10] = 14; m->best_match[12][11] =  8; m->best_match[12][12] =  9; m->best_match[12][13] =  6; m->best_match[12][14] = 1; m->best_match[12][15] = 4;
  m->best_match[13][ 0] = 15; m->best_match[13][ 1] =  1; m->best_match[13][ 2] =  3; m->best_match[13][ 3] =  6; m->best_match[13][ 4] =  8; m->best_match[13][ 5] = 10; m->best_match[13][ 6] = 13; m->best_match[13][ 7] =  0; m->best_match[13][ 8] = 14; m->best_match[13][ 9] = 11; m->best_match[13][10] = 12; m->best_match[13][11] =  7; m->best_match[13][12] =  9; m->best_match[13][13] =  5; m->best_match[13][14] = 4; m->best_match[13][15] = 2;
  m->best_match[14][ 0] =  0; m->best_match[14][ 1] =  1; m->best_match[14][ 2] =  2; m->best_match[14][ 3] =  6; m->best_match[14][ 4] =  7; m->best_match[14][ 5] =  9; m->best_match[14][ 6] = 14; m->best_match[14][ 7] = 15; m->best_match[14][ 8] = 12; m->best_match[14][ 9] = 11; m->best_match[14][10] = 13; m->best_match[14][11] =  8; m->best_match[14][12] = 10; m->best_match[14][13] =  5; m->best_match[14][14] = 4; m->best_match[14][15] = 3;
  m->best_match[15][ 0] = 15; m->best_match[15][ 1] = 10; m->best_match[15][ 2] = 11; m->best_match[15][ 3] = 12; m->best_match[15][ 4] = 13; m->best_match[15][ 5] = 14; m->best_match[15][ 6] =  5; m->best_match[15][ 7] =  6; m->best_match[15][ 8] =  7; m->best_match[15][ 9] =  8; m->best_match[15][10] =  9; m->best_match[15][11] =  1; m->best_match[15][12] =  2; m->best_match[15][13] =  3; m->best_match[15][14] = 0; m->best_match[15][15] = 4;
  m->gap_penalty_set[0][0]=m->gap_penalty_set[1][0]=16; /*SAVE PENALTIES*/
  m->gap_penalty_set[0][1]=m->gap_penalty_set[1][1]=6;
  m->gap_penalty_set[0][2]=m->gap_penalty_set[1][2]=1;
  m->max_gap_length = 30;
  m->trunc_gap_length = 20;
  m->decay_gap_length = 10;
  CALLOC (m->gap_penalty_x, m->max_gap_length+2, LPOScore_T);
  CALLOC (m->gap_penalty_y, m->max_gap_length+2, LPOScore_T);
    /*** GAP OPENING PENALTY @ L=0->1 */
  m->gap_penalty_x[0] = m->gap_penalty_set[0][0];
  m->gap_penalty_y[0] = m->gap_penalty_set[1][0];
  
  /*** 1st AFFINE EXTENSION PENALTY (A1) @ L=1->2,2->3,...T-1->T */
  for (i=1;i<m->trunc_gap_length;i++) {
    m->gap_penalty_x[i] = m->gap_penalty_set[0][1];
    m->gap_penalty_y[i] = m->gap_penalty_set[1][1];
  }
  
  /*** DECAYING EXTENSION PENALTY (A1-->A2; skipped if D=0) @ L=T->T+1,...T+D-1->T+D */
  for (i=0;i<m->decay_gap_length;i++) {
    double dec_x = (m->gap_penalty_set[0][1] - m->gap_penalty_set[0][2]) / ((double)(m->decay_gap_length + 1));
    double dec_y = (m->gap_penalty_set[1][1] - m->gap_penalty_set[1][2]) / ((double)(m->decay_gap_length + 1));
    m->gap_penalty_x[i+m->trunc_gap_length] = m->gap_penalty_set[0][1] - (i+1) * dec_x;
    m->gap_penalty_y[i+m->trunc_gap_length] = m->gap_penalty_set[1][1] - (i+1) * dec_y;
  }
  
  /*** 2nd AFFINE EXTENSION PENALTY (A2) @ L>=T+D */
  m->gap_penalty_x[m->max_gap_length] = m->gap_penalty_set[0][2];
  m->gap_penalty_y[m->max_gap_length] = m->gap_penalty_set[1][2];
  
  m->gap_penalty_x[m->max_gap_length+1] = 0;  /* DON'T REMOVE THIS!... SPECIAL STATE USED IN align_lpo. */
  m->gap_penalty_y[m->max_gap_length+1] = 0;  /* DON'T REMOVE THIS!... SPECIAL STATE USED IN align_lpo. */

}

/** prints a scoring matrix, only including those symbols in subset[] */
void print_score_matrix(FILE *ifile,ResidueScoreMatrix_T *m,char subset[])
{
  int i,i_m,j,j_m,nsubset;

  nsubset=strlen(subset);

  printf(" ");
  LOOPF (i,nsubset)
    printf("  %c",subset[i]);
  printf("\n");

  LOOPF (i,nsubset) {
    LOOP (i_m,m->nsymbol)
      if (m->symbol[i_m]==subset[i])
	break;
    printf("%c",subset[i]);
    LOOPF (j,nsubset) {
      LOOP (j_m,m->nsymbol)
	if (m->symbol[j_m]==subset[j])
	  break;
      printf("%3d",m->score[i_m][j_m]);
    }
    printf("\n");
  }
  return;
}



/** restricts seq[] to the set of allowed characters given by symbol[];
 other characters will be replaced by the default symbol[0] */
int limit_residues(char seq[],char symbol[])
{
  int i,len,nreplace=0;

  len=strlen(seq);
  for (i=strspn(seq,symbol);i<len;i=i+1+strspn(seq+i+1,symbol)) {
    seq[i]=symbol[0]; /* FORCE IT TO BE A VALID SYMBOL */
    nreplace++; /* COUNT TOTAL REPLACEMENTS */
  }
  return nreplace;
}





/** RETURNS THE COMPLEMENTARY BASE *IN LOWER CASE* */
char complementary_base(char base) 
{
  switch (base) {
  case 'A': case 'a':                     return 't';
  case 'T': case 't': case 'U': case 'u': return 'a';
  case 'G': case 'g':                     return 'c';
  case 'C': case 'c':                     return 'g';
  default: return base;
  }
}

/** REVERSE COMPLEMENTS seq[] IN PLACE, AND RETURNS POINTER TO seq  
---------------------------------------------------------------
-----------------------------------------------------------
*/
char *reverse_complement(char seq[])
{
  int i,j;
  char c;
  for ((i=0),(j=strlen(seq)-1);i<=j;i++,j--) {/* SWAP FROM ENDS TO CENTER*/
    c=complementary_base(seq[i]); /* SWAP ENDS AND COMPLEMENT */
    seq[i]=complementary_base(seq[j]);
    seq[j]=c;
  }
  return seq;
}

