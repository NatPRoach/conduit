
type
  PFile* {.importc: "FILE*", header: "<stdio.h>".} = distinct pointer
    # import C's FILE* type; Nim will treat it as a new pointer type

const
  REV_COMP_STRING* : cstring = "/rev_comp"

## * THE NULL LETTER-REFERENCE

const
  INVALID_LETTER_POSITION* = (-1)

type
  LPOLetterRefT* = cint
  LPOScoreT* = cint

## * NEEDED FOR seq_util.h

# type
#   ResidueScoreT* = LPOScoreT

const
  RESIDUE_SCORE_DEFINED* = true

## * linked list for storing source origin (sequence position)
## * from which this letter was derived

type
  LPOLetterSourceS* {.bycopy,extern:"LPOLetterSource_S".} = object
    ## * index of the sequence, referencing the sourceSeq[] array
    iseq*: cint
    ## * index of the corresponding position in that sequence
    ipos*: LPOLetterRefT
    ## * next node in the linked list
    more*: ptr LPOLetterSourceS

  LPOLetterSourceT* = LPOLetterSourceS

## * linked list for connecting an LPOLetter to either right
## * or left

when not defined(USE_WEIGHTED_LINKS):
  type
    LPOLetterLinkS* {.bycopy,extern:"LPOLetterLink_S".} = object
      ## * ADJACENT LETTER LINKED TO THIS LETTER
      ipos*: LPOLetterRefT    
      ## * next node in the linked list
      more*: ptr LPOLetterLinkS

when defined(USE_WEIGHTED_LINKS):
  type
    LPOLetterLinkS* {.bycopy,extern:"LPOLetterLink_S".} = object
      ## * ADJACENT LETTER LINKED TO THIS LETTER
      ipos*: LPOLetterRef_T
      ## * transition cost for traversing this link
      score*: LPOScore_T
      ## * next node in the linked list
      more*: ptr LPOLetterLink_S

type
  LPOLetterLinkT* = LPOLetterLinkS

## * the chunk size for allocating additional
## *    letters in an LPOLetter_T array

const
  LPO_LETTER_BUFFER_CHUNK* = 64

## * Structure for storing individual LPO Letters

type
  LPOLetterS* {.bycopy,extern:"LPOLetter_S".} = object
    left*: LPOLetterLinkT       ## * ADJACENT LETTER(S) TO THE LEFT
    right*: LPOLetterLinkT      ## * ADJACENT LETTER(S) TO THE RIGHT
    source*: LPOLetterSourceT   ## * SOURCE SEQ POSITION(S)
    alignRing*: LPOLetterRefT   ## * CIRCULAR LIST OF ALIGNED POSITIONS
    ringId*: LPOLetterRefT      ## * MINIMUM INDEX OF ALL POSITIONS ON THE RING
    score*: cfloat              ## * SCORE FOR BALANCING PARTIAL ORDER EFFECTS 
                                ## * ON MATRIX NEUTRALITY
    letter*: char               ## * THE ACTUAL RESIDUE CODE!

  LPOLetterT* = LPOLetterS

## * maximum length of a sequence name

const
  SEQUENCE_NAME_MAX* = 128

## * buffer chunk size for expanding a block of seq storage

const
  SEQUENCE_BUFFER_CHUNK* = 8

## * buffer chunk size for expanding a sourceSeq[] array

const
  SOURCE_SEQ_BUFFER_CHUNK* = 16
  NUMDATA_BUFFER_CHUNK* = 4

## * storage for quantitative data attached to a sequence

type
  LPONumericDataS* {.bycopy,extern:"LPONumericData_S".} = object
    name*: array[SEQUENCE_NAME_MAX, char] ## *
    ## *
    title*: cstring            ## *
    data*: ptr cdouble

  LPONumericDataT* {.bycopy,extern:"LPONumericData_T".} = LPONumericDataS

## * Structure for storing individual source sequence information,
##  stuff like name, title etc.

type
  LPOSourceInfoS* {.bycopy,extern:"LPOSourceInfo_S".} = object
    name*: array[SEQUENCE_NAME_MAX, char] ## *
    ## *
    title*: cstring            ## *
    sequence*: cstring         ## *
    seqToPo*: ptr cint        ## *
    poToSeq*: ptr cint        ## *
    data*: ptr LPONumericDataT ## *
    ndata*: cint               ## *
    length*: cint              ## *
    istart*: cint              ## * FOR PURPOSES OF HEAVIEST BUNDLE CALCULATION
    weight*: cint              ## * WHAT BUNDLE IS THIS A MEMBER OF?
    bundleId*: cint

  LPOSourceInfoT* {.bycopy,extern:"LPOSourceInfo_T".} = LPOSourceInfoS

## * the NULL bundle-reference

const
  NO_BUNDLE* = (-1)

## * bundle-reference meaning "include all bundles"

const
  ALL_BUNDLES* = (-1)

## * holder for an LPO sequence, its letters,
##   and associated information

type
  LPOSequenceS* {.bycopy,extern:"LPOSequence_S".} = object
    length*: cint              ## *
    ## *
    letter*: ptr LPOLetterT    ## *
    title*: cstring            ## *
    sequence*: cstring         ## *
    name*: array[SEQUENCE_NAME_MAX, char] ## *
    nsourceSeq*: cint         ## *
    sourceSeq*: ptr LPOSourceInfoT

  LPOSequenceT*{.bycopy,extern:"LPOSequence_T".} = LPOSequenceS
  SequenceT*{.bycopy,extern:"Sequence_T".} = LPOSequenceT

## *@memo GENERAL FORM IS seq_y[j].left.ipos

template SEQ_Y_LEFT*(j: untyped): untyped =
  (j - 1)

template SEQ_Y_RIGHT*(j: untyped): untyped =
  (j + 1)

## *@memo Data structure for analyzing sequence differences in MSA

type
  LPOLetterCountS* {.bycopy,extern:"LPOLetterCount_S".} = object
    isError* {.bitsize: 2.}: cuint
    meetsCriteria* {.bitsize: 1.}: cuint
    seqCount* {.bitsize: 29.}: cuint

  LPOLetterCountT* {.bycopy,extern:"LPOLetterCount_T".} = LPOLetterCountS

## * classification of sequence differences

# const
#   no_error* = 0
#   substitution_error* = 1
#   insertion_error* = 2
#   deletion_error* = 3
#   max_error_states* = 4

## * DON'T ALLOCATE MORE THAN THIS TOTAL AMOUNT OF MEMORY
## ---------------------------------------------------------------
## ---------------------------------------------------------------
##
const
  POA_MAX_ALLOC* = 300000000


type
  ResidueScoreT* = cint
    ##  DEFAULT: USE int FOR SCORING
# const
#   dont_switch_case* = 0
#   switch_case_to_lower* = 1
#   switch_case_to_upper* = 2

const
  MATRIX_SYMBOL_MAX* = 16

type
  ResidueScoreMatrixT* {.bycopy,extern:"ResidueScoreMatrix_T".} = object
    nsymbol*: cint
    # Need the extra char to add null terminator to c
    symbol*: array[MATRIX_SYMBOL_MAX + 1, char]
    score*: array[MATRIX_SYMBOL_MAX, array[MATRIX_SYMBOL_MAX, ResidueScoreT]]
    bestMatch*: array[MATRIX_SYMBOL_MAX, array[MATRIX_SYMBOL_MAX, cint]]
    gapPenaltySet*: array[2, array[3, ResidueScoreT]]
    truncGapLength*: cint
    decayGapLength*: cint
    gapPenaltyX*: ptr array[32,ResidueScoreT]
    gapPenaltyY*: ptr array[32,ResidueScoreT]
    maxGapLength*: cint
    nfreq*: cint    ##  STORE FREQUENCIES OF AMINO ACIDS FOR BALANCING MATRIX...
    freqSymbol*: array[MATRIX_SYMBOL_MAX, char]
    freq*: array[MATRIX_SYMBOL_MAX, cfloat]

# type
#   FILE* {.bycopy.} = object

type
  LPOReturnResultS* {.bycopy,extern:"LPOReturnResult_S".} = object
    sequence : ptr LPOSequenceT
    lpoSeqs* : ptr LPOSequenceT
    inputSeqs : ptr ptr LPOSequenceT
    matrix* : ptr ResidueScoreMatrixT
    nInputSeqs : cint
    nseq : cint

  LPOReturnResultT*{.bycopy,extern:"LPOReturnResult_T".} = LPOReturnResultS