
type
  PFile* {.importc: "FILE*", header: "<stdio.h>".} = distinct pointer
    # import C's FILE* type; Nim will treat it as a new pointer type

const
  REV_COMP_STRING* : cstring = "/rev_comp"

## * THE NULL LETTER-REFERENCE

const
  INVALID_LETTER_POSITION* = (-1)

type
  LPOLetterRef_T* = cint
  LPOScore_T* = cint

## * NEEDED FOR seq_util.h

# type
#   ResidueScore_T* = LPOScore_T

const
  RESIDUE_SCORE_DEFINED* = true

## * linked list for storing source origin (sequence position)
##  from which this letter was derived

type
  LPOLetterSource_S* {.bycopy.} = object
    ## * index of the sequence, referencing the source_seq[] array
    iseq*: cint
    ## * index of the corresponding position in that sequence
    ipos*: LPOLetterRef_T
    ## * next node in the linked list
    more*: ptr LPOLetterSource_S

  LPOLetterSource_T* = LPOLetterSource_S

## * linked list for connecting an LPOLetter to either right
##  or left

when not defined(USE_WEIGHTED_LINKS):
  type
    LPOLetterLink_S* {.bycopy.} = object
      ## * ADJACENT LETTER LINKED TO THIS LETTER
      ipos*: LPOLetterRef_T    
      ## * next node in the linked list
      more*: ptr LPOLetterLink_S

when defined(USE_WEIGHTED_LINKS):
  type
    LPOLetterLink_S* {.bycopy.} = object
      ## * ADJACENT LETTER LINKED TO THIS LETTER
      ipos*: LPOLetterRef_T
      ## * transition cost for traversing this link
      score*: LPOScore_T
      ## * next node in the linked list
      more*: ptr LPOLetterLink_S

type
  LPOLetterLink_T* = LPOLetterLink_S

## * the chunk size for allocating additional
##     letters in an LPOLetter_T array

const
  LPO_LETTER_BUFFER_CHUNK* = 64

## * Structure for storing individual LPO Letters

type
  LPOLetter_S* {.bycopy.} = object
    left*: LPOLetterLink_T     ## * ADJACENT LETTER(S) TO THE LEFT
    ## * ADJACENT LETTER(S) TO THE RIGHT
    right*: LPOLetterLink_T    ## * SOURCE SEQ POSITION(S)
    source*: LPOLetterSource_T ## * CIRCULAR LIST OF ALIGNED POSITIONS
    align_ring*: LPOLetterRef_T ## * MINIMUM INDEX OF ALL POSITIONS ON THE RING
    ring_id*: LPOLetterRef_T   ## * SCORE FOR BALANCING PARTIAL ORDER EFFECTS ON MATRIX NEUTRALITY
    score*: cfloat             ## * THE ACTUAL RESIDUE CODE!
    letter*: char

  LPOLetter_T* = LPOLetter_S

## * maximum length of a sequence name

const
  SEQUENCE_NAME_MAX* = 128

## * buffer chunk size for expanding a block of seq storage

const
  SEQUENCE_BUFFER_CHUNK* = 8

## * buffer chunk size for expanding a source_seq[] array

const
  SOURCE_SEQ_BUFFER_CHUNK* = 16
  NUMDATA_BUFFER_CHUNK* = 4

## * storage for quantitative data attached to a sequence

type
  LPONumericData_S* {.bycopy.} = object
    name*: array[SEQUENCE_NAME_MAX, char] ## *
    ## *
    title*: cstring            ## *
    data*: ptr cdouble

  LPONumericData_T* = LPONumericData_S

## * Structure for storing individual source sequence information,
##  stuff like name, title etc.

type
  LPOSourceInfo_S* {.bycopy.} = object
    name*: array[SEQUENCE_NAME_MAX, char] ## *
    ## *
    title*: cstring            ## *
    sequence*: cstring         ## *
    seq_to_po*: ptr cint        ## *
    po_to_seq*: ptr cint        ## *
    data*: ptr LPONumericData_T ## *
    ndata*: cint               ## *
    length*: cint              ## *
    istart*: cint              ## * FOR PURPOSES OF HEAVIEST BUNDLE CALCULATION
    weight*: cint              ## * WHAT BUNDLE IS THIS A MEMBER OF?
    bundle_id*: cint

  LPOSourceInfo_T* = LPOSourceInfo_S

## * the NULL bundle-reference

const
  NO_BUNDLE* = (-1)

## * bundle-reference meaning "include all bundles"

const
  ALL_BUNDLES* = (-1)

## * holder for an LPO sequence, its letters,
##   and associated information

type
  LPOSequence_S* {.bycopy.} = object
    length*: cint              ## *
    ## *
    letter*: ptr LPOLetter_T    ## *
    title*: cstring            ## *
    sequence*: cstring         ## *
    name*: array[SEQUENCE_NAME_MAX, char] ## *
    nsource_seq*: cint         ## *
    source_seq*: ptr LPOSourceInfo_T

  LPOSequence_T*{.bycopy.} = LPOSequence_S
  Sequence_T*{.bycopy.} = LPOSequence_T

## *@memo GENERAL FORM IS seq_y[j].left.ipos

template SEQ_Y_LEFT*(j: untyped): untyped =
  (j - 1)

template SEQ_Y_RIGHT*(j: untyped): untyped =
  (j + 1)

## *@memo Data structure for analyzing sequence differences in MSA

type
  LPOLetterCount_S* {.bycopy.} = object
    is_error* {.bitsize: 2.}: cuint
    meets_criteria* {.bitsize: 1.}: cuint
    seq_count* {.bitsize: 29.}: cuint

  LPOLetterCount_T* = LPOLetterCount_S

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
  ResidueScore_T* = cint
    ##  DEFAULT: USE int FOR SCORING
# const
#   dont_switch_case* = 0
#   switch_case_to_lower* = 1
#   switch_case_to_upper* = 2

const
  MATRIX_SYMBOL_MAX* = 16

type
  ResidueScoreMatrix_T* {.bycopy.} = object
    nsymbol*: cint
    symbol*: array[MATRIX_SYMBOL_MAX, char]
    score*: array[MATRIX_SYMBOL_MAX, array[MATRIX_SYMBOL_MAX, ResidueScore_T]]
    best_match*: array[MATRIX_SYMBOL_MAX, array[MATRIX_SYMBOL_MAX, cint]]
    gap_penalty_set*: array[2, array[3, ResidueScore_T]]
    trunc_gap_length*: cint
    decay_gap_length*: cint
    gap_penalty_x*: ptr array[32,ResidueScore_T]
    gap_penalty_y*: ptr array[32,ResidueScore_T]
    max_gap_length*: cint
    nfreq*: cint               ##  STORE FREQUENCIES OF AMINO ACIDS FOR BALANCING MATRIX...
    freq_symbol*: array[MATRIX_SYMBOL_MAX, char]
    freq*: array[MATRIX_SYMBOL_MAX, cfloat]

# type
#   FILE* {.bycopy.} = object

type
  LPOReturnResult_S* {.bycopy.} = object
    sequence : ptr LPOSequence_T
    lpo_seqs* : ptr LPOSequence_T
    input_seqs : ptr ptr LPOSequence_T
    matrix* : ptr ResidueScoreMatrix_T
    n_input_seqs : cint
    nseq : cint

  LPOReturnResult_T*{.bycopy.} = LPOReturnResult_S