import sets
import fasta
import hts
import strutils
import strformat

type
  FastqRecord* = object
    readId* : string
    sequence* : string
    qualities* : string


proc writeFASTQrecordToFile*( outfile : File,
                              record : FastqRecord) = 
  outfile.write(&"@{record.readId}\n")
  outfile.write(&"{record.sequence}\n")
  outfile.write("+\n")
  outfile.write(&"{record.qualities}\n")

proc writeFASTQrecordsToFile*( outfile : File,
                               records : seq[FastqRecord]) = 
  for record in records:
    writeFASTQrecordToFile( outfile, record)


proc convertQualityIntsToString( intQuals: seq[uint8]) : string =
  var charQuals : seq[char]
  for qual in intQuals:
    charQuals.add(char(qual + 33'u8))
  result = charQuals.join("")


proc writeBamRecordToFASTQfile*( outfile : File,
                                 record : Record) = 
  var sequence : string
  discard record.sequence(sequence)
  var qualities : seq[uint8]
  discard record.base_qualities(qualities)
  let qualStr = convertQualityIntsToString(qualities)
  # echo record.qname
  writeFASTQrecordToFile( outfile, FastqRecord(readId : record.qname,
                                               sequence : sequence,
                                               qualities : qualStr))


proc convertFASTQfileToFASTAfile*(infilepath,outfilepath : string ) =
  var infile,outfile : File
  discard open(infile,infilepath,fmRead)
  discard open(outfile,outfilepath,fmWrite)
  var readIds : HashSet[string]
  while true:
    try:
      let line1 = infile.readLine()
      try:
        assert line1[0] == '@'
      except AssertionDefect:
        echo "File not in FASTQ format"
        raise
      if line1 notin readIds:
        readIds.incl(line1)
        outfile.write(&">{line1[1..^1]}\n")
        outfile.write(&"{infile.readLine().replace(sub='U',by='T')}\n")
        discard infile.readLine()
        discard infile.readLine()
      else:
        discard infile.readLine()
        discard infile.readLine()
        discard infile.readLine()
    except EOFError:
      break
  infile.close()
  outfile.close()


proc convertFASTQtoFASTA*(record : FastqRecord) : FastaRecord = 
  result.readId = record.readId
  result.sequence = record.sequence


proc convertFASTQtoFASTA*(records : openArray[FastqRecord]) : seq[FastaRecord] =
  for record in records:
    result.add(convertFASTQtoFASTA(record))


proc parseFASTQ*(infile : File) : seq[FastqRecord] = 
  while true:
    try:
      let line1 = infile.readLine()
      try:
        assert line1[0] == '@'
      except AssertionDefect:
        echo "File not in FASTQ format"
        raise
      let readId = line1[1..^1]
      let sequence = infile.readLine().replace(sub='U',by='T')
      discard infile.readLine()
      let quals = infile.readLine()
      result.add(FastqRecord( readId   : readId,
                              sequence  : sequence,
                              qualities : quals))
    except EOFError:
      break


iterator iterFASTQ*(infile : File) : FastqRecord =
  while true:
    try:
      let line1 = infile.readLine()
      try:
        assert line1[0] == '@'
      except AssertionDefect:
        echo "File not in FASTQ format"
        raise
      let readId = line1[1..^1]
      let sequence = infile.readLine().replace(sub='U',by='T')
      discard infile.readLine()
      let quals = infile.readLine()
      yield FastqRecord( readId   : readId,
                         sequence  : sequence,
                         qualities : quals)
    except EOFError:
      break