import fasta

type
  FastqRecord* = object
    read_id* : string
    sequence* : string
    qualities* : string


proc convertFASTQfileToFASTAfile*(infilepath,outfilepath:string) = 
  var infile,outfile : File
  discard open(infile,infilepath,fmRead)
  discard open(outfile,outfilepath,fmWrite)
  var read_ids : HashSet[string]
  while true:
    try:
      let line1 = infile.readLine()
      try:
        assert line1[0] == '@'
      except AssertionError:
        echo "File not in FASTQ format"
        raise
      if line1 notin read_ids:
        read_ids.incl(line1)
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
  result.read_id = record.read_id
  result.sequence = record.sequence


proc convertFASTQtoFASTA*(records : openArray[FastqRecord]) : seq[FastaRecord] =
  for record in records:
    result.add(convertFASTQtoFASTA(record))


#TODO - Need to do this in a way where it's an iterator instead. This is loading the entire file in to memory, which is fine for small cluster files but will use huge memory at scale
proc parseFASTQ*(infile : File) : seq[FastqRecord] = 
  while true:
    try:
      let line1 = infile.readLine()
      try:
        assert line1[0] == '@'
      except AssertionError:
        echo "File not in FASTQ format"
        raise
      let read_id = line1[1..^1]
      let sequence = infile.readLine().replace(sub='U',by='T')
      discard infile.readLine()
      let quals = infile.readLine()
      result.add(FastqRecord( read_id   : read_id,
                              sequence  : sequence,
                              qualities : quals))
    except EOFError:
      break
