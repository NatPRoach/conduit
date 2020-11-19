import hts
import fastq

type
  FastaRecord* = object
    read_id* : string
    sequence* : string

proc writeFASTArecordToFile*(outfile : File, record : Record, wrap_len : int = 75) = 
  var sequence : string
  record.sequence(sequence)
  outfile.write(&">{record.qname}\n")
  for i in 0..<(sequence.len div wrap_len):
    outfile.write(&"{sequence[i*wrap_len..(i+1)*wrap_len - 1]}\n")
  if sequence.len mod wrap_len != 0 :
    outfile.write(&"{sequence[wrap_len*(sequence.len div wrap_len)..^1]}\n")


proc convertUtoTinFASTA(infilepath,outfilepath:string) =
  var infile,outfile : File
  discard open(infile,infilepath,fmRead)
  discard open(outfile,outfilepath,fmWrite)
  while true:
    try:
      let line = infile.readLine()
      if line[0] == '>':
        outfile.write(&"{line}\n")
      else:
        outfile.write(&"{line.replace(sub='U',by='T')}\n")
    except EOFError:
      break
  infile.close()
  outfile.close()