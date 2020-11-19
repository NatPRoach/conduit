import hts
import strutils
import strformat

type
  FastaRecord* = object
    read_id* : string
    sequence* : string


proc writeFASTArecordToFile*( outfile : File, record : FastaRecord, wrap_len : int = 75) = 
  outfile.write(&">{record.read_id}\n")
  for i in 0..<(record.sequence.len div wrap_len):
    outfile.write(&"{record.sequence[i*wrap_len..(i+1)*wrap_len - 1]}\n")
  if record.sequence.len mod wrap_len != 0 :
    outfile.write(&"{record.sequence[wrap_len*(record.sequence.len div wrap_len)..^1]}\n")

# proc writeCorrectedReads*( outfile : File, records : seq[FastaRecord], wrap_len : int = 75) = 
proc writeFASTArecordsToFile*( outfile : File, records : seq[FastaRecord], wrap_len : int = 75) = 
  for record in records:
    writeFASTArecordToFile( outfile, record, wrap_len)


proc writeBamRecordToFASTAfile*(outfile : File, record : Record, wrap_len : int = 75) = 
  var sequence : string
  record.sequence(sequence)
  writeFASTArecordToFile( outfile, FastaRecord(read_id : record.qname, sequence : sequence), wrap_len)


proc convertUtoTinFASTA*(infilepath,outfilepath:string) =
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


#TODO - Need to do this in a way where it's an iterator instead. This is loading the entire file in to memory, which is fine for small cluster files but will use huge memory at scale
proc parseFasta*(file : File) : seq[FastaRecord] = 
  var records : seq[FastaRecord]
  var read_id : string
  var sequence : string
  var count = 0
  while true:
    try:
      let line = file.readLine()
      if line[0] == '>':
        if count != 0:
          records.add(FastaRecord( read_id : read_id, sequence : sequence))
          sequence = ""
        count += 1
        read_id = line.strip(leading=true,trailing=false,chars = {'>'}).strip(leading=false,trailing=true,chars = {'\n'})
      else:
        sequence = sequence & line.strip(leading=false,trailing=true,chars = {'\n'})
    except EOFError:
      records.add(FastaRecord( read_id : read_id, sequence : sequence))
      break
  return records

proc splitFASTA*(infilepath,outfilepath_prefix : string, split_num : int = 200) : (int,int) =
  # Takes in a fasta file with more reads than some arbitrary number split_num
  # produces fasta files split into sizes of that size or smaller. Similar to RATTLE implementation in that to avoid size biases the sampling is performed with an offset.
  var infile : File
  discard open(infile,infilepath,fmRead)
  var records : seq[FastaRecord]
  var read_id : string
  var sequence : string
  var count = 0
  while true:
    try:
      let line = infile.readLine()
      if line[0] == '>':
        if count != 0:
          records.add(FastaRecord( read_id : read_id, sequence : sequence))
          sequence = ""
        count += 1
        read_id = line.strip(leading=true,trailing=false,chars = {'>'}).strip(leading=false,trailing=true,chars = {'\n'})
      else:
        sequence = sequence & line.strip(leading=false,trailing=true,chars = {'\n'})
    except EOFError:
      records.add(FastaRecord( read_id : read_id, sequence : sequence))
      break
  infile.close()
  let num_outfiles = int(records.len mod split_num != 0) + (records.len div split_num)
  if num_outfiles > 1:
    for i in 0..<num_outfiles:
      var outfile : File
      discard open(outfile,&"{outfilepath_prefix}_subfasta{i}.fa",fmWrite)
      for j in 0..<split_num:
        let idx = i + (j * num_outfiles)
        if idx < records.len:
          outfile.write(">",records[idx].read_id,"\n")
          outfile.writeLine(records[idx].sequence)
      outfile.close()
  return (num_outfiles,records.len)


proc splitFASTA2*(infilepath,outfilepath_prefix : string, split_num : int = 200) : (int,int) =
  # Takes in a fasta file with more reads than some arbitrary number split_num
  # produces fasta files split into sizes of that size or smaller. Different from RATTLE implementation in that size bias is not considered and linear blocks of sequences are extracted
  var infile : File
  discard open(infile,infilepath,fmRead)
  var records : seq[FastaRecord]
  var read_id : string
  var sequence : string
  var count = 0
  while true:
    try:
      let line = infile.readLine()
      if line[0] == '>':
        if count != 0:
          records.add(FastaRecord( read_id : read_id, sequence : sequence))
          sequence = ""
        count += 1
        read_id = line.strip(leading=true,trailing=false,chars = {'>'}).strip(leading=false,trailing=true,chars = {'\n'})
      else:
        sequence = sequence & line.strip(leading=false,trailing=true,chars = {'\n'})
    except EOFError:
      records.add(FastaRecord( read_id : read_id, sequence : sequence))
      break
  infile.close()
  let num_outfiles = int(records.len mod split_num != 0) + (records.len div split_num)
  if num_outfiles > 1:
    for i in 0..<num_outfiles:
      var outfile : File
      discard open(outfile,&"{outfilepath_prefix}_subfasta{i}.fa",fmWrite)
      for j in 0..<split_num:
        let idx = (i*split_num) + j
        if idx < records.len:
          outfile.write(">",records[idx].read_id,"\n")
          outfile.writeLine(records[idx].sequence)
      outfile.close()
  return (num_outfiles,records.len)