import hts
import strutils
import strformat

type
  FastaRecord* = object
    readId* : string
    sequence* : string


proc writeFASTArecordToFile*( outfile : File,
                              record : FastaRecord,
                              wrap_len : int = 75) = 
  outfile.write(&">{record.readId}\n")
  for i in 0..<(record.sequence.len div wrap_len):
    let line = record.sequence[i *
                               wrap_len..(i + 1) * wrap_len - 1]
    outfile.write(&"{line}\n")
  if record.sequence.len mod wrap_len != 0 :
    let line = record.sequence[wrap_len *
                               (record.sequence.len div wrap_len)..^1]
    outfile.write(&"{line}\n")

# proc writeCorrectedReads*( outfile : File,
#                            records : seq[FastaRecord],
#                            wrap_len : int = 75) = 
proc writeFASTArecordsToFile*( outfile : File,
                               records : seq[FastaRecord],
                               wrap_len : int = 75) = 
  for record in records:
    writeFASTArecordToFile( outfile, record, wrap_len)


proc writeBamRecordToFASTAfile*( outfile : File,
                                 record : Record,
                                 wrap_len : int = 75) = 
  var sequence : string
  record.sequence(sequence)
  writeFASTArecordToFile( outfile, FastaRecord(readId : record.qname,
                                               sequence : sequence),
                                               wrap_len)


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


proc parseFasta*(file : File) : seq[FastaRecord] = 
  var records : seq[FastaRecord]
  var readId : string
  var sequence : string
  var start = true
  while true:
    try:
      let line = file.readLine()
      if line[0] == '>':
        if not start:
          records.add(FastaRecord( readId : readId, sequence : sequence))
          sequence = ""
        else:
          start = false
        readId = line.strip(leading=true,
                            trailing=false,
                            chars = {'>'}).strip(leading=false,
                                                 trailing=true,
                                                 chars = {'\n'})
      else:
        sequence = sequence &
          line.strip(leading=false,
                     trailing=true,
                     chars = {'\n'})
    except EOFError:
      records.add(FastaRecord( readId : readId, sequence : sequence))
      break
  return records


iterator iterFasta*(file : File) : FastaRecord =
  var readId : string
  var sequence : string
  var start = true
  while true:
    try:
      let line = file.readLine()
      if line[0] == '>':
        if not start:
          yield FastaRecord( readId : readId, sequence : sequence)
          sequence = ""
        else:
          start = false
        readId = line.strip(leading=true,
                            trailing=false,
                            chars = {'>'}).strip(leading=false,
                                                 trailing=true,
                                                 chars = {'\n'})
      else:
        sequence = sequence &
          line.strip(leading=false,
                     trailing=true,
                     chars = {'\n'})
    except EOFError:
      yield FastaRecord( readId : readId, sequence : sequence)
      break


proc splitFASTA*(infilepath,outfilepath_prefix : string,
                 split_num : int = 200) : (int,int) =
  # Takes in a fasta file with more reads than some arbitrary number split_num
  # produces fasta files split into sizes of that size or smaller.
  # Similar to RATTLE implementation in that to avoid size biases
  # the sampling is performed with an offset.
  var infile : File
  discard open(infile,infilepath,fmRead)
  var records : seq[FastaRecord]
  var readId : string
  var sequence : string
  var count = 0
  while true:
    try:
      let line = infile.readLine()
      if line[0] == '>':
        if count != 0:
          records.add(FastaRecord( readId : readId, sequence : sequence))
          sequence = ""
        count += 1
        readId = line.strip(leading=true,
                            trailing=false,
                            chars = {'>'}).strip(leading=false,
                                                 trailing=true,
                                                 chars = {'\n'})
      else:
        sequence = sequence & line.strip(leading=false,
                                         trailing=true,
                                         chars = {'\n'})
    except EOFError:
      records.add(FastaRecord( readId : readId, sequence : sequence))
      break
  infile.close()
  let numOutfiles = int(records.len mod split_num != 0) +
                     (records.len div split_num)
  if numOutfiles > 1:
    for i in 0..<numOutfiles:
      var outfile : File
      discard open(outfile,&"{outfilepath_prefix}_subfasta{i}.fa",fmWrite)
      for j in 0..<split_num:
        let idx = i + (j * numOutfiles)
        if idx < records.len:
          outfile.write(">",records[idx].readId,"\n")
          outfile.writeLine(records[idx].sequence)
      outfile.close()
  return (numOutfiles,records.len)


proc splitFASTA2*(infilepath,outfilepath_prefix : string,
                  split_num : int = 200) : (int,int) =
  # Takes in a fasta file with more reads than some arbitrary number split_num
  # produces fasta files split into sizes of that size or smaller.
  # Different from RATTLE implementation in that size bias is not
  # considered and linear blocks of sequences are extracted
  var infile : File
  discard open(infile,infilepath,fmRead)
  var records : seq[FastaRecord]
  var readId : string
  var sequence : string
  var count = 0
  while true:
    try:
      let line = infile.readLine()
      if line[0] == '>':
        if count != 0:
          records.add(FastaRecord( readId : readId, sequence : sequence))
          sequence = ""
        count += 1
        readId = line.strip(leading=true,
                            trailing=false,
                            chars = {'>'}).strip(leading=false,
                                                 trailing=true,
                                                 chars = {'\n'})
      else:
        sequence = sequence & line.strip(leading=false,
                                         trailing=true,
                                         chars = {'\n'})
    except EOFError:
      records.add(FastaRecord( readId : readId, sequence : sequence))
      break
  infile.close()
  let numOutfiles = int(records.len mod split_num != 0) +
                    (records.len div split_num)
  if numOutfiles > 1:
    for i in 0..<numOutfiles:
      var outfile : File
      discard open(outfile,&"{outfilepath_prefix}_subfasta{i}.fa",fmWrite)
      for j in 0..<split_num:
        let idx = (i*split_num) + j
        if idx < records.len:
          outfile.write(">",records[idx].readId,"\n")
          outfile.writeLine(records[idx].sequence)
      outfile.close()
  return (numOutfiles,records.len)