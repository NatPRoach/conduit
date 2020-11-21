#small library for dealing with nucleic acid sequences.
import strutils
import tables

proc revComp*(nts : string) : string =
  let wc_pairs = {'A' : 'T',
                  'a' : 'T',
                  'C' : 'G',
                  'c' : 'G',
                  'T' : 'A',
                  't' : 'A',
                  'U' : 'A',
                  'u' : 'A',
                  'G' : 'C',
                  'g' : 'C'}.toTable()
  var revcomp : seq[char]
  for i in 1..nts.len:
    revcomp.add(wc_pairs[nts[^i]])
  result = revcomp.join("")