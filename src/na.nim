#small library for dealing with nucleic acid sequences.
import strutils
import tables

proc revComp*(nts : string) : string =
  let wcPairs = {'A' : 'T',
                  'a' : 'T',
                  'C' : 'G',
                  'c' : 'G',
                  'T' : 'A',
                  't' : 'A',
                  'U' : 'A',
                  'u' : 'A',
                  'G' : 'C',
                  'g' : 'C'}.toTable()
  var revComp : seq[char]
  for i in 1..nts.len:
    revComp.add(wcPairs[nts[^i]])
  result = revComp.join("")