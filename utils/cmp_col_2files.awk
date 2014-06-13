BEGIN {
  FS = OFS = "\t"
  }
NR == FNR {
  # while reading the 1st file
  # store its records in the array f
  f[$1] = $0
  next
  }
$1 in f {
  # when match is found
  # print all values
  print f[$1], $0
  }
