BEGIN {
  FS = OFS = "\t"
  }
NR == FNR {
  # while reading the 1st file
  # store its records in the array f
  f[$1] = $0
  next
  }


{
#print "Dolar 1 = " $1;

if ($1 in f) 
{
  # when match is found
  # print all values
  print $0, f[$1]
}
else
{
  print $0, $1, 0
}

}

