
FILES=$(find tests/ -name "error.txt")
for file in $FILES
do
  echo $file
  echo -n "      "
  cat $file
  echo ""
done