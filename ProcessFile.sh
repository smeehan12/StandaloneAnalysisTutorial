echo STARTING

echo $1 
echo $2

./ExRootSTDHEPConverter $1 temp.root
./NTuple temp.root $2

rm temp.root

echo ENDING