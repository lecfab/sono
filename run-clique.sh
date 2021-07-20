#!/bin/bash
# ./run-clique.sh ../datasets/net_friendster.txt friendster 3

order=tmp-$2.txt
if [ -f "$order" ]; then
    echo "Optimised dpm already exists."
else
    echo "Optimising dpm"
    ./ord $1 trianglesdpm -u -o $order
fi
echo
echo "Test with dpm"
echo "-------------------------"
gcc ../kCliqueListing/Degen/kClist.c -O9 -o ../kCliqueListing/Degen/kClist
../kCliqueListing/Degen/kClist $3 $1 $order

echo
echo "Test with core"
echo "-------------------------"
../kCliqueListing/Degen/kClist $3 $1

echo
echo "Test with core-colour"
echo "-------------------------"
gcc ../kCliqueListing/DegenCol/DegenCol.c -O9 -o ../kCliqueListing/DegenCol/DegenCol
../kCliqueListing/DegenCol/DegenCol $3 $1

echo
echo "Test with core-colour-improved"
echo "-------------------------"
gcc ../kCliqueListing/DDegCol/DDegCol.c -O9 -o ../kCliqueListing/DDegCol/DDegCol
../kCliqueListing/DDegCol/DDegCol $3 $1

echo
echo "Test with degree"
echo "-------------------------"
gcc ../kCliqueListing/Degree/Degree.c -O9 -o ../kCliqueListing/Degree/Degree
../kCliqueListing/Degree/Degree $3 $1

echo
echo "Test with degree-colour"
echo "-------------------------"
gcc ../kCliqueListing/DegCol/DegCol.c -O9 -o ../kCliqueListing/DegCol/DegCol
../kCliqueListing/DegCol/DegCol $3 $1

echo
echo "Test with degree-colour-improved"
echo "-------------------------"
gcc ../kCliqueListing/DDegree/DDegree.c -O9 -o ../kCliqueListing/DDegree/DDegree
../kCliqueListing/DDegree/DDegree $3 $1

echo
echo "Test with arboricity: eliminated because it is too slow"
echo "-------------------------"
# gcc ../kCliqueListing/Arboricity/Arboricity.c -O9 -o ../kCliqueListing/Arboricity/Arboricity
# ../kCliqueListing/Arboricity/Arboricity $3 $1

echo

if [ $3 == 3 ]; then
    echo
    echo "Test with LDegree"
    echo "-------------------------"
    g++ ../kCliqueListing/LD/LwithOrdering.cpp -std=c++11 -O3 -o ../kCliqueListing/LD/LD
    ../kCliqueListing/LD/LD $3 $1 1
    echo
    echo "Test with LDegen"
    echo "-------------------------"
    ../kCliqueListing/LD/LD $3 $1 2
    # echo
    # echo "Test with LDPM"
    # echo "------------------------- TODO -------------------------"


    echo "Test with my algos"
    echo "-------------------------"
    orderedges=tmp-edges-$2.txt
    if [ -f "$orderedges" ]; then
        echo "Optimised dpm ORDER already exists."
    else
        echo "Creating dpm edgelist"
        ./rankedges $1 $order $orderedges -u
    fi
    echo "-------------------------"
    ./alg $orderedges -a triangles -u
fi
