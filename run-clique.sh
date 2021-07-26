#!/bin/bash
# ./run-clique.sh ../datasets/net_friendster.txt friendster 3

if [ -z "${1}" ] || [ -z "${2}" ] || [ -z "${3}" ]; then
    echo "Use format ./run-clique.sh DATASET NAME K"
    exit
fi

orderdpm=../tmp/$2-dpm-order.txt
edgesdpm=../tmp/$2-dpm-edges.txt
orderdpm2=../tmp/$2-dpm2-order.txt
edgesdpm2=../tmp/$2-dpm2-edges.txt

echo "- - - - - - - - - - - - -"
if [ -f "$orderdpm" ]; then
    echo "Order dpm already exists."
else
    echo "Optimising dpm"
    ./ord $1 trianglesdpm -u -o $orderdpm
fi
if [ -f "$edgesdpm" ]; then
    echo "Edges dpm already exists."
else
    echo "Creating dpm edgelist"
    ./rankedges $1 $orderdpm $edgesdpm -u
fi

if [ -f "$orderdpm2" ]; then
    echo "Order dpm already exists."
else
    echo "Optimising dpm"
    ./ord $1 trianglesdpm2 -u -o $orderdpm2
fi
if [ -f "$edgesdpm2" ]; then
    echo "Edges dpm already exists."
else
    echo "Creating dpm edgelist"
    ./rankedges $1 $orderdpm2 $edgesdpm2 -u
fi
echo "- - - - - - - - - - - - -"

echo
echo "Test with dpm"
echo "-------------------------"
gcc ../kCliqueListing/Degen/kClist.c -O9 -o ../kCliqueListing/Degen/kClist
../kCliqueListing/Degen/kClist $3 $1 $orderdpm

echo
echo "Test with dpm2"
echo "-------------------------"
gcc ../kCliqueListing/Degen/kClist.c -O9 -o ../kCliqueListing/Degen/kClist
../kCliqueListing/Degen/kClist $3 $1 $orderdpm2

# echo
# echo "Test with core"
# echo "-------------------------"
# ../kCliqueListing/Degen/kClist $3 $1

# echo
# echo "Test with core-colour"
# echo "-------------------------"
# gcc ../kCliqueListing/DegenCol/DegenCol.c -O9 -o ../kCliqueListing/DegenCol/DegenCol
# ../kCliqueListing/DegenCol/DegenCol $3 $1
#
# echo
# echo "Test with core-colour-improved"
# echo "-------------------------"
# gcc ../kCliqueListing/DDegCol/DDegCol.c -O9 -o ../kCliqueListing/DDegCol/DDegCol
# ../kCliqueListing/DDegCol/DDegCol $3 $1

# echo
# echo "Test with degree"
# echo "-------------------------"
# gcc ../kCliqueListing/Degree/Degree.c -O9 -o ../kCliqueListing/Degree/Degree
# ../kCliqueListing/Degree/Degree $3 $1
#
# echo
# echo "Test with degree-colour"
# echo "-------------------------"
# gcc ../kCliqueListing/DegCol/DegCol.c -O9 -o ../kCliqueListing/DegCol/DegCol
# ../kCliqueListing/DegCol/DegCol $3 $1
#
# echo
# echo "Test with degree-colour-improved"
# echo "-------------------------"
# gcc ../kCliqueListing/DDegree/DDegree.c -O9 -o ../kCliqueListing/DDegree/DDegree
# ../kCliqueListing/DDegree/DDegree $3 $1

# echo
# echo "Test with arboricity: eliminated because it is too slow"
# echo "-------------------------"
# gcc ../kCliqueListing/Arboricity/Arboricity.c -O9 -o ../kCliqueListing/Arboricity/Arboricity
# ../kCliqueListing/Arboricity/Arboricity $3 $1



echo
echo "Test with kiteratif (original)"
echo "-------------------------"
./alg $1 -a cliques -u -n $3

echo
echo "Test with kiteratif (dpm)"
echo "-------------------------"
./alg $edgesdpm -a cliques -u -n $3

echo
echo "Test with kiteratif (dpm2)"
echo "-------------------------"
./alg $edgesdpm2 -a cliques -u -n $3

echo

if [ $3 == 3 ]; then
    # echo
    # echo "Test with LDegree"
    # echo "-------------------------"
    # g++ ../kCliqueListing/LD/LwithOrdering.cpp -std=c++11 -O3 -o ../kCliqueListing/LD/LD
    # ../kCliqueListing/LD/LD $3 $1 1
    # echo
    # echo "Test with LDegen"
    # echo "-------------------------"
    # ../kCliqueListing/LD/LD $3 $1 2
    # echo
    # echo "Test with LDPM"
    # echo "------------------------- TODO -------------------------"


    echo "Test with my algos"
    echo "-------------------------"
    ./alg $edgesdpm -a triangles -u
fi
