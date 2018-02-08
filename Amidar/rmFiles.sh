cd log
shopt -s extglob
rm !(gps*)
cd ..
./graph.sh
cd log
shopt -s extglob
rm !(*pdf)
