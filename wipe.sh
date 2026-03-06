rm -rf ~/miniconda3/envs/topiary-hell/lib/python3.14/site-packages/topiary* 
rm -rf build/ 
for x in `find . -iname "__pycache__"`; do rm -rf $x; done
