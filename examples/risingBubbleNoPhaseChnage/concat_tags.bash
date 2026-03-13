source ~/FluTAS/run/env/bin/activate
cp log.out data
mkdir data/post/concatenate
python3 concatenate_tags.py  > ./data/post/concatenate/concat.dat
deactivate
