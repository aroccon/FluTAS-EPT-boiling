grep 'Max T_wall:    ' log.out > temps.txt
source ../env/bin/activate
python3 tempEvolution.py
deactivate
