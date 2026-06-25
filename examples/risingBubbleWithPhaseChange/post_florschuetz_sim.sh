echo 'Warning: Create a virtual enviornment for Python and link it in this script before execution!'
source ../env/bin/activate
cp log.out data
python3 post_florschuetz_sim.py
deactivate
