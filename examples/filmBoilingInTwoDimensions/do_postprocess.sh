cp log.out data
source /home/weber/FluTAS/run/env/bin/activate
python3 extract_Nusselt_numbers.py
python3 extract_VTK_data.py
python3 plotTheData.py
deactivate
