echo "Installing Python package numpy..."
pip install numpy
echo "Installing Python package cython..."
pip install cython
echo "Installing Python package h5py..."
pip install h5py
echo "Cloning pbcore from https://github.com/PacificBiosciences/pbcore.git"
git clone https://github.com/PacificBiosciences/pbcore.git
cd pbcore
echo "Installing pbcore..."
python setup.py install
cd ../
echo "Installing SMALR..."
python setup.py install
