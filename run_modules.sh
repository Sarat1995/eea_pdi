cd VibroNISE
make all
cd runs/n-ph-pdi
../../bin/VibroNISE InputDirect1D.dat
cd ../tp-pdi
../../bin/VibroNISE InputDirect1D.dat
cd ../../..

cd eea_theory
python bands.py 0
cd ..

cd trpl
python fit_trpl.py 030.csv -D -a -f 1
cd ..
