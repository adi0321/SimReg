# comment
# (note: the <tab> in the command line is necessary for make to work) 
#target:  dependency1 dependency2 ...
#      <tab> command

all: simreg.cpp
	g++ -fopenmp -o simreg simreg.cpp
	g++ -std=c++0x ./lib/compute-freq-42t.cpp ./include/current_time.cpp ./include/norm.cpp ./include/print.cpp ./include/prepro_extract_read_classes_cc.cpp ./include/load_obsCounts.cpp -o ./lib/compute-freq-42t -I /usr/local/include/ -L /usr/local/lib/
	g++ -std=c++0x ./lib/compute-freq.cpp ./include/current_time.cpp ./include/norm.cpp ./include/print.cpp ./include/prepro_extract_read_classes_cc.cpp ./include/load_obsCounts.cpp -o ./lib/compute-freq -I /usr/local/include/ -L /usr/local/lib/
	g++ -std=c++0x ./lib/compute_obsRC_CC.cpp ./include/current_time.cpp ./include/norm.cpp ./include/print.cpp ./include/extract_obsRC.cpp -o ./lib/compute_obsRC_CC -I /usr/local/include/ -L /usr/local/lib/
	g++ -std=c++0x ./lib/compute_sRC_d.cpp ./include/current_time.cpp ./include/norm.cpp ./include/print.cpp ./include/prepro_extract_read_classes_cc.cpp ./include/extract_obsCounts.cpp -o ./lib/compute_sRC_d -I /usr/local/include/ -L /usr/local/lib/
	g++ -std=c++0x ./utils/subSampling.cpp ./include/current_time.cpp ./include/loadSam.cpp -o ./utils/subSampling

clean:
	$(RM) simreg

