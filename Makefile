simulator:
	g++ -o simulator 3dvibe_main.cpp memory.cpp ThreeD.cpp random.cpp 

rdf:
	g++ -o rdf rdf_main.cpp pair.cpp memory.cpp random.cpp 
