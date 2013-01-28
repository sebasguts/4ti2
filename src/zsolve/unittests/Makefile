test_sources = \
	ValueTree.cpp \
	Vector.cpp \
	VectorArray.cpp \
	testTest.cpp

testMain: testMain.cpp $(test_sources)
	g++ -pthread -lgtest -I"../.." -std=c++0x $(test_sources) $< -o $@

clean:
	rm -f testMain

