test_sources = \
	ValueTree.cpp \
	testTest.cpp

testMain: testMain.cpp $(test_sources)
	g++ -pthread -lgtest -std=c++0x $(test_sources) $< -o $@

clean:
	rm -f testMain

