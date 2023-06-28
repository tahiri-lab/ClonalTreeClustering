CXX = g++
CXXFLAGS = -std=c++11

TARGET = readLineage

$(TARGET): readLineage.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

.PHONY: clean
clean:
	rm -f $(TARGET)