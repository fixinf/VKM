CXXFLAGS =	-O2 -g -Wall -fmessage-length=0 `gsl-config --cflags` 

OBJS =		VKM.o EoS.o set_const.o fitting.o

LIBS =		`gsl-config --libs` 

TARGET =	VKM

$(TARGET):	$(OBJS)
	$(CXX) -static -o $(TARGET) $(OBJS) $(LIBS)

all:	$(TARGET)

EoS.o: EoS.cpp
		$(CXX) -c EoS.cpp -o EoS.o
		
set_const.o: set_const.cpp
		$(CXX) -c set_const.cpp -o set_const.o

fitting.o : fitting.cpp
		$(CXX) -c fitting.cpp -o fitting.o
clean:
	rm -f $(OBJS) $(TARGET)
