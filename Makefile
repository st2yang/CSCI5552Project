CC = g++
CFLAGS = -Wall -fPIC -I/usr/local/Aria/include -L/usr/local/Aria/lib -lAria -lpthread -ldl -lrt -std=c++11

OUTPATH = ./build
MAINOUT = slam_dunk

MAININ = main.cpp
ODOMIN = kalman-filter.cpp
MOVEIN = movementcontroller.cpp
FEATIN = LineExtraction_RANSAC.cpp line-match.cpp

all: slam

slam: $(MAININ) $(ODOMIN) $(MOVEIN) $(FEATIN)
	mkdir -p $(OUTPATH)
	$(CC) -o $(OUTPATH)/$(MAINOUT) $(MOVEIN) $(FEATIN) $(ODOMIN) $(MAININ) $(CFLAGS)

clean:
	rm -rf $(OUTPATH)/$(MAINOUT)
	rm -rf $(OUTPATH)
