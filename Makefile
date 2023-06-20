NODES_FILE = nodes
CC = upcxx
CFLAGS = -std=c++17 -O2
LDFLAGS = -lstdc++ -lupcxx
TARGET = ising.out
SRC = main.cpp utils.cpp

$(TARGET): $(SRC)
	UPCXX_GASNET_CONDUIT=udp $(CC) $(CFLAGS) $(LDFLAGS) -o $(TARGET) $(SRC)

# -----------------    Run    -----------------
run: $(TARGET)
	UPCXX_GASNET_CONDUIT=udp UPCXX_NODES=$(NODES_FILE) upcxx-run -shared-heap 256M -n 16 ./$(TARGET)

# -----------------    Clean    -----------------
clean:
	rm -f $(TARGET)