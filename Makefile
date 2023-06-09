CC = upcxx
CFLAGS = -O2 -std=c++17
LDFLAGS = -L/opt/nfs/mpich-3.2/lib -L/opt/nfs/mpe2-2.4.9b/lib -lmpe -lX11 -L/opt/nfs/mpich-3.2/bin/lib -lmpicxx -Wl,-rpath -Wl,/opt/nfs/mpich-3.2/lib -Wl,--enable-new-dtags -lmpi
TARGET = ising.out
SRC = main.cpp utils.cpp

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $(TARGET) $(SRC)

run: $(TARGET)
	source /opt/nfs/config/source_upcxx_2023.3.sh && \
		UPCXX_GASNET_CONDUIT=udp upcxx-run -shared-heap 256M -n 4 $(shell /opt/nfs/config/station204_name_list.sh 1 16) ./$(TARGET)

clean:
	rm -f $(TARGET)