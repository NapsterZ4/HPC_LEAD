all: shear

shear:
	g++ -openmp shear.cpp -o shear
	g++ -w shear.cpp -o shear_serial

clean:
	rm -f shear shear_serial
