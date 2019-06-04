lib.so: lib.cpp
	g++ -O3 -fPIC -shared -Wall -Wpedantic lib.cpp -o lib.so

synchrotron.so: synchrotron.cpp
	g++ -O3 -fPIC -shared -Wall -Wpedantic synchrotron.cpp -o synchrotron.so

readme.html: README.md
	markdown README.md > readme.html

circletestlib.so: circletestlib.cpp
	g++ -O3 -fPIC -shared -Wall -Wpedantic circletestlib.cpp -o circletestlib.so
