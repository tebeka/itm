all:
	python setup.py build_ext -i

clean:
	rm -f itm_wrap.c* itm.py _itm.cpython-36m-x86_64-linux-gnu.so
	rm -rf build __pycache__

fresh: clean all

.PHONY: all clean fresh
