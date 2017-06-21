all:
	python setup.py build_ext -i

clean:
	rm -f itm_wrap.c* itm.cpp itm*.so
	rm -rf build __pycache__

fresh: clean all

.PHONY: all clean fresh
