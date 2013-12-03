test:
	python3.3 tests/test_all.py
package:
	python3.3 setup.py sdist

html_doc:
	cd docs && make html && cd ..

show_html_docs:
	firefox docs/build/html/index.html &

