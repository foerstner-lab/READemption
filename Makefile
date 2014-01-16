test:
	python3.3 tests/test_all.py

coverage:
	coverage3 run tests/test_all.py
	coverage3 report

package:
	python3.3 setup.py sdist

html_doc:
	cd docs && make html && cd ..

show_html_docs:
	firefox docs/build/html/index.html &

readme_txt:
	pandoc --from=markdown --to=plain README.md -o README.tex

readme_html:
	pandoc --from=markdown --to=html README.md -o README.html

readme_rst:
	pandoc --from=markdown --to=rst README.md -o README.rst

readme_clean:
	rm -f README.tex README.html README.rst
