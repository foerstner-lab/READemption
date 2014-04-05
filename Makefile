test:
	python3.3 tests/test_all.py

coverage:
	coverage3 run tests/test_all.py
	coverage3 report

package:
	python3 setup.py sdist

package_to_pypi:
	python3 setup.py sdist upload
	@echo "Go to https://pypi.python.org/pypi/READemption/"

html_doc:
	cd docs && make html && cd ..

upload_doc:
	cd docs/build/html/ && zip -r READemption_docs.zip * && cd ../../.. && mv docs/build/html/READemption_docs.zip .
	@echo "Upload at https://pypi.python.org/pypi?%3Aaction=pkg_edit&name=READemption"

show_html_docs:
	firefox docs/build/html/index.html &

readme_txt:
	pandoc --from=markdown --to=plain README.md -o README.txt

readme_html:
	pandoc --from=markdown --to=html README.md -o README.html

readme_rst:
	grep -v "^\[!" README.md | sed -e "1d" > README.md.tmp
	pandoc --from=markdown --to=rst README.md.tmp -o README.rst
	rm README.md.tmp

readme_clean:
	rm -f README.tex README.html README.rst
	rm -f README.tex README.html README.txt

pylint:
	pylint bin/reademption reademptionlib/* tests/*

new_release:
	@echo "* Please do this manually:"
	@echo "* ------------------------"
	@echo "* Create/checkout a release branch"
	@echo "* Change bin/reademption"
	@echo "* Change setup.py"
	@echo "* Change docs/source/conf.py"
	@echo "* Change CHANGELOG.txt"
	@echo "* Create new docs"
	@echo "* Test package creation"
	@echo "* Test doc creation"
	@echo "* Commit changes e.g. 'git commit -m 'Set version to 0.2.X'"	
	@echo "* Tag the commit e.g. 'git tag -a v0.2.X -m 'version v0.2.X''"
	@echo "* Merge release into dev and master"
	@echo "* After pushing generate a new release based on this tag at"
	@echo "  https://github.com/konrad/READemption/releases/new"
